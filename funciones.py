# Disponible en GitHub: https://github.com/vilcagamarracf/Inv_arroz_METRIC

# Importe de librerías

import math
import os

import ee

import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

from geemap import cartoee

# Inicio de Earth Engine
ee.Initialize()

# Funciones
def convert_RAW_to_TOA(img_raw):

    """Conversión de imágenes RAW a Radiancia y Reflectancia al tope de la atmósfera (TOA)
    Ref: https://www.usgs.gov/landsat-missions/using-usgs-landsat-level-1-data-product

    Reflectancia TOA: B, G, R, NIR, SWIR1, SWIR2
    Radiancia TOA: B10

    Parameters
    ----------
    img_raw : ee.Image 
        Imagen Landsat Raw (no es necesario aplicar factor de escala 0.00001)

    Returns
    -------
    img_toa : ee.Image
        Imagen procesada a TOA con bandas 2 al 7 y 10 
    """

    # Conversion to TOA Reflectance for Bands R, G, B, Nir

    # Constante 'Local sun elevation angle' convertido a radianes (ee.Image.sin necesita radianes)
    sun_elev = ee.Number(img_raw.get('SUN_ELEVATION')).multiply(math.pi /180)
    
    # Procesar las bandas seleccionadas
    for i in range(2,7+1): # B, G, R, NIR, SWIR1, SWIR2
        reflectance_mult_band_x = ee.Number(img_raw.get(f'REFLECTANCE_MULT_BAND_{i}'))
        reflectance_add_band_x =  ee.Number(img_raw.get(f'REFLECTANCE_ADD_BAND_{i}'))
        band = (
            img_raw.select(f'B{i}')
            .multiply(reflectance_mult_band_x)
            .add(reflectance_add_band_x)
            .divide(sun_elev.sin()) # Corrección por ángulo solar de elev
            )
        img_raw = img_raw.addBands(band, None, True)

    # Conversion to TOA Radiance for Thermal Band B10 
    # RADIANCE_MULT_BAND_10 & 11 = 3.3420E-04
    # RADIANCE_ADD_BAND_10  & 11 = 0.1
    b10 = img_raw.select(['B10']).multiply(3.3420E-04).add(0.1) # TOA spectral radiance

    # Agregar bandas a imagen final
    bandas = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10']
    img_toa = img_raw.addBands(b10, None, True).select(bandas)

    return img_toa

# ========================================================
# Índices de vegetación aplicables a Landsat 8 y 9 TOA
# - NDVI, SAVI, LAI y EVI
# ========================================================

def get_ndvi_L8(img_toa):

    """ Indice de vegetación normalizada NDVI
    Ref: https://www.usgs.gov/landsat-missions/landsat-normalized-difference-vegetation-index
    """

    img_ndvi = img_toa.expression(
        '(nir - red) / (nir + red)',
        {'red': img_toa.select('B4'),
         'nir': img_toa.select('B5')}
        ).rename('NDVI')

    return img_ndvi

def get_savi_L8(img_toa, L=None):

    """ Indice de vegetación ajustada al suelo
    Ref: https://www.usgs.gov/landsat-missions/landsat-soil-adjusted-vegetation-index
    """

    if L is None:
        L = 0.5
        
    img_savi = img_toa.expression(
        '(nir - red) / (nir + red + L) * (1+L)',
        {'red': img_toa.select('B4'),
         'nir': img_toa.select('B5'),
         'L'  : ee.Number(L)}
        ).rename('SAVI')

    return img_savi

def get_lai_L8(img_savi):

    """ Índice de área foliar
    Variable útil para caracterizar la dinámica, productividad y requerimientos hídricos de 
    cultivos.
    """

    img_lai = img_savi.expression(
        '- (log((0.69 - savi)/0.59))/(0.91)',
        {'savi': img_savi.where(img_savi.gt(0.689), 0.689)} # Evitar huecos
        ).rename('LAI')

    return img_lai

def get_evi_L8(img_toa):

  img_evi = img_toa.expression(
      '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)',
      {'red' : img_toa.select('B4'),
       'nir' : img_toa.select('B5'),
       'blue': img_toa.select('B2')}
       ).rename('EVI')

  return img_evi


def get_emisividades(img_ndvi, img_lai):

    """Obtener Emisividad y Transmisividad (Tasumi, 2003)
    
    - e_nb: Narrow band transmissivity
    - e_0 : broad-band surface emissivity
    
    Argumentos:
    - img_ndvi, img_lai (ee.Image): Indices a partir de una imagen 
        calibrada a Reflectancia a tope de la atmósfera (TOA)

    Retorna:
    - e_nb, e_0 (ee.Image) : dos imágenes de emisividad 
    """

    # e_0
    e_0 = img_lai.expression('0.95 + 0.01*LAI', 
                             {'LAI': img_lai}
                             ).rename('e_0') # LAI <= 3 y NDVI > 0

    e_0 = e_0.where(img_ndvi.lt(0), 0.985) # NDVI < 0
    e_0 = e_0.where(img_lai.gt(3), 0.98).rename('e_0') # LAI > 3; NDVI > 0

    # e_NB
    e_NB = img_lai.expression('0.97 + 0.0033*LAI',
                              {'LAI':img_lai}).rename('e_NB') # LAI <= 3 y NDVI > 0

    e_NB = e_NB.where(img_ndvi.lt(0), 0.985) # NDVI < 0
    e_NB = e_NB.where(img_lai.gt(3).And(img_ndvi.gt(0)), 0.98) # LAI > 3; NDVI > 0

    return e_NB, e_0


def get_decl_lat_hra(img, roi, doy): 

    """Obtener ángulos necesarios para estimar theta

    Inputs:
    - img : ee.Image
    - roi : ee.Geometry
    - doy : número o ee.Number

    Retorna (en radianes):
    - angl_decl : Ángulo de declinación terrestre (ee.Number)
    - latitud   : Mapa de latitudes (ee.Image)
    - angle_hra : Mapa de ángulos horarios (ee.Image)
    
    """

    # Datos a partir de imagen satelital
    # doy = img.date().getRelative('day', 'year') # ee.Number 
    proj = img.select([0]).projection() # EPSG:32617 -> WGS 84 / UTM zone 17N

    img_lonlat = ee.Image.pixelLonLat().reproject(proj).clip(roi) # ee.Image con bandas 'Longitude' y 'Latitude'

    factor_rad = ee.Number(math.pi/180) # Factor de conversión a radianes


    # Declinacion solar δ [rad] : ee.Number 
    angl_decl = ee.Number.expression('-23.45*cos( 360*(d+10)/365 * factor )', 
                                     {'d':doy, 'factor':factor_rad}).multiply(factor_rad) 

    # Latitud ϕ [rad] : ee.Image
    latitud = img_lonlat.select('latitude').multiply(factor_rad)


    # Angulo horario [rad] : ee.Image
    B   = ee.Number.expression('360/365*(doy-81)', 
                                {'doy': doy}).multiply(factor_rad) # radianes

    EoT = ee.Number.expression('9.87*sin(2*B) - 7.53*cos(B) - 1.5*sin(B)', 
                                {'B': B}) # grados 
    
    longitud = img_lonlat.select('longitude') # grados 

    # LSTM = -75 = 15° * (Diferencia entre hora local y GMT = -5) 
    tc  = longitud.expression('4*(long-(-75)) + EoT',
                              {'EoT': EoT, 'long':longitud}) # radianes 

    lt  = img.date().getFraction('day', 'America/Lima').multiply(24)#.getInfo() # 10.4712 aprox

    # Hora solar local
    lst = tc.divide(60).add(lt)

    angle_hra = lst.subtract(12).multiply(15).multiply(factor_rad).rename('angle_hra') # radianes


    return angl_decl, latitud, angle_hra




def get_cos_theta(angl_decl, latitud, angle_hra, slope_rad, aspect_rad):

    """Obtener el coseno del ángulo theta
    
    Requiere: 
    - Función `get_decl_lat_hra(img, roi, doy)`

    Args:
    - angl_decl, latitud, angle_hra (de la función)
    - slope_rad  (ee.Image) : Mapa de pendientes en radianes 
    - aspect_rad (ee.Image) : Mapa de aspecto en radianes

    Retorna:
    - cos_theta (ee.Image) : Mapa de Cosenos del ángulo incidencia solar 
    
    """
    # 1. Mountain Model
    # Duffie and Beckman (1991)
    cos_theta = latitud.expression(
        """
        sin(delta)*sin(phi)*cos(s) 
        - sin(delta)*cos(phi)*sin(s)*cos(gamma)
        + cos(delta)*cos(phi)*cos(s)*cos(omega)
        + cos(delta)*sin(phi)*sin(s)*cos(gamma)*cos(omega)
        + cos(delta)*sin(gamma)*sin(s)*sin(omega)
        """,
        {
            'delta': angl_decl,
            'phi': latitud,
            's': slope_rad,
            'gamma': aspect_rad,
            'omega': angle_hra
        }
    ).rename('cos_theta') # ee.Image

    return cos_theta


def get_surface_temp(img_toa, e_nb):

    """Obtener mapa de temperatura de superficie en Kelvin

    Inputs:
    - img_toa (ee.Image): Imagen calibrada al tope de la atmósfera TOA
    - e_nb    (ee.Image): Mapa de emisividades de banda estrecha 

    Retorna: 
    - ee.Image con banda 'Ts_k'
    
    """

    # Factor de correción Rc dónde: Rp = 0.91, Tnb = 0.866, Rsky = 1.32
    rc = img_toa.expression('(Lt10 - 0.91)/0.866 - (1-e_nb)*1.32',
                            {'Lt10':img_toa.select('B10'), 'e_nb':e_nb})
    
    # Termperatura de superficie: b10_k1 = 774.89, b10_k2 = 1321.08
    ts_k = rc.expression('1321.08/( log( (e_nb*774.89/rc) + 1 ) )',
                         {'e_nb':e_nb, 'rc':rc}).rename('Ts_k')
    
    return ts_k



def convert_TOA_to_SR(img_toa, P_air, w, cos_theta_hor):

    """Conversión de imagen TOA a Reflectancia de Superficie
    
    Descripción
    -----------
    En METRIC se deriva la reflectancia de superficie bidireccional 
    utilizando funciones calibradas de transmitancia atmosférica 
    y reflectancia de la trayectoria por Tasumi et al. (2007).

    Se corrige los valores de las bandas de la imagen TOA para la 
    dispersión y absorción de la radiación solar entrante y 
    reflejada desde la superficie, basándose en una función de 
    corrección atmosférica simplificada que sólo requiere mediciones
    puntuales o estimaciones de la presión de vapor cerca
    de la superficie. 

    Parametros
    ----------
    img_toa : ee.Image
    P_air : ee.Image 
        Presión del aire
    w : ee.Image
        Mapa de probabilidad de lluvias
    cos_theta_hor : ee.Image
        Mapa del coseno del ángulo theta

    Retorna
    -------
    img_sr : ee.Image
        ee.Image con bandas de Reflectancia de Superficie con bandas
        ['B2_SR','B3_SR','B4_SR','B5_SR','B6_SR','B7_SR']

    """

    bands = ['B2','B3','B4','B5','B6','B7']

    img_toa = img_toa.select(bands)
    
    for band in bands:
        if band == 'B2':
            C1 = 0.987
            C2 = -0.00071
            C3 = 0.000036
            C4 = 0.0880
            C5 = 0.0789
            Cb = 0.640

        elif band == 'B3':
            C1 = 2.319
            C2 = -0.00016
            C3 = 0.000105
            C4 = 0.0437
            C5 = -1.2697
            Cb = 0.310

        elif band == 'B4':
            C1 = 0.951
            C2 = -0.00033
            C3 = 0.00028
            C4 = 0.0875
            C5 = 0.1014
            Cb = 0.286

        elif band == 'B5':
            C1 = 0.375
            C2 = -0.00048
            C3 = 0.005018
            C4 = 0.1355
            C5 = 0.6621
            Cb = 0.189

        elif band == 'B6':
            C1 = 0.234
            C2 = -0.00101
            C3 = 0.004336
            C4 = 0.0560
            C5 = 0.7757
            Cb = 0.274

        elif band == 'B7':
            C1 = 0.365
            C2 = -0.00097
            C3 = 0.004296
            C4 = 0.0155
            C5 = 0.639
            Cb = -0.186

        t_in_b = img_toa.expression(
            """C1*exp( C2*P_air/(K_t*cos(theta_hor))
            - (C3*W+C4)/(cos(theta_hor)) )+C5""",
            {'C1': C1, 'C2': C2, 'C3': C3, 'C4': C4, 'C5': C5, 
             'P_air' : P_air,
             'K_t': 1,
             'W': w,
             'theta_hor': cos_theta_hor}
            ).rename(f't_in_{band}')

        t_out_b = img_toa.expression(
            """C1*exp( C2*P_air/(K_t*cos(theta_hor))
            - (C3*W+C4)/(cos(theta_hor)) )+C5""",
            {'C1': C1, 'C2': C2, 'C3': C3, 'C4': C4, 'C5': C5, 
             'P_air' : P_air,
             'K_t': 1,
             'W': w,
             'theta_hor': 1}
            ).rename(f't_out_{band}')

        R_in_s_b = t_in_b.multiply(t_out_b)

        R_out_s_b = img_toa.expression(
            'band - Cb*(1-t_in_b)', 
            {'Cb':Cb, 't_in_b': t_in_b, 'band':img_toa.select(band)}
        )

        p_s_b = R_out_s_b.divide(R_in_s_b)

        img_toa = img_toa.addBands(p_s_b, None, True)
    
    bands_sr = ['SR_B2', 
                'SR_B3', 
                'SR_B4', 
                'SR_B5', 
                'SR_B6', 
                'SR_B7']

    img_sr = img_toa.select(bands).rename(bands_sr)

    return img_sr


def get_albedo(img_sr, albedo_method):

    """Obtener el albedo de una imagen mediante una calibración con coeficientes de 
    ponderación por banda (Tasumi et al., 2008) https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51)

    Inputs:
    - img_sr : ee.Image
    - albedo_method: str

    Retorna:
    - albedo : ee.Image
    """

    # '0.300*B2 + 0.277*B3 + 0.233*B4 + 0.143*B5 + 0.036*B6 + 0.012*B7',  # Silva et al. (2016)
    # '0.254*B2 + 0.149*B3 + 0.147*B4 + 0.311*B5 + 0.103*B6 + 0.036*B7',  # Tasumi et al. (2007)
    # '0.3561*B2 + 0.3972*B3 + 0.3904*B4 + 0.6966*B5 + 0.2286*B6 + 0.1596*B7', # Javier
    # '(0.356*B2 + 0.130*B4 + 0.373*B5 + 0.085*B6 + 0.072*B7 - 0.018)/1.016',  # Javier 2
    
    albedo = img_sr.expression(
        albedo_method,
        {
            'B2' : img_sr.select('SR_B2'),
            'B3' : img_sr.select('SR_B3'),
            'B4' : img_sr.select('SR_B4'),
            'B5' : img_sr.select('SR_B5'),
            'B6' : img_sr.select('SR_B6'),
            'B7' : img_sr.select('SR_B7')        
        }
    ).rename('albedo')

    return albedo



def getRadiacionNeta(img_ee, roi, dem, lai_method, albedo_method, HR):

    """Obtener mapa de Radiación Neta

    Parametros
    ----------
    img_ee : ee.Image
    roi : ee.Geometry
    dem : ee.Image

    Retorna
    -------
    R_n : ee.Image
    img_sr_tasumi : ee.Image
    img_productos : ee.Image
    d2 : ee.Number
    doy : ee.Number
    """

    # =================================================================
    # Constantes
    # =================================================================
    
    # Factor de conversión de grados a radianes
    factor_rad = ee.Number(math.pi/180)

    # Fechas
    img_date = img_ee.date() # ee.Date
    doy = img_date.getRelative('day', 'year').add(1) # ee.Number
    fecha = img_date.format('YYYY-MM-dd').getInfo() # string

    # =================================================================
    # Procesamiento
    # =================================================================

    # Procesar imagen RAW a TOA
    img_toa = convert_RAW_to_TOA(img_ee).clip(roi) # ee.Image
    
    # Índices de vegetación (eq. 23, 19 y 18)
    img_ndvi = get_ndvi_L8(img_toa)
    
    if lai_method == 0:
        img_savi = get_savi_L8(img_toa, L=0.1)
        img_lai  = get_lai_L8(img_savi)
        
    if lai_method == 1:
        img_savi = get_savi_L8(img_toa, L=0.5)
        img_lai  = get_lai_L8(img_savi)
    
    # Obtener LAI mediante relación NDVI - IAF
    if lai_method == 2:
        img_savi = get_savi_L8(img_toa) 
        img_lai = img_ndvi.multiply(2.1362).add(0.0869).rename('LAI') # Usando 10 datos (1 fecha - excel de Martin)
        
    if lai_method == 3:
        img_savi = get_savi_L8(img_toa) 
        img_lai = img_ndvi.expression('(2.3523*img_ndvi)**2 - 1.9013*img_ndvi + 1.7714', 
                                      {'img_ndvi': img_ndvi.select('NDVI')}
                                      ).rename('LAI') # Usando 30 datos


    # img_lai = img_lai.where(img_lai.lte(0), 0)
    
    # A partir del DEM: Pendiente y Aspect [rad]
    img_dem_clippped = dem.clip(roi)
    img_slopes = ee.Terrain.slope(img_dem_clippped)  # grados
    img_slopes_rad = img_slopes.multiply(factor_rad) # ee.Image, radianes 
        
    img_aspect = ee.Terrain.aspect(img_dem_clippped) # grados
    img_aspect_rad = img_aspect.multiply(factor_rad) # ee.Image, radianes


    # Ángulos Declinación, Latitud y Horario [rad]
    angle_decl, latitud, angle_hra = get_decl_lat_hra(img_ee, roi, doy)


    # Parámetro: Emisividad e_0 (eq. 17)
    e_nb, e_0 = get_emisividades(img_ndvi, img_lai)

    # =================================================================
    # Parámetro: Radiación de onda corta entrante 
    # =================================================================
    
    # Requiere: t_sw, d2, cos_theta_rel 
    # t_sw requiere: P_air, w, cos_theta_hor
    # w requiere: e_s

    # 1. d2 INVERSE RELATIVE DISTANCE EARTH-SUN (eq. 9)
    d2 = ee.Number.expression('1/( 1+0.033*cos(doy*2*pi/365) )',
                              {'doy': doy, 'pi': math.pi})

    # 2. cos_theta_rel [rad] (eq. 7)
    cos_theta_rel = get_cos_theta(angle_decl, 
                                  latitud, 
                                  angle_hra, 
                                  img_slopes_rad, 
                                  img_aspect_rad)

    # 3. t_sw (eq. 4)
    # 3.1. P_air Atmospheric Pressure [kPa] (eq. 5)
    P_air = dem.expression('101.3*( (293-0.0065*z)/293 )**5.26', 
                           {'z': dem.select(0)}) # ee.Image

    # 3.2. cos_theta_hor (eq. 8)
    cos_theta_hor = get_cos_theta(angle_decl, 
                                  latitud, 
                                  angle_hra, 
                                  slope_rad=0, 
                                  aspect_rad=0)
    
    # 3.3. w Water in the atmosphere [mm]
    # 3.3.1. Temperatura de superficie (eq. 20)
    ts = get_surface_temp(img_toa, e_nb)      # [K]
    ts_c = ts.subtract(273.15).rename('Ts_c') # [°C]

    # 3.3.2. Near-surface vapor pressure ea [kPa]
    # Ojo temperatura se requiere en °C
    # e0_eq = '6.112*exp(17.67*t / (t + 243.5))' # https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
    e0_eq = '0.6108*exp(17.27*t / (t + 237.3))' # FAO 56, Allen (2012) https://www.fao.org/3/x0490s/x0490s01.pdf
    e0 = ts_c.expression(e0_eq, {'t': ts_c.select('Ts_c')})
    ea = e0.multiply(HR/100).rename('vapor_pressure')

    # 3.3.3. Agua precipitable en la atmósfera w [mm] Garrison & Adler (1990)
    w = ea.expression('0.14*ea*P + 2.1', {'ea':ea, 'P':P_air}).rename('w')
    
    # 3.3.4. t_sw: broad-band atmospheric transmissivity (eq. 4), ASCE-EWRI (2005)
    t_sw = P_air.expression('''
                            0.35 + 0.627*exp(
                            - 0.00146*P_air/(Kt*cos_theta_hor) 
                            - 0.075*(W/cos_theta_hor)**0.4
                            )
                            ''',
                            {'P_air':P_air, 'Kt':1, 'cos_theta_hor':cos_theta_hor, 'W':w}
                            ).rename('t_sw')
    
    # Finalmente: R_s_incoming (eq. 3)
    R_s_in = cos_theta_rel.expression('1367*cos_theta_rel*t_sw/d2',
                                      {'cos_theta_rel':cos_theta_rel, 't_sw':t_sw, 'd2':d2}
                                      ).rename('R_s_in')

    # =================================================================
    # Parámetro: Radiación de onda larga entrante
    # =================================================================

    # ea: effective atmospheric emissivity (eq. 25)
    atm_emissivity_ea = t_sw.expression('0.85*(- log(t_sw))**0.09',
                                        {'t_sw': t_sw}
                                        ).rename('atm_emissivity_ea')

    # Finalmente: R_l_in (eq. 24)
    R_l_in = atm_emissivity_ea.multiply(5.67E-08).multiply(ts.pow(4)).rename('R_l_in')

    # =================================================================
    # Parámetro: Radiación de onda larga saliente (eq.  16)
    # =================================================================

    R_l_out = e_0.multiply(5.67E-08).multiply(ts.pow(4)).rename('R_l_out')

    # =================================================================
    # Parámetro: Albedo
    # =================================================================

    # Corrección atmosférica - Tasumi et al. (2007) (eqs. 10 - 14)
    img_sr_tasumi = convert_TOA_to_SR(img_toa, P_air, w, cos_theta_hor)

    # Albedo (eq. 15)
    img_albedo = get_albedo(img_sr_tasumi, albedo_method)

    # Radiación Neta (eq. 2)
    Rn = img_albedo.expression(
        '(1-albedo)*R_s_in + (R_l_in - R_l_out) - (1-e_0)*R_l_in',
        {'albedo':img_albedo,
         'R_s_in':R_s_in,
         'R_l_in':R_l_in,
         'R_l_out':R_l_out,
         'e_0':e_0}
         ).rename('R_n')

    img_Rn = Rn.addBands([R_s_in, R_l_in, R_l_out])

    # Juntar los parámetros obtenidos en una sola imagen
    img_productos = ee.Image([
        img_ndvi.select('NDVI'),
        img_savi.select('SAVI'),
        img_lai.select('LAI'),
        img_albedo,
        ts,
        ts_c,
        t_sw,
        e_0,
        e_nb,
        cos_theta_rel,
        cos_theta_hor,
        img_slopes,
        dem
        ])

    img_productos_dict = {
        'img_Rn': img_Rn,
        'img_SR': img_sr_tasumi,
        'img_toa': img_toa,
        'img_productos': img_productos,
        'd2': d2,
        'doy': doy,
        'fecha': fecha
    }
        
    return img_productos_dict


# =====================================================
# Funciones de ayuda: Obtener estadísticas y Gráficas
# =====================================================

def get_stats(img, geometry, scale):

    ''' Obtener valores estadísticos de la imagen

    Descripción
    -----------
    Uso de ee.Image.reduceRegion para obtener estadísticas
    con reductores `.unweighted()`
    
    Parametros
    ----------
    img : ee.Image
    geometry : ee.Geometry
    scale : int
        tamaño de pixel
        
    Retorna
    -------
    stats : dict
        Diccionario que contiene las siguientes estadísticas:  
        Media, Mediana, Moda, Desviación estándar, Mínimo y Máximo (unweighted)
    '''
    
    values = ee.List([
        img.reduceRegion(ee.Reducer.mean().unweighted()  , geometry=geometry, scale=scale),
        img.reduceRegion(ee.Reducer.median().unweighted(), geometry=geometry, scale=scale),
        img.reduceRegion(ee.Reducer.mode().unweighted()  , geometry=geometry, scale=scale),
        img.reduceRegion(ee.Reducer.stdDev(), geometry=geometry, scale=scale),
        img.reduceRegion(ee.Reducer.min(), geometry=geometry, scale=scale),
        img.reduceRegion(ee.Reducer.max(), geometry=geometry, scale=scale)
    ])

    columns = ['mean', 'median', 'mode', 'stdDev', 'min', 'max']
                        
    stats = ee.Dictionary.fromLists(columns, values).getInfo()
    # stats_df = pd.DataFrame.from_dict(stats, orient='index')

    return stats


def get_grafica_cartoee_color(image, 
                              vis_params,                              
                              text=None, 
                              title_map=None,
                              label=None, 
                              save_fig=None,
                              nogrid=None):
    
    """Obtener gráficas con cartoee
    La variable zoom_region debe asignarse según la zona de estudio.
    Para modificaciones seguir: https://geemap.org/cartoee/#cartoee-module
    """

    # Establecer área de plot usando las coordenadas de predios_bound.coordinates().getInfo()
    extent = 0.0005
    zoom_region = [
        -79.77332525015065+extent, -6.605665317455976-extent,
        -79.78719338866794-extent, -6.594549+extent
    ]
    
    # Establecer figsize para plot
    # fig = plt.figure(figsize=(5,5)) # Para juntar en forma de mosaicos
    fig = plt.figure(figsize=(8,6)) # Para analizar
    # fig = plt.figure(figsize=(16,12)) # Para recortar la barra
    
    # ee.Image a plotear
    ax = cartoee.get_map(image, region=zoom_region, vis_params=vis_params)

    # Añadir grillas
    if nogrid is None:
        cartoee.add_gridlines(ax, 
                              interval=0.005, 
                              ytick_rotation=90, 
                              linestyle=":", 
                              linewidth=0 # Grillas invisibles
                              ) # xtick_rotation=45

    # Añadir barra de color
    if label is not None:
        cartoee.add_colorbar(ax, 
                             vis_params=vis_params, 
                             loc='right', 
                             # label=label, 
                             #posOpts=[0.86, 0.25, 0.02, 0.5],
                             tick_font_size=12
                             # ticks=[15, 20, 25, 30, 35], # LAI
                             # drawedges=True, 
                             # extend='both', # Genera flechas hacia los extremos
                             )
        ax.text(-79.77332525015065+2*extent, -6.5956, label, fontsize=12)

    # Añadir texto
    if title_map is not None:
        ax.set_title(title_map) # , fontsize=11
        
    if text is not None:
        ax.text(-79.7872, -6.6056, text, fontsize=12) 
        # fontsize=18 para mejor visibilidad en mosaicos
        
    # add scale bar
    scale_bar_dict = {
          "length": 100, 
          "xy": (0.9, 0.05), 
          "linewidth": 2,
          "fontsize": 12,
          "color": "black",
          "unit": "m",
          "ha": "center",
          "va": "bottom",    
    }
    cartoee.add_scale_bar_lite(ax, **scale_bar_dict)
    
    ax.tick_params(axis = 'both') # , labelsize = 9
    
    # Guardar graficas
    if save_fig is not None:

        ruta = r'C:/Users/usuario/Documents/00-notebooks-2022/output'
        ruta_img = os.path.join(ruta, save_fig + '.jpg')
        plt.savefig(ruta_img, bbox_inches = 'tight', pad_inches = .1)#, dpi=400)
        # Recortar márgenes con ayuda de bbox_inches y pad_inches

    plt.show()

# Para probar la función de elaboración de gráficos 

# image = predios_vis#.blend(predios_bound_vis)
# vis_params = {'color':'black'}
# title_map = 'Área de interés'

# get_grafica_cartoee_color(
#     image, 
#     vis_params, 
#     title_map=title_map, 
#     text='Hola mundo', 
#     save_fig='img_prueba_pad')



# ======================================================
# Flujo de calor sensible (H)
# - Funciones para la selección de pixeles Frío y Caliente
# ======================================================

def recorte_por_percentiles(n1, n2, img_ndvi, img_ts, geometry, filtrado=None):
    
    """ Recortar mapas de temperatura en función de mapa NDVI para obtener pixeles candidatos.

    Parametros
    ----------
    n1, n2 : int
    img_ndvi : ee.Image
    img_ts : ee.Image
    geometry : ee.Geometry
    filtrado : str

    Retorna
    -------
    pixeles : ee.Image
        Mapa recortado (por updateMask) de pixeles candidatos.
    """

    # 1. Procesar n1
    # Para NDVI > 0, retorna el percentil n1
    perc_n1 = ee.Number(img_ndvi
                        .updateMask(img_ndvi.gte(0))
                        .reduceRegion(ee.Reducer.percentile([n1]), geometry=geometry, scale=30)
                        .values()
        )
    
    if filtrado == 'cold':
        # Pixel frío
        # Para los valores de ndvi mayores al percentil perc_n1, retornar sus temperaturas
        ts_recortado_p_n1 = img_ts.updateMask(img_ndvi.gte(perc_n1)) # Temperatura
    else: 
        # Pixel caliente
        img_ndvi = img_ndvi.clip(geometry)
        img_ts = img_ts.clip(geometry)
        # Para los valores de ndvi menores al percentil perc_n1, retornar sus temperaturas
        ts_recortado_p_n1 = img_ts.updateMask(img_ndvi.lte(perc_n1))
        
    # 2. Procesar n2
    perc_n2 = ee.Number(ts_recortado_p_n1
                        .reduceRegion(ee.Reducer.percentile([n2]), geometry=geometry, scale=30)
                        .values()
        )
    
    if filtrado == 'cold':
        # Los valores más bajos de temperatura
        pixeles = ts_recortado_p_n1.updateMask(ts_recortado_p_n1.lte(perc_n2))
    else:
        # Los valores más altos de temperatura
        pixeles = ts_recortado_p_n1.updateMask(ts_recortado_p_n1.gte(perc_n2))
            
    return pixeles

def recorte_por_percentiles_alb(n1, n2, img_ndvi, img_albedo, geometry, filtrado=None):
    
    """ Recortar mapas de temperatura en función de mapa NDVI para obtener pixeles candidatos.

    Parametros
    ----------
    n1, n2 : int
    img_ndvi : ee.Image
    img_ts : ee.Image
    geometry : ee.Geometry
    filtrado : str

    Retorna
    -------
    pixeles : ee.Image
        Mapa recortado (por updateMask) de pixeles candidatos.
    """

    # 1. Procesar n1
    # Para NDVI > 0, retorna el percentil n1
    perc_n1 = ee.Number(img_ndvi
                        .updateMask(img_ndvi.gte(0))
                        .reduceRegion(ee.Reducer.percentile([n1]), geometry=geometry, scale=30)
                        .values()
        )
    
    if filtrado == 'cold':
        # Pixel frío
        # Para los valores de ndvi mayores al percentil perc_n1, retornar sus temperaturas
        alb_recortado_p_n1 = img_albedo.updateMask(img_ndvi.gte(perc_n1)) # Albedo
    else: 
        # Pixel caliente
        img_ndvi = img_ndvi.clip(geometry)
        img_albedo = img_albedo.clip(geometry)
        # Para los valores de ndvi menores al percentil perc_n1, retornar sus temperaturas
        alb_recortado_p_n1 = img_albedo.updateMask(img_ndvi.lte(perc_n1))
        
    # 2. Procesar n2
    perc_n2 = ee.Number(alb_recortado_p_n1
                        .reduceRegion(ee.Reducer.percentile([n2]), geometry=geometry, scale=30)
                        .values()
        )
    
    if filtrado == 'cold':
        # Los valores más bajos de temperatura
        pixeles = alb_recortado_p_n1.updateMask(alb_recortado_p_n1.lte(perc_n2))
    else:
        # Los valores más altos de temperatura
        pixeles = alb_recortado_p_n1.updateMask(alb_recortado_p_n1.gte(perc_n2))
            
    return pixeles


def get_pixel_values(punto, d2, img_productos):

    """ Obtener datos de un pixel mediante sus coordenadas.

    Parametros
    ----------
    punto : ee.Geomety.Point
        Punto seleccionado para extraer valores
    pixeles : ee.Image
        Mapa de pixeles candidatos
    d2 : float
        Parámetro en función del DOY necesario para Rn 
    img_productos : ee.Image
        Conjunto de imágenes en forma de bandas producto de función getRadiacionNeta
    dem : ee.Image
        Mapa de DEM SRTM

    Retorna
    -------
    pix_values : dict
        Diccionario con valores de 
        ['R_n', 'G', 'NDVI', 'LAI', 'Albedo', 'Ts', 'T_sw', 'e_0', 'Elev_m', 'Slope_d']
    """
    
    # Extraer datos de pixeles mediante reductor
    pix_datos_extraidos = img_productos.reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=punto,
        scale=30
    )
    
    # Asignar los datos extraidos a variables
    pix_ndvi  = ee.Number(pix_datos_extraidos.get('NDVI'))
    pix_savi  = ee.Number(pix_datos_extraidos.get('SAVI'))
    pix_lai   = ee.Number(pix_datos_extraidos.get('LAI'))
    
    pix_slope = ee.Number(pix_datos_extraidos.get('slope')) # Grados
    pix_elev  = ee.Number(pix_datos_extraidos.get('elevation')) # m
    
    pix_alb   = ee.Number(pix_datos_extraidos.get('albedo'))
    pix_ts_K  = ee.Number(pix_datos_extraidos.get('Ts_k'))
    pix_ts_C  = ee.Number(pix_datos_extraidos.get('Ts_c'))
    
    pix_theta = ee.Number(pix_datos_extraidos.get('cos_theta'))
    pix_t_sw  = ee.Number(pix_datos_extraidos.get('t_sw'))
    pix_e_0   = ee.Number(pix_datos_extraidos.get('e_0'))                
    
    # Procesar Radiación Neta Rn
    Rn_calculated = ee.Number.expression(
        """
        (1-pix_alb)*1367*pix_theta/d2*pix_t_sw 
        + pix_e_0*5.67*10**(-8)*( 0.85*(-pix_t_sw_log)**0.09 - 1 )*pix_ts_K**4
        """,
        {
            'pix_alb': pix_alb,
            'pix_theta': pix_theta,
            'd2': d2,
            'pix_t_sw': pix_t_sw,
            'pix_t_sw_log': pix_t_sw.log(),
            'pix_e_0': pix_e_0,
            'pix_ts_K': pix_ts_K 
        }
    )
    
    # Procesar flujo de calor del suelo G
    
    # Forma 1: Bastiananssen (1995) representing values near midday.
    # G_pix = ee.Number.expression(
    #     '(Ts_K - 273.15)*(0.0038 + 0.0074*albedo)*(1-0.98*pix_ndvi**4)*Rn_calculated',
    #     {
    #         'Ts_K': pix_ts_K,
    #         'albedo': pix_alb,
    #         'pix_ndvi': pix_ndvi,
    #         'Rn_calculated': Rn_calculated
    #     }
    # )
    
    # Forma 2: Tasumi (2003) using soil heat flux data collected by Wright (1982), 
    # USDA-ARS for irrigated crops near Kimberly, Idaho
    G_pix = ee.Number.expression(
        '''
        (pix_lai < 0.5) 
            ? 1.8*(pix_ts_K - 273.15) + 0.084*Rn_calculated 
            : (0.05 + 0.18*exp(-0.521*pix_lai))*Rn_calculated
        ''', 
        {
            'pix_lai': pix_lai,
            'pix_ts_K': pix_ts_K,
            'Rn_calculated': Rn_calculated
        })
    
    # Producir un diccionario para exportar
    pix_values_list = ee.List([Rn_calculated, 
                               G_pix, 
                               pix_ndvi, 
                               pix_savi,
                               pix_lai, 
                               pix_alb, 
                               pix_ts_K,
                               pix_ts_C,
                               pix_t_sw,
                               pix_e_0,
                               pix_elev,
                               pix_slope,
                               pix_theta])
    
    columns = ['R_n', 'G', 'NDVI', 'SAVI', 'LAI', 
               'Albedo', 'Ts_k', 'Ts_c',
               'T_sw', 'e_0', 'Elev_m', 'Slope_d', 'cos_theta'] 
                        
    pix_values = ee.Dictionary.fromLists(columns, pix_values_list).getInfo()
    
    return pix_values


# =========================================
# Parte Iterativa
# =========================================

# Funciones a emplear
# 0. Velocidad de fricción
def get_u_star(u_200, Z_om, correccion=None):
    
    """Velocidad de fricción [m s-1]

    Inputs:
    - u_200
    - Z_om
    - correccion (opcional): Aplicado durante iteraciones con valores de 
        correción por estabilidad y200m.

    Retorna:
    - u_star
    """
    
    k_const = 0.41 # von Karman's constant
    if correccion == None:
        u_star = k_const*u_200*1/(math.log(200/Z_om)) # blending height: 200 m 
    else: 
        y200m = correccion['y200m']
        u_star = k_const*u_200*1/(math.log(200/Z_om) - y200m)
    
    return u_star


# 1. Resistencia Aerodinámica
def get_rah(u_star, correccion=None):
    
    '''Resistencia aerodinámica (eq 30) [s/m]

    Inputs:
    - u_star
    - correccion (opcional): Aplicado durante iteraciones con valores de correción por 
      estabilidad yhz2 e yhz1

    Retorna:
    - r_ah 
    '''
    
    k_const = 0.41
    z2 = 2
    z1 = 0.1
    
    if correccion == None:
        r_ah = math.log(z2/z1)/(u_star*k_const)
    else:
        yhz1 = correccion['yhz1']
        yhz2 = correccion['yhz2']
        r_ah = (math.log(z2/z1) + yhz1 - yhz2)/(u_star*k_const)
    
    return r_ah


# 2. Densidad del aire
def get_air_dens(pix_elev, air_temp_K):
    
    '''Densidad del aire [kg m-3]

    Standard equations for mean atmospheric pressure and the universal gas law 
    and simplifying for the effect of vapor pressure 

    Inputs:
    - pix_elev
    - air_temp_K [K]

    Retorna:
    - air_dens_p
    '''
    
    # (eq 37) For 20 °C: 293.15° = 273.15° + 20°
    # P estimated using a digital elevation map
    air_pressure_P = 101.3*( (293.15-0.0065*pix_elev)/293.15)**5.26  # [kPa]
    air_dens_p = 1000*air_pressure_P/(1.01*(air_temp_K)*287) # [kg m-3]
                                      
    return air_dens_p

# 3. dT
def get_dT(pix_values, Inst_ETr, Kc, r_ah, air_dens):
    
    '''Gradiente de temperatura cercana a superficie
    
    Inputs
    - Inst_ETr : Instantaneous ET [mm/hr]
    - Cp : Air specific heat [J/kg/K]
    - lambda_LE : Latent heat of vaporization [J/kg]
    - Kc : Crop coefficient [-]
    - Ts_c : Surface Temperature [°C]

    Retorna diccionario
    - R_n
    - G
    - LE
    - H0    [W/m2]
    - dT0
    '''
    
    # 1. LE Latent Heat Flux
    # Kc = 1.05 para pixel frio
    # Kc = 0 para pixel caliente 
    Cp = 1004 
    Ts = pix_values['Ts_k']
    
    lambda_LE = ( 2.501-0.00236*(Ts-273.15) ) * 10**6 # Eq 53
    LE = Inst_ETr * Kc * lambda_LE / 3600  # [W m-2]
    
    # 2. H a partir de balance de energía (Eq 1)
    # H = R_n - G - LE  , todos en [W m-2] 
    R_n = pix_values['R_n']
    G = pix_values['G']
    if G < 0:
        G = 0
    H = R_n - G - LE 

    # 3. dT para pixel frío y caliente
    dT = r_ah * H /(air_dens * Cp)
    
    resultados = {'R_n': R_n, 
                  'G' : G,
                  'LE': LE,
                  'H0' : H,
                  'dT0': dT}
    
    # dT0 es producto de rah, p y H. con dT0 se obtiene los coeficientes a y b.
    # dT1 se obtiene a partir de a, b y tsdem

    return resultados

# 4. H final
def get_H(pixel_c, pixel_h):
    
    """ Obtener H a partir de r_ah, dT y air_dens_p 
    
    Inputs
    - pixel_c : Pixel Frío (dict)
    - pixel_h : Pixel Caliente (dict)
    
    Retorna:
    - Coeficientes a y b
    - dT1 a partir de coeficientes a y b
    - H1 a partir de r_ah, dT1 y air_dens_p 
    """

    # 1. Coeficientes a y b a partir de los dT0 obtenidos por pixel
    pix_c_dt = pixel_c['dT0'] 
    pix_h_dt = pixel_h['dT0'] 
    pix_c_Ts_dem = pixel_c['Ts_dem']
    pix_h_Ts_dem = pixel_h['Ts_dem']

    a = (pix_h_dt - pix_c_dt)/(pix_h_Ts_dem - pix_c_Ts_dem) # Eq 48
    b = pix_c_dt - a * pix_c_Ts_dem  # a*Ts_dem + b = dT, Eq 49

    # 2. dT1 a partir de los coeficientes
    pix_c_dt1 = a*pix_c_Ts_dem + b
    pix_h_dt1 = a*pix_h_Ts_dem + b

    # 3. H
    pix_c_H = pixel_c['air_dens_p']*1004*pix_c_dt1/pixel_c['r_ah']
    pix_h_H = pixel_h['air_dens_p']*1004*pix_h_dt1/pixel_h['r_ah']
    
    coeficientes = {'a_coef': a, 'b_coef': b}
    
    pixel_c = {**pixel_c, **coeficientes, **{'dT1': pix_c_dt1, 'H1': pix_c_H}}
    pixel_h = {**pixel_h, **coeficientes, **{'dT1': pix_h_dt1, 'H1': pix_h_H}}
    
    return pixel_c, pixel_h



# El procesado inicial genera los valores que se ingresarán a la sección iterativa

def procesado_inicial(pixel_values, elev_station, filtrado_ws, filtrado_ETr, Kc):
    
    '''Procesa datos de pixeles escogidos y de estación meteorológica.
    
    Retorna dict
    - Ts: Temperatura del pixel [K]
    - Ts_dem: Temperatura corregida por elevación [K]
    - Z_om
    - u_200
    - u*
    - r_ah
    - air_dens_p
    '''

    # Selección de parámetros
    pixel_elev = pixel_values['Elev_m']
    pixel_ts = pixel_values['Ts_k'] # K

    # 1. Ts_dem
    Tlapse_rate = 6.5   # °C/km
    Ts_dem = pixel_ts + Tlapse_rate/1000*(pixel_elev-elev_station) # K


    # 2. Emisividad e broadband

    # 2.1. Z_om : Momentum roughness length 
    pixel_lai = pixel_values['LAI']
    pixel_slope = pixel_values['Slope_d'] # grados

    # Flat model momemtum roughness length
    if pixel_lai < 0.2778:
        Z_om_flat = 0.005 # m
        # represents roughness typical of bare agricultural soils
    else:
        Z_om_flat = 0.018*pixel_lai
    
    # Adjusting momentm roughness length for slopes (Z_om_mtn)
    if pixel_slope < 5:
        Z_om = Z_om_flat
    else:
        Z_om = Z_om_flat * ( 1+(pixel_slope-5)/20 )
    
    # Resultado: Z_om
    
    
    # 3. u_star: Friction velocity
    # z_om_w: Roughness length for the weather station surface (Brutsaert, 1982) eq 40 en manual
    # z_om_w = 0.12*h     # h: Average vegetation height around the weather station 
    # siendo h la altura del grass de la estación met (relativo a cada estación)
    h = 0.15
    z_om_w = 0.12*h

    # 3.1. u_200: Wind speed at a blending height assumed to be 200 m (eq 32)
    # filtrado_ws: Velocidad del viento medido por estación meteorológica a las 10:30 [Kh/hr]
    blending_h = 200 
    u_w = filtrado_ws*1000/3600 # [m/s]
    z_x = elev_station  # altura del anemómetro, varía de 2 a 2.80m (relativo a cada estación)
    u_200 = u_w*math.log(blending_h/z_om_w)/math.log(z_x/z_om_w) # [m/s]    
        
    # 3.2. u*1: Friction velocity (eq 31)
    u_star = get_u_star(u_200, Z_om)

    
    # 4. r_ah1: Aerodynamic resistance (eq 30) [s/m]
    r_ah = get_rah(u_star)
    
    # 5. Air density p (eq 37) for 20 °C: 293° = 273° + 20°
    air_dens_p = get_air_dens(pixel_elev, air_temp_K=20+273.15)
    
    # 6. dT (eq 28)
    resultado1 = get_dT(pixel_values, filtrado_ETr, Kc, r_ah, air_dens_p)
    
    # resultado1
    # {'R_n': R_n, 
    #  'G' : G,
    #  'LE': LE,
    #  'H0' : H,
    #  'dT0': dT}

    resultado2 = {
        'Ts_k': pixel_ts,   # Temperatura del pixel 
        'Ts_dem': Ts_dem, # Temperatura corregida por elevación
        'Z_om': Z_om,
        'u_200': u_200,
        'u*': u_star,
        'r_ah': r_ah,
        'air_dens_p': air_dens_p
    }
    
    # Combinar diccionarios 1 y 2
    resultado = {**resultado2, **resultado1}
    
    return resultado



# Stability Correction Functions 

def stability_corr(pixel_procesado):
    
    '''Monin-Obukhov theory

    Retorna diccionario
    - y200m, yhz2, yhz1
    '''
    
    # Inputs
    H1 = pixel_procesado['H1']
    air_dens = pixel_procesado['air_dens_p']
    u_star = pixel_procesado['u*']
    Ts_k = pixel_procesado['Ts_k']
    Cp = 1004
    
    # The Monin-Obukhov length (L) defines 
    # the stability conditions of the atmosphere in the iterative process
    if H1 == 0:
        L = -1000
    else: 
        L = -air_dens*Cp*u_star**3*Ts_k/(0.41*9.81*H1)
    
   
    # When L<0, the lower atmospheric boundary layer is unstable
    # and when L>0, the boundary layer is stable
   
    # Stable conditions
    if L > 0:  
        y200m = -5*2/L
        yhz2  = -5*2/L
        yhz1  = -5*0.1/L
    
    # Unstable conditions
    else:      
        # y200m
        x200m = (1-16*200/L)**0.25
        y200m = 2*math.log((1+x200m)/2) + math.log((1+x200m**2)/2) - 2*math.atan(x200m)+0.5*math.pi
        
        # yhz2
        x2m  = (1-16*2/L)**0.25
        yhz2 = 2*math.log((1+x2m**2)/2)
        
        # yhz1
        x01m = (1-16*0.1/L)**0.25
        yhz1 = 2*math.log((1+x01m**2)/2)
    
    # Resultados
    valores_L = {
        'L': L,
        'y200m': y200m,
        'yhz2': yhz2,
        'yhz1': yhz1
    }
    
    return valores_L



# N iteraciones

def iteracion(u_200, pixel_procesado, pixel_valores, H0):

    """ Iterar parámetros de entrada para obtener r_ah corregido por estabilidad.

    Parámetros
    ----------
    u_star: float
        Requiere u_200, Z_om, correccion
    r_ah: float
        Requiere u_star, correccion
    p: float
        Requiere pix_elev, air_temp_K= Ts - dT0 + 273.15
    dT0

    Retorna
    -------
    output: dict
        Resultados de iteración.
    """
    # Valores iniciales a partir de pixel procesado
    z_om = pixel_procesado['Z_om']
    ts_k = pixel_procesado['Ts_k']
    dt0  = pixel_procesado['dT0']

    # Procesar parametros para obtener dT
    u_star_adj = get_u_star(u_200, z_om, correccion=pixel_procesado)
    
    r_ah_adj = get_rah(u_star_adj, correccion=pixel_procesado)
    
    # Densidad del aire (depende del dT0 de la iteracion anterior)
    p_adj = get_air_dens(pixel_valores['Elev_m'], 
                         air_temp_K=ts_k-dt0)
    
    # Procesar dT
    dT = r_ah_adj * H0 /(p_adj * 1004)

    # Output
    output = {
        'Ts_k': pixel_procesado['Ts_k'],
        'Ts_dem': pixel_procesado['Ts_dem'], # Añadido para obtener H
        'Z_om': pixel_procesado['Z_om'],
        'u*': u_star_adj,
        'r_ah': r_ah_adj,
        'dT0': dT,
        'air_dens_p': p_adj
    }

    return output


def parte_iterativa(n_iteraciones, 
                    elev_station,
                    pix_f_values, pix_c_values, 
                    dato_ws, dato_et):
    
    '''
    Realizar primera y n iteraciones que se indique. 

    Retorna diccionario con información de iteraciones. 
    - Temperatura y Temperatura corregida por dem
    - Z_om, u_200, u*, r_ah, air_dens
    - Coeficientes a y b
    - Valores de corrección por estabilidad L, y200m, yhz2, yhz1
    '''

    # Establecer listas de inicio
    resultados_f = []
    resultados_c = []

    # 1er Iteración: Condición estable de la atmóstfera

    pix_f_procesadoinicial = procesado_inicial(pix_f_values, elev_station, dato_ws, dato_et, Kc=1.05)
    pix_c_procesadoinicial = procesado_inicial(pix_c_values, elev_station, dato_ws, dato_et, Kc=0) 

    u_200 = pix_f_procesadoinicial['u_200']

    pix_f_procesado, pix_c_procesado = get_H(pix_f_procesadoinicial, pix_c_procesadoinicial)

    pix_f_stability_cor = stability_corr(pix_f_procesado)
    pix_c_stability_cor = stability_corr(pix_c_procesado)
    
    resultados_f.append({**{'iter':1}, **pix_f_procesado, **pix_f_stability_cor})
    resultados_c.append({**{'iter':1}, **pix_c_procesado, **pix_c_stability_cor})

    # Parte iterativa: N Iteraciones 
    x = 0

    while x < n_iteraciones-1:

        pix_f_iter = iteracion(u_200, resultados_f[-1], pix_f_values, pix_f_procesadoinicial['H0'])
        pix_c_iter = iteracion(u_200, resultados_c[-1], pix_c_values, pix_c_procesadoinicial['H0'])

        pix_f_post_iter, pix_c_post_iter = get_H(pix_f_iter, pix_c_iter)

        pix_f_stability_cor = stability_corr(pix_f_post_iter)
        pix_c_stability_cor = stability_corr(pix_c_post_iter)

        resultados_f.append({**{'iter':x+2}, **pix_f_post_iter, **pix_f_stability_cor})
        resultados_c.append({**{'iter':x+2}, **pix_c_post_iter, **pix_c_stability_cor})

        x+=1

    return resultados_f, resultados_c


def get_grafica_iteracion(fecha, 
                          resultados_dict, 
                          title_name, 
                          color_l, 
                          fin_iteraciones,
                          save_fig=None):

    """Obtener gráficas de sección iterativa para pixeles frío y caliente.
    
    Modificar la ruta para resultados.
    """
    
    list_variables = ['r_ah', 'air_dens_p', 'dT0']
    ylabels = ['Rah', 'Air density', 'dT']
    valores = np.arange(1, fin_iteraciones+2+1)
    
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), layout='constrained', sharex=True)
        
    axs[0].set_title(f'{fecha} - Resultados {title_name}')

    for i, variable in enumerate(list_variables):

        # Gráfica
        axs[i].plot(valores, resultados_dict[variable], color=color_l, marker='.')
        axs[i].set(xticks=np.arange(0, 2+fin_iteraciones+1, step=2), ylabel=ylabels[i])
        axs[i].grid(alpha=0.2)

    axs[2].set_xlabel('N° Iteraciones')

    # Guardar gráfica
    if save_fig is True:
        ruta = f'output/{fecha}/'
        nombre_img = f'{fecha} - {title_name}.png'
        fig.savefig(ruta+nombre_img) # ,dpi=400)


def get_H_corregido(img_productos_clipped, iteracion_output, dem, elev_station):
    
    """ Función aerodinámica para estimar el parámetro H
    
    Retorna:
    img_H: ee.Image
    """
    
    # =====================================================
    # Inputs 

    # Mapas elaborados al inicio
    img_tempK = img_productos_clipped.select('Ts_k')
    img_tempC = img_productos_clipped.select('Ts_c')
    img_lai = img_productos_clipped.select('LAI')

    # Resultados de proceso iterativo
    u_200 = iteracion_output[0]['u_200']
    
    # =====================================================
    # Sensible Heat Flux H

    # 1. dT
    # 1.1. Ts_dem
    Ts_dem = img_tempK.expression('img_tempK + Tlapse_rate/1000*(dem_elev-elev_station)',
                                  {'img_tempK': img_tempK,
                                   'Tlapse_rate': 6.5,
                                   'dem_elev': dem,
                                   'elev_station': elev_station})

    # Observar como actúa la corrección por temperatura usando el DEM
    # display(get_stats(img_tempC, predios_agricolas, 30))
    # display(get_stats(Ts_dem, predios_agricolas, 30))


    # 1.2. dT (para temperaturas en K)
    coef_a = iteracion_output[-1]['a_coef'] 
    coef_b = iteracion_output[-1]['b_coef'] 

    # img_dT = b + a*Ts_dem
    img_dT = Ts_dem.expression('b + a*Ts_dem',
                               {'a': coef_a, 'b': coef_b, 'Ts_dem': Ts_dem}).rename('dT')

    # get_stats(img_dT, predios_agricolas, 30)

    # 2. Air density
    air_pressureP = dem.expression(
        '101.3*( (293.15-0.0065*dem_elev)/293.15)**5.26',
        {'dem_elev': dem})

    air_dens_p = air_pressureP.expression(
        '1000*air_pressureP/(1.01*(air_tempK)*287)',
        {'air_pressureP': air_pressureP,
         'air_tempK': img_tempK.subtract(img_dT)}).rename('air_dens')

    # get_stats(air_dens_p, predios_agricolas, 30)


    # 4. rah
    # 4.1. Z_om
    # 2. Momentum roughness length Z_om
    # Flat model momemtum roughnes length
    # img_lai < 0.2778 -> Z_om = 0.005 # represents roughness typical of bare agricultural soils
    # img_lai >= 0.2778 -> Z_om = 0.018*img_lai
    # img_zom1 = img_lai.multiply(0.018)
    img_zom1 = img_lai.where(img_lai.gte(0.2778), img_lai.multiply(0.018)).where(img_lai.lt(0.2778), 0.005).rename('Z_om')
    # get_stats(img_zom1, predios_agricolas, 30)

    # Adjusting momentum roughness length for slopes
    # img_zom2 = img_zom1.where(img_slope.gte(5), img_slope.subtract(5).divide(20).add(1).multiply(img_zom1))
    # get_stats(img_zom2, predios_agricolas, 30)

    # rah
    img_ustar = img_zom1.expression('0.41*u_200/( log(200/zom) - y200m )',
                                    {'u_200': u_200,
                                     'zom': img_zom1,
                                     'y200m': iteracion_output[-1]['y200m']}).rename('u_star')

    img_rah = img_ustar.expression('( log(2/0.1) - yhz2 + yhz1 )/(img_ustar*0.41)',
                                   {'yhz2': iteracion_output[-1]['yhz2'],
                                    'yhz1': iteracion_output[-1]['yhz1'],
                                    'img_ustar': img_ustar}).rename('r_ah')

    # get_stats(img_rah, predios_agricolas, 30)

    # Air Specific Heat Cp: 1004
    img_H = air_dens_p.multiply(1004).multiply(img_dT).divide(img_rah).rename('H')

    return img_H


# Utilitarios

def grafica_coefs(df_resultados, fecha, coef_a, coef_b, save_files=None):

    # Gráfica 1x2
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

    axs[0].plot(df_resultados['iter'], df_resultados['a_coef'])
    axs[0].axhline(0, color='black', alpha=.2)
    axs[0].set(title=f'{fecha}\nCoeficiente a: {coef_a:.3f}', ylabel='Valores', xlabel='Iteraciones')
    axs[0].grid(alpha=0.2)
    
    axs[1].plot(df_resultados['iter'], df_resultados['b_coef'])
    axs[1].axhline(0, color='black', alpha=.2)
    axs[1].set(title=f'{fecha}\nCoeficiente b: {coef_b:.3f}', xlabel='Iteraciones')
    axs[1].grid(alpha=0.2)

    if save_files != None:
        fig.savefig(f'output/{fecha}_coefs.png')

    plt.show()



def tabla_coefs(lista_imgs, n_imgs, lista_nombres, save_files=None):
    
    lista_pixeles = []
    for index in range(n_imgs):
        pixeles = lista_imgs[index]['pixeles']
        lista_pixeles.append(pixeles)

    df_pixeles = pd.DataFrame(lista_pixeles, index=lista_nombres)

    if save_files != None:
        df_pixeles.to_csv(save_files)
        
    return df_pixeles