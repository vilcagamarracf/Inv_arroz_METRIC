import math
import pandas as pd
import matplotlib.pyplot as plt

# gee
import ee
import geemap
import geemap.colormaps as cmp # Paletas para visualización de imágenes (cm se confunde con matplotlib)import ee
import math
import pandas as pd

ee.Initialize()

# @markdown Función `convert_RAW_to_TOA(img)`

# @markdown **Retorna:**
# @markdown - `ee.Image`: Imagen TOA (no es necesario aplicar factor de escala 0.00001)

def convert_RAW_to_TOA(img_raw):
    """
    Conversión de imágenes RAW a Radiancia y Reflectancia
    al tope de la atmósfera (TOA)

    Reflectancia TOA: B, G, R, NIR, SWIR1, SWIR2 
    Radiancia TOA: B10
    """
    # TOA Reflectance for Bands R, G, B, Nir

    sun_elev = ee.Number(img_raw.get('SUN_ELEVATION')).multiply(math.pi /180)

    for i in range(2,8): # B, G, R, NIR, SWIR1, SWIR2
        reflectance_mult_band_x = ee.Number(img_raw.get(f'REFLECTANCE_MULT_BAND_{i}'))
        reflectance_add_band_x =  ee.Number(img_raw.get(f'REFLECTANCE_ADD_BAND_{i}'))
        band = (
            img_raw.select(f'B{i}')
            .multiply(reflectance_mult_band_x)
            .add(reflectance_add_band_x)
            .divide(sun_elev.sin()) # Corrección por ángulo solar de elev
            )
        img_raw = img_raw.addBands(band, None, True)

    # TOA Radiance for Thermal Band B10 
    # RADIANCE_MULT_BAND_10 & 11 = 3.3420E-04
    # RADIANCE_ADD_BAND_10  & 11 = 0.1

    b10 = img_raw.select(['B10']).multiply(3.3420E-04).add(0.1)

    bandas = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10']

    img_toa = img_raw.addBands(b10, None, True).select(bandas)

    return img_toa



# @markdown Índices de vegetación aplicables a Landsat 8 y 9 TOA
# @markdown - NDVI: Función `get_ndvi_L8(image)`
# @markdown - SAVI: Función `get_savi_L8(image)` **(Asumiendo L = 0.5)**
# @markdown - LAI : Función `get_lai_L8(image)`
# # @markdown - EVI : Función `get_evi_L8(image)`

# @markdown **Retorna:**
# @markdown - `ee.Image`: Se agregará el resultado como banda a la imagen input.

# https://www.usgs.gov/landsat-missions/landsat-normalized-difference-vegetation-index
def get_ndvi_L8(image_toa):

    ndvi = image_toa.expression(
        '(nir - red) / (nir + red)',
        {'red' : image_toa.select('B4'),
        'nir' : image_toa.select('B5')}
        ).rename('NDVI')
    return ndvi

# https://www.usgs.gov/landsat-missions/landsat-soil-adjusted-vegetation-index
def get_savi_L8(image_toa):
    savi = image_toa.expression(
        '(nir - red) / (nir + red + L) * (1+L)',
        {'red' : image_toa.select('B4'),
        'nir' : image_toa.select('B5'),
        'L'   : ee.Number(0.5)}
        ).rename('SAVI')
    return savi

def get_lai_L8(image_savi):
    lai = image_savi.expression(
        '- (log((0.69 - savi)/0.59))/(0.91)',
        {'savi' : image_savi.where(image_savi.gt(0.689), 0.689)} # Evitar huecos
        ).rename('LAI')
    return lai
  
# def get_evi_L8(image):
#   evi = image.expression(
#       '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)',
#       {'red' : image.select('B4'),
#        'nir' : image.select('B5'),
#        'blue': image.select('B2')}
#        ).rename('EVI')
#   return image.addBands(evi)



# @markdown Función `get_emisividades(img_ndvi, img_savi, img_lai)` 

# Emisividad de banda estrecha / Narrow band emissivity

# @markdown **Retorna:**
# @markdown - `ee.Image`: Transmisividad de banda estrecha `e_nb` (Tasumi, 2003)
# @markdown - `ee.Image`: Emisividad de superficie de banda ancha `e_0` (Tasumi, 2003)

def get_emisividades(img_ndvi, img_lai):

    """Obtener Emisividad de la superficie
    
    e_nb: Narrow band transmissivity
    e_0: broad-band surface emissivity
    
    Argumentos:
    - img_ndvi, img_savi, img_lai : Indices a partir de una imagen 
        calibrada a Reflectancia a tope de la atmósfera (TOA)
    Retorna:
    - e_nb, e_0 : dos imágenes de emisividad 
    """
    # e_0
    e_0 = img_lai.expression('0.95 + 0.01*LAI',
                            {'LAI':img_lai}).rename('e_0') # LAI <= 3 y NDVI > 0

    e_0 = e_0.where(img_ndvi.lt(0), 0.985) # NDVI < 0
    e_0 = e_0.where(img_lai.gt(3), 0.98).rename('e_0') # LAI > 3; NDVI > 0

    # e_NB
    e_NB = img_lai.expression('0.97 + 0.0033*LAI',
                                {'LAI':img_lai}).rename('e_NB') # LAI <= 3 y NDVI > 0

    e_NB = e_NB.where(img_ndvi.lt(0), 0.985) # NDVI < 0
    e_NB = e_NB.where(img_lai.gt(3).And(img_ndvi.gt(0)), 0.98) # LAI > 3; NDVI > 0

    return e_NB, e_0


# @markdown Función `get_decl_lat_hra(img, roi, doy)` 

# @markdown **Retorna**: 
# @markdown - `ee.Number`: Ángulo de declinación terrestre `angl_decl`
# @markdown - `ee.Image` : Mapa de latitudes `latitud`  
# @markdown - `ee.Image` : Mapa de ángulos horario `angle_hra`

# @markdown Valores en radianes

def get_decl_lat_hra(img, roi, doy): 

    # Datos a partir de imagen satelital
    # doy = img.date().getRelative('day', 'year') # ee.Number 
    proj = img.select([0]).projection() # EPSG:32617 -> WGS 84 / UTM zone 17N

    img_lonlat = ee.Image.pixelLonLat().reproject(proj).clip(roi) # ee.Image con bandas 'Longitude' y 'Latitude'

    factor_rad = ee.Number(math.pi/180) # Factor de conversión a radianes


    # Declinacion solar δ [rad] : ee.Number 
    angl_decl = ee.Number.expression(
        '-23.45*cos( 360*(d+10)/365 * factor )', 
        {'d':doy, 'factor':factor_rad}).multiply(factor_rad) 


    # Latitud ϕ [rad] : ee.Number

    latitud = img_lonlat.select('latitude').multiply(factor_rad)



    # Angulo horario [rad] : ee.Number

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

# @markdown Función `get_cos_theta(angl_decl, latitud, angle_hra, slope_rad, aspect_rad)`
# @markdown - Requiere Función `get_decl_lat_hra(img, roi, doy)`

# @markdown **Retorna**: 
# @markdown - `ee.Image`: Mapa de Cosenos del ángulo incidencia solar `cos_theta`

def get_cos_theta(angl_decl, latitud, angle_hra, slope_rad, aspect_rad):

    cos_theta = latitud.expression(
            """
            sin(delta)*sin(phi)*cos(s) 
            - sin(delta)*cos(phi)*sin(s)*cos(gamma)
            + cos(delta)*cos(phi)*cos(s)*cos(omega)
            + cos(delta)*sin(phi)*sin(s)*cos(gamma)*cos(omega)
            + cos(delta)*sin(gamma)*sin(s)*sin(omega)
            """,
            {'delta': angl_decl,
            'phi': latitud,
            's': slope_rad,
            'gamma': aspect_rad,
            'omega': angle_hra}
            ).rename('cos_theta') # ee.Image

    return cos_theta

# @markdown Función `get_surface_temp(img_toa, e_nb)`

# @markdown **Retorna**:
# @markdown - `ee.Image`: Temperatura de superficie en Kelvin

def get_surface_temp(img_toa, e_nb):

    # Factor de correción Rc dónde: Rp = 0.91, Tnb = 0.866, Rsky = 1.32
    rc = img_toa.expression(
        '(Lt10 - 0.91)/0.866 - (1-e_nb)*1.32',
        {'Lt10':img_toa.select('B10'), 'e_nb':e_nb}
        )
    
    # Termperatura de superficie: b10_k1 = 774.89, b10_k2 = 1321.08
    ts_k = rc.expression(
        '1321.08/( log( (e_nb*774.89/rc) + 1 ) )',
        {'e_nb':e_nb, 'rc':rc}
        ).rename('Ts_k')
    
    return ts_k



# @markdown Función: `convert_TOA_to_SR(img_toa, P_air
# @markdown , w, cos_theta_hor)`

# @markdown **Retorna:**
# @markdown  - `ee.Image`: Imagen Reflectancia de superficie con bandas 
# @markdown  `['B2_SR','B3_SR','B4_SR','B5_SR','B6_SR','B7_SR']`
 
def convert_TOA_to_SR(img_toa, P_air, w, cos_theta_hor):

    """Conversión de imagen TOA a Reflectancia de Superficie
    En METRIC se deriva la reflectancia de superficie bidireccional 
    utilizando funciones calibradas de transmitancia atmosférica 
    y reflectancia de la trayectoria por Tasumi et al. (2007).

    Se corrige los valores de las bandas de la imagen TOA para la 
    dispersión y absorción de la radiación solar entrante y 
    reflejada desde la superficie, basándose en una función de 
    corrección atmosférica simplificada que sólo requiere mediciones
    puntuales o estimaciones de la presión de vapor cerca
    de la superficie. 
    """
    bands = ['B2','B3','B4','B5','B6','B7']

    img_toa = img_toa.select(bands)
    
    for b in bands:
        if b == 'B2':
            C1 = 0.987
            C2 = -0.00071
            C3 = 0.000036
            C4 = 0.0880
            C5 = 0.0789
            Cb = 0.640

        elif b == 'B3':
            C1 = 2.319
            C2 = -0.00016
            C3 = 0.000105
            C4 = 0.0437
            C5 = -1.2697
            Cb = 0.310

        elif b == 'B4':
            C1 = 0.951
            C2 = -0.00033
            C3 = 0.00028
            C4 = 0.0875
            C5 = 0.1014
            Cb = 0.286

        elif b == 'B5':
            C1 = 0.375
            C2 = -0.00048
            C3 = 0.005018
            C4 = 0.1355
            C5 = 0.6621
            Cb = 0.189

        elif b == 'B6':
            C1 = 0.234
            C2 = -0.00101
            C3 = 0.004336
            C4 = 0.0560
            C5 = 0.7757
            Cb = 0.274

        elif b == 'B7':
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
            ).rename(f't_in_{b}')

        t_out_b = img_toa.expression(
            """C1*exp( C2*P_air/(K_t*cos(theta_hor))
            - (C3*W+C4)/(cos(theta_hor)) )+C5""",
            {'C1': C1, 'C2': C2, 'C3': C3, 'C4': C4, 'C5': C5, 
            'P_air' : P_air,
            'K_t': 1,
            'W': w,
            'theta_hor': 1}
            ).rename(f't_out_{b}')

        R_in_s_b = t_in_b.multiply(t_out_b)

        R_out_s_b = img_toa.expression(
            'b - Cb*(1-t_in_b)', 
            {'Cb':Cb, 't_in_b': t_in_b, 'b':img_toa.select(b)}
        )

        p_s_b = R_out_s_b.divide(R_in_s_b)

        img_toa = img_toa.addBands(p_s_b, None, True)
    
    bands_sr = ['B2_SR','B3_SR','B4_SR','B5_SR','B6_SR','B7_SR']

    img_sr = img_toa.select(bands).rename(bands_sr)

    return img_sr

# @markdown Función: `get_albedo(img_sr)`

# @markdown Se usa coeficientes de ponderación propuestos por Tasumi et al. (2008) 

# @markdown **Retorna**:
# @markdown - `ee.Image`: Mapa de albedo

def get_albedo(img_sr):

    """
    Obtener el albedo de una imagen mediante una calibración con 
    coeficientes de ponderación por banda (Tasumi et al., 2008)
    https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51)

    Argumentos:
        image (ee.Image) : Imagen satelital SR.

    Retorna: 
        albedo (ee.Image) : Albedo de una imagen satelital. 
    """

    albedo = img_sr.expression(
        '0.254*B2 + 0.149*B3 + 0.147*B4 + 0.311*B5 + 0.103*B6 + 0.036*B7',
        {
            'B2' : img_sr.select('B2_SR'),
            'B3' : img_sr.select('B3_SR'),
            'B4' : img_sr.select('B4_SR'),
            'B5' : img_sr.select('B5_SR'),
            'B6' : img_sr.select('B6_SR'),
            'B7' : img_sr.select('B7_SR')
        }
    ).rename('albedo')

    return albedo


# @markdown Función `getRadiacionNeta(img_ee, roi, dem)`

# @markdown **Retorna**:
# @markdown - `R_n`, `img_sr_tasumi`, 
# @markdown - `img_albedo`, `img_ndvi`, `img_lai`, 
# @markdown - `cos_theta_rel`, `t_sw`, `e_0`, `d2`, `ts`, `doy`

def getRadiacionNeta(img_ee, roi, dem):

    doy = img_ee.date().getRelative('day', 'year') # ee.Number
    
    img_toa = convert_RAW_to_TOA(img_ee) # ee.Image

    # Índices de vegetación (eq. 23, 19 y 18)
    img_ndvi = get_ndvi_L8(img_toa)
    img_savi = get_savi_L8(img_toa)
    img_lai  = get_lai_L8(img_savi)
    
    # Factor de conversión de grados a radianes
    factor_rad = ee.Number(math.pi/180)

    # A partir del DEM: Pendiente y Aspect [rad]
    slope_rad  = ee.Terrain.slope(dem).multiply(factor_rad)   # ee.Image 
    aspect_rad = ee.Terrain.aspect(dem).multiply(factor_rad)  # ee.Image 

    # Ángulos Declinación, Latitud y Horario [rad]
    angle_decl, latitud, angle_hra = get_decl_lat_hra(img_ee, roi, doy)


    # Parámetro: Emisividad e_0 (eq. 17)
    e_nb, e_0 = get_emisividades(img_ndvi, img_lai)


    # Parámetro: Radiación de onda corta entrante 
    # Requiere: t_sw, d2, cos_theta_rel 
    # t_sw requiere: P_air, w, cos_theta_hor
    # w requiere: e_s

    # 1. d2 INVERSE RELATIVE DISTANCE EARTH-SUN (eq. 9)
    d2 = ee.Number.expression(
        '1/( 1+0.033*cos(doy*2*pi/365) )',
        {'doy':doy, 'pi':math.pi}
        )

    # 2. cos_theta_rel [rad] (eq. 7)
    cos_theta_rel = get_cos_theta(
        angle_decl, latitud, angle_hra, slope_rad, aspect_rad
        )

    # 3. t_sw (eq. 4)
    # 3.1. P_air Atmospheric Pressure [kPa] (eq. 5)
    P_air = dem.expression(
        '101.3*( (293-0.0065*z)/293 )**(5.26)', 
        {'z': dem.select(0)}
        ) # ee.Image

    # 3.2. cos_theta_hor (eq. 8)
    cos_theta_hor = get_cos_theta(
        angle_decl, latitud, angle_hra, slope_rad=0, aspect_rad=0)
    
    # 3.3. w Water in the atmosphere [mm]
    # 3.3.1. Temperatura de superficie (eq. 20)
    ts = get_surface_temp(img_toa, e_nb)      # [K]
    ts_c = ts.subtract(273.15).rename('Ts_c') # [°C]

    # 3.3.2. Near-surface vapor pressure ea [kPa]
    # Ojo temperatura se requiere en °C

    ea = ts_c.expression(
        '6.112*exp(17.67 / (t + 243.5))', 
        {'t' : ts_c.select('Ts_c')}
        ).rename('vapor_pressure')

    # 3.3.3. Agua precipitable en la atmósfera w [mm]
    w = ea.multiply(P_air).multiply(0.14).add(2.1).rename('w')

    # 3.3.4. t_sw: broad-band atmospheric transmissivity (eq. 4)
    t_sw = P_air.expression(
        '''
        0.35 
        + 0.627*exp( 
            - 0.00146*P_air/(Kt*cos_theta_hor) 
            - 0.075*(W/cos_theta_hor)**0.4
            ) 
        ''',
        {'P_air':P_air, 'Kt':1, 'cos_theta_hor':cos_theta_hor, 'W':w}
        ).rename('t_sw')
    
    # Finalmente: R_s_incoming (eq. 3)
    R_s_in = cos_theta_rel.expression(
        '1367*cos_theta_rel*t_sw/d2',
        {'cos_theta_rel':cos_theta_rel, 't_sw':t_sw, 'd2':d2}
        ).rename('R_s_in')


    # Parámetro: Radiación de onda larga entrante
    # ea: effective atmospheric emissivity (eq. 25)
    atm_emissivity_ea = t_sw.expression(
        '0.85*(- log(t_sw))**0.09',{'t_sw': t_sw}).rename('atm_emissivity_ea')

    # Finalmente: R_l_in (eq. 24)
    R_l_in = (atm_emissivity_ea
                .multiply(5.67E-08)
                .multiply(ts.pow(4))
                .rename('R_l_in')
    )

    # Parámetro: Radiación de onda larga saliente (eq.  16)
    R_l_out = (e_0
                .multiply(5.67E-08)
                .multiply(ts.pow(4))
                .rename('R_l_out')
    )


    # Parámetro: Albedo

    # Corrección Tasumi (eqs. 10 - 14)
    img_sr_tasumi = convert_TOA_to_SR(img_toa, P_air, w, cos_theta_hor)
    # Bandas: ['B2_SR', 'B3_SR', 'B4_SR', 'B5_SR', 'B6_SR', 'B7_SR']

    # Albedo (eq. 15)
    img_albedo = get_albedo(img_sr_tasumi)

    # Radiación Neta (eq. 2)
    R_n = img_albedo.expression('(1-albedo)*R_s_in + (R_l_in - R_l_out) - (1-e_0)*R_l_in',
        {'albedo':img_albedo,
        'R_s_in':R_s_in,
        'R_l_in':R_l_in,
        'R_l_out':R_l_out,
        'e_0':e_0}
        ).rename('R_n')

    R_n = R_n.addBands([R_s_in, R_l_in, R_l_out])

    # Juntar los parámetros obtenidos en una sola imagen
    img_productos = ee.Image(
        [img_ndvi.select('NDVI'),
        img_savi.select('SAVI'),
        img_lai.select('LAI'),
        img_albedo,
        ts,
        ts_c,
        t_sw,
        e_0,
        e_nb,
        cos_theta_rel
        ]
    )

    return R_n, img_sr_tasumi, img_productos, d2, doy



# @markdown Función `get_stats(img, geometry, scale)`
# @markdown Los reductores son `.unweighted()` 

def get_stats(img, geometry, scale):

    ''' Obtener valores estadísticos de la imagen
    Uso de ee.Image.reduceRegion para obtener estadísticas
    Inputs: 
        ee.Image, ee.Geometry, scale (valor de pixel)
    Retorna: 
        Media, Mediana, Moda, Desviación estándar, Mínimo y Máximo
        (unweighted)
    '''
    
    mean   = img.reduceRegion(ee.Reducer.mean().unweighted()  , geometry=geometry, scale=scale)
    median = img.reduceRegion(ee.Reducer.median().unweighted(), geometry=geometry, scale=scale)
    mode   = img.reduceRegion(ee.Reducer.mode().unweighted()  , geometry=geometry, scale=scale)
    stdDev = img.reduceRegion(ee.Reducer.stdDev(), geometry=geometry, scale=scale)
    min_   = img.reduceRegion(ee.Reducer.min(), geometry=geometry, scale=scale)
    max_   = img.reduceRegion(ee.Reducer.max(), geometry=geometry, scale=scale)

    values = ee.List([mean, median, mode, stdDev, min_, max_])
    columns = ['mean','median','mode','stdDev','min','max']
                        
    stats = ee.Dictionary.fromLists(columns, values).getInfo()
    # stats_df = pd.DataFrame.from_dict(stats, orient='index')
    return stats