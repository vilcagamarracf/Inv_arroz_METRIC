import ee

import geemap
from geemap import cartoee

import os
import pandas as pd
import matplotlib.pyplot as plt


# Función `generar_reporte(icol)`

def generar_reporte(icol):

    """Generar reportes para ee.ImageCollection 's
    
    A partir de un ee.ImageCollection devolver una tabla 
    con información de ID, Fechas (Año, Mes, Día, Día Juliano) y 
    propiedades como Ángulo Zenital y Porcentaje de Nubosidad.

    Argumentos:
    - icol : ee.ImageCollection
    Retorna:
    - pandas.Dataframe
    """

    ## Generando campos para la tabla
    ID_snippet_name = icol.get('system:id').getInfo()

    # Lista con fechas UTC+00 (en milisegundos)
    lista_fechas = icol.aggregate_array('system:time_start').getInfo()

    # Lista con ID's
    imgCol_ids = icol.aggregate_array('system:id').getInfo() # ID

    # Armando dataframe
    df = pd.DataFrame(lista_fechas, columns = ['millis'])

    df["ID GEE"] = imgCol_ids
    df['Timestamp'] = pd.to_datetime(df['millis'], unit='ms') - pd.DateOffset(hours=5) # UTC-5
    df['Fecha Precisa'] = pd.DatetimeIndex(df['Timestamp']) # Con hora
    df['Fecha Corta'] = pd.DatetimeIndex(df['Timestamp']).date
    df['Año'] = pd.DatetimeIndex(df['Timestamp']).year
    df['Mes'] = pd.DatetimeIndex(df['Timestamp']).month
    df['Día'] = pd.DatetimeIndex(df['Timestamp']).day
    # df['Hora'] = pd.DatetimeIndex(df['Timestamp']).hour
    df['Día Juliano'] = pd.DatetimeIndex(df['Timestamp']).dayofyear

    l_icols_id = [
        'LANDSAT/LC08/C01/T1', 'LANDSAT/LC08/C02/T1', 'LANDSAT/LC08/C02/T1_TOA', 
        'LANDSAT/LC09/C02/T1','LANDSAT/LC09/C02/T1_L2', 'LANDSAT/LC09/C02/T1_TOA',
        'LANDSAT/LC09/C02/T1_RT_TOA',
        'LANDSAT/LE07/C02/T1', 'LANDSAT/LE07/C02/T1_TOA',
        'LANDSAT/LT05/C01/T1_TOA', 'LANDSAT/LM05/C01/T1'
        ]

    l_icols_sr_id = [
        'LANDSAT/LT05/C02/T1_L2',
        'LANDSAT/LC08/C02/T1_L2',
        'LANDSAT/LE07/C02/T1_L2'
        ]

    icol_s2 = [
        'COPERNICUS/S2',
        'COPERNICUS/S2_SR'
        ]

    # Agregando propiedades
    if ID_snippet_name == 'LANDSAT/LC08/C01/T1_SR':
        df['CLOUD_COVER'] = icol.aggregate_array('CLOUD_COVER').getInfo()
        df['SOLAR_ZENITH_ANGLE'] = icol.aggregate_array('SOLAR_ZENITH_ANGLE').getInfo()

    elif ID_snippet_name in l_icols_id:
        df['CLOUD_COVER'] = icol.aggregate_array('CLOUD_COVER').getInfo()
        df['CLOUD_COVER'] = df['CLOUD_COVER'].round(decimals=3)
        df['SUN_ELEVATION'] = icol.aggregate_array('SUN_ELEVATION').getInfo()

    elif ID_snippet_name in l_icols_sr_id:
        df['CLOUD_COVER'] = icol.aggregate_array('CLOUD_COVER').getInfo()
        df['SUN_ELEVATION'] = icol.aggregate_array('SUN_ELEVATION').getInfo()
        df['ZENITH ANGLE'] = 90. - df['SUN_ELEVATION']
        df['ZENITH ANGLE'] = df['ZENITH ANGLE'].round(decimals=7)
        df = df.drop(columns=['SUN_ELEVATION'])

    elif ID_snippet_name in icol_s2:
        df['MEAN SOLAR ZENITH ANGLE'] = icol.aggregate_array('MEAN_SOLAR_ZENITH_ANGLE').getInfo()
        df['CLOUDY PIXEL PERCENTAGE'] = icol.aggregate_array('CLOUDY_PIXEL_PERCENTAGE').getInfo()
    
    df = df.drop(columns=['millis', 'Timestamp'])
    
    return df


# Función `ver_imgs_mensual(mes, df, snippet_name, roi)`

def ver_imgs_mensual(mes, df, snippet_name, roi):

    """Visualizar imágenes a partir del reporte generado con la función generar_reporte
    
    Argumentos:
    - mes : lista de valores [1,2,3] del 1 al 12 que representen los meses
    - df : pandas.DataFrame con columna 'ID GEE' de ID's 
    - roi : ee.geometry o ee.FeatureCollection 
    Retorna:
    - ee.Image

    Usarse con interact
    """
    
    # Del dataframe obtenemos las imagenes de acuerdo al mes que queramos
    lista_imagenes = df[df['Mes'] == mes]['ID GEE'].tolist()

    # Parametros de visualización RGB
    bands = ['B4', 'B3', 'B2']
    vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': bands}# 'gamma':1.4
    
    # Visualización
    # import geemap.eefolium as geemap
    Map = geemap.Map(layer_ctrl=True) # basemap='OpenStreetMap.Mapnik', 
    Map.centerObject(roi, 11)  # Map.setCenter(-79.809, -6.746, 9)

    icol_t1 = [
        'LANDSAT/LC08/C01/T1', 
        'LANDSAT/LC08/C02/T1', 
        'LANDSAT/LC09/C02/T1'
        ]
                    
    icol_t1_l2 = [
        'LANDSAT/LC08/C02/T1_L2', 
        'LANDSAT/LC09/C02/T1_L2'
        ]

    icol_le07_lt05 = [
        'LANDSAT/LT05/C02/T1_L2', 
        'LANDSAT/LE07/C02/T1_L2'
        ]

    icol_toa = [
        'LANDSAT/LC08/C02/T1_TOA',
        'LANDSAT/LC09/C02/T1_TOA'
        ]

    icol_s2 = ['COPERNICUS/S2','COPERNICUS/S2_SR']

    # Sentinel-2
    if snippet_name in icol_s2:
        for i in lista_imagenes:
            image = ee.Image(i).multiply(0.0001)                 
            Map.addLayer(image, vis_rgb, f'Imagen {i}')

    # Landsat 5, 7, 8 y 9
    elif snippet_name == 'LANDSAT/LC08/C01/T1_SR':
        for i in lista_imagenes:
            image = ee.Image(i).multiply(0.0001)
            Map.addLayer(image, vis_rgb, f'Imagen {i[-4:-2]}/{i[-2:]}')

    elif snippet_name in icol_t1:
        for i in lista_imagenes:
            image = ee.Image(i)
            vis_rgb = {'min': 0.0, 'max': 30000, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'{i[:12]} - {i[-4:-2]}/{i[-2:]}')

    elif snippet_name in icol_toa:
        for i in lista_imagenes:
            image = ee.Image(i)
            vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'Imagen {i[-4:-2]}/{i[-2:]}')

    elif snippet_name in icol_t1_l2:
        for i in lista_imagenes:
            image = ee.Image(i).multiply(0.0000275).add(-0.2)
            bands = ['SR_B4', 'SR_B3', 'SR_B2']
            vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': bands} 
            Map.addLayer(image, vis_rgb, f'Imagen {i[-4:-2]}/{i[-2:]}')
    
    elif snippet_name in icol_le07_lt05:
        for i in lista_imagenes:
            image = ee.Image(i).multiply(0.0000275).add(-0.2)
            bands = ['SR_B3', 'SR_B2', 'SR_B1']
            vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'Imagen {i[-4:-2]}/{i[-2:]}')

    elif snippet_name == 'LANDSAT/LE07/C02/T1_TOA':
        for i in lista_imagenes:
            image = ee.Image(i)
            bands = ['B3', 'B2', 'B1']
            vis_rgb = {'min': 0.0, 'max': 0.4, 'gamma': 1.2, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'Imagen {i[-4:-2]}/{i[-2:]}')

    elif snippet_name == 'LANDSAT/LE07/C02/T1':
        for i in lista_imagenes:
            image = ee.Image(i)
            bands = ['B3', 'B2', 'B1']
            vis_rgb = {'min': 0.0, 'max': 300, 'gamma': 1.2, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'Imagen {i[-4:-2]}/{i[-2:]}')
    
    # ROI
    empty = ee.Image().paint(roi, 3,5) # https://github.com/google/earthengine-api/blob/6445cae4c371a8244f70ae08c01a6da05dbc4c7d/python/examples/py/FeatureCollection/from_polygons.py
    Map.addLayer(empty, {}, 'ROI')

    return Map

# Función `get_grafica_cartoee(image, vis_params, title_map, label, img_date, save_fig=None)`

def get_grafica_cartoee(image, vis_params, title_map, label, save_fig=None):
    
    zoom_region = [-79.787761, -6.607186, -79.771781, -6.593991]

    fig = plt.figure(figsize=(8,6))

    # ee.Image a plotear
    ax = cartoee.get_map(image, region=zoom_region, vis_params=vis_params)

    # Añadir barra de color
    # https://matplotlib.org/stable/api/colorbar_api.html#module-matplotlib.colorbar   
    cb = cartoee.add_colorbar(
        ax, 
        vis_params=vis_params, 
        loc='right', 
        label=label, 
        posOpts=[0.86, 0.25, 0.02, 0.5]
    )# drawedges=True, extend='both')

    # Añadir grillas
    cartoee.add_gridlines(
        ax, 
        interval=0.005, 
        ytick_rotation=90, 
        linestyle=":", 
        linewidth=0
    ) # xtick_rotation=45

    # Añadir texto
    ax.text(-79.7828, -6.5962, title_map, fontsize=14)

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
    
    ax.tick_params(axis = 'both', labelsize = 11)
    
    # Guardar graficas
    if save_fig != None:
        # plt.savefig(f'/Users/usuario/Downloads/{title_map}_{img_date}.jpg')
        # plt.savefig(f'/content/{title_map}_{img_date}.jpg')
        ruta = r'C:/Users/usuario/Documents/00-notebooks-2022/images/cartoee/'
        ruta_img = ruta + save_fig + '.jpg'
        plt.savefig(ruta_img)#, dpi=300)
    else:
        pass

    plt.show()