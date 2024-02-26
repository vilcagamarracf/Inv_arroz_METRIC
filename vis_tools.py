import ee
import geemap
from geemap import cartoee

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

    # Generando campos para la tabla
    icol_id = icol.get('system:id').getInfo()

    # Lista con fechas UTC+00 (en milisegundos)
    list_icol_dates = icol.aggregate_array('system:time_start').getInfo()

    # Lista con ID's
    list_icol_ids = icol.aggregate_array('system:id').getInfo()  # ID
  
    # Armando dataframe
    icol_df = pd.DataFrame(list_icol_dates, columns=['millis'])

    icol_df["ID GEE"] = list_icol_ids
    icol_df['Timestamp'] = pd.to_datetime(icol_df['millis'], unit='ms') - pd.DateOffset(hours=5)  # UTC-5
    icol_df['Fecha Precisa'] = pd.DatetimeIndex(icol_df['Timestamp'])  # Con hora
    icol_df['Fecha Corta'] = pd.DatetimeIndex(icol_df['Timestamp']).date
    icol_df['Año'] = pd.DatetimeIndex(icol_df['Timestamp']).year
    icol_df['Mes'] = pd.DatetimeIndex(icol_df['Timestamp']).month
    icol_df['Día'] = pd.DatetimeIndex(icol_df['Timestamp']).day
    # df['Hora'] = pd.DatetimeIndex(df['Timestamp']).hour
    icol_df['Día Juliano'] = pd.DatetimeIndex(icol_df['Timestamp']).dayofyear

    list_icol_landsat = [
        'LANDSAT/LC08/C01/T1', 'LANDSAT/LC08/C02/T1', 'LANDSAT/LC08/C02/T1_TOA',
        'LANDSAT/LC09/C02/T1', 'LANDSAT/LC09/C02/T1_L2', 'LANDSAT/LC09/C02/T1_TOA',
        'LANDSAT/LC09/C02/T1_RT_TOA',
        'LANDSAT/LE07/C02/T1', 'LANDSAT/LE07/C02/T1_TOA',
        'LANDSAT/LT05/C01/T1_TOA', 'LANDSAT/LM05/C01/T1'
    ]

    list_icol_landsat_sr = [
        'LANDSAT/LT05/C02/T1_L2',
        'LANDSAT/LC08/C02/T1_L2',
        'LANDSAT/LE07/C02/T1_L2'
    ]

    list_icol_sentinel2 = [
        'COPERNICUS/S2',
        'COPERNICUS/S2_SR'
    ]

    # Agregando propiedades
    if icol_id == 'LANDSAT/LC08/C01/T1_SR':
        icol_df['CLOUD_COVER'] = icol.aggregate_array('CLOUD_COVER').getInfo()
        icol_df['SOLAR_ZENITH_ANGLE'] = icol.aggregate_array('SOLAR_ZENITH_ANGLE').getInfo()

    elif icol_id in list_icol_landsat:
        icol_df['CLOUD_COVER']   = icol.aggregate_array('CLOUD_COVER').getInfo()
        icol_df['CLOUD_COVER']   = icol_df['CLOUD_COVER'].round(decimals=3)
        icol_df['SUN_ELEVATION'] = icol.aggregate_array('SUN_ELEVATION').getInfo()

    elif icol_id in list_icol_landsat_sr:
        icol_df['CLOUD_COVER']   = icol.aggregate_array('CLOUD_COVER').getInfo()
        icol_df['SUN_ELEVATION'] = icol.aggregate_array('SUN_ELEVATION').getInfo()
        icol_df['ZENITH ANGLE']  = 90. - icol_df['SUN_ELEVATION']
        icol_df['ZENITH ANGLE']  = icol_df['ZENITH ANGLE'].round(decimals=7)
        icol_df = icol_df.drop(columns=['SUN_ELEVATION'])

    elif icol_id in list_icol_sentinel2:
        icol_df['MEAN SOLAR ZENITH ANGLE'] = icol.aggregate_array('MEAN_SOLAR_ZENITH_ANGLE').getInfo()
        icol_df['CLOUDY PIXEL PERCENTAGE'] = icol.aggregate_array('CLOUDY_PIXEL_PERCENTAGE').getInfo()

    icol_df = icol_df.drop(columns=['millis', 'Timestamp'])

    return icol_df


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
    vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': bands}  # 'gamma':1.4

    # Visualización
    # import geemap.eefolium as geemap
    Map = geemap.Map(layer_ctrl=True)  # basemap='OpenStreetMap.Mapnik',
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

    icol_s2 = ['COPERNICUS/S2', 'COPERNICUS/S2_SR']

    # Sentinel-2
    if snippet_name in icol_s2:
        for i in lista_imagenes:
            image = ee.Image(i).multiply(0.0001)
            Map.addLayer(image, vis_rgb, f'{i}')

    # Landsat 5, 7, 8 y 9
    elif snippet_name == 'LANDSAT/LC08/C01/T1_SR':
        for i in lista_imagenes:
            image = ee.Image(i).multiply(0.0001)
            Map.addLayer(image, vis_rgb, f'{i[-4:-2]}/{i[-2:]}')

    elif snippet_name in icol_t1:
        for i in lista_imagenes:
            image = ee.Image(i)
            vis_rgb = {'min': 0.0, 'max': 30000, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'{i[:12]} - {i[-4:-2]}/{i[-2:]}')

    elif snippet_name in icol_toa:
        for i in lista_imagenes:
            image = ee.Image(i)
            vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'{i[-4:-2]}/{i[-2:]} L{i[11:12]}')

    elif snippet_name in icol_t1_l2:
        for i in lista_imagenes:
            image = ee.Image(i).multiply(0.0000275).add(-0.2)
            bands = ['SR_B4', 'SR_B3', 'SR_B2']
            vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'{i[-4:-2]}/{i[-2:]}')

    elif snippet_name in icol_le07_lt05:
        for i in lista_imagenes:
            image = ee.Image(i).multiply(0.0000275).add(-0.2)
            bands = ['SR_B3', 'SR_B2', 'SR_B1']
            vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'{i[-4:-2]}/{i[-2:]}')

    elif snippet_name == 'LANDSAT/LE07/C02/T1_TOA':
        for i in lista_imagenes:
            image = ee.Image(i)
            bands = ['B3', 'B2', 'B1']
            vis_rgb = {'min': 0.0, 'max': 0.4, 'gamma': 1.2, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'{i[-4:-2]}/{i[-2:]}')

    elif snippet_name == 'LANDSAT/LE07/C02/T1':
        for i in lista_imagenes:
            image = ee.Image(i)
            bands = ['B3', 'B2', 'B1']
            vis_rgb = {'min': 0.0, 'max': 300, 'gamma': 1.2, 'bands': bands}
            Map.addLayer(image, vis_rgb, f'{i[-4:-2]}/{i[-2:]}')

    # ROI
    # https://github.com/google/earthengine-api/blob/6445cae4c371a8244f70ae08c01a6da05dbc4c7d/python/examples/py/FeatureCollection/from_polygons.py
    empty = ee.Image().paint(roi, 3, 5)
    Map.addLayer(empty, {}, 'ROI')

    return Map


