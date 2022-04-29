import ee
import geemap
import pandas as pd

# @markdown Función `generar_reporte(icol)`

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

    l_icols_id = ['LANDSAT/LC08/C01/T1','LANDSAT/LC08/C02/T1', 'LANDSAT/LC08/C02/T1_TOA', 
                    'LANDSAT/LC09/C02/T1','LANDSAT/LC09/C02/T1_L2', 'LANDSAT/LC09/C02/T1_TOA',
                    'LANDSAT/LC09/C02/T1_RT_TOA',
                    'LANDSAT/LE07/C02/T1', 'LANDSAT/LE07/C02/T1_TOA',
                    'LANDSAT/LT05/C01/T1_TOA', 'LANDSAT/LM05/C01/T1']

    l_icols_sr_id = ['LANDSAT/LT05/C02/T1_L2','LANDSAT/LC08/C02/T1_L2',
                    'LANDSAT/LE07/C02/T1_L2']

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

    elif ID_snippet_name == 'COPERNICUS/S2_SR':
        df['MEAN_SOLAR_ZENITH_ANGLE'] = icol.aggregate_array('MEAN_SOLAR_ZENITH_ANGLE').getInfo()
        df['CLOUDY_PIXEL_PERCENTAGE'] = icol.aggregate_array('CLOUDY_PIXEL_PERCENTAGE').getInfo()
    
    df = df.drop(columns=['millis', 'Timestamp'])
    # df.to_csv('datos_2020_L8_SR.csv')
    
    return df

# @markdown Función `ver_imgs_mensual(mes, df, snippet_name, roi)`

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

    icol_t1 = ['LANDSAT/LC08/C01/T1', 'LANDSAT/LC08/C02/T1', 
                'LANDSAT/LC09/C02/T1']
                    
    icol_t1_l2 = ['LANDSAT/LC08/C02/T1_L2', 
                    'LANDSAT/LC09/C02/T1_L2']

    icol_le07_lt05 = ['LANDSAT/LT05/C02/T1_L2', 
                        'LANDSAT/LE07/C02/T1_L2']

    icol_toa = ['LANDSAT/LC08/C02/T1_TOA',
                'LANDSAT/LC09/C02/T1_TOA']

    # Sentinel-2
    if snippet_name == 'COPERNICUS/S2_SR':
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