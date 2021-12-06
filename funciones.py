# Funciones usadas en Inv-Arroz-METRIC

import ee
import geemap
import pandas as pd
from ipywidgets import interact, fixed

# Función generar_reporte(icol)

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
  ID_snippet_name = icol.getInfo()['id']

  # Lista con fechas (en milisegundos)
  lista_fechas = icol.aggregate_array('system:time_start').getInfo() # Fechas
  # Lista con ID's
  imgCol_ids = [f'{ID_snippet_name}/{i}' for i in icol.aggregate_array("system:index").getInfo()] # ID

  # Armando dataframe
  df = pd.DataFrame(lista_fechas, columns = ['millis'])

  df["ID GEE"] = imgCol_ids
  df['Timestamp'] = pd.to_datetime(df['millis'], unit='ms')
  df['Fecha Precisa'] = pd.DatetimeIndex(df['Timestamp']) # Con hora
  # df['Fecha Corta'] = pd.DatetimeIndex(df['Timestamp']).date
  df['Año'] = pd.DatetimeIndex(df['Timestamp']).year
  df['Mes'] = pd.DatetimeIndex(df['Timestamp']).month
  df['Día'] = pd.DatetimeIndex(df['Timestamp']).day
  # df['Hora'] = pd.DatetimeIndex(df['Timestamp']).hour
  df['Día Juliano'] = pd.DatetimeIndex(df['Timestamp']).dayofyear

  # Agregando propiedades
  if ID_snippet_name == 'LANDSAT/LC08/C01/T1_SR':
    df["Zenith Angle"] = icol.aggregate_array('SOLAR_ZENITH_ANGLE').getInfo()
    df["Porcentaje Nubes (%)"] = icol.aggregate_array('CLOUD_COVER').getInfo()
  
  elif ID_snippet_name == 'COPERNICUS/S2_SR':
    df["Zenith Angle"] = icol.aggregate_array('MEAN_SOLAR_ZENITH_ANGLE').getInfo()
    df["Porcentaje Nubes (%)"] = icol.aggregate_array('CLOUDY_PIXEL_PERCENTAGE').getInfo()
  
  df = df.drop(columns=['millis', 'Timestamp'])
  # df.to_csv('datos_2020_L8_SR.csv')
  
  return df


# Función ver_imgs_mensual(mes, df, roi)

def ver_imgs_mensual(mes, df, roi):
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
  vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': ['B4', 'B3', 'B2']}
  # 'gamma':1.4
  
  # geemap para la visualización
  Map = geemap.Map(basemap='OpenStreetMap.Mapnik', layer_ctrl=True)
  Map.centerObject(roi, 12)  # Map.setCenter(-79.809, -6.746, 9)
  if lista_imagenes[0][:4] == 'COPE':
    for i in lista_imagenes:                   
      Map.addLayer(ee.Image(i).multiply(0.0001), vis_rgb, f'Imagen {i[21:23]}/{i[23:25]}') # .clip(roiDep)
  else:
    for i in lista_imagenes:
      Map.addLayer(ee.Image(i).multiply(0.0001), vis_rgb, f'Imagen {i[-4:-2]}/{i[-2:]}') # .clip(roiDep)
    
  # Visualizar al último el distrito de Chongoyape
  Map.addLayer(roi, {'color':'00FF00'}, 'Chongoyape') # roiChongoyapeDraw
  
  return Map


# Función `maskS2clouds(image)` para enmascarar nubes en imágenes Sentinel-2

## Operadores de bit a bit
# Operador bit a bit desplazamiento a la izquierda <<
# Little indian (los bits se cuentan de izq a derecha) y Big indian (los bits se cuentan de derecha a izq)
def maskS2clouds(image):
  qa = image.select('QA60')
  opaque_cloud = 1 << 10
  cirrus_cloud = 1 << 11
  mask = qa.bitwiseAnd(opaque_cloud).eq(0)\
           .And(qa.bitwiseAnd(cirrus_cloud).eq(0))
  clean_image = image.updateMask(mask)
  return clean_image