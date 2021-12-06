# Funciones usadas en Inv-Arroz-METRIC

#------------- LANDSAT 8 ------------- 
def generar_reporte_landsat8_T1_SR(path, row, year):

  ID_snippet_name = "LANDSAT/LC08/C01/T1_SR"
  ## Filtrado de ee.imageCollection
  icol_sr = (
      ee.ImageCollection(ID_snippet_name)\
      .filterDate(str(year), str(year+1))\
      .filterMetadata('WRS_PATH', 'equals', path)\
      .filterMetadata('WRS_ROW', 'equals', row)
  )

  ## Generando campos para la tabla
  # Lista con fechas (en milisegundos)
  lista_fechas = icol_sr.aggregate_array('system:time_start').getInfo()

  # Lista con ID's
  imgCol_ids = [f'{ID_snippet_name}/{i}' 
                for i in icol_sr.aggregate_array("system:index").getInfo()]

  # Tabla con pandas
  # import pandas as pd

  df = pd.DataFrame(lista_fechas, columns = ['millis'])

  df["Landsat ID"] = imgCol_ids
  df['Timestamp'] = pd.to_datetime(df['millis'], unit='ms')
  df['Fecha Precisa'] = pd.DatetimeIndex(df['Timestamp']) # Con hora
  df['Fecha Corta'] = pd.DatetimeIndex(df['Timestamp']).date
  df['Año'] = pd.DatetimeIndex(df['Timestamp']).year
  df['Mes'] = pd.DatetimeIndex(df['Timestamp']).month
  df['Día'] = pd.DatetimeIndex(df['Timestamp']).day
  df['Hora'] = pd.DatetimeIndex(df['Timestamp']).hour
  df['Día Juliano'] = pd.DatetimeIndex(df['Timestamp']).dayofyear
  # df['Sensor'] = 'OLI/TIRS'
  
  df["Zenith Angle"] = icol_sr.aggregate_array('SOLAR_ZENITH_ANGLE').getInfo()
  df["Porcentaje Nubes (%)"] = icol_sr.aggregate_array('CLOUD_COVER').getInfo()
  
  df = df.drop(columns=['millis', 'Timestamp'])
  # df.to_csv('datos_2020_L8_SR.csv')
  
  return df

# @markdown Función `visualizar_mes_L8_SR(mes, df, roi)`

def visualizar_mes_L8_SR(mes, df, roi):
  
  # Del dataframe obtenemos las imagenes de acuerdo al mes que queramos
  lista_imagenes = df[df['Mes'] == mes]['Landsat ID'].tolist()

  # Parametros de visualización RGB
  vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': ['B4', 'B3', 'B2']}
  # 'gamma':1.4
  
  # geemap para la visualización
  Map = geemap.Map(basemap='OpenStreetMap.Mapnik', layer_ctrl=True)

  Map.centerObject(roi, 12)  # Map.setCenter(-79.809, -6.746, 9)
  for i in lista_imagenes:
    Map.addLayer(ee.Image(i).multiply(0.0001), vis_rgb, f'Imagen {i[-2:]}') # .clip(roiDep)

  # Visualizar al último el distrito de Chongoyape
  Map.addLayer(chongoyape, {'color':'00FF00'}, 'Chongoyape') # roiChongoyapeDraw
  
  return Map

# ------------- SENTINEL-2 -------------
# Función: `generar_reporte_S2_SR(roi, year)`
def generar_reporte_S2_SR(roi, year):

  ID_snippet_name = 'COPERNICUS/S2_SR'

  ## Filtrado de ee.imageCollection
  icol_sr = (
      ee.ImageCollection(ID_snippet_name)\
      .filterDate(str(year), str(year+1))\
      .filterBounds(roi)
  )

  ## Generando campos para la tabla
  # Lista con fechas (en milisegundos)
  lista_fechas = icol_sr.aggregate_array('system:time_start').getInfo()

  # Lista con ID's
  imgCol_ids = [f'{ID_snippet_name}/{i}' for i in icol_sr.aggregate_array("system:index").getInfo()]

  # Tabla con pandas
  # import pandas as pd

  df = pd.DataFrame(lista_fechas, columns = ['millis'])

  df['Sentinel-2 ID'] = imgCol_ids
  df['Timestamp'] = pd.to_datetime(df['millis'], unit='ms')
  df['Fecha Precisa'] = pd.DatetimeIndex(df['Timestamp']) # Con hora
  df['Fecha Corta'] = pd.DatetimeIndex(df['Timestamp']).date
  df['Año'] = pd.DatetimeIndex(df['Timestamp']).year
  df['Mes'] = pd.DatetimeIndex(df['Timestamp']).month
  df['Día'] = pd.DatetimeIndex(df['Timestamp']).day
  df['Hora'] = pd.DatetimeIndex(df['Timestamp']).hour
  df['Día Juliano'] = pd.DatetimeIndex(df['Timestamp']).dayofyear
  # df['Sensor'] = 'OLI/TIRS'
  
  df["Zenith Angle"] = icol_sr.aggregate_array('MEAN_SOLAR_ZENITH_ANGLE').getInfo()
  df["Porcentaje Nubes (%)"] = icol_sr.aggregate_array('CLOUDY_PIXEL_PERCENTAGE').getInfo()
  
  df = df.drop(columns=['millis', 'Timestamp'])
  # df.to_csv('datos_2020_L8_SR.csv')

  return df # ,icol_sr

# @markdown Función `maskS2clouds(image)` para enmascarar nubes

# Operadores de bit a bit
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

# @markdown Función `visualizar_mes_S2_SR(mes, df, roi)`

def visualizar_mes_S2_SR(mes, df, roi):
  
  # Del dataframe obtenemos las imagenes de acuerdo al mes que queramos
  lista_imagenes = df[df['Mes'] == mes]['Sentinel-2 ID'].tolist()

  # Parametros de visualización RGB
  vis_rgb = {'min': 0.0, 'max': 0.3, 'bands': ['B4', 'B3', 'B2']} # 'gamma':1.4

  # geemap para la visualización 
  Map = geemap.Map(layer_ctrl=True) # basemap='HYBRID',OpenStreetMap.Mapnik
  Map.centerObject(roi, 12)  #11, Map.setCenter(-79.809, -6.746, 9)

  # Agregar imágenes por mes
  # Visualizar al último el distrito de Chongoyape
  Map.addLayer(roi, {'color':'00FF00'}, 'ROI') # roiChongoyapeDraw

  for i in lista_imagenes:
    Map.addLayer(maskS2clouds(ee.Image(i)).multiply(0.0001).clip(roi), 
                 vis_rgb, f'Imagen {i[21:23]}/{i[23:25]} con máscara básica') # {i} .clip(roiChongoyape)
    Map.addLayer(ee.Image(i).multiply(0.0001).clip(roi), vis_rgb, f'Imagen {i[21:23]}/{i[23:25]}') # {i} .clip(roiChongoyape)
  
  return Map