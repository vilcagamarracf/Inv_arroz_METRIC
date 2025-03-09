<p align='center'>
   <img src="figures/Github Cover - METRIC.png" alt="Github Cover - METRIC"/>
</p>

<h4 align="center"> Aplicación de Hidrología usando Python :ear_of_rice: </h4>

<p align='center'>
   <img src="https://img.shields.io/github/last-commit/vilcagamarracf/Inv_arroz_METRIC?style=flat-square" alt="GitHub last commit"/>
</p>

## Descripción 
El presente repositorio alberga el código usado para generar mapas de  evapotranspiración, en este caso para cultivos de arroz ubicados en el distrito de Ferreñafe, usando el modelo METRIC ([Allen et al., 2007](https://www.researchgate.net/publication/228615269_Satellite-Based_Energy_Balance_for_Mapping_Evapotranspiration_With_Internalized_Calibration_METRIC_-_Model)), el cual se basa en el modelo SEBAL ([Bastiaanssen et al., 1998](https://www.sciencedirect.com/science/article/abs/pii/S0022169498002534)), cuya formulación posee como fundamento la ecuación de balance energía superficial. 

<img src="figures/Mapa Ubicacion.png" alt="Zona de trabajo" width='80%'/>

Ubicación de la región Lambayeque en la costa norte del Perú (a), ubicación de la zona de estudio en la provincia y distrito Ferreñafe, región Lambayeque (b), y su delimitación en el Fundo Zapote Figueroa (c). Ciclo de cultivo enero-junio 2022. 

<img src="figures/ET_mosaico.png" alt="ET_mosaico" width='65%'/>

Resultados obtenidos del modelo METRIC para las fechas evaluadas. 

## Proceso para generar mapas de ET
Para generar mapas de evapotranspiración (ET) con el modelo METRIC es necesario obtener los siguientes inputs:
1. **Lista de ID's de imágenes Landsat**: Para el rango temporal de evaluación se filtran las imágenes de la colección de imágenes Landsat con el objetivo de obtener una lista de ID's mediante la libreta `1_Explorar_imagenes_local.ipynb`.
2. **Datos meteorológicos**: La libreta `2_Explorar_datos_estacion.ipynb` se usa para obtener los datos meteorológicos para cada fecha y hora de evaluación (11 am. aproximadamente para imágenes Landsat en UTC-5) a partir de un archivo `.csv` de estación meteorológica.
3. **Coordenadas de pixeles fríos y pixeles calientes**: Para las fechas que dispongan de imagen satelital Landsat se requiere la selección de pixeles, proceso realizado en la libreta `3_Seleccion_pixeles.ipynb`.

Los inputs mencionados son incorporados en la libreta final de nombre `5_Modelo_Metric.ipynb`. La libreta `4_Proceso_Iterativo.ipynb` representa la replicación del archivo Excel del proceso iterativo para H.

### Código de Python

En GitHub:
- Carpeta [`notebooks`](https://github.com/vilcagamarracf/Inv_arroz_METRIC/tree/main/notebooks)

Jupyter Notebook Viewer:
- [![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.org/github/vilcagamarracf/Inv_arroz_METRIC/tree/main/notebooks/)

## Referencias

Artículos
- Allen, R. G., Tasumi, M., & Trezza, R. (2007). Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC)—Model. *Journal of irrigation and drainage engineering*, 133(4), 380-394. [http://dx.doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380)](https://www.researchgate.net/publication/228615269_Satellite-Based_Energy_Balance_for_Mapping_Evapotranspiration_With_Internalized_Calibration_METRIC_-_Model)
- Bhattarai, N., Quackenbush, L. J., Im, J., & Shaw, S. B. (2017). A new optimized algorithm for automating endmember pixel selection in the SEBAL and METRIC models. *Remote Sensing of Environment*, 196, 178-192. [https://doi.org/10.1016/j.rse.2017.05.009](https://doi.org/10.1016/j.rse.2017.05.009)
- Van der Tol, C., & Parodi, G. N. (2012). Guidelines for remote sensing of evapotranspiration. *Evapotranspiration — Remote sensing and modeling*, 227, 250. [https://doi.org/10.5772/18582](https://doi.org/10.5772/18582)

Curso
- Archibald, J., Gutierrez, H., Xu, S. (2020) Evapotranspiration. HydroLearn. [https://edx.hydrolearn.org/courses/course-v1:HumboldtState+ENGR440+2020_Fall/about](https://edx.hydrolearn.org/courses/course-v1:HumboldtState+ENGR440+2020_Fall/about)

Libros
- Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines for computing crop water requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300(9), D05109.
   - [Español](https://www.fao.org/3/x0490s/x0490s00.htm) - Guías para la determinación de los requerimientos de agua de los cultivos
   - [English](https://www.fao.org/3/X0490E/x0490e00.htm) - Guidelines for computing crop water requirements 
- Dingman, S. L. (2015). Physical hydrology. Waveland press. 
   - Información sobre el intercambio de agua y energía en superficie-atmósfera


<!-- Código de terceros
- [Thyago Anthony Soares Lima en ResearchGate](https://www.researchgate.net/publication/348906683_SEBAL_for_LANDSAT_8_Python) : Procesamiento de imágenes Landsat 8 en Python usando la metodología SEBAL para estimar evapotranspiración. -->

## Licencia/License
MIT License - Copyright (c) 2021 Cesar Francisco Vilca Gamarra \
Mayor información [aquí](https://github.com/vilcagamarracf/Inv_arroz_METRIC/blob/main/LICENSE)


## Contacto
En estas redes sociales comparto una variedad de contenidos sobre Agricultura de Precisión, programación con Earth Engine y más.

> &nbsp;&middot;&nbsp; Website [vilcagamarracf.github.io](https://vilcagamarracf.github.io/) &nbsp;&middot;&nbsp;
> LinkedIn [@cesarvilca](https://www.linkedin.com/in/cesarvilca/) &nbsp;&middot;&nbsp;
> Twitter [@vilcagamarracf](https://twitter.com/vilcagamarracf)
