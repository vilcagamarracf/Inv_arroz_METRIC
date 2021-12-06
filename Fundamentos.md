Al parecer los archivos .md en github no renderiza codigo latex.

**Problemática**

- Estimación de ET común:
  $$
  \text{ET} = \text{ET}_{ref} * K_c
  $$
  Dónde:
  - $\text{ET}$ : Evapotranspiración
  - $\text{ET}_{ref}$ : Evapotranspiración de referencia basada en el clima
  - $K_c$: Coeficiente de cultivo (de acuerdo al tipo y etapa de crecimiento)

- Dificultades de la estimación común:
  - Confiabilidad de usar valores idealizados de $K_c$ ya que se generaron para ciertas condiciones reales de crecimiento y vegetación.
  Identificación de etapas vegetativas y condiciones de crecimiento para comparación con valores de $K_c$
  - Predicción correcta de fechas según etapa de crecimiento para cultivos de grandes extensiones (ya que 

- Una solución:
  - **Los satélites son capaces de obtener información espacial de evapotranspiración de numerosas extensiones a partir de técnicas de balance de energía.**  

> Estimación de la evapotranspiración de los cultivos bajo riego mediante imágenes de sensores remotos (tipo de sensor y bandas del espectro electromagnético)

**¿Por qué medir/estimar la evapotranspiración?**

- La evapotranspiración (ET) es un componente esencial en modelos hidrológicos y de circulación general.
- Es usado para inferir la humedad del suelo, como dato de entrada para pronósticos climáticos y de inundaciones.

**Métodos de estimación**

En Operational Remote Sensing of ET and Challenges [(2011)](https://www.intechopen.com/books/evapotranspiration-remote-sensing-and-modeling/operational-remote-sensing-of-et-and-challenges) se menciona a Kalma et al. (2008), donde reunen las **metodologías existentes** para estimar la evapotranspiración de los cultivos usando Teledetección, tales como:
  - **Balance Energético de la Superficie**
  - Métodos estadísticos que utilizan diferencias entre la temperatura de la superficie y el aire
  - Correlaciones simplificadas o relaciones entre extremos de temperatura superficial de una imagen y puntos finales del ET anticipado
  - ET relativa basada en la vegetación que se multiplica por una ET de referencia basada en el tiempo

En la presente investigación, profundizaremos el método de Balance Energético de la Superficie.

**Métodos: Balance Energético de la Superficie**

Se subdivide en:
- Balance energético completo para la imagen satelital: 

  $$ \lambda{E}=R_n - G -H $$

  dónde:
  - $\lambda{E}$ es la densidad latente del flujo térmico, representando la energía "consumida" por la evaporación del agua, 
  - $R_n$ es la densidad neta del flujo de radiación, 
  - $G$ es la densidad del flujo de calor del suelo, y 
  - $H$ es la densidad sensible del flujo de calor al aire.
- Índice de estrés hídrico basado en la temperatura superficial y las cantidades de vegetación. 
- Aplicación de un Modelo Continuo de Superficie Terrestre (MSL) que se inicializa parcialmente y avanzado, en el tiempo, utilizando imágenes satelitales

Todas las metodologías anteriormente mencionadas solamente trabajan sobre imágenes disponibles y dependiendo de la revisita del mismo, quedan vacíos de información entre imágenes.
