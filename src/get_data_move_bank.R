rm(list=ls())
library(move)
library(tidyr)

# Set credential
loginStored <- movebankLogin(username="", password="")

gallopavo_all <- getMovebankData(study="Eastern Wild Turkey_Margadant", login=loginStored)
gallopavo_all

# Get track GPS data
gallopavo <- getMovebankLocationData(study="Eastern Wild Turkey_Margadant" , sensorID="GPS", login=loginStored,removeDuplicatedTimestamps=TRUE)
# Remove NAN

gallopavo <- gallopavo %>% drop_na()


# ********************************************
# Área de actividad utilizando el método AKDE.
# Autocorrelated Kernel Density Estimator (AKDE)
# *********************************************
library(ctmm)

# Convertimos estos datos en un objeto MoveStack ya que incluye a varios individuos
gallo_move <- move(x=gallopavo$location.long, y=gallopavo$location.lat,
                       time=as.POSIXct(gallopavo$timestamp, format= "%Y-%m-%d %H:%M:S", tz="UTC"),
                       proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                       data=gallopavo, animal=gallopavo$individual.local.identifier)

gallo_tel <- as.telemetry(gallo_move)

# 3.5 Con lo siguiente podemos hacer un plot de todos los individuos.
plot(gallo_tel, col=rainbow(length(gallo_tel)))

gallo_1<-variogram(gallo_tel$Sky_462)
var_gallo_1<-variogram.fit(gallo_1)

# Modelo sin autocorrelaci�n, que bas�camente equivale a ajustar o estimar el Kernel Density Estimator (KDE)
M.IID_gallo_1 <- ctmm.fit(gallo_tel$Sky_462)

# Estimamos automaticamente los mejores parametros para el modelo utilizando ctmm.guess
m.ouf_gallo_1 <- ctmm.guess(gallo_tel$Sky_462,interactive=FALSE)
A_gallo_1 <- ctmm.select(gallo_tel$Sky_462, m.ouf_gallo_1, verbose=TRUE,level=1)

# 5.3 Podemos implementar la selecci�n del modelo m�s conveniente utilizando la funci�n
# ctmm.select. Esta funci�n considera la estimaci�n inicial y posteriomente hace iteraciones de este modelo para
# seleccionar el mejor modelo basado en criterio de informaci�n (AIC)

M.OUF_gallo_1 <- ctmm.fit(gallo_tel$Sky_462, m.ouf_gallo_1)
M.OUF_gallo_1

# 6.3 Ahora lo hacemos para el mejor modelo pero utilizando el m�todo de "pesos optimos" (optimal weighting).
# Esta opci�n es adecuada cuando tenemos datos donde las localizaciones no fueron tomadas secuencialmente o hay
# varios datos p�rdidos. Para esto se utiliza el argumento weights=TRUE. El m�todo de pesos optimos estima un error
# menor, una resoluci�n m�s fina y ayuda a mitigar el efecto de alg�n sesgo en el muestreo. Puede ayudar a tener
# estimaciones m�s uniformes en sitios donde hay mayor intensidad de muestreo (un menor ancho de banda)
UD2w_gallo_1 <- akde(gallo_tel$Sky_462,M.OUF_gallo_1, weights=TRUE)

# 8. Finalmente exportamos los shapefiles y los rasters
# 8.1 Primero exportamos el shapefile del �rea de actividad calculada con el m�todo M.OUF con pesos optimos
writeShapefile (UD2w_gallo_1, "akde_Mek_Fish_", overwrite = TRUE)

# 8.2 Creamos el raster con UD del weighted OUF AKDE
gallo_1_raster<-raster(UD2w_gallo_1)
writeRaster(gallo_1_raster, filename="Mek_Fish_AKDE_w.tif", overwrite=TRUE)


build_akde <- function(animal){
  print(paste0("Ejecución para el animal ",animal@info$identity))
  print("Estimamos el modelo sin considerar la autocorrelacion")
  M.IID <- ctmm.fit(animal)
  print("Estimamos automaticamente los mejores parametros para el modelo utilizando ctmm.guess")
  m.ouf <- ctmm.guess(animal,interactive=FALSE)
  print("Seleccionamos mejor modelo de acuerdo al AIC")
  M.OUF <- ctmm.fit(animal, m.ouf)
  print("Estimamos el akde del modelo iid")
  UD0 <- akde(animal,M.IID)
  print("Estimamos el mejor modelo sin el método de pesos optimos (optimal weighting)")
  UD2 <- akde(animal,M.OUF)
  print("Estimamos el mejor modelo pero utilizando el método de pesos optimos (optimal weighting)")
  UD2w <- akde(animal,M.OUF, weights=TRUE)
  # calculate one extent for all UDs
  EXT <- extent(list(UD0,UD2,UD2w),level=0.95)

  setwd("/home/milo/PCIC/Maestría/4toSemestre/ecol_mov/proyecto/output/akde/plots")
  print("Creamos plot IID AKDE")
  cairo_ps(file = paste0("iid_kde_",animal@info$identity,".eps"),onefile = FALSE,fallback_resolution=600)
  plot(animal,UD=UD0,xlim=EXT$x,ylim=EXT$y)
  title(paste0("IID KDE-",animal@info$identity))
  dev.off()

  print("Creamos plot con UD del OUF AKDE")
  cairo_ps(file = paste0("ouf_kde_",animal@info$identity,".eps"),onefile = FALSE,fallback_resolution=600)
  plot(animal,UD=UD2,xlim=EXT$x,ylim=EXT$y)
  title(paste0("OUF AKDE-",animal@info$identity))
  dev.off()

  print("Creamos plot con UD del weighted OUF AKDE")
  cairo_ps(file = paste0("w_ouf_kde_",animal@info$identity,".eps"),onefile = FALSE,fallback_resolution=600)
  plot(animal,UD=UD2w,xlim=EXT$x,ylim=EXT$y)
  title(paste0("weighted OUF AKDE-",animal@info$identity))
  dev.off()

  #print("Exportamos el shapefile del área de actividad calculada con el método M.OUF con pesos optimos")
  #setwd("/home/milo/PCIC/Maestría/4toSemestre/ecol_mov/proyecto/output/akde/shp")
  #writeShapefile(UD2w_gallo,paste0("akde_",gallo@info$identity,"_") , overwrite = TRUE)
  setwd("/home/milo/PCIC/Maestría/4toSemestre/ecol_mov/proyecto/output/akde/raster")

  print("Creamos el raster con IID AKDE")
  animal_ud0 <-raster(UD0)
  writeRaster(animal_ud0, filename=paste0("iid_kde_",animal@info$identity,".tif"), overwrite=TRUE)

  print("Creamos el raster con UD del OUF AKDE")
  animal_ud2 <-raster(UD2)
  writeRaster(animal_ud2, filename=paste0("ouf_kde_",animal@info$identity,".tif"), overwrite=TRUE)

  print("Creamos el raster con UD del weighted OUF AKDE")
  animal_ud2w <-raster(UD2w)
  writeRaster(animal_ud2w, filename=paste0("w_ouf_kde_",animal@info$identity,".tif"), overwrite=TRUE)

}

for (gallo in gallo_tel) {
  build_akde(gallo)
}


### EVI DATA: https://developers.google.com/earth-engine/datasets/catalog/MODIS_MOD09GA_006_EVI
### Global surface water : https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_3_Metadata
### Wetness DAT : https://developers.google.com/earth-engine/datasets/catalog/Oxford_MAP_TCW_5km_Monthly#description

## recursos
# https://www.cartoscience.com/dtm
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12559
# https://ctmm-initiative.github.io/ctmm/articles/akde.html
# https://onlinelibrary.wiley.com/doi/epdf/10.1002/ece3.4823
# https://mbonnema.github.io/GoogleEarthEngine/07-SAR-Water-Classification/
