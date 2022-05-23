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
