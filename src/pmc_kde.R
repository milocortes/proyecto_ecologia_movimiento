
# ****************************
# Poligono Minimo Convexo
# ****************************

library(adehabitatHR)


# 1.2  Primero voy a covertir los datos en "SpatialPolygonsDataFrame" que es el formato que utiliza adehabitatHR para trabajar los datos y poder calcular las �reas de actividad
coords <- cbind(gallopavo$location.long,gallopavo$location.lat)

# 1.3  Ahora separo la columna con los nombre y creo un nuevo objeto pero como dataframe
id <- cbind(as.data.frame(gallopavo$individual.local.identifier))

# 1.4 Ahora creo un objeto "SpatialPointsDataFrame" utilizando las coordenadas y los nombres
gallopavo_points<-SpatialPointsDataFrame(coords,id)
ele_points

# 1.5 Ya podemos calcular los poligonos minimos convexos. Esto lo vamos a hacer al 90, 95 y 100 %. Esto
# se hace con el siguiente c�digo. Llamando al objeto que creamos nos da el �rea para cada uno de los
# individuos. En este caso en km2:
mcp_100 <- mcp(gallopavo_points[,1], unout="km2", percent=100)
mcp_100

# ****************************
# Kernel Density Estimator
# ****************************

kde_gallo <- kernelUD(gallopavo_points[,1], h ="href", kern="bivnorm")
# Con el comando image vemos los kerneles Ahora vemos los kerneles con el comando image
image(kde_gallo)
