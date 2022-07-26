# Script para reproducir los análisis del trabajo:
## Análisis de patrones espaciales de puntos para el estudio de tendencias 
## locacionales en distribuciones de yacimientos arqueológicos
# Autor: Miguel Carrero Pazos
# Email: miguel.carrero.pazos@gmail.com

# Definir el espacio de trabajo
setwd("~/Desktop/barbanza_rutas/datos") #MacOSX

# Cargar los paquetes de R
spatpack<-c("raster","spatstat","rgdal","maptools", "sp", "sf")
lapply(spatpack, require, character.only=TRUE)

# Parte 1. Estudio del patrón monovariante
# Cargar datos de base y definir el área de estudio
lcps <- raster("grid/densidad_lcps_grid_1000.tiff")
area_estudio <-readOGR(dsn="shp/area_estudio.shp", layer="area_estudio")
area_w <- as.owin((area_estudio))
sites <-readOGR(dsn="shp/megalitos.shp", layer="megalitos")

## Cortar el raster de tránsito con la máscara del área de estudio
lcps_crop <- crop(lcps, extent(area_estudio))
lcps.m <- mask(lcps_crop, area_estudio)

# Crear el patrón espacial de puntos
sppp <- ppp(x=sites$UMTX, y=sites$UMTY, window=area_w)

# Preparar el raster para spatstat y crear la tendencia de primer orden
lcps_im <- as.im(as(lcps.m,"SpatialGridDataFrame"))
covlist <- list(lcps_im)
names(covlist) <- c("lcps_im")
fotrend <- ~ lcps_im

# Estudiar el patrón monovariante (rhohat)
transito.rh <- rhohat(sppp, lcps_im, confidence=0.95)

# Definir el modelo de primer orden (modelo de regresión logística)
mod1 <- step(ppm(sppp, trend = fotrend, 
                 interaction = NULL, 
                 covariates = covlist, 
                 method = "logi"))

# Crear la superficie de intensidad de primer orden
summary(mod1)
logodds <- -17.52146+(lcps_im*26245.14551)

# Generar el gráfico rhohat y la intensidad de predicción de megalitos (regresión logística) 
png(filename = "figuras/Figura 5.png",
    width = 7, height = 5, units="in",res=300)
par(mfrow=c(1,2))
plot(transito.rh, main="", xlab="Densidad", ylab="", legend=FALSE, cex.axis=0.7)
legend("topleft", legend="Tránsito", cex=0.7, bty='n', text.font=2)
plot(logodds, main="")
dev.off()

# Definir el modelo nulo (intensidad homogénea), con objetivos comparativos
mod0 <- ppm(sppp, ~1)

# Comprobar los AICS de los modelos
AIC(mod0) # Modelo nulo
AIC(mod1) # Modelo de primer orden (tránsito)

# Crear las funciones de correlación par para cada modelo
## Definir primero las simulaciones y el rango de significación (95%)
sims <- 999
nrank <- round((sims + 1) / 100 * 2.5, 0)

Pcfinhom_mod0 <- envelope(mod0, fun=pcfinhom, correction="best", rank=nrank, nsim=sims)
Pcfinhom_mod1 <- envelope(mod1, fun=pcfinhom,correction="best", rank=nrank, nsim=sims)

# Guardar las funciones de correlación par
save(Pcfinhom_mod0, file = "~/Desktop/barbanza_rutas/datos/Rdata/Pcfinhom_mod0.RData")
save(Pcfinhom_mod1, file = "~/Desktop/barbanza_rutas/datos/Rdata/Pcfinhom_mod1.RData")

# Para cargar los resultados de las funciones de correlación par, ejecutar:
load(file = "~/Desktop/barbanza_rutas/datos/Rdata/Pcfinhom_mod0.RData")
load(file = "~/Desktop/barbanza_rutas/datos/Rdata/Pcfinhom_mod1.RData")

# Crear los gráficos de correlación par
png(file="figuras/Figura 6.png",
    width=11, height=5.5, units="in",res=300)
par(mfrow=c(1,2))
par(mar=c(3,3,3,3))
plot(Pcfinhom_mod0,ylim=c(0,50), xlim=c(0, 2000),legend=FALSE,xlab="",ylab="",bty="l",main="a. Modelo nulo (CSR)")
title(ylab=expression(italic("g(r)")), line=2,cex.lab=1)
title(xlab=expression(italic("Distancia entre puntos (metros)")), line=2,cex.lab=1)
legend(500,35, legend=c("Patrón de puntos observado","Estadístico aleatorio esperado"),lty=c(1,2),col=c("black","red"),bty="n",cex=0.9)
legend(570,30.5,legend=c("Simulación de Monte Carlo"),pch=0,fill=c("grey"),col=NA,border=NA,bty="n",cex=0.9)
plot(Pcfinhom_mod1,ylim=c(0,50), xlim=c(0, 2000),legend=FALSE,xlab="",ylab="",bty="l",main="b. Modelo de primer orden\n(tránsito a través de la sierra)")
title(ylab=expression(italic("g(r)")), line=2,cex.lab=1)
title(xlab=expression(italic("Distancia entre puntos (metros)")), line=2,cex.lab=1)
legend(500,35, legend=c("Patrón de puntos observado","Estadístico aleatorio esperado"),lty=c(1,2),col=c("black","red"),bty="n",cex=0.9)
legend(570,30.5,legend=c("Simulación de Monte Carlo"),pch=0,fill=c("grey"),col=NA,border=NA,bty="n",cex=0.9)
dev.off()

# Parte 2. Simulaciones poisson no homogéneas
# Cargar datos de base y definir la nueva área de estudio reducida a las zonas altas de la sierra
lcps <- raster("grid/densidad_lcps_grid_1000.tiff")
area_estudio <-readOGR(dsn="shp/area_red.shp", layer="area_red")
area_w <- as.owin((area_estudio))
sites <-readOGR(dsn="shp/megalitos.shp", layer="megalitos")

## Cortar el raster de tránsito con la máscara del área de estudio
lcps_crop <- crop(lcps, extent(area_estudio))
lcps.m <- mask(lcps_crop, area_estudio)
lcps_im <- as.im(as(lcps.m,"SpatialGridDataFrame"))

# Generar un patrón de puntos Poisson no homogéneo influenciado por el tránsito
sppp_ipoiss <- rpoispp(lambda = predict(rhohat(sppp,
                                               lcps_im, 
                                               method="ratio",
                                               confidence=0.95, 
                                               eps=50)), nsim=999)

# Ejecutar la simulación (función de correlación par)
sim_PCF_prim_orden <- envelope(sppp,
                               W = area_w,
                               pcfinhom,
                               simulate = sppp_ipoiss,
                               nsim = 999,
                               rank = nrank,
                               correction = "best")

# Guardar el resultado de la función
save(sim_PCF_prim_orden, file = "~/Desktop/barbanza_rutas/datos/Rdata/sim_PCF_prim_orden.RData")

# Cargar el resultado de la función
load(file="~/Desktop/barbanza_rutas/datos/Rdata/sim_PCF_prim_orden.RData")

# Generar el gráfico
png(file="figuras/Figura 7.png",
    width=7, height=5, units="in",res=300)
l <- layout(matrix(c(1, 2, 3, 
                     4, 4, 4),
                   nrow = 2,
                   ncol = 3,
                   byrow = TRUE))
layout.show(l)
par(mar=c(2,2,2,2))
plot(lcps_im, main="Monumentos megalíticos", cex.main = 1.5)
plot(sppp, add=T, pch = 21, col = "white", bg="black")
plot(lcps_im, main="Simulación no. 41", cex.main = 1.5)
plot(sppp_ipoiss$`Simulation 41`, add=T, pch = 20, col = "white")
plot(lcps_im, main="Simulación no. 947", cex.main = 1.5)
plot(sppp_ipoiss$`Simulation 947`, add=T, pch = 20, col = "white")
par(mar=c(4,10,2,10))
plot(sim_PCF_prim_orden,ylim=c(0,50), xlim=c(0, 2000),legend=FALSE,xlab="",ylab="",bty="l",main="")
title(ylab=expression(italic("g(r)")), line=2,cex.lab=1)
title(xlab=expression(italic("Distancia entre puntos (metros)")), line=2,cex.lab=1)
legend(500,35, legend=c("Patrón de puntos observado","Estadístico aleatorio esperado"),lty=c(1,2), col=c("black","red"),bty="n",cex=1)
legend(555,27,legend=c("Simulación de Monte Carlo"),pch=0,fill=c("grey"),col=NA,border=NA,bty="n",cex=1)
dev.off()

