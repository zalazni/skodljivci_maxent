setwd("C:/Users/ZalaZn/OneDrive - Univerza v Ljubljani/DR/Ecoclimatic niche modeling")
# setwd("C:/Users/zalaz/Work Folders/ZALA BF delo/MR/DR/Ecoclimatic niche modeling")

# install.packages("maxnet")
# install.packages("terra")
# library(terra)

library(maxnet)
library(dismo) 
library(maps)
library(mapdata)
library(rJava) 
# library(maptools)
library(ggplot2)
library(jsonlite)
library(pracma)
library(spData)
library(sf)
library(dplyr)
library(rasterVis)
library(ggspatial)
library(raster)
library(ncdf4)
library(sp)
library(caret)
library(rgdal)
library(pROC)
library(verification)
library(scales)
library(RColorBrewer)

zuzelka1 = 'Spodoptera frugiperda'
# zuzelka1 = 'Trichogramma pretiosum'
# zuzelka1 = 'Cotesia marginiventris'
# zuzelka1 = 'Telenomus remus'
# zuzelka1 = 'Eiphosoma laphygmae'

rattler<-gbif(zuzelka1)
rattler<-subset(rattler, !is.na(lon) & !is.na(lat))
rattlerdups=duplicated(rattler[, c("lon", "lat")])
rattler <-rattler[!rattlerdups, ]

# path <- "C:/Users/zalaz/Work Folders/ZALA BF delo/MR/DR/podatki_hist/OPSI/"
path <- "C:/Users/ZalaZn/OneDrive - Univerza v Ljubljani/DR/podatki_hist/OPSI/" # sluzba

# nuts3_mapdata <- st_read("C:/Users/zalaz/Work Folders/ZALA BF delo/MR/DR/NUTS3_ID.gpkg")
nuts3_mapdata <- st_read("C:/Users/ZalaZn/OneDrive - Univerza v Ljubljani/DR/NUTS3_ID.gpkg")
slovenia_nuts3_mapdata <- filter(nuts3_mapdata, cntr_code=="SI") # filtriram samo slovenske NUTS3 regije

# # narisemo, kje se pojavlja nasa zuzelka
data(wrld_simpl)
plot(wrld_simpl, xlim=c(min(rattler$lon)-1,max(rattler$lon)+1), ylim=c(min(rattler$lat)-1,max(rattler$lat)+1), axes=TRUE, col="light yellow")
points(rattler$lon, rattler$lat, col="lightblue", pch=20, cex=0.75)
save.image(paste0("occurrence_",zuzelka1,"_GBIF.png", sep=""))

# klimatski podatki

currentEnv=getData("worldclim", var="bio", res=2.5)
futureEnv=getData('CMIP5', var='bio', res=2.5, rcp=85, model='HE', year=70)
names(futureEnv)=names(currentEnv)
# 4km resolution (0.04166667x0.04166667 decimal degrees)

# https://pjbartlein.github.io/REarthSysSci/raster_intro.html
currentEnv=dropLayer(currentEnv, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))
futureEnv=dropLayer(futureEnv, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))
# Napovedne klimatske spremenljivke:
# BIO1 = Annual Mean Temperature
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO12 = Annual Precipitation
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

# naredimo ovojnico 10 stopinj okoli nasih podatkov, za ta obmocja potem
# "izrezemo" tudi klimatske podatke v naslednji vrstici
model.extent<-extent(min(rattler$lon)-10,max(rattler$lon)+10,min(rattler$lat)-10,max(rattler$lat)+10)
modelEnv=crop(currentEnv,model.extent)
modelFutureEnv=crop(futureEnv, model.extent)

plot(modelEnv[["bio1"]]/10, main="Annual Mean Temperature")
map('worldHires',xlim=c(min(rattler$lon)-10,max(rattler$lon)+10), ylim=c(min(rattler$lat)-10,max(rattler$lat)+10), fill=FALSE, add=TRUE)
points(rattler$lon, rattler$lat, pch=20, cex=0.75)

# Naredimo model na 80 % podatkov in s CV preverimo na 20 % podatkov, kako dobro napove
rattlerocc=cbind.data.frame(rattler$lon,rattler$lat) #first, just make a data frame of latitudes and longitudes for the model
fold <- kfold(rattlerocc, k=5) # add an index that makes five random groups of observations
# rattlertest <- rattlerocc[fold == 1, ] # hold out one fifth as test data
# rattlertrain <- rattlerocc[fold != 1, ] # the other four fifths are training data
# # In real applications, since the particulars of the model depend on the data used to fit it, we would actually fit the model multiple times, withholding each fifth of the data separately,
# # then average the results. This is called k-fold cross-validation (in our case 5-fold). However, for our purposes here, we will just fit the model once.
#
# # Now we can fit the SDM using the Maximum Entropy (Maxent) algorithm,
# # which tries to define the combination of environmental responses that best predicts the occurrence of the species.
# rattler.me <- maxent(modelEnv, rattlertrain) #note we just using the training data
# plot(rattler.me) # vidimo, da je dalec najpomembnejsa bio12
# # ter delno bio1 = letna T, bio16 = padavine najvlazn. cetrt. in bio19 = padavine najhlad. cetrt.
# # kako je verjetnost pojava zuzelke odvisna od klim. spremenljivk?
# response(rattler.me)# pri visji letni T (bio1 > 300 K) se mocno poveca verjetnost, podobno z letnimi padavinami (bio12) > 2000
# # Napovemo verjetnost za pojav zuzelke
# rattler.pred <- predict(rattler.me, modelEnv)
# 
# # The area under this curve, which varies from zero to one, provides an assessment of the model. An AUC value of 0.5 is the same as random guessing of presence/absence, while values towards one mean our predictions are more reliable. To generate and evaluate the AUC for our model,
# # we first generate background points for pseudoabsences. The ?randomPoints()? function even makes
# # sure that the points occur only in areas where the predictor variables exist. That is, none of our points will occur in the ocean.
bg <- randomPoints(modelEnv, 1000)
# # Then we can use ?evaluate()? to generate several diagnostics as well as the AUC, using our test data as our presences against our pseudoabsences.
# # e1 <- evaluate(rattler.me, p=rattlertest, a=bg, x=modelEnv)
# # plot(e1, 'ROC') # precej stran od 0.5, mnogo bolje napove kot random
# 
# current <- plot(rattler.pred, main="Predicted Current Suitability", xlim=c(13.3,16.7), ylim=c(45.3,47))
# points(rattler$lon, rattler$lat, pch=20, cex=0.75)
# # plot(rattler.pred, main="Predicted Current Suitability", xlim=c(130,150), ylim=c(-10,-8))
# # rattler.2070 = predict(rattler.me, modelFutureEnv)# napoved za prihodnosti
# # saveRDS(rattler.2070, file = "rattler_2070.rds")
# rattler.2070 <- readRDS(file = "rattler_2070.rds")
# 
# plot(rattler.2070, main="Predicted Future Suitability", xlim=c(13.3,16.7), ylim=c(45.3,47))
# # plot(rattler.2070, main="Predicted Future Suitability", xlim=c(130,150), ylim=c(-10,-8))
# # map('worldHires', fill=FALSE, add=TRUE)
# # points(rattler$lon, rattler$lat, pch=20, cex=0.75)
# 
# rattler.2070 <- data.frame(rasterToPoints(rattler.2070))
# 
# ggplot() +
#   geom_raster(data=rattler.2070,aes(x = x, y = y, fill = layer)) +
#   geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
#   # scale_fill_distiller("ET", palette = "Spectral") +
#   labs(x="Longitude",y="Latitude",title="Spodoptera predicted future suitability, RCP8.5, 2070") +
#   coord_sf(crs = st_crs(4326)) + xlim(c(13.3,16.7)) + ylim(c(45.3,47)) +
#   theme_light() + theme(legend.position = c(.9, .2), 
#                         legend.title = element_text(face = "bold"),
#                         text = element_text(size=18))

# za vsako od petih sort moram narediti ?e validacijo!!!
auc_c <- c()
p_c <- c()
k <- 5
# Perform k-fold cross-validation
for (i in 1:k) {
  rattlertest <- rattlerocc[fold == i, ]
  rattlertrain <- rattlerocc[fold != i, ]
  rattler.me <- maxent(modelEnv, rattlertrain)
  # rattler.pred <- predict(rattler.me, modelEnv)
  # rattler.test_pred <- predict(rattler.me, rattlertest, type = "response")
  
  e1 <- evaluate(rattler.me, p=rattlertest, a=bg, x=modelEnv) #https://rpubs.com/mlibxmda/GEOG70922_Week5
  # roc_result <- roc(e1@presence, rattler.pred@data@values)
  # # Perform the hypothesis test using roc.test
  # roc_test_result <- roc.test(roc_result, method = "bootstrap", R = 1000) 
  # # You can adjust the number of bootstrap samples (R)
  # 
  # # Print the results
  # print(roc_test_result)
  # # Specify the theoretical AUC value
  # theoretical_auc <- 0.5  # Replace with your desired theoretical AUC value
  # # Create a reference ROC curve with the theoretical AUC value
  # reference_roc <- roc(c(0, 1), c(0, 1), auc = theoretical_auc)
  # # Perform the roc.test to compare the observed ROC curve to the reference ROC curve
  # roc_test_result <- roc.test(e1, reference_roc, method = "delong")
  # # Extract the p-value
  # p_value <- roc_test_result$p.value
  
  plot(e1, 'ROC')
  auc_c <- c(auc_c,e1@auc)
  # p_c <- c(p_c,p_value)
  print(i)
}

# SEDAJ pa naredimo polni model, na vseh podatkih!
rattler.me <- maxent(modelEnv, rattlerocc)
e1 <- evaluate(rattler.me, p=rattlerocc, a=bg, x=modelEnv)
plot(e1, 'ROC')

saveRDS(rattler.me, file = paste0("rattler_me",zuzelka1,".rds"))
rattler.me <- readRDS(file = paste0("rattler_me",zuzelka1,".rds"))
plot(rattler.me) 
save.image(paste0("slike_clanek/var_contr_",zuzelka1,".png", sep=""))

response(rattler.me)
save.image(paste0("response_",zuzelka1,"_GBIF.png", sep=""))



plot(modelEnv[["bio19"]], main="bio19 Coldest qrt prec") # okoli 500 pas
map('worldHires',xlim=c(min(rattler$lon)-10,max(rattler$lon)+10), ylim=c(min(rattler$lat)-10,max(rattler$lat)+10), fill=FALSE, add=TRUE)
points(rattler$lon, rattler$lat, pch=20, cex=0.75)

plot(modelEnv[["bio12"]], main="bio12 annual prec") # okoli 4000 step func
xlim(min(rattler$lon)-2,max(rattler$lon)+2)
map('worldHires',xlim=c(min(rattler$lon)-2,max(rattler$lon)+2), ylim=c(min(rattler$lat)-2,max(rattler$lat)+2), fill=FALSE, add=TRUE)

# rattler.pred <- predict(rattler.me, modelEnv)
# saveRDS(rattler.pred, file = "rattler_pred.rds")
# rattler.pred <- readRDS(file = "rattler_pred.rds")








# plot(rattler.2070, main="Predicted Future Suitability", xlim=c(13.3,16.7), ylim=c(45.3,47))
# map('worldHires', fill=FALSE, add=TRUE)
# points(rattler$lon, rattler$lat, pch=20, cex=0.75)

############

# # poskusimo z modelom, ki ima le glavnih 6 napovednih spremenljvik
# currentEnv_6=dropLayer(currentEnv, c("bio9", "bio8", "bio6", "bio18", "bio12"))
# futureEnv_6=dropLayer(futureEnv, c("bio9", "bio8", "bio6", "bio18", "bio12"))
# 
# modelEnv_6=crop(currentEnv_6,model.extent)
# modelFutureEnv_6=crop(futureEnv_6, model.extent)
# 
# rattler.me_6 <- maxent(modelEnv_6, rattlertrain) 
# plot(rattler.me_6)
# 
# e6 <- evaluate(rattler.me_6, p=rattlertest, a=bg, x=modelEnv_6)
# plot(e6, 'ROC')
# 
# rattler.2070_6 = predict(rattler.me_6, modelFutureEnv_6)
# plot(rattler.2070_6, main="Predicted Future Suitability - 6 vars")
# map('worldHires', fill=FALSE, add=TRUE)
# points(rattler$lon, rattler$lat, pch="+", cex=0.2)
# 
# # narediti model na na?ih podatkih, na vi?ji resoluciji!!!
# path <- "C:/Users/zalaz/Work Folders/ZALA BF delo/MR/DR/Ecoclimatic niche modeling/wc2.1_30s_bio/"
# current_hires <- raster(paste0(path,"wc2.1_30s_bio_1.tif"))
# 
# prec_futureEnv=getData('worldclim', var='prec', res=0.5, lon=5, lat=45)
# tmin_futureEnv=getData('worldclim', var='tmin', res=0.5, lon=5, lat=45)
# tmax_futureEnv=getData('worldclim', var='tax', res=0.5, lon=5, lat=45)
# b <- biovars(prec, tmin, tmax, ...)
# as.matrix(b)


########################  Historicni podatki za Slovenijo

nc_data <- nc_open(paste0(path,"evspsblpot","_12km_ARSO_v5_day_19810101_20101231.nc"))
#nc_data <- nc_open("dataset-insitu-gridded-observations-europe/tn_ens_mean_0.1deg_reg_v27.0e.nc")
names(nc_data$var)
print(nc_data)
lon <- ncvar_get(nc_data, "X") # drugi parameter je ime spremenljivke v datoteki
lat <- ncvar_get(nc_data, "Y")
start_lon <- 13.2 # meje https://sl.wikipedia.org/wiki/Geografija_Slovenije
end_lon <- 16.8

start_lat<-45.3 #+ 25/60 +18.34/3600
end_lat<-47 #+ 54/60 + 37.52/3600

time <- ncvar_get(nc_data, "time")
nx = length(lon)
ny = length(lat)
start1<-c(1,1,1)
count1<-c(nx,ny, 10957)
# povpre?na dnevna T
datumi0 <- seq(0,count1[3]-1,by = 1)
pretvorba_sek_v_dan1 = 60*60*24 #?eprav imam dnevne podatke moram to dat, ker as.positxct meri v sekundah
datumi <- as.POSIXct((start1[3]+datumi0)*pretvorba_sek_v_dan1,origin="1980-12-31 00:00:00")
leta <- format(datumi,format="%Y")
datum2 = format(as.POSIXct(datumi, format="%Y-%m-%d"), "%m/%d/%Y")
nc_close(nc_data) # konec branja

# minimalna dnevna T
nc_data_min <- nc_open(paste0(path,"tasmin","_12km_ARSO_v5_day_19810101_20101231.nc"))
data_min <- ncvar_get(nc_data_min, "tasmin", start=start1, count=count1)
# dataMIN<-structure(data_min, .Names = datum2)
nc_close(nc_data_min) # konec branja

# maksimalna dnevna T
nc_data_max<- nc_open(paste0(path,"tasmax","_12km_ARSO_v5_day_19810101_20101231.nc"))
data_max <- ncvar_get(nc_data_max, "tasmax", start=start1, count=count1)
# dataMAX<-structure(data_max, .Names = datum2)
nc_close(nc_data_max) # konec branja

# padavine dnevne
nc_data_rr<- nc_open(paste0(path,"pr","_12km_ARSO_v5_day_19810101_20101231.nc"))
data_rr <- na.omit(ncvar_get(nc_data_rr, "pr", start=start1, count=count1))
# dataRR<-structure(data_rr, .Names = datum2)
# data_rr <- data.frame(datum = datumi, rr = data_rr)
nc_close(nc_data_rr) # konec branja

lats1 <- seq(1,24,by=1)
lons1 <- seq(1,40,by=1)

# Define the dimensions of the 3D matrix
n_rows <- 24
n_cols <- 40
n_layers <- 19

# Create a 3D matrix with zeros (you can change the fill value if desired)
bios <- array(0, dim = c(n_rows, n_cols, n_layers))
leto = format(datumi, format = "%Y")
leta <- seq(min(leto),max(leto),by=1)

for(j in lats1){
  for(i in lons1){
    print(i)
    print(j)
    data_rr1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), rr = data_rr[i,j,1:length(time)]*24*60*60) # ARSO podatki so v mm/s zato za dnevne pomno?im s ?tevilom s v dnevu
    data_min1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tn = (data_min[i,j,1:length(time)]-273)) # ARSO podatki so v Kelvinih
    data_max1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tx = (data_max[i,j,1:length(time)]-273))

    # mese?ne vsote padavin in mese?na povpre?ja temperature
    dataRR <- aggregate(rr ~ mes + leto, data_rr1 , FUN = sum , na.rm=TRUE, na.action=na.pass)
    dataMIN <- aggregate(tn ~ mes + leto, data_min1 , mean , na.rm=TRUE, na.action=na.pass)
    dataMAX <- aggregate(tx ~ mes + leto, data_max1, mean , na.rm=TRUE, na.action=na.pass)
    
    # povpre?je let obdobja 2071-2100
    dataRR <- aggregate(rr ~ leto, dataRR , mean , na.rm=TRUE, na.action=na.pass)
    dataMIN <- aggregate(tn ~ leto, dataMIN , mean , na.rm=TRUE, na.action=na.pass)
    dataMAX <- aggregate(tx ~ leto, dataMAX, mean , na.rm=TRUE, na.action=na.pass)

    # tukaj mi mora vrniti NA za biohistSI, ?e so vshodni trije NA!!!!!!!!!!!!
    bio_histSI <- biovars(prec = dataRR$rr, tmin = dataMIN$tn, tmax = dataMAX$tx)
    # tukaj naredim za vsako to?ko vrstico 19-ih spremenljivk bio1 - bio19

    bios[j,i,1:19] = bio_histSI

  }}

bios <- bios[1:24, 1:40, c(1,5:9,12,16:19)]

# tukaj odpiram plasti datoteke bios
# rabim pa RasterStack z 19 plastmi za bio1 - bio19
# nc_file <-  nc_open(paste0(path,"pr","_12km_ARSO_v5_day_19810101_20101231.nc"))
# nc_data <- ncvar_get(nc_file, varid = "pr")
# nrows <- length(ncvar_get(nc_file, "Y"))
# ncols <- length(ncvar_get(nc_file, "X"))
nlayers <- dim(bios)[3]  # nlayers mora biti lat_dim x lon_dim x BIO1-BIO19

# Create an empty RasterStack object
raster_stack <- stack()

# Loop through the layers and add them to the RasterStack
for (i in 1:nlayers) {
  layer <- matrix(bios[, , i], nrow = n_rows, ncol = n_cols)
  raster_layer <- raster(layer)
  extent_of_raster <- extent(min(lon), max(lon), min(lat), max(lat))
  # Assign the extent to the raster layer
  extent(raster_layer) <- extent_of_raster
  # Define the projection of the raster layer (WGS84 in this example)
  projection(raster_layer) <- CRS("+proj=longlat +datum=WGS84")
  raster_stack <- addLayer(raster_stack, raster_layer)
}

currentEnvSI <- raster_stack
layer_names <- c("bio1", "bio5", "bio6", "bio7", "bio8", "bio9", "bio12", "bio16", "bio17", "bio18", "bio19")
names(currentEnvSI) <- layer_names
modelEnvSI=crop(currentEnvSI,extent_of_raster)

plot(modelEnvSI[["bio1"]]/10, main="Annual Mean Temperature")
plot(modelEnv[["bio1"]]/10, main="Annual Mean Temperature", xlim=c(13.3,16.7), ylim=c(45.3,47))
map('worldHires',xlim=c(min(rattler$lon)-10,max(rattler$lon)+10), ylim=c(min(rattler$lat)-10,max(rattler$lat)+10), fill=FALSE, add=TRUE)
points(rattler$lon, rattler$lat, pch=20, cex=0.75)

plot(modelEnvSI[["bio5"]]/10, main="Max Temperature of Warmest Month")
plot(modelEnv[["bio5"]]/10, main="Max Temperature of Warmest Month", xlim=c(13.3,16.7), ylim=c(45.3,47))

plot(modelEnvSI[["bio6"]]/10, main="Min Temperature of Coldest Month")
plot(modelEnv[["bio6"]]/10, main="Min Temperature of Coldest Month", xlim=c(13.3,16.7), ylim=c(45.3,47))

plot(modelEnvSI[["bio7"]]/10, main="Temperature Annual Range (BIO5-BIO6)")
plot(modelEnv[["bio7"]]/10, main="Temperature Annual Range (BIO5-BIO6)", xlim=c(13.3,16.7), ylim=c(45.3,47))

plot(modelEnvSI[["bio8"]]/10, main="Mean Temperature of Wettest Quarter")
plot(modelEnv[["bio8"]]/10, main="Mean Temperature of Wettest Quarter", xlim=c(13.3,16.7), ylim=c(45.3,47))

plot(modelEnvSI[["bio9"]]/10, main="Mean Temperature of Driest Quarter")
plot(modelEnv[["bio9"]]/10, main="Mean Temperature of Driest Quarter", xlim=c(13.3,16.7), ylim=c(45.3,47))

plot(modelEnvSI[["bio12"]], main="Annual Precipitation")
plot(modelEnv[["bio12"]], main="Annual Precipitation", xlim=c(10,19), ylim=c(40,49))
plot(modelEnv[["bio12"]], main="Annual Precipitation", xlim=c(130,150), ylim=c(-10,-8))

plot(modelEnvSI[["bio16"]], main="Precipitation of Wettest Quarter")
plot(modelEnv[["bio16"]], main="Precipitation of Wettest Quarter", xlim=c(-80,-70), ylim=c(0,10))

plot(modelEnvSI[["bio17"]], main="Precipitation of Driest Quarter")
plot(modelEnv[["bio17"]], main="Precipitation of Driest Quarter", xlim=c(13.3,16.7), ylim=c(45.3,47))

plot(modelEnvSI[["bio18"]], main="Precipitation of Warmest Quarter")
plot(modelEnv[["bio18"]], main="Precipitation of Warmest Quarter", xlim=c(13.3,16.7), ylim=c(45.3,47))

plot(modelEnvSI[["bio19"]], main="Precipitation of Coldest Quarter")
plot(modelEnv[["bio19"]], main="Precipitation of Coldest Quarter", xlim=c(13.3,16.7), ylim=c(45.3,47))

rattler.me <- readRDS(file = paste0("rattler_me",zuzelka1,".rds"))
rattler.predSI <- predict(rattler.me, modelEnvSI)
plot(rattler.predSI, main="Predicted Suitability")

rattler.pred_gg <- data.frame(rasterToPoints(rattler.predSI))
start_lon <- 13.2 # meje https://sl.wikipedia.org/wiki/Geografija_Slovenije
end_lon <- 16.8

start_lat<-45.3 #+ 25/60 +18.34/3600
end_lat<-47 #+ 54/60 + 37.52/3600

ggplot() +
  geom_raster(data=rattler.pred_gg,aes(x = x, y = y, fill = layer)) +
  geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
  scale_fill_distiller("Probability", palette = "Spectral", limits = c(0,1)) +
  labs(x="Longitude",y="Latitude",title=paste0("Predicted suitability for ",zuzelka1,", 1981-2010")) +
  coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
  theme_light() + theme(legend.position = c(.9, .2), 
                        legend.title = element_text(face = "bold"),
                        text = element_text(size=18))
ggsave(paste0("hist/Predicted_suitability_1981-2010_",zuzelka1,".png"),width = 10, height = 7)
 

plot(shape)# # #nc_data <- nc_open("dataset-insitu-gridded-observations-europe/tn_ens_mean_0.1deg_reg_v27.0e.nc")
names(nc_data$var)
print(nc_data)
lonlat <- ncvar_get(nc_data, "crs") # drugi parameter je ime spremenljivke v datoteki
# lonlat
nc_data$dim$longitude$vals
nc_data$dim$latitude$vals











# 
# 
# ################# PROJEKCIJE
# 
# model <- c("CNRM-CERFACS-CNRM-CM5", 
#             "ICHEC-EC-EARTH",
#             "IPSL-IPSL-CM5A-MR",
#             "MOHC-HadGEM2-ES", # samo do 30. 11. 2099 pri padavinah, zato tudi pri Tmin in Tmax .nc fajl spremenim, da gre le do 30.11.
#             "MPI-M-MPI-ESM-LR",
#             "MPI-M-MPI-ESM-LR")
# rcp <- "rcp85"
# drugo <- c("_r1i1p1_CLMcom-CCLM4-8-17_v1",
#            "_r3i1p1_DMI-HIRHAM5_v1",
#            "_r1i1p1_IPSL-INERIS-WRF331F_v1",
#            "_r1i1p1_KNMI-RACMO22E_v2",
#            "_r1i1p1_CLMcom-CCLM4-8-17_v1",
#            "_r1i1p1_SMHI-RCA4_v1a")
# leto_zac <- "2071"
# leto_kon <- c("21001231","21001231","21001231","20991130","21001231","21001231")
# z=4
# leto_zac1 <- "2070"
# 
# # padavine dnevne
# nc_data_rr<- nc_open(paste0(path,"pr","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
# data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
# time <- ncvar_get(nc_data_rr, "time")
# start1<-c(1,1,1)
# count1<-c(nx,ny, length(time))
# # povpre?na dnevna T
# datumi0 <- seq(0,count1[3]-1,by = 1)
# pretvorba_sek_v_dan1 = 60*60*24 #?eprav imam dnevne podatke moram to dat, ker as.positxct meri v sekundah
# datumi <- as.POSIXct((start1[3]+datumi0)*pretvorba_sek_v_dan1,origin=paste0(leto_zac1,"-12-31 00:00:00"))
# leta <- format(datumi,format="%Y")
# datum2 = format(as.POSIXct(datumi, format="%Y-%m-%d"), "%m/%d/%Y")
# nc_close(nc_data_rr) # konec branja
# # maksimalna dnevna T
# nc_data_max<- nc_open(paste0(path,"tasmax","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
# data_max <- ncvar_get(nc_data_max, "tasmax", start=start1, count=count1)
# nc_close(nc_data_max) # konec branja
# # minimalna dnevna T
# nc_data_min <- nc_open(paste0(path,"tasmin","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
# data_min <- ncvar_get(nc_data_min, "tasmin", start=start1, count=count1)
# nc_close(nc_data_min) # konec branja
# 
# leto = format(datumi, format = "%Y")
# leta1 <- seq(min(leto),max(leto),by=1)
# lats1 <- seq(1,24,by=1)
# lons1 <- seq(1,40,by=1)
# n_rows <- 24
# n_cols <- 40
# n_layers <- 19
# bios2070 <- array(0, dim = c(n_rows, n_cols, n_layers))
# 
# # IZBRANO LETO
# l = 2075
# 
# for(j in lats1){
#   for(i in lons1){
#     print(i)
#     print(j) # podatki za T so v worldclim mno?eni z 10, zato da nimajo decimalnih mest, zato tudi jaz svoje mno?im, da bom lahko uporabila model
#     data_rr1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), rr = data_rr[i,j,1:length(time)]*24*60*60) # ARSO podatki so v mm/s zato za dnevne pomno?im s ?tevilom s v dnevu
#     data_min1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tn = 10*(data_min[i,j,1:length(time)]-273)) # ARSO podatki so v Kelvinih
#     data_max1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tx = (data_max[i,j,1:length(time)]-273)*10)
#     
#     dataRR <- aggregate(rr ~ mes + leto, data_rr1 , FUN = sum , na.rm=TRUE, na.action=na.pass)
#     dataMIN <- aggregate(tn ~ mes + leto, data_min1 , mean , na.rm=TRUE, na.action=na.pass)
#     dataMAX <- aggregate(tx ~ mes + leto, data_max1, mean , na.rm=TRUE, na.action=na.pass)
#     
#     dataRR <- dataRR[dataRR$leto == l, ]
#     dataMIN <- dataMIN[dataMIN$leto == l, ]
#     dataMAX <- dataMAX[dataMAX$leto == l, ]
#     bio_histSI <- biovars(prec = dataRR$rr, tmin = dataMIN$tn, tmax = dataMAX$tx)
#     bios2070[j,i,1:19] = bio_histSI
#   }}
# 
# bios2070 <- bios2070[1:24, 1:40, c(1,5:9,12,16:19)]
# nlayers <- dim(bios2070)[3]  # nlayers mora biti lat_dim x lon_dim x BIO1-BIO19
# 
# raster_stack2070 <- stack(
#   lapply(1:nlayers, function(i) {
#     layer <- matrix(bios2070[, , i], nrow = n_rows, ncol = n_cols)
#     layer <- flipdim(layer, 1)
#     raster_layer <- raster(layer)  # Define raster layer inside the loop
#     extent_of_raster <- extent(min(lon), max(lon), min(lat), max(lat))
#     # Assign the extent to the raster layer
#     extent(raster_layer) <- extent_of_raster
#     # Define the projection of the raster layer (WGS84 in this example)
#     projection(raster_layer) <- CRS("+proj=longlat +datum=WGS84")
#     return(raster_layer)
#   })
# )
# 
# 
# futureEnvSI <- raster_stack2070
# layer_names <- c("bio1", "bio5", "bio6", "bio7", "bio8", "bio9", "bio12", "bio16", "bio17", "bio18", "bio19")
# names(futureEnvSI) <- layer_names
# modelFutureEnvSI=crop(futureEnvSI, model.extent)
# # plot(modelFutureEnvSI[["bio1"]]/10, main="Annual Mean Temperature")
# # plot(modelFutureEnvSI[["bio12"]], main="Annual Precipitation")
# # plot(modelFutureEnvSI[["bio16"]], main="Precipitation of Wettest Quarter")
# 
# rattler.2070SI <- predict(rattler.me, modelFutureEnvSI)
# # plot(rattler.2070SI, main="Predicted Future Suitability ARSO")
# 
# rattler.2070SI_gg <- data.frame(rasterToPoints(rattler.2070SI))
# xi <- seq(min(rattler.2070SI_gg$x),max(rattler.2070SI_gg$x),length.out = 10)
# yi <- seq(min(rattler.2070SI_gg$y),max(rattler.2070SI_gg$y),length.out = 6)
# lon <- round(seq(min(lon),max(lon),length.out = 10), digits = 2)
# lat <- round(seq(min(lat),max(lat),length.out = 6), digits = 2)
# 
# ggplot(data=rattler.2070SI_gg) +
#   geom_raster(aes(x = x, y = y, fill = layer)) +
#   scale_fill_distiller("Probability", palette = "Spectral", limits=c(0,1)) +
#   labs(x="Longitude",y="Latitude",title=paste0("Predicted suitability for ",zuzelka,", ",l,", model ",model[z]," and ",rcp)) +
#   #scale_x_continuous(breaks = xi, labels = lon) + scale_y_continuous(breaks = yi,labels=lat) +
#   theme_light() + theme(legend.position = c(.9, .2), 
  # legend.title = element_text(face = "bold"),
  # text = element_text(size=18))
# ggsave(paste0("proj/Predicted_suitability_",l,"_",model[z],"_",rcp,".png"),width = 10, height = 7)
# 
# ########## RPOJEKCIJE FUNKCIJA ZA RISANJE VEC LET, MODELOV in SCENARIJ RCP4.5
# leta <- c(2075, 2080, 2085, 2090, 2095)
# 
# # for(l in leta){ 
#   model <- c("CNRM-CERFACS-CNRM-CM5", 
#              "ICHEC-EC-EARTH",
#              "IPSL-IPSL-CM5A-MR",
#              "MOHC-HadGEM2-ES", # samo do 30. 11. 2099 pri padavinah, zato tudi pri Tmin in Tmax .nc fajl spremenim, da gre le do 30.11.
#              "MPI-M-MPI-ESM-LR",
#              "MPI-M-MPI-ESM-LR")
#   rcp <- "rcp45"
#   drugo <- c("_r1i1p1_CLMcom-CCLM4-8-17_v1",
#              "_r3i1p1_DMI-HIRHAM5_v1",
#              "_r1i1p1_IPSL-INERIS-WRF331F_v1",
#              "_r1i1p1_KNMI-RACMO22E_v2",
#              "_r1i1p1_CLMcom-CCLM4-8-17_v1",
#              "_r1i1p1_SMHI-RCA4_v1a")
#   leto_zac <- "2071"
#   leto_kon <- c("21001231","21001231","21001231","20991130","21001231","21001231")
#   z=4
#   leto_zac1 <- "2070"
#   
#   # padavine dnevne
#   nc_data_rr<- nc_open(paste0(path,"pr","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
#   data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
#   time <- ncvar_get(nc_data_rr, "time")
#   start1<-c(1,1,1)
#   count1<-c(nx,ny, length(time))
#   # povpre?na dnevna T
#   datumi0 <- seq(0,count1[3]-1,by = 1)
#   pretvorba_sek_v_dan1 = 60*60*24 #?eprav imam dnevne podatke moram to dat, ker as.positxct meri v sekundah
#   datumi <- as.POSIXct((start1[3]+datumi0)*pretvorba_sek_v_dan1,origin=paste0(leto_zac1,"-12-31 00:00:00"))
#   leta <- format(datumi,format="%Y")
#   datum2 = format(as.POSIXct(datumi, format="%Y-%m-%d"), "%m/%d/%Y")
#   nc_close(nc_data_rr) # konec branja
#   # maksimalna dnevna T
#   nc_data_max<- nc_open(paste0(path,"tasmax","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
#   data_max <- ncvar_get(nc_data_max, "tasmax", start=start1, count=count1)
#   nc_close(nc_data_max) # konec branja
#   # minimalna dnevna T
#   nc_data_min <- nc_open(paste0(path,"tasmin","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
#   data_min <- ncvar_get(nc_data_min, "tasmin", start=start1, count=count1)
#   nc_close(nc_data_min) # konec branja
#   
#   leto = format(datumi, format = "%Y")
#   # leta <- seq(min(leto),max(leto),by=1)
#   lats1 <- seq(1,24,by=1)
#   lons1 <- seq(1,40,by=1)
#   n_rows <- 24
#   n_cols <- 40
#   n_layers <- 19
#   bios2070 <- array(0, dim = c(n_rows, n_cols, n_layers))
#     
# 
#   # print(l)
#   bios <- array(0, dim = c(n_rows, n_cols, n_layers))
#   for(j in lats1){
#     for(i in lons1){
#       data_rr1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), rr = data_rr[i,j,1:length(time)]*24*60*60) # ARSO podatki so v mm/s zato za dnevne pomno?im s ?tevilom s v dnevu
#       data_min1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tn = 10*(data_min[i,j,1:length(time)]-273)) # ARSO podatki so v Kelvinih
#       data_max1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tx = (data_max[i,j,1:length(time)]-273)*10)
#       
#       dataRR <- aggregate(rr ~ mes, data_rr1 , FUN = sum , na.rm=TRUE, na.action=na.pass)
#       dataMIN <- aggregate(tn ~ mes, data_min1 , mean , na.rm=TRUE, na.action=na.pass)
#       dataMAX <- aggregate(tx ~ mes, data_max1, mean , na.rm=TRUE, na.action=na.pass)
#       
#       # dataRR <- aggregate(rr ~ mes + leto, data_rr1 , FUN = sum , na.rm=TRUE, na.action=na.pass)
#       # dataMIN <- aggregate(tn ~ mes + leto, data_min1 , mean , na.rm=TRUE, na.action=na.pass)
#       # dataMAX <- aggregate(tx ~ mes + leto, data_max1, mean , na.rm=TRUE, na.action=na.pass)
#       # 
#       # dataRR <- dataRR[dataRR$leto == l, ]
#       # dataMIN <- dataMIN[dataMIN$leto == l, ]
#       # dataMAX <- dataMAX[dataMAX$leto == l, ]
#       bio_histSI <- biovars(prec = dataRR$rr, tmin = dataMIN$tn, tmax = dataMAX$tx)
#       bios2070[j,i,1:19] = bio_histSI
#     }}
#   
#   bios2070 <- bios2070[1:24, 1:40, c(1,5:9,12,16:19)]
#   nlayers <- dim(bios2070)[3]  # nlayers mora biti lat_dim x lon_dim x BIO1-BIO19
# 
#   raster_stack2070 <- stack(
#     lapply(1:nlayers, function(i) {
#       layer <- matrix(bios2070[, , i], nrow = n_rows, ncol = n_cols)
#       layer <- flipdim(layer, 1)
#       raster_layer <- raster(layer)  # Define raster layer inside the loop
#       extent_of_raster <- extent(min(lon), max(lon), min(lat), max(lat))
#       # Assign the extent to the raster layer
#       extent(raster_layer) <- extent_of_raster
#       # Define the projection of the raster layer (WGS84 in this example)
#       projection(raster_layer) <- CRS("+proj=longlat +datum=WGS84")
#       return(raster_layer)
#     })
#   )
#   
#   
#   futureEnvSI <- raster_stack2070
#   layer_names <- c("bio1", "bio5", "bio6", "bio7", "bio8", "bio9", "bio12", "bio16", "bio17", "bio18", "bio19")
#   names(futureEnvSI) <- layer_names
#   modelFutureEnvSI=crop(futureEnvSI, model.extent)
# 
#   rattler.2070SI <- predict(rattler.me, modelFutureEnvSI)
#   rattler.2070SI_gg <- data.frame(rasterToPoints(rattler.2070SI))
#   xi <- seq(min(rattler.2070SI_gg$x),max(rattler.2070SI_gg$x),length.out = 10)
#   yi <- seq(min(rattler.2070SI_gg$y),max(rattler.2070SI_gg$y),length.out = 6)
#   lon <- round(seq(min(lon),max(lon),length.out = 10), digits = 2)
#   lat <- round(seq(min(lat),max(lat),length.out = 6), digits = 2)
#   
#   ggplot(data=rattler.2070SI_gg) +
#     geom_raster(aes(x = x, y = y, fill = layer)) + geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
#     # geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
#     scale_fill_distiller("Probability", palette = "Spectral", limits = c(0,1)) +
#     labs(x="Longitude",y="Latitude",title=paste0("Predicted suitability for ",zuzelka,", 2071-2100, model ",model[z]," and ",rcp)) +
#     theme_light() 
#   ggsave(paste0("proj/2071-2100/Predicted_suitability_2071-2100_",model[z],"_",rcp,".png"),width = 10, height = 7)
# # }


######################################### PROJEKCIJE MEDIANA MODELOV 2071 - 2100

model <- c("CNRM-CERFACS-CNRM-CM5", 
           "ICHEC-EC-EARTH",
           "IPSL-IPSL-CM5A-MR",
           "MOHC-HadGEM2-ES", # samo do 30. 11. 2099 pri padavinah, zato tudi pri Tmin in Tmax .nc fajl spremenim, da gre le do 30.11.
           "MPI-M-MPI-ESM-LR",
           "MPI-M-MPI-ESM-LR")
drugo <- c("_r1i1p1_CLMcom-CCLM4-8-17_v1",
           "_r3i1p1_DMI-HIRHAM5_v1",
           "_r1i1p1_IPSL-INERIS-WRF331F_v1",
           "_r1i1p1_KNMI-RACMO22E_v2",
           "_r1i1p1_CLMcom-CCLM4-8-17_v1",
           "_r1i1p1_SMHI-RCA4_v1a")
leto_kon <- c("21001231","21001231","21001231","20991130","21001231","21001231")

risi0 <- function(rcp, leto_zac, leto_zac1, obd, zuzelka){

  nc_data_rr<- nc_open(paste0(path,"pr","_12km_MOHC-HadGEM2-ES_",rcp,drugo[4],"_day_",leto_zac,"0101_",leto_kon[4],".nc"))
  time <- ncvar_get(nc_data_rr, "time")
  start1<-c(1,1,1)
  count1<-c(nx,ny, length(time))
  data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
  # povpre?na dnevna T
  datumi0 <- seq(0,count1[3]-1,by = 1)
  pretvorba_sek_v_dan1 = 60*60*24 #?eprav imam dnevne podatke moram to dat, ker as.positxct meri v sekundah
  datumi <- as.POSIXct((start1[3]+datumi0)*pretvorba_sek_v_dan1,origin=paste0(leto_zac1,"-12-31 00:00:00"))
  leta <- format(datumi,format="%Y")
  datum2 = format(as.POSIXct(datumi, format="%Y-%m-%d"), "%m/%d/%Y")
  nc_close(nc_data_rr) # konec branja
  
  rattler.vsi <- data.frame()
  
  for(z in seq(1,6,by=1)){
    model <- c("CNRM-CERFACS-CNRM-CM5", 
               "ICHEC-EC-EARTH",
               "IPSL-IPSL-CM5A-MR",
               "MOHC-HadGEM2-ES", # samo do 30. 11. 2099 pri padavinah, zato tudi pri Tmin in Tmax .nc fajl spremenim, da gre le do 30.11.
               "MPI-M-MPI-ESM-LR",
               "MPI-M-MPI-ESM-LR")
    drugo <- c("_r1i1p1_CLMcom-CCLM4-8-17_v1",
               "_r3i1p1_DMI-HIRHAM5_v1",
               "_r1i1p1_IPSL-INERIS-WRF331F_v1",
               "_r1i1p1_KNMI-RACMO22E_v2",
               "_r1i1p1_CLMcom-CCLM4-8-17_v1",
               "_r1i1p1_SMHI-RCA4_v1a")
    leto_zac <- "2071"
    leto_kon <- c("21001231","21001231","21001231","20991130","21001231","21001231")
    leto_zac1 <- "2070"
    
    # padavine dnevne
    nc_data_rr<- nc_open(paste0(path,"pr","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
    data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
    nc_close(nc_data_rr) # konec branja
    # maksimalna dnevna T
    nc_data_max<- nc_open(paste0(path,"tasmax","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
    data_max <- ncvar_get(nc_data_max, "tasmax", start=start1, count=count1)
    nc_close(nc_data_max) # konec branja
    # minimalna dnevna T
    nc_data_min <- nc_open(paste0(path,"tasmin","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
    data_min <- ncvar_get(nc_data_min, "tasmin", start=start1, count=count1)
    nc_close(nc_data_min) # konec branja
    
    leto = format(datumi, format = "%Y")
    leta1 <- seq(min(leto),max(leto),by=1)
    lats1 <- seq(1,24,by=1)
    lons1 <- seq(1,40,by=1)
    n_rows <- 24
    n_cols <- 40
    n_layers <- 19
    bios2070 <- array(0, dim = c(n_rows, n_cols, n_layers))
  
    for(j in lats1){
      for(i in lons1){
        print(i)
        print(j) # podatki za T so v worldclim mno?eni z 10, zato da nimajo decimalnih mest, zato tudi jaz svoje mno?im, da bom lahko uporabila model
        data_rr1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), rr = data_rr[i,j,1:length(time)]*24*60*60) # ARSO podatki so v mm/s zato za dnevne pomno?im s ?tevilom s v dnevu
        data_min1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tn = 10*(data_min[i,j,1:length(time)]-273)) # ARSO podatki so v Kelvinih
        data_max1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tx = (data_max[i,j,1:length(time)]-273)*10)
        
        # mese?ne vsote padavin in mese?na povpre?ja temperature
        dataRR <- aggregate(rr ~ mes + leto, data_rr1 , FUN = sum , na.rm=TRUE, na.action=na.pass)
        dataMIN <- aggregate(tn ~ mes + leto, data_min1 , mean , na.rm=TRUE, na.action=na.pass)
        dataMAX <- aggregate(tx ~ mes + leto, data_max1, mean , na.rm=TRUE, na.action=na.pass)
        
        # povpre?je let obdobja 2071-2100
        dataRR <- aggregate(rr ~ leto, dataRR , mean , na.rm=TRUE, na.action=na.pass)
        dataMIN <- aggregate(tn ~ leto, dataMIN , mean , na.rm=TRUE, na.action=na.pass)
        dataMAX <- aggregate(tx ~ leto, dataMAX, mean , na.rm=TRUE, na.action=na.pass)
        
        bio_histSI <- biovars(prec = dataRR$rr, tmin = dataMIN$tn, tmax = dataMAX$tx)
        bios2070[j,i,1:19] = bio_histSI
      }}
    
    bios2070 <- bios2070[1:24, 1:40, c(1,5:9,12,16:19)]
    nlayers <- dim(bios2070)[3]  # nlayers mora biti lat_dim x lon_dim x BIO1-BIO19

    raster_stack2070 <- stack(
      lapply(1:nlayers, function(i) {
        layer <- matrix(bios2070[, , i], nrow = n_rows, ncol = n_cols)
        layer <- flipdim(layer, 1)
        raster_layer <- raster(layer)  # Define raster layer inside the loop
        extent_of_raster <- extent(min(lon), max(lon), min(lat), max(lat))
        # Assign the extent to the raster layer
        extent(raster_layer) <- extent_of_raster
        # Define the projection of the raster layer (WGS84 in this example)
        projection(raster_layer) <- CRS("+proj=longlat +datum=WGS84")
        return(raster_layer)
      })
    )
    
    
    futureEnvSI <- raster_stack2070
    layer_names <- c("bio1", "bio5", "bio6", "bio7", "bio8", "bio9", "bio12", "bio16", "bio17", "bio18", "bio19")
    names(futureEnvSI) <- layer_names
    modelFutureEnvSI=crop(futureEnvSI, extent_of_raster)
  
    rattler.me <- readRDS(file = paste0("rattler_me",zuzelka,".rds"))
    rattler.2070SI <- predict(rattler.me, modelFutureEnvSI)
  
    rattler.2070SI_gg <- data.frame(rasterToPoints(rattler.2070SI))
    rattler.2070SI_gg$model_ime <- rep(z,length(rattler.2070SI_gg$x))
    
    rattler.vsi <- rbind(rattler.vsi, rattler.2070SI_gg)
  }
  
  
  saveRDS(rattler.vsi, file = paste0("rattler_",zuzelka,"_",rcp,"_",obd,".rds"))
  return(rattler.vsi)
}


rattler.vsi0 <- risi0(rcp = "rcp45", leto_zac = "2071", leto_zac1 = "2070", obd = "2071-2100", zuzelka = zuzelka1)
rattler.vsi1 <- risi0(rcp = "rcp85", leto_zac = "2071", leto_zac1 = "2070", obd = "2071-2100", zuzelka = zuzelka1)


grafi <- function(rattler.vsi, rcp, obd, zuzelka){
  rattler.mediana <- aggregate(layer ~  x + y, rattler.vsi, FUN = median , na.rm=TRUE, na.action=na.pass)
  rattler.min <- aggregate(layer ~  x + y, rattler.vsi, FUN = "min" , na.rm=TRUE, na.action=na.pass)
  rattler.max <- aggregate(layer ~  x + y, rattler.vsi, FUN = "max" , na.rm=TRUE, na.action=na.pass)
  
  koruza <- read_sf("koruza/GRID_Koruza_ha.shp")
  koruza <- st_transform(koruza, crs = 4326)
  koruza$mid <- st_centroid(koruza$geometry) 
  koruza$Area <- koruza$Koruza_ha
  labels <- c("Low", "Moderate", "High", "Very High")
  # labels <- c("Unsuitable", "High") 
  
  ggplot(data=rattler.mediana) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    # geom_sf(data = koruza, aes(geometry = mid, size = Area), color = "black", show.legend = "point") +
    scale_fill_gradientn("Suitability",  
                         colours = rev(brewer.pal(6, "Spectral")),
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),                         
                         labels = c(0, 0.25, 0.5, 0.75, 1)) +
    geom_sf(data = slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
    labs(x="Longitude",y="Latitude")+#,title=paste0("Predicted suitability for ",zuzelka,", ",obd,", Median and ",rcp)) +
    coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
    theme_light() + theme(legend.position = c(.9, .3), legend.title = element_text(face = "bold", size=16),
                          text = element_text(size=16))#,
                          # plot.title = element_text(size=16))
# 
#   ggplot(data=rattler.mediana) +
#     geom_raster(aes(x = x, y = y, fill = layer)) +
#     scale_fill_gradientn("Suitability", 
#                          colours = rev(brewer.pal(6, "Spectral")),
#                          limits = c(0, 1),
#                          breaks = c(0, 1),
#                          labels = labels) +
#     geom_sf(data = slovenia_nuts3_mapdata, color=alpha("black",0.4), fill = NA) +
#     labs(x="Longitude", y="Latitude", title=paste0("Predicted suitability for ", zuzelka, ", ", obd, ", Median and ", rcp)) +
#     coord_sf(crs = st_crs(4326)) +
#     xlim(c(start_lon, end_lon)) + ylim(c(start_lat, end_lat)) +
#     theme_light() + 
#     theme(
#       legend.position = c(.9, .3), 
#       legend.title = element_text(face = "bold", size=16),
#       text = element_text(size=16),
#       plot.title = element_text(size=16)
#     )
  ggsave(paste0("proj/",obd,"/Predicted_suitability_",zuzelka,"_",obd,"_median_",rcp,".png"),width = 10, height = 7)
  
  ggplot(data=rattler.min) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    scale_fill_gradientn("Suitability",  
                         colours = rev(brewer.pal(6, "Spectral")),
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),                         
                         labels = c(0, 0.25, 0.5, 0.75, 1)) +
    geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
    # geom_sf(data = koruza, aes(geometry = mid, size = Area), color = "black", show.legend = "point") +
    labs(x="Longitude",y="Latitude")+#,title=paste0("Predicted suitability for ",zuzelka,", ",obd,", Min and ",rcp)) +
    coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
    theme_light() + theme(legend.position = c(.9, .3), legend.title = element_text(face = "bold", size=16),
                          text = element_text(size=16))#,
                          # plot.title = element_text(size=16))
  ggsave(paste0("proj/",obd,"/Predicted_suitability_",zuzelka,"_",obd,"_min_",rcp,".png"),width = 10, height = 7)
  
  ggplot(data=rattler.max) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    scale_fill_gradientn("Suitability",  
                         colours = rev(brewer.pal(6, "Spectral")),
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),                         
                         labels = c(0, 0.25, 0.5, 0.75, 1)) +
    # geom_sf(data = koruza, aes(geometry = mid, size = Area), color = "black", show.legend = "point") +
    geom_sf(data = slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
    labs(x="Longitude",y="Latitude")+#,title=paste0("Predicted suitability for ",zuzelka,", ",obd,", Max and ",rcp)) +
    coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
    theme_light() + theme(legend.position = c(.9, .3), legend.title = element_text(face = "bold", size=16),
                          text = element_text(size=16))+#,
                          # plot.title = element_text(size=16))
  ggsave(paste0("proj/",obd,"/Predicted_suitability_",zuzelka,"_",obd,"_max_",rcp,".png"),width = 10, height = 7)
}


rattler.vsi0 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp45_2071-2100.rds"))
rattler.vsi1 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp85_2071-2100.rds"))

grafi(rattler.vsi = rattler.vsi0, rcp = "rcp45", obd = "2071-2100", zuzelka = zuzelka1)
grafi(rattler.vsi = rattler.vsi1, rcp = "rcp85", obd = "2071-2100", zuzelka = zuzelka1)

######################################### PROJEKCIJE MEDIANA MODELOV 2041 - 2070

model <- c("CNRM-CERFACS-CNRM-CM5", 
           "ICHEC-EC-EARTH",
           "IPSL-IPSL-CM5A-MR",
           "MOHC-HadGEM2-ES", 
           "MPI-M-MPI-ESM-LR",
           "MPI-M-MPI-ESM-LR")
drugo <- c("_r1i1p1_CLMcom-CCLM4-8-17_v1",
           "_r3i1p1_DMI-HIRHAM5_v1",
           "_r1i1p1_IPSL-INERIS-WRF331F_v1",
           "_r1i1p1_KNMI-RACMO22E_v2",
           "_r1i1p1_CLMcom-CCLM4-8-17_v1",
           "_r1i1p1_SMHI-RCA4_v1a")

risi <- function(rcp, leto_zac, leto_kon, leto_zac1, obd, zuzelka){
  nc_data_rr<- nc_open(paste0(path,"pr","_12km_MOHC-HadGEM2-ES_",rcp,drugo[4],"_day_",leto_zac,"0101_",leto_kon,".nc"))
  time <- ncvar_get(nc_data_rr, "time")
  start1<-c(1,1,1)
  count1<-c(nx,ny, length(time))
  data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
  # povpre?na dnevna T
  datumi0 <- seq(0,count1[3]-1,by = 1)
  pretvorba_sek_v_dan1 = 60*60*24 #?eprav imam dnevne podatke moram to dat, ker as.positxct meri v sekundah
  datumi <- as.POSIXct((start1[3]+datumi0)*pretvorba_sek_v_dan1,origin=paste0(leto_zac1,"-12-31 00:00:00"))
  leta <- format(datumi,format="%Y")
  datum2 = format(as.POSIXct(datumi, format="%Y-%m-%d"), "%m/%d/%Y")
  nc_close(nc_data_rr) # konec branja
  rattler.vsi <- data.frame()
  
  for(z in seq(1,6,by=1)){
    # padavine dnevne
    nc_data_rr<- nc_open(paste0(path,"pr","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon,".nc"))
    data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
    nc_close(nc_data_rr) # konec branja
    # maksimalna dnevna T
    nc_data_max<- nc_open(paste0(path,"tasmax","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon,".nc"))
    data_max <- ncvar_get(nc_data_max, "tasmax", start=start1, count=count1)
    nc_close(nc_data_max) # konec branja
    # minimalna dnevna T
    nc_data_min <- nc_open(paste0(path,"tasmin","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon,".nc"))
    data_min <- ncvar_get(nc_data_min, "tasmin", start=start1, count=count1)
    nc_close(nc_data_min) # konec branja
    
    leto = format(datumi, format = "%Y")
    leta1 <- seq(min(leto),max(leto),by=1)
    lats1 <- seq(1,24,by=1)
    lons1 <- seq(1,40,by=1)
    n_rows <- 24
    n_cols <- 40
    n_layers <- 19
    bios2070 <- array(0, dim = c(n_rows, n_cols, n_layers))
    
    for(j in lats1){
      for(i in lons1){
        print(i)
        print(j) # podatki za T so v worldclim mno?eni z 10, zato da nimajo decimalnih mest, zato tudi jaz svoje mno?im, da bom lahko uporabila model
        data_rr1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), rr = data_rr[i,j,1:length(time)]*24*60*60) # ARSO podatki so v mm/s zato za dnevne pomno?im s ?tevilom s v dnevu
        data_min1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tn = 10*(data_min[i,j,1:length(time)]-273)) # ARSO podatki so v Kelvinih
        data_max1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tx = (data_max[i,j,1:length(time)]-273)*10)
        
        dataRR <- aggregate(rr ~ mes + leto, data_rr1 , FUN = sum , na.rm=TRUE, na.action=na.pass)
        dataMIN <- aggregate(tn ~ mes + leto, data_min1 , mean , na.rm=TRUE, na.action=na.pass)
        dataMAX <- aggregate(tx ~ mes + leto, data_max1, mean , na.rm=TRUE, na.action=na.pass)
        
        dataRR <- aggregate(rr ~ leto, dataRR , mean , na.rm=TRUE, na.action=na.pass)
        dataMIN <- aggregate(tn ~ leto, dataMIN , mean , na.rm=TRUE, na.action=na.pass)
        dataMAX <- aggregate(tx ~ leto, dataMAX, mean , na.rm=TRUE, na.action=na.pass)
        
        bio_histSI <- biovars(prec = dataRR$rr, tmin = dataMIN$tn, tmax = dataMAX$tx)
        bios2070[j,i,1:19] = bio_histSI
      }}
    
    bios2070 <- bios2070[1:24, 1:40, c(1,5:9,12,16:19)]
    nlayers <- dim(bios2070)[3]  # nlayers mora biti lat_dim x lon_dim x BIO1-BIO19
    raster_stack2070 <- stack(
      lapply(1:nlayers, function(i) {
        layer <- matrix(bios2070[, , i], nrow = n_rows, ncol = n_cols)
        layer <- flipdim(layer, 1)
        raster_layer <- raster(layer)  # Define raster layer inside the loop
        extent_of_raster <- extent(min(lon), max(lon), min(lat), max(lat))
        # Assign the extent to the raster layer
        extent(raster_layer) <- extent_of_raster
        # Define the projection of the raster layer (WGS84 in this example)
        projection(raster_layer) <- CRS("+proj=longlat +datum=WGS84")
        return(raster_layer)
      })
    )
    
    futureEnvSI <- raster_stack2070
    layer_names <- c("bio1", "bio5", "bio6", "bio7", "bio8", "bio9", "bio12", "bio16", "bio17", "bio18", "bio19")
    names(futureEnvSI) <- layer_names
    modelFutureEnvSI=crop(futureEnvSI, extent_of_raster)
    rattler.me <- readRDS(file = paste0("rattler_me",zuzelka,".rds"))
    rattler.2070SI <- predict(rattler.me, modelFutureEnvSI)
    rattler.2070SI_gg <- data.frame(rasterToPoints(rattler.2070SI))
    rattler.2070SI_gg$model_ime <- rep(z,length(rattler.2070SI_gg$x))
    rattler.vsi <- rbind(rattler.vsi, rattler.2070SI_gg)
  }
  
  saveRDS(rattler.vsi, file = paste0("rattler_",zuzelka,"_",rcp,"_",obd,".rds"))
  return(rattler.vsi)
}


# rattler.vsi2 <- risi(rcp = "rcp45", leto_zac = "2041", leto_kon = "20701231", leto_zac1 = "2040", obd = "2041-2070", zuzelka = zuzelka1)
# rattler.vsi3 <- risi(rcp = "rcp85", leto_zac = "2041", leto_kon = "20701231", leto_zac1 = "2040", obd = "2041-2070", zuzelka = zuzelka1)
rattler.vsi2 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp45_2041-2070.rds"))
rattler.vsi3 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp85_2041-2070.rds"))

# rattler.vsi4 <- risi(rcp = "rcp45", leto_zac = "2011", leto_kon = "20401231", leto_zac1 = "2010", obd = "2011-2040", zuzelka = zuzelka1)
# rattler.vsi5 <- risi(rcp = "rcp85", leto_zac = "2011", leto_kon = "20401231", leto_zac1 = "2010", obd = "2011-2040", zuzelka = zuzelka1)
rattler.vsi4 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp45_2011-2040.rds"))
rattler.vsi5 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp85_2011-2040.rds"))

grafiK <- function(rattler.vsi, rcp, obd, zuzelka){
  rattler.mediana <- aggregate(layer ~  x + y, rattler.vsi, FUN = median , na.rm=TRUE, na.action=na.pass)
  rattler.min <- aggregate(layer ~  x + y, rattler.vsi, FUN = "min" , na.rm=TRUE, na.action=na.pass)
  rattler.max <- aggregate(layer ~  x + y, rattler.vsi, FUN = "max" , na.rm=TRUE, na.action=na.pass)
  
  koruza <- read_sf("koruza/GRID_Koruza_ha.shp")
  koruza <- st_transform(koruza, crs = 4326)
  koruza$mid <- st_centroid(koruza$geometry) 
  koruza$Koruza_ha[koruza$Koruza_ha == -9] <- NA
  koruza$Area <- koruza$Koruza_ha
  labels <- c("Low", "Moderate", "High", "Very high") 
  
  ggplot(data=rattler.mediana) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    geom_sf(data = koruza, aes(geometry = mid, size = Area), color = "black", show.legend = "point") +
    scale_fill_gradientn("Suitability",  
                         colours = rev(brewer.pal(6, "Spectral")),
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),                         
                         labels = c(0, 0.25, 0.5, 0.75, 1)) +
    geom_sf(data = slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
    labs(x="Longitude",y="Latitude")+#,title=paste0("Predicted suitability for ",zuzelka,", ",obd,", Median and ",rcp)) +
    coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
    theme_light() + theme(legend.position = c(.9, .3), legend.title = element_text(face = "bold", size=16),
                          text = element_text(size=16))#,
                          # plot.title = element_text(size=16))
  
  
  ggsave(paste0("proj/",obd,"/Predicted_suitability_",zuzelka,"_",obd,"_median_",rcp,".png"),width = 10, height = 7)
  
  ggplot(data=rattler.min) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    scale_fill_gradientn("Suitability",  
                         colours = rev(brewer.pal(6, "Spectral")),
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),                         
                         labels = c(0, 0.25, 0.5, 0.75, 1)) +
    geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
    geom_sf(data = koruza, aes(geometry = mid, size = Area), color = "black", show.legend = "point") +
    labs(x="Longitude",y="Latitude")+#,title=paste0("Predicted suitability for ",zuzelka,", ",obd,", Min and ",rcp)) +
    coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
    theme_light() + theme(legend.position = c(.9, .3), legend.title = element_text(face = "bold", size=16),
                          text = element_text(size=16))#,
#                          plot.title = element_text(size=16))
  ggsave(paste0("proj/",obd,"/Predicted_suitability_",zuzelka,"_",obd,"_min_",rcp,".png"),width = 10, height = 7)
  
  ggplot(data=rattler.max) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    scale_fill_gradientn("Suitability",  
                         colours = rev(brewer.pal(6, "Spectral")),
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),                         
                         labels = c(0, 0.25, 0.5, 0.75, 1)) +
    geom_sf(data = koruza, aes(geometry = mid, size = Area), color = "black", show.legend = "point") +
    geom_sf(data = slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
    labs(x="Longitude",y="Latitude")+#,title=paste0("Predicted suitability for ",zuzelka,", ",obd,", Max and ",rcp)) +
    coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
    theme_light() + theme(legend.position = c(.9, .3), legend.title = element_text(face = "bold", size=16),
                          text = element_text(size=16))#,
#                          plot.title = element_text(size=16))
  ggsave(paste0("proj/",obd,"/Predicted_suitability_",zuzelka,"_",obd,"_max_",rcp,".png"),width = 10, height = 7)
}

grafiK(rattler.vsi = rattler.vsi2, rcp = "rcp45", obd = "2041-2070", zuzelka = zuzelka1)
grafiK(rattler.vsi = rattler.vsi3, rcp = "rcp85", obd = "2041-2070", zuzelka = zuzelka1)

grafiK(rattler.vsi = rattler.vsi4, rcp = "rcp45", obd = "2011-2040", zuzelka = zuzelka1)
grafiK(rattler.vsi = rattler.vsi5, rcp = "rcp85", obd = "2011-2040", zuzelka = zuzelka1)


##################################################################################################################
#
#                     DRUGE ZUZELKE Trichogramma pretiosum, 
#         Cotesia marginiventris, Telenomus remus in Eiphosoma laphygmae
#
##################################################################################################################
start_lon <- 13.2
end_lon <- 16.8

start_lat<-45.3 #+ 25/60 +18.34/3600
end_lat<-47

# zuzelka1 = 'Spodoptera frugiperda'
# zuzelka1 = 'Trichogramma pretiosum'
# zuzelka1 = 'Cotesia marginiventris'
# zuzelka1 = 'Telenomus remus'
zuzelka1 = 'Eiphosoma laphygmae'
#
# r1=gbif(zuzelka1)
# r1<-subset(r1, !is.na(lon) & !is.na(lat))
# rattlerdups=duplicated(r1[, c("lon", "lat")])
# r1 <-r1[!rattlerdups, ]
# 
# # data(wrld_simpl)
# # plot(wrld_simpl, xlim=c(min(r1$lon)-1,max(r1$lon)+1), ylim=c(min(r1$lat)-1,max(r1$lat)+1), axes=TRUE, col="light yellow")
# # points(r1$lon, r1$lat, col="lightblue", pch=20, cex=0.75)
# # save.image(paste0("occurrence_",zuzelka1,"_GBIF.png"))
# currentEnv=getData("worldclim", var="bio", res=2.5)
# currentEnv=dropLayer(currentEnv, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio13", "bio14", "bio15"))
# model.extent<-extent(min(r1$lon)-10,max(r1$lon)+10,min(r1$lat)-10,max(r1$lat)+10)
# modelEnv=crop(currentEnv,model.extent)
# 
# rattlerocc=cbind.data.frame(r1$lon,r1$lat) # a data frame of latitudes and longitudes for the model
# # fold <- kfold(rattlerocc, k=5) # add an index that makes five random groups of observations
# # rattlertest <- rattlerocc[fold == 1, ]
# # rattlertrain <- rattlerocc[fold != 1, ]
# # rattler.me1 <- maxent(modelEnv, rattlertrain)
# # saveRDS(rattler.me1, file = paste0("rattler_me",zuzelka1,".rds"))
# # plot(rattler.me1)
# # response(rattler.me1)
# 
# rattler.me1FULL <- maxent(modelEnv, rattlerocc)
# saveRDS(rattler.me1FULL, file = paste0("rattler_me",zuzelka1,".rds"))
# plot(rattler.me1FULL)
# response(rattler.me1FULL)
# # rattler.pred1 <- predict(rattler.me1, modelEnv)
# # saveRDS(rattler.pred1, file = paste0("rattler_pred",zuzelka1,".rds"))
# bg <- randomPoints(modelEnv, 1000)
# e1 <- evaluate(rattler.me1FULL, p=rattlerocc, a=bg, x=modelEnv)
# plot(e1, 'ROC') # precej stran od 0.5, mnogo bolje napove kot random

######## HISTORICNI DRUGE ZUZELKE

rattler.me1 <- readRDS(file = paste0("rattler_me",zuzelka1,".rds"))

e1 <- evaluate(rattler.me1, p=rattlerocc, a=bg, x=modelEnv)
plot(e1, 'ROC')
e_data <- data.frame(TPR = e1@TPR, FPR = e1@FPR)
ee<-ggplot(e_data, aes(x=FPR, y=TPR)) +
  geom_point(size=2, colour="darkred")+
  labs(title = paste("AUC = ",round(e1@auc,digits = 3),sep="")) +  ylab("dele? resni?nih pozitivnih primerov") + xlab("dele? la?nih pozitivnih primerov") + 
  geom_abline(intercept = 0, slope = 1) +
  theme_light(base_size = 18)
ggsave(paste("hist/ROC_plot_",zuzelka1,".png",sep=""),ee, dpi=500, height=6, width=9, units="in")
print(ee)

plot(rattler.me1)
response_data <- response(rattler.me1)
library(tidyverse)
var_contr<-rattler.me1@results[7:17,1]
var_imena<-c("bio1","bio12","bio16","bio17","bio18","bio19","bio5","bio6","bio7","bio8","bio9")
var_c <- data.frame(var_contr,var_imena)
var_c1<- var_c %>%
  arrange(desc(var_contr)) 
aa<-ggplot(var_c1, aes(x=reorder(var_imena, -var_contr), y=var_contr)) +
  geom_point(size=2, colour="darkred")+
  labs(title = paste("Prispevki spremenljivk - ",zuzelka1,sep="")) +  ylab("odstotni dele? [%]") + xlab("ime spremenljivke") + 
  theme_light(base_size = 18)
ggsave(paste("hist/variable_contribution",zuzelka1,".png",sep=""),aa, dpi=500, height=6, width=9, units="in")
print(aa)

# First, choose the environmental variable you want to plot (e.g., bio1)
bio1_range <- seq(min(modelEnvSI$bio1[], na.rm=TRUE), max(modelEnvSI$bio1[], na.rm=TRUE), length.out=100)
# Create a data frame holding the range of bio1 and constant values for the other variables
# You can use the mean of the other variables for simplicity
predictors_df <- data.frame(
  bio1 = bio1_range,
  bio12 = seq(min(modelEnvSI$bio12[], na.rm=TRUE), max(modelEnvSI$bio12[], na.rm=TRUE), length.out=100),  # Adjust for your other variables
  bio17 = seq(min(modelEnvSI$bio17[], na.rm=TRUE), max(modelEnvSI$bio17[], na.rm=TRUE), length.out=100),  # Adjust for your other variables
  bio18 = seq(min(modelEnvSI$bio18[], na.rm=TRUE), max(modelEnvSI$bio18[], na.rm=TRUE), length.out=100)
  # Add other necessary variables depending on your model
)

# Make predictions using the MaxEnt model
response_predictions <- predict(rattler.me1, predictors_df)


nc_data <- nc_open(paste0(path,"evspsblpot","_12km_ARSO_v5_day_19810101_20101231.nc"))
lon <- ncvar_get(nc_data, "X") # drugi parameter je ime spremenljivke v datoteki
lat <- ncvar_get(nc_data, "Y")
time <- ncvar_get(nc_data, "time")
nx = length(lon)
ny = length(lat)
start1<-c(1,1,1)
count1<-c(nx,ny, 10957)
datumi0 <- seq(0,count1[3]-1,by = 1)
pretvorba_sek_v_dan1 = 60*60*24 #?eprav imam dnevne podatke moram to dat, ker as.positxct meri v sekundah
datumi <- as.POSIXct((start1[3]+datumi0)*pretvorba_sek_v_dan1,origin="1980-12-31 00:00:00")
leta <- format(datumi,format="%Y")
datum2 = format(as.POSIXct(datumi, format="%Y-%m-%d"), "%m/%d/%Y")
nc_close(nc_data) # konec branja

# minimalna dnevna T
nc_data_min <- nc_open(paste0(path,"tasmin","_12km_ARSO_v5_day_19810101_20101231.nc"))
data_min <- ncvar_get(nc_data_min, "tasmin", start=start1, count=count1)
nc_close(nc_data_min) # konec branja
# maksimalna dnevna T
nc_data_max<- nc_open(paste0(path,"tasmax","_12km_ARSO_v5_day_19810101_20101231.nc"))
data_max <- ncvar_get(nc_data_max, "tasmax", start=start1, count=count1)
nc_close(nc_data_max) # konec branja
# padavine dnevne
nc_data_rr<- nc_open(paste0(path,"pr","_12km_ARSO_v5_day_19810101_20101231.nc"))
data_rr <- na.omit(ncvar_get(nc_data_rr, "pr", start=start1, count=count1))
nc_close(nc_data_rr) # konec branja
lats1 <- seq(1,24,by=1)
lons1 <- seq(1,40,by=1)

n_rows <- 24
n_cols <- 40
n_layers <- 19
bios <- array(0, dim = c(n_rows, n_cols, n_layers))
leto = format(datumi, format = "%Y")
leta <- seq(min(leto),max(leto),by=1)

for(j in lats1){
  for(i in lons1){
    print(i)
    print(j)
    data_rr1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), rr = data_rr[i,j,1:length(time)]*24*60*60) # ARSO podatki so v mm/s zato za dnevne pomno?im s ?tevilom s v dnevu
    data_min1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tn = 10*(data_min[i,j,1:length(time)]-273)) # ARSO podatki so v Kelvinih
    data_max1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tx = (data_max[i,j,1:length(time)]-273)*10)
    # mese?ne vsote padavin in mese?na povpre?ja temperature
    dataRR <- aggregate(rr ~ mes + leto, data_rr1 , FUN = sum , na.rm=TRUE, na.action=na.pass)
    dataMIN <- aggregate(tn ~ mes + leto, data_min1 , FUN = mean , na.rm=TRUE, na.action=na.pass)
    dataMAX <- aggregate(tx ~ mes + leto, data_max1, FUN = mean , na.rm=TRUE, na.action=na.pass)
    # povpre?je let obdobja 1981-2010
    dataRR <- aggregate(rr ~ leto, dataRR , FUN = mean , na.rm=TRUE, na.action=na.pass)
    dataMIN <- aggregate(tn ~ leto, dataMIN , FUN = mean , na.rm=TRUE, na.action=na.pass)
    dataMAX <- aggregate(tx ~ leto, dataMAX, FUN = mean , na.rm=TRUE, na.action=na.pass)
    
    bio_histSI <- biovars(prec = dataRR$rr, tmin = dataMIN$tn, tmax = dataMAX$tx)
    bios[j,i,1:19] = bio_histSI
  }}

bios <- bios[1:24, 1:40, c(1,5:9,12,16:19)]
nlayers <- dim(bios)[3]  # nlayers mora biti lat_dim x lon_dim x BIO1-BIO19

raster_stack <- stack()

# Loop through the layers and add them to the RasterStack
for (i in 1:nlayers){
  layer <- matrix(bios[, , i], nrow = n_rows, ncol = n_cols)
  raster_layer <- raster(layer)
  extent_of_raster <- extent(min(lon), max(lon), min(lat), max(lat))
  # Assign the extent to the raster layer
  extent(raster_layer) <- extent_of_raster
  # Define the projection of the raster layer (WGS84 in this example)
  projection(raster_layer) <- CRS("+proj=longlat +datum=WGS84")
  raster_stack <- addLayer(raster_stack, raster_layer)
}

currentEnvSI <- raster_stack
layer_names <- c("bio1", "bio5", "bio6", "bio7", "bio8", "bio9", "bio12", "bio16", "bio17", "bio18", "bio19")
names(currentEnvSI) <- layer_names
modelEnvSI=crop(currentEnvSI,extent_of_raster)


rattler.predSI <- predict(rattler.me1, modelEnvSI)
rattler.pred_gg <- data.frame(rasterToPoints(rattler.predSI))
xi <- seq(min(rattler.pred_gg$x),max(rattler.pred_gg$x),length.out = 10)
yi <- seq(min(rattler.pred_gg$y),max(rattler.pred_gg$y),length.out = 6)
lon <- round(seq(min(lon),max(lon),length.out = 10), digits = 2)
lat <- round(seq(min(lat),max(lat),length.out = 6), digits = 2)


labels <- c("Low", "Moderate", "High") 

ggplot(data=rattler.pred_gg) +
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
  scale_fill_gradientn("Suitability",  
                       colours = rev(brewer.pal(6, "Spectral")),
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),                         
                       labels = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x="Longitude",y="Latitude")+#,title=paste0("Predicted suitability for ",zuzelka1,", 1981-2010")) +
  coord_sf(crs = st_crs(4326)) + #xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
  theme_light() + theme(legend.position = c(.9, .2), 
                       legend.title = element_text(face = "bold"),
                       text = element_text(size=18))
ggsave(paste0("hist/Predicted_suitability_1981-2010",zuzelka1,".png"),width = 10, height = 7)

bb<-ggplot(data=rattler.pred_gg) +
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
  scale_fill_gradientn("Ustreznost", 
                       colours = rev(brewer.pal(6, "Spectral")),
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),                         
                       labels = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x="G. dol?ina",y="G. ?irina",title=paste0("Verjetnost pojava ?u?elke ",zuzelka1,", 1981-2010")) +
  coord_sf(crs = st_crs(4326)) + #xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
  theme_light() + theme(legend.position = c(.9, .2), 
                        legend.title = element_text(face = "bold"),
                        text = element_text(size=18))
print(bb)
ggsave(paste0("hist/SLO_Predicted_suitability_1981-2010",zuzelka1,".png"),bb,width = 10, height = 7)

######### PROJEKCIJE MEDIANA MODELOV DRUGE ZUZELKE

zuzelka1 = 'Spodoptera frugiperda'
# zuzelka1 = 'Trichogramma pretiosum'
# zuzelka1 = 'Cotesia marginiventris'
# zuzelka1 = 'Telenomus remus'
# zuzelka1 = 'Eiphosoma laphygmae'
# 
rattler.me1 <- readRDS(file = paste0("rattler_me",zuzelka1,".rds"))
# rattler.pred1 <- readRDS(file = paste0("rattler_pred",zuzelka1,".rds"))

model <- c("CNRM-CERFACS-CNRM-CM5", 
           "ICHEC-EC-EARTH",
           "IPSL-IPSL-CM5A-MR",
           "MOHC-HadGEM2-ES", 
           "MPI-M-MPI-ESM-LR",
           "MPI-M-MPI-ESM-LR")
drugo <- c("_r1i1p1_CLMcom-CCLM4-8-17_v1",
           "_r3i1p1_DMI-HIRHAM5_v1",
           "_r1i1p1_IPSL-INERIS-WRF331F_v1",
           "_r1i1p1_KNMI-RACMO22E_v2",
           "_r1i1p1_CLMcom-CCLM4-8-17_v1",
           "_r1i1p1_SMHI-RCA4_v1a")
leto_kon <- c("21001231","21001231","21001231","20991130","21001231","21001231")


risi0_druge <- function(rcp, leto_zac, leto_zac1, obd, zuzelka){
  nc_data_rr<- nc_open(paste0(path,"pr","_12km_MOHC-HadGEM2-ES_",rcp,drugo[4],"_day_",leto_zac,"0101_",leto_kon[4],".nc"))
  time <- ncvar_get(nc_data_rr, "time")
  start1<-c(1,1,1)
  count1<-c(nx,ny, length(time))
  data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
  datumi0 <- seq(0,count1[3]-1,by = 1)
  pretvorba_sek_v_dan1 = 60*60*24 #?eprav imam dnevne podatke moram to dat, ker as.positxct meri v sekundah
  datumi <- as.POSIXct((start1[3]+datumi0)*pretvorba_sek_v_dan1,origin=paste0(leto_zac1,"-12-31 00:00:00"))
  leta <- format(datumi,format="%Y")
  datum2 = format(as.POSIXct(datumi, format="%Y-%m-%d"), "%m/%d/%Y")
  nc_close(nc_data_rr) # konec branja
  rattler.vsi <- data.frame()
  
  for(z in seq(1,6,by=1)){
    leto_zac <- "2071"
    leto_kon <- c("21001231","21001231","21001231","20991130","21001231","21001231")
    leto_zac1 <- "2070"
    # padavine dnevne
    nc_data_rr<- nc_open(paste0(path,"pr","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
    data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
    nc_close(nc_data_rr) # konec branja
    # maksimalna dnevna T
    nc_data_max<- nc_open(paste0(path,"tasmax","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
    data_max <- ncvar_get(nc_data_max, "tasmax", start=start1, count=count1)
    nc_close(nc_data_max) # konec branja
    # minimalna dnevna T
    nc_data_min <- nc_open(paste0(path,"tasmin","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon[z],".nc"))
    data_min <- ncvar_get(nc_data_min, "tasmin", start=start1, count=count1)
    nc_close(nc_data_min) # konec branja
    
    leto = format(datumi, format = "%Y")
    leta1 <- seq(min(leto),max(leto),by=1)
    lats1 <- seq(1,24,by=1)
    lons1 <- seq(1,40,by=1)
    n_rows <- 24
    n_cols <- 40
    n_layers <- 19
    bios2070 <- array(0, dim = c(n_rows, n_cols, n_layers))
    
    for(j in lats1){
      for(i in lons1){
        print(i)
        print(j) # podatki za T so v worldclim mnozeni z 10, zato da nimajo decimalnih mest, zato tudi jaz svoje mnozim, da bom lahko uporabila model
        data_rr1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), rr = data_rr[i,j,1:length(time)]*24*60*60) # ARSO podatki so v mm/s zato za dnevne pomno?im s ?tevilom s v dnevu
        data_min1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tn = 10*(data_min[i,j,1:length(time)]-273)) # ARSO podatki so v Kelvinih
        data_max1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tx = (data_max[i,j,1:length(time)]-273)*10)
        # mesecne vsote padavin in mesecna povprecja temperature
        dataRR <- aggregate(rr ~ mes + leto, data_rr1 , FUN = sum , na.rm=TRUE, na.action=na.pass)
        dataMIN <- aggregate(tn ~ mes + leto, data_min1 , mean , na.rm=TRUE, na.action=na.pass)
        dataMAX <- aggregate(tx ~ mes + leto, data_max1, mean , na.rm=TRUE, na.action=na.pass)
        # povprecje let obdobja 2071-2100
        dataRR <- aggregate(rr ~ leto, dataRR , mean , na.rm=TRUE, na.action=na.pass)
        dataMIN <- aggregate(tn ~ leto, dataMIN , mean , na.rm=TRUE, na.action=na.pass)
        dataMAX <- aggregate(tx ~ leto, dataMAX, mean , na.rm=TRUE, na.action=na.pass)
        bio_histSI <- biovars(prec = dataRR$rr, tmin = dataMIN$tn, tmax = dataMAX$tx)
        bios2070[j,i,1:19] = bio_histSI
      }}
    
    bios2070 <- bios2070[1:24, 1:40, c(1,5:9,12,16:19)]
    nlayers <- dim(bios2070)[3]  # nlayers mora biti lat_dim x lon_dim x BIO1-BIO19
    
    raster_stack2070 <- stack(
      lapply(1:nlayers, function(i) {
        layer <- matrix(bios2070[, , i], nrow = n_rows, ncol = n_cols)
        layer <- flipdim(layer, 1)
        raster_layer <- raster(layer)  # Define raster layer inside the loop
        extent_of_raster <- extent(min(lon), max(lon), min(lat), max(lat))
        extent(raster_layer) <- extent_of_raster
        projection(raster_layer) <- CRS("+proj=longlat +datum=WGS84")
        return(raster_layer)
      })
    )
    
    futureEnvSI <- raster_stack2070
    layer_names <- c("bio1", "bio5", "bio6", "bio7", "bio8", "bio9", "bio12", "bio16", "bio17", "bio18", "bio19")
    names(futureEnvSI) <- layer_names
    modelFutureEnvSI=crop(futureEnvSI, extent_of_raster)
    rattler.me <- readRDS(file = paste0("rattler_me",zuzelka,".rds"))
    rattler.2070SI <- predict(rattler.me, modelFutureEnvSI)
    rattler.2070SI_gg <- data.frame(rasterToPoints(rattler.2070SI))
    rattler.2070SI_gg$model_ime <- rep(z,length(rattler.2070SI_gg$x))
    rattler.vsi <- rbind(rattler.vsi, rattler.2070SI_gg)
  }
  saveRDS(rattler.vsi, file = paste0("rattler_",zuzelka,"_",rcp,"_",obd,".rds"))
  return(rattler.vsi)
}

risi_druge <- function(rcp, leto_zac, leto_kon, leto_zac1, obd, zuzelka){
  
  nc_data_rr<- nc_open(paste0(path,"pr","_12km_MOHC-HadGEM2-ES_",rcp,drugo[4],"_day_",leto_zac,"0101_",leto_kon,".nc"))
  time <- ncvar_get(nc_data_rr, "time")
  start1<-c(1,1,1)
  count1<-c(nx,ny, length(time))
  data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
  # povpre?na dnevna T
  datumi0 <- seq(0,count1[3]-1,by = 1)
  pretvorba_sek_v_dan1 = 60*60*24 #?eprav imam dnevne podatke moram to dat, ker as.positxct meri v sekundah
  datumi <- as.POSIXct((start1[3]+datumi0)*pretvorba_sek_v_dan1,origin=paste0(leto_zac1,"-12-31 00:00:00"))
  leta <- format(datumi,format="%Y")
  datum2 = format(as.POSIXct(datumi, format="%Y-%m-%d"), "%m/%d/%Y")
  nc_close(nc_data_rr) # konec branja
  rattler.vsi <- data.frame()
  
  for(z in seq(1,6,by=1)){
    # padavine dnevne
    nc_data_rr<- nc_open(paste0(path,"pr","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon,".nc"))
    data_rr <- ncvar_get(nc_data_rr, "pr", start=start1, count=count1)
    nc_close(nc_data_rr) # konec branja
    # maksimalna dnevna T
    nc_data_max<- nc_open(paste0(path,"tasmax","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon,".nc"))
    data_max <- ncvar_get(nc_data_max, "tasmax", start=start1, count=count1)
    nc_close(nc_data_max) # konec branja
    # minimalna dnevna T
    nc_data_min <- nc_open(paste0(path,"tasmin","_12km_",model[z],"_",rcp,drugo[z],"_day_",leto_zac,"0101_",leto_kon,".nc"))
    data_min <- ncvar_get(nc_data_min, "tasmin", start=start1, count=count1)
    nc_close(nc_data_min) # konec branja
    
    leto = format(datumi, format = "%Y")
    leta1 <- seq(min(leto),max(leto),by=1)
    lats1 <- seq(1,24,by=1)
    lons1 <- seq(1,40,by=1)
    n_rows <- 24
    n_cols <- 40
    n_layers <- 19
    bios2070 <- array(0, dim = c(n_rows, n_cols, n_layers))
    
    for(j in lats1){
      for(i in lons1){
        print(i)
        print(j) # podatki za T so v worldclim mno?eni z 10, zato da nimajo decimalnih mest, zato tudi jaz svoje mno?im, da bom lahko uporabila model
        data_rr1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), rr = data_rr[i,j,1:length(time)]*24*60*60) # ARSO podatki so v mm/s zato za dnevne pomno?im s ?tevilom s v dnevu
        data_min1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tn = 10*(data_min[i,j,1:length(time)]-273)) # ARSO podatki so v Kelvinih
        data_max1 <- data.frame(datum = datumi, mes = months(datumi), dan = format(datumi, format = "%d"), leto = format(datumi, format = "%Y"), tx = (data_max[i,j,1:length(time)]-273)*10)
        
        dataRR <- aggregate(rr ~ mes + leto, data_rr1 , FUN = sum , na.rm=TRUE, na.action=na.pass)
        dataMIN <- aggregate(tn ~ mes + leto, data_min1 , mean , na.rm=TRUE, na.action=na.pass)
        dataMAX <- aggregate(tx ~ mes + leto, data_max1, mean , na.rm=TRUE, na.action=na.pass)
        
        dataRR <- aggregate(rr ~ leto, dataRR , mean , na.rm=TRUE, na.action=na.pass)
        dataMIN <- aggregate(tn ~ leto, dataMIN , mean , na.rm=TRUE, na.action=na.pass)
        dataMAX <- aggregate(tx ~ leto, dataMAX, mean , na.rm=TRUE, na.action=na.pass)
        
        bio_histSI <- biovars(prec = dataRR$rr, tmin = dataMIN$tn, tmax = dataMAX$tx)
        bios2070[j,i,1:19] = bio_histSI
      }}
    
    bios2070 <- bios2070[1:24, 1:40, c(1,5:9,12,16:19)]
    nlayers <- dim(bios2070)[3]  # nlayers mora biti lat_dim x lon_dim x BIO1-BIO19
    raster_stack2070 <- stack(
      lapply(1:nlayers, function(i) {
        layer <- matrix(bios2070[, , i], nrow = n_rows, ncol = n_cols)
        layer <- flipdim(layer, 1)
        raster_layer <- raster(layer)  # Define raster layer inside the loop
        extent_of_raster <- extent(min(lon), max(lon), min(lat), max(lat))
        # Assign the extent to the raster layer
        extent(raster_layer) <- extent_of_raster
        # Define the projection of the raster layer (WGS84 in this example)
        projection(raster_layer) <- CRS("+proj=longlat +datum=WGS84")
        return(raster_layer)
      })
    )
    
    futureEnvSI <- raster_stack2070
    layer_names <- c("bio1", "bio5", "bio6", "bio7", "bio8", "bio9", "bio12", "bio16", "bio17", "bio18", "bio19")
    names(futureEnvSI) <- layer_names
    modelFutureEnvSI=crop(futureEnvSI, extent_of_raster)
    rattler.me <- readRDS(file = paste0("rattler_me",zuzelka,".rds"))
    rattler.2070SI <- predict(rattler.me, modelFutureEnvSI)
    rattler.2070SI_gg <- data.frame(rasterToPoints(rattler.2070SI))
    rattler.2070SI_gg$model_ime <- rep(z,length(rattler.2070SI_gg$x))
    rattler.vsi <- rbind(rattler.vsi, rattler.2070SI_gg)
  }

  saveRDS(rattler.vsi, file = paste0("rattler_",zuzelka,"_",rcp,"_",obd,".rds"))
  return(rattler.vsi)
}
  
# 
# rattler.vsi0 <- risi0_druge(rcp = "rcp45", leto_zac = "2071", leto_zac1 = "2070", obd = "2071-2100", zuzelka = zuzelka1)
# rattler.vsi1 <- risi0_druge(rcp = "rcp85", leto_zac = "2071", leto_zac1 = "2070", obd = "2071-2100", zuzelka = zuzelka1)
rattler.vsi0 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp45_2071-2100.rds"))
rattler.vsi1 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp85_2071-2100.rds"))

grafi(rattler.vsi = rattler.vsi0, rcp = "rcp45", obd = "2071-2100", zuzelka = zuzelka1)
grafi(rattler.vsi = rattler.vsi1, rcp = "rcp85", obd = "2071-2100", zuzelka = zuzelka1)

# rattler.vsi2 <- risi_druge(rcp = "rcp45", leto_zac = "2041", leto_kon = "20701231", leto_zac1 = "2040", obd = "2041-2070", zuzelka = zuzelka1)
# rattler.vsi3 <- risi_druge(rcp = "rcp85", leto_zac = "2041", leto_kon = "20701231", leto_zac1 = "2040", obd = "2041-2070", zuzelka = zuzelka1)
rattler.vsi2 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp45_2041-2070.rds"))
rattler.vsi3 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp85_2041-2070.rds"))

# rattler.vsi4 <- risi_druge(rcp = "rcp45", leto_zac = "2011", leto_kon = "20401231", leto_zac1 = "2010", obd = "2011-2040", zuzelka = zuzelka1)
# rattler.vsi5 <- risi_druge(rcp = "rcp85", leto_zac = "2011", leto_kon = "20401231", leto_zac1 = "2010", obd = "2011-2040", zuzelka = zuzelka1)
rattler.vsi4 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp45_2011-2040.rds"))
rattler.vsi5 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp85_2011-2040.rds"))

grafiK(rattler.vsi = rattler.vsi2, rcp = "rcp45", obd = "2041-2070", zuzelka = zuzelka1)
grafiK(rattler.vsi = rattler.vsi3, rcp = "rcp85", obd = "2041-2070", zuzelka = zuzelka1)

grafiK(rattler.vsi = rattler.vsi4, rcp = "rcp45", obd = "2011-2040", zuzelka = zuzelka1)
grafiK(rattler.vsi = rattler.vsi5, rcp = "rcp85", obd = "2011-2040", zuzelka = zuzelka1)





zuzelka1 = 'Eiphosoma laphygmae'
rattler<-gbif(zuzelka1)
a <- data.frame(ime = rep(zuzelka1,length(rattler$datasetKey)),rattler$datasetKey)

zuzelka1 = 'Spodoptera frugiperda'
rattler<-gbif(zuzelka1)
b <- data.frame(ime = rep(zuzelka1,length(rattler$datasetKey)),rattler$datasetKey)

zuzelka1 = 'Trichogramma pretiosum'
rattler<-gbif(zuzelka1)
c <- data.frame(ime = rep(zuzelka1,length(rattler$datasetKey)),rattler$datasetKey)

zuzelka1 = 'Cotesia marginiventris'
rattler<-gbif(zuzelka1)
d <- data.frame(ime = rep(zuzelka1,length(rattler$datasetKey)),rattler$datasetKey)

zuzelka1 = 'Telenomus remus'
rattler<-gbif(zuzelka1)
e <- data.frame(ime = rep(zuzelka1,length(rattler$datasetKey)),rattler$datasetKey)

skupaj <- rbind(a,b,c,d,e)
saveRDS(skupaj, file = paste0("Occurrence_datasetKey.rds"))
write.csv(skupaj, file = paste0("Occurrence_datasetKey.csv"))

# zuzelka1 = 'Spodoptera frugiperda'
# zuzelka1 = 'Trichogramma pretiosum'
# zuzelka1 = 'Cotesia marginiventris'
# zuzelka1 = 'Telenomus remus'
zuzelka1 = 'Eiphosoma laphygmae'
rattler.me1 <- readRDS(file = paste0("rattler_me",zuzelka1,".rds"))
rattler.vsi2 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp45_2041-2070.rds"))
rattler.vsi3 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp85_2041-2070.rds"))

# rattler.vsi4 <- risi_druge(rcp = "rcp45", leto_zac = "2011", leto_kon = "20401231", leto_zac1 = "2010", obd = "2011-2040", zuzelka = zuzelka1)
# rattler.vsi5 <- risi_druge(rcp = "rcp85", leto_zac = "2011", leto_kon = "20401231", leto_zac1 = "2010", obd = "2011-2040", zuzelka = zuzelka1)
rattler.vsi4 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp45_2011-2040.rds"))
rattler.vsi5 <- readRDS(file = paste0("rattler_",zuzelka1,"_rcp85_2011-2040.rds"))


grafi_SLO <- function(rattler.vsi, rcp, obd, zuzelka){
  rattler.mediana <- aggregate(layer ~  x + y, rattler.vsi, FUN = median , na.rm=TRUE, na.action=na.pass)
  rattler.min <- aggregate(layer ~  x + y, rattler.vsi, FUN = "min" , na.rm=TRUE, na.action=na.pass)
  rattler.max <- aggregate(layer ~  x + y, rattler.vsi, FUN = "max" , na.rm=TRUE, na.action=na.pass)
  
  koruza <- read_sf("koruza/GRID_Koruza_ha.shp")
  koruza <- st_transform(koruza, crs = 4326)
  koruza$mid <- st_centroid(koruza$geometry) 
  koruza$Koruza <- koruza$Koruza_ha
  
  ggplot(data=rattler.mediana) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    geom_sf(data = koruza, aes(geometry = mid, size = Koruza), color = "black", show.legend = "point") +
    scale_fill_distiller("Ustreznost", palette = "Spectral", limits=c(0,1)) +
    geom_sf(data = slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
    labs(x="G. dol?ina",y="G. ?irina",title=paste0("Verjetnost pojava ?u?elke ",zuzelka,", ",obd,", Mediana, ",rcp)) +
    coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
    theme_light() + theme(legend.position = c(.9, .3), legend.title = element_text(face = "bold", size=16),
                          text = element_text(size=18),
                          plot.title = element_text(size=18))
  ggsave(paste0("proj/",obd,"/SLOPredicted_suitability_",zuzelka,"_",obd,"_median_",rcp,".png"),width = 10, height = 7)
  
  ggplot(data=rattler.min) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    scale_fill_distiller("Ustreznost", palette = "Spectral", limits=c(0,1)) +
    geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
    geom_sf(data = koruza, aes(geometry = mid, size = Koruza), color = "black", show.legend = "point") +
    labs(x="G. dol?ina",y="G. ?irina",title=paste0("Verjetnost pojava ?u?elke ",zuzelka,", ",obd,", Min, ",rcp)) +
    coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
    theme_light() + theme(legend.position = c(.9, .3), legend.title = element_text(face = "bold", size=16),
                          text = element_text(size=18),
                          plot.title = element_text(size=18))
  ggsave(paste0("proj/",obd,"/SLOPredicted_suitability_",zuzelka,"_",obd,"_min_",rcp,".png"),width = 10, height = 7)
  
  ggplot(data=rattler.max) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    scale_fill_distiller("Ustreznost", palette = "Spectral", limits=c(0,1)) +
    geom_sf(data = koruza, aes(geometry = mid, size = Koruza), color = "black", show.legend = "point") +
    geom_sf(data = slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
    labs(x="G. dol?ina",y="G. ?irina",title=paste0("Verjetnost pojava ?u?elke ",zuzelka,", ",obd,", Max, ",rcp)) +
    coord_sf(crs = st_crs(4326)) + xlim(c(start_lon,end_lon)) + ylim(c(start_lat,end_lat)) +
    theme_light() + theme(legend.position = c(.9, .3), legend.title = element_text(face = "bold", size=16),
                          text = element_text(size=18),
                          plot.title = element_text(size=18))
  ggsave(paste0("proj/",obd,"/SLOPredicted_suitability_",zuzelka,"_",obd,"_max_",rcp,".png"),width = 10, height = 7)
}
grafi_SLO(rattler.vsi = rattler.vsi2, rcp = "rcp45", obd = "2041-2070", zuzelka = zuzelka1)
grafi_SLO(rattler.vsi = rattler.vsi3, rcp = "rcp85", obd = "2041-2070", zuzelka = zuzelka1)

grafi_SLO(rattler.vsi = rattler.vsi4, rcp = "rcp45", obd = "2011-2040", zuzelka = zuzelka1)
grafi_SLO(rattler.vsi = rattler.vsi5, rcp = "rcp85", obd = "2011-2040", zuzelka = zuzelka1)


######### KORUZA povr?ina njiv s koruzo za leto 2022

# Load shapefile
koruza <- read_sf("koruza/GRID_Koruza_ha.shp")
ggplot() +
  geom_sf(data = koruza, aes(fill = Koruza_ha), colour=alpha("black",0.1)) +
  geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
  labs(x="Longitude",y="Latitude",title=paste0("Total area of corn fields, 2022")) +
  scale_fill_distiller(name = "Area (ha)", palette = "YlOrBr", na.value="white") +
  geom_point() + theme_light(base_size = 17) 
ggsave(paste0("koruza_2022.png"),width = 10, height = 7)


ggplot() +
  geom_sf(data = koruza, aes(fill = Koruza_ha), colour=alpha("black",0.1)) +
  geom_sf(data=slovenia_nuts3_mapdata, color=alpha("black",0.4),fill = NA) +
  labs(x="Longitude",y="Latitude",title=paste0("Total area of corn fields, 2022")) +
  scale_fill_distiller(name = "Area (ha)", palette = "YlOrBr", na.value="white") +
  geom_point() + theme_light(base_size = 17) 



rattler.mediana <- aggregate(layer ~  x + y, rattler.vsi, FUN = median , na.rm=TRUE, na.action=na.pass)
rattler.min <- aggregate(layer ~  x + y, rattler.vsi, FUN = "min" , na.rm=TRUE, na.action=na.pass)
rattler.max <- aggregate(layer ~  x + y, rattler.vsi, FUN = "max" , na.rm=TRUE, na.action=na.pass)

koruza <- read_sf("koruza/GRID_Koruza_ha.shp")
koruza <- st_transform(koruza, crs = 4326)
koruza$mid <- st_centroid(koruza$geometry) 
koruza$Area <- koruza$Koruza_ha

# Plot the data
ggplot() +
  geom_raster(data = rattler.mediana, aes(x = x, y = y, fill = layer)) +
  scale_fill_distiller("Probability", palette = "Spectral", limits = c(0, 1)) +
  geom_sf(data = slovenia_nuts3_mapdata, color = alpha("black", 0.4), fill = NA) +
  geom_sf(data = koruza, aes(geometry = mid, size = Area), color = "black", show.legend = "point") +  # Adjusted line
  labs(x = "Longitude", y = "Latitude", title = paste0("Predicted suitability for ", zuzelka, ", ", obd, ", Median and ", rcp)) +
  coord_sf(crs = st_crs(4326)) + xlim(c(start_lon, end_lon)) + ylim(c(start_lat, end_lat)) +
  theme_light() + theme(legend.position = c(.9, .2),
    legend.title = element_text(face = "bold"),
    text = element_text(size = 17),
    plot.title = element_text(size = 17))




### IZ HELPA - kako narisati na sredini poligona to?ko velikosti vrednosti za tisto mre?no to?ko
if (requireNamespace("sf", quietly = TRUE)) {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  ggplot(nc) +
    geom_sf(aes(fill = AREA))
  
  # If not supplied, coord_sf() will take the CRS from the first layer
  # and automatically transform all other layers to use that CRS. This
  # ensures that all data will correctly line up
  nc_3857 <- sf::st_transform(nc, 3857)
  ggplot() +
    geom_sf(data = nc) +
    geom_sf(data = nc_3857, colour = "red", fill = NA)
  
  # Unfortunately if you plot other types of feature you'll need to use
  # show.legend to tell ggplot2 what type of legend to use
  nc_3857$mid <- sf::st_centroid(nc_3857$geometry)
  ggplot(nc_3857) +
    geom_sf(colour = "white") +
    geom_sf(aes(geometry = mid, size = AREA), show.legend = "point")
}







