################################################################################
### Satellite Data
################################################################################

# Package names
packages <- c("fields","raster","terra","luna", "dplyr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

##########
## Load MODIS data
##########

# Domain
aoi <- c(-95.91153, -91.28381, 34.29519, 37.06811)
# Path
datadir <- getwd()
dir.create(datadir, showWarnings=FALSE)
mf <- getNASA("MOD11A1", "2023-08-20","2023-08-20", aoi=aoi, download = T,
              path=datadir, username="tgyger", password="Julia4$Marianne")

# Raster
r <- raster(mf[1])
# Raster to data
data <- rasterToPoints(r)
# Scaled Kelvin to degrees Celcius
data[,3] <- data[,3]*0.02-273.15

data <- as.data.frame(data)
colnames(data) <- c("east","north","temp")
quilt.plot(data[,1],data[,2],data[,3],nx = 200)

# Define Test coordinates
mf <- getNASA("MOD11A1", "2023-09-06","2023-09-06", aoi=aoi, download = T,
              path=datadir, username="tgyger", password="Julia4$Marianne")

# Raster
r <- raster(mf[1])
# Raster to data
data2 <- rasterToPoints(r)
quilt.plot(data2[,1],data2[,2],data2[,3],nx = 200)
data2 <- as.data.frame(data2)
colnames(data2) <- c("east","north","temp2")

# Combine data frames
data <- data %>% dplyr::mutate(eastnorth = paste0(data$east,data$north))
data2 <- data2 %>% dplyr::mutate(eastnorth = paste0(data2$east,data2$north))
data2 <- data2 %>% dplyr::select("eastnorth","temp2")
data_combined <- data %>% left_join(data2,by = join_by(eastnorth)) %>% dplyr::select(east,north,temp,temp2)

data_combined <- data_combined %>% dplyr::mutate(temp2 = ifelse(!is.na(temp2),temp,NA))  %>% dplyr::mutate(temp = ifelse(is.na(temp2),temp,NA)) %>% dplyr::mutate(temp3 = ifelse(is.na(temp2),temp,temp2)) %>%
  sample_n(800000, replace = F)

# Order (not necessary)
normed1 <- (data_combined[,1]-min(data_combined[,1]))/max(data_combined[,1]-min(data_combined[,1]))
normed2 <- (data_combined[,2]-min(data_combined[,2]))/max(data_combined[,2]-min(data_combined[,2]))
data_combined <- data_combined[order((normed1)^2+(normed2)^2),]
# Train and Test data
set.seed(10)
data_train <- data_combined %>% tidyr::drop_na(temp2) %>% dplyr::select(east,north,temp2) %>% rename(temp = temp2) %>% sample_n(400000)
data_test <- data_combined %>% tidyr::drop_na(temp) %>% dplyr::select(east,north,temp) %>% sample_n(200000)

# Plots
par(mfrow = c(1, 2),mai = c(0.5, 0.5, 0.5, 0.05))
quilt.plot(data_combined[,1],data_combined[,2],data_combined[,5],nx = 400,ny = 400,legend.width = 0.5,xlab = "Easting (m)",ylab = "Northing (m)",main="Total (degrees Celcius)",cex.axis=1.8,axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)
quilt.plot(data_train[,1],data_train[,2],data_train[,3],nx = 400,ny = 400,legend.width = 0.5,xlab = "Easting (m)",ylab = "Northing (m)",cex.axis=1.8,main="Training (degrees Celcius)",axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)