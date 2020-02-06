rm(list = ls())

library(raster)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(colorRamps)
library(readr)
library(Metrics)

# df = stack(paste0("/Users/", Sys.info()[7], "/Desktop/HadISST_sst.nc"), varname = "sst")
load("~/CC_LME_SST/Hadley_SST.RData")

e = extent(-140, -100, 22.50, 47.50) #California Current LME lat-lon range

cc = crop(df, e); rm(e)

# change to data.frame
df <- cc %>% rasterToPoints() %>% data.frame()

# 1955-2012 climatology
cc_1955_2012 <- cc[[1021:1716]] #january 1955 - december 2012
plot(mean(cc_1955_2012), col = matlab.like(100), axes = F); maps::map(add = T, fill = T, col = "gray"); box(); degAxis(1); degAxis(2, las = 1)

# trim Hadley SST to 1955-2018
cc_1870_2019 <- cc[[1:1798]] #Jan 1870 - Oct 2019
plot(mean(cc_1870_2019), col = matlab.like(100), axes = F); maps::map(add = T, fill = T, col = "gray"); box(); degAxis(1); degAxis(2, las = 1)

# calcualte SST anomalies to mean(1955-2012)
cc_anomaly = cc_1870_2019 - mean(cc_1955_2012, na.rm = T)
plot(mean(cc_anomaly), col = matlab.like(100), axes = F); maps::map(add = T, fill = T, col = "gray"); box(); degAxis(1); degAxis(2, las = 1)
names(cc_anomaly) = names(cc_1870_2019)
plot(cc_anomaly, col = matlab.like(100), axes = F)

# change to data.frame
df <- cc_anomaly %>% rasterToPoints() %>% data.frame()

#add lme
names(df)
latlon = df[,c(2,1)]
coordinates(latlon)=~x+y
lme<-rgdal::readOGR("/Users/ktanaka/Google Drive/Research/GIS/LME66/LMEs66.shp")
CRS.new<-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  #EPSG:102003
proj4string(latlon) <- CRS.new 
proj4string(lme) <- CRS.new
area <- over(latlon,lme)
colnames(area)[1] = "lme"
df = cbind(df, area)
df = df[which(df$LME_NUMBER == 3),]
rm(area, latlon, lme, CRS.new)

#aggregate SST mean by year, skip last 9 rows
sst_bin = df[,c(3:(length(df)-9))]
group = names(sst_bin)
year = substr(group, 2,5)
month = substr(group, 7,8)
group = cbind(year, month)
colnames(group) = c("Year", "Month")
sst_bin = as.data.frame(t(sst_bin))
sst_bin = cbind(group, sst_bin)
sst_bin$Mean = rowMeans(sst_bin[,3:dim(sst_bin)[2]], na.rm = T)
sst_bin = sst_bin[c("Year", "Month", "Mean")]
# sst_bin = sst_bin[which(sst_bin$Year %in% c(1986:2017)),]
rm(group, month, year)

# calcuate trend
d = summarySE(sst_bin, measurevar = "Mean", groupvars = c("Year"))
d$Year = as.numeric(as.character(d$Year))
summary(lm(Mean~Year, data = d))

# plot trend
ggplot(d, aes(x=Year, y=Mean, col = Mean)) +
  geom_errorbar(aes(ymin = Mean-ci, ymax = Mean+ci), width = 0.5) +
  geom_line() +
  geom_point() +
  theme_pubr() + 
  geom_smooth(method='lm', se = F, col = "gray", size = 0.5) +
  scale_color_gradientn(colours = matlab.like(100), "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous("", breaks = round(seq(min(d$Year), max(d$Year), by = 10),1)) +
  scale_y_continuous("SST anomaly relatuve to 1955-2012 (deg C)", limit = c(-3,3), breaks = round(seq(-4, 4, by = 1),1)) + 
  theme(legend.position = c(0.05,0.9))

rm(d)

###################################################
### use NWA climatology and add CM2.1 anomalies ###
###################################################

#read Calicornia Current climatology 1955-2012
setwd("~/Desktop/")
# clim = read_csv("nep_all_t00an10.csv", skip = 1) #1955-2012 0.1 * 0.1
# clim = read_csv("nep_5564_t00an10.csv", skip = 1) #1955-1964, 0.1 dec deg
# clim = read_csv("nep_6574_t00an10.csv", skip = 1) #1965-1974, 0.1 dec deg
# clim = read_csv("nep_7584_t00an10.csv", skip = 1) #1975-1984, 0.1 dec deg
clim = read_csv("nep_8594_t00an10.csv", skip = 1) #1985-1994, 0.1 dec deg
# clim = read_csv("nep_95A4_t00an10.csv", skip = 1) #1995-2004, 0.1 dec deg
# clim = read_csv("nep_A5B2_t00an10.csv", skip = 1) #2005-2012, 0.1 dec deg


clim = clim[,1:3]
colnames(clim) = c("lat", "lon", "sst")

#add lme
names(clim)
latlon = clim[,c(2,1)]
coordinates(latlon)=~lon+lat
lme<-rgdal::readOGR("/Users/ktanaka/Google Drive/Research/GIS/LME66/LMEs66.shp")
CRS.new<-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  #EPSG:102003
proj4string(latlon) <- CRS.new 
proj4string(lme) <- CRS.new
area <- over(latlon,lme)
colnames(area)[1] = "lme"
clim = cbind(clim, area)
clim = clim[which(clim$LME_NUMBER == 3),]
clim = clim[,c("lon", "lat", "LME_NAME", "SUM_GIS_KM", "sst")]
rm(area, latlon, lme, CRS.new)

cc = vector("list")

for (i in 1:dim(sst_bin)[1]) {
  
  sst = clim$sst + sst_bin[i,3]
  cc[[i]] = sst
  print(i)
  
}

cc = bind_cols(cc)

cc = as.data.frame(cc)
colnames(cc) = rownames(sst_bin)

names(cc) <- substring(names(cc), 2, 8)
names(cc) <- gsub(x = names(cc), pattern = "\\.", replacement = "-")
names(cc)

cc = cbind(clim[,c("lon", "lat", "LME_NAME", "SUM_GIS_KM")], cc)

d = cc
d$mean = rowMeans(cc[,c(5:1802)])

world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf") #worldwide country polygon

ggplot(d) +
  geom_raster(aes(lon, lat, fill = mean), show.legend = T) +
  geom_sf(data = world, fill = "gray", size = .1, color = "black") +
  scale_fill_gradientn(colours = matlab.like(100), "mean SST 1870-2019") +
  coord_sf(xlim = range(cc$lon), ylim = range(cc$lat)) +
  scale_x_continuous(expand = c(-0, 0)) +
  scale_y_continuous(expand = c(-0, 0)) +
  theme_pubr() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "right")

###############################
### copmare with Hadley_SST ###
###############################
hadley = stack(paste0("/Users/", Sys.info()[7], "/Desktop/HadISST_sst.nc"), varname = "sst")

e = extent(-131.5, -109.5, 22.50, 47.50) #CC

hadley = crop(hadley, e); rm(e)

# change to data.frame
df <- hadley %>% rasterToPoints() %>% data.frame()

#add lme
names(df)
latlon = df[,c(1,2)]; plot(latlon)
coordinates(latlon)=~x+y
lme<-rgdal::readOGR("/Users/ktanaka/Google Drive/Research/GIS/LME66/LMEs66.shp")
CRS.new<-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  #EPSG:102003
proj4string(latlon) <- CRS.new 
proj4string(lme) <- CRS.new
area <- over(latlon,lme)
colnames(area)[1] = "lme"
df = cbind(df, area)
df = df[which(df$LME_NUMBER == 3),]
rm(area, latlon, lme, CRS.new)
plot(df$x, df$y)

#aggregate SST mean by year, skip last 9 rows
sst_bin = df[,c(3:(length(df)-9))]
group = names(sst_bin)
year = substr(group, 2,5)
month = substr(group, 7,8)
group = cbind(year, month)
colnames(group) = c("Year", "Month")
sst_bin = as.data.frame(t(sst_bin))
sst_bin = cbind(group, sst_bin)
sst_bin$Mean = rowMeans(sst_bin[,3:dim(sst_bin)[2]], na.rm = T)
sst_bin = sst_bin[c("Year", "Month", "Mean")]
# sst_bin = sst_bin[which(sst_bin$Year %in% c(1982:2017)),]
rm(group, month, year)

# calcuate trend
d = summarySE(sst_bin, measurevar = "Mean", groupvars = c("Year"))
d$Year = as.numeric(as.character(d$Year))
summary(lm(Mean~Year, data = d))

hadley = d

# hadley sst trend 1955-2018
cc = cc[,c(5:dim(cc)[2])]
group = names(cc)
group = stringr::str_split_fixed(group, "-", 2)
colnames(group) = c("Year", "Month")
cc = as.data.frame(t(cc))
cc = cbind(group, cc)
cc$mean = rowMeans(cc[,c(3:dim(cc)[2])], na.rm = T)
cc = cc[c("Year", "mean")]

d = summarySE(cc, measurevar = "mean", groupvars = c("Year"))
d$Year = as.numeric(as.character(d$Year))

cc = d

plot(hadley$Year, hadley$Mean, 
     type = 'o', col = 'blue', pch = 20, axes = F, lwd = 2, 
     ylab = "deg C", xlab = "Year")
box(bty = "l"); axis(1); axis(2, las = 2)
points(cc$Year, cc$mean, type = 'o', 
       col = 'red', pch = 20, lwd = 2,  bty = "l")

