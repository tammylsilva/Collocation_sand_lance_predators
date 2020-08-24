#Calculate indicies of coolocation and make maps of metrics and sand lance sites
#1/9/2019 Tammy Silva

#refs: Bez & Riovoiranard 2000, ICES Geostats in R

#This script processes data, analyses and create figures in the manuscript:
#High collocation between sand lance and top predators: implications for conservation and management
#Silva et al 2020 Conservation Science and Practice 

#this script uses visual observations (count data) of humpback whales and great 
#shearwaters and count data of sand lance from USGS SeaBOSS grab samples to calculate
#Center of Gravity, Inertia and Global Index of Collocation by cruise and for all cruises combined and
#creates figures

library(RGeostats)
library(tidyverse)
library(rgdal)
library(raster)
library(RColorBrewer)
library(maps) #for scale bars on maps
library(GISTools) #for north arrow on maps
library(png) #for bringing in png and embedding in maps


###Load in and clean data

#read in subset of visual observations with just humpback and great shearwater and sand lance data
vis_obs <- read.csv(file="./data/seaboss/vis_obs_subset_coloc_analysis.csv", header=TRUE, sep=",")

#pull out just sl, great shearwater and humpback whale observations
vis_obs<-dplyr::select(vis_obs, local_date, local_time, site, sl, contains("GRSH"), contains("HUWH"))

#add up columns for total shearwater and whale sightings
vis_obs$b_total <- rowSums(dplyr::select(vis_obs, contains("GRSH")))
vis_obs$w_total <- rowSums(dplyr::select(vis_obs, contains("HUWH")))

#filter rows with no data 'nd' for sand lance counts
vis_obs <- filter(vis_obs, sl !="nd")

#for some reason, sl is read in as a factor (bc of 'nd') so need to convert to character and the integer
vis_obs$sl <- as.integer(as.character(vis_obs$sl))

#format local_date as dates
vis_obs$local_date<-as.Date(vis_obs$local_date, format="%m/%d/%Y")

#create a column with just the year
vis_obs$year<-format(vis_obs$local_date,format="%Y")
vis_obs$year<- as.factor(vis_obs$year)

#create a column with just the month
vis_obs$month<-format(vis_obs$local_date,format="%m")


#read in site data
sites<-read.csv(file="./data/seaboss/sites.csv", header=TRUE)

#remove rows with NA (no positions for stations)
sites<-na.omit(sites)

#get rid of LSB and R sites
sites <- filter(sites, !grepl('LSB|R', site))

#split site into 2 columns-letters vs numbers
sites <- separate(sites, site, into = c("text", "num"), sep = "(?<=[A-Za-z])(?=[0-9])")
sites$num <- as.numeric(sites$num)

#only keep the standard 44 sites and remove all others
sites <- filter(sites, num < 17)

#put columns back together to make sites normal again
sites <- unite(sites, site, c("text", "num"),sep = "")

#sites to factor
sites$site <- as.factor(sites$site)

#remove extra S and C sites
sites <- filter(sites, site !="S15")
sites <- filter(sites, site !="S16")
sites <- filter(sites, site !="C15")
sites <- filter(sites, site !="C16")

#make sites lon minus
sites$lon<- -sites$lon

#merge sites and vis obs data frames to pair lat/lon with site name for plotting / mapping
data <- inner_join(vis_obs, sites, by="site")

#make the data frame smaller with just cols I need so RGeostats db.dataframe won't crash (RGeostats only takes numeric data)
data <- dplyr::select(data, local_date, site, sl, w_total, b_total, lat, lon, year, month)

#fill in season
data <- mutate(data, season = ifelse(month=='05'| month=='06', 'spring', ifelse(month=='11' | month=='09'| month=='10','fall','summer')))

#remove 2019 since not included in analysis
data <- data %>% filter(year !='2019')

#remove fall 2018 since not a complete cruise
data <- data %>% filter(local_date !='2018-11-30')
data <- data %>% filter(local_date !='2018-11-9')

#assign a cruise number (or season and year)
data <- data %>% 
  #unite("cruise",season, year, sep="", remove=FALSE)
  mutate(cruise = group_indices_(data, .dots=c("season", "year")))

data$cruise<-as.factor(data$cruise)

#########################################################

###Calculations


###Tally the number of sl, whales and birds for each cruise and number of sites surveyed
#(Table 1)
counts<-data %>% group_by(cruise) %>% summarise(sl=sum(sl), whales=sum(w_total),
                                                birds=sum(b_total),n_distinct(site))


###Calculate Global Index of Collocation by cruise
#(Table 3)

#remove non numeric columns - RGeostats db's cannot handle anything but numeric columns and it crashes R
cruise_data<-dplyr::select(data,cruise, sl,w_total,b_total,lat,lon)

#select cruise
cruise<-filter(cruise_data, cruise=='11')

#create db data 
db.data = db.create(cruise)

#locators allow the user to define the role attributed to each attribute
#define the attribtues (variables) that will serve as coordinates (lat/lon)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)

#designate the attribute(variable) that is of interest-no limit to number of attributes that can be attached to this locator
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

#compute global index of colocation and automatically generate a generic plot
SI.gic(db1=db.data,db2=db.data,name1="sl",name2="w_total",flag.plot=T,flag.inertia=T,
       flag.ellipse=T,ylim=c(42.1,42.45))



###Calculate GIC for all cruises combined

# #sum counts over all years
data_total <- data %>% group_by(site) %>% summarise(sl=sum(sl), w=sum(w_total), b=sum(b_total))

cruise<-inner_join(data_total, sites, by='site')

#remove non numeric columns - RGeostats db's cannot handle anything but numeric columns and it crashes R
cruise_data<-dplyr::select(cruise, sl,w,b,lat,lon)

#create db data 
db.data = db.create(cruise_data)

#locators allow the user to define the role attributed to each attribute
#define the attribtues (variables) that will serve as coordinates (lat/lon)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)

#designate the attribute(variable) that is of interest-no limit to number of attributes that can be attached to this locator
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w","z",2)
db.data = db.locate(db.data,"b","z",3)

#compute global index of colocation and plot
SI.gic(db1=db.data,db2=db.data,name1="sl",name2="w",flag.plot=T,flag.inertia=T,
       flag.ellipse=T,ylim=c(42.1,42.45))

########################################################

###Figures

#read in shapefiles

#read in sbnms shapefile
sbnms<-rgdal::readOGR(dsn="./data/environmental", layer="sbnms")#with rgdal package

#load in sediment data
str_name2<-'./data/environmental/mb_backpc30m.tif'
sub_raster=raster(str_name2)
newproj<-proj4string(sbnms)
sub_raster<-projectRaster(sub_raster, crs=newproj)
sub_crop<-crop(sub_raster, sbnms)
sub_mask<-mask(sub_crop, sbnms)
plot(sub_mask, breaks=c(75,165,255),col=c("gray78","gray54","gray27"))
plot(sbnms, add=T)

#load in a 40m isobath shapefile
bathy40<-rgdal::readOGR(dsn="./data/environmental", layer="bathy40m")#with rgdal package
bathy40<-spTransform(bathy40, CRS(proj4string(sbnms)))


##Figure 1 - map of study area

tiff( "./figures/Collocation_manuscript/20190711_sites_gray_substrate.tif",width =7, height = 10, units = 'in', res = 300) #create a tiff device, sized and 300 dpi res
plot(sub_crop,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),xlab="", ylab="",cex.lab=1.3,cex.axis=1.3, axis.args=list(cex.axis=1.3))
points(sites$lon, sites$lat, pch=21,cex=1.4, col="black",bg="black")
plot(sbnms,lwd=4, add=T) 
plot(bathy40, add=T)
maps::map.scale(x=-70.57, y=42.75, cex=1.4,metric=TRUE, ratio=FALSE)
north.arrow(xb=-70.57, y=42.67,len=0.008, lab="N")
dev.off() 



###Figure 2 - center of gravity and inertia for all species by cruise

tiff( "./figures/Collocation_manuscript/revisions_conscipractice/colloc_all_spp_by_cruise.tif",width = 10, height = 8, units = 'in', res = 300)
par(mfrow=c(2,6))
par(mar=c(0,0.5,0,0.5),oma = c(4, 0, 3, 0))

#fall 2013 no sand lance, no bird data collected, 1 whale, no spatial stats
cruise<-filter(data, cruise=='1')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
points(cruise[4,5], cruise[4,4], col="blue", pch=16,cex=1.2)
mtext("A) Fall 2013", side=3, line=-1.8, cex=1.3,adj=0.1)
mtext("No sand lance", side=3, line=-3.4, cex=1,adj=0)
mtext("No bird data", side=3, line=-4.9, cex=1,adj=0)

#spring 2014 #whales, 1 bird, no sand lance
cruise<-filter(data, cruise=='6')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="w_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="blue",lwd=2.5,add=TRUE)
points(cruise[3,5], cruise[3,4], col="gold", pch=16,cex=1.2)
mtext("B) Spring 2014", side=3, line=-1.8, cex=1.3,adj=0.1)
mtext("No sand lance", side=3, line=-3.4, cex=1,adj=0)

#fall 2014
cruise<-filter(data, cruise=='2')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=3,add=TRUE)
SI.cgi(db.data, name="w_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="blue",lwd=3,add=TRUE)
SI.cgi(db.data, name="b_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="gold",lwd=3,add=TRUE)
mtext("C) Fall 2014", side=3, line=-1.8, cex=1.3,adj=0.1)

#spring 2015
cruise<-filter(data, cruise=='7')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
cruise <- cruise[-c(6),] #removing duplicated site
cruise <- cruise[-c(21),]#removing duplicated site
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=2.5,add=TRUE)
SI.cgi(db.data, name="w_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="blue",lwd=2.5,add=TRUE)
mtext("D) Spring 2015", side=3, line=-1.8, cex=1.3,adj=0.1)
mtext("No birds", side=3, line=-3.4, cex=1,adj=0)

#Fall 2015
cruise<-filter(data, cruise=='3')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=2.5,add=TRUE)
SI.cgi(db.data, name="w_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="blue",lwd=2.5,add=TRUE)
SI.cgi(db.data, name="b_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="gold",lwd=2.5,add=TRUE)
mtext("E) Fall 2015", side=3, line=-1.8, cex=1.3,adj=0.1)

#Spring 2016
cruise<-filter(data, cruise=='8')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=2.5,add=TRUE)
mtext("F) Spring 2016", side=3, line=-1.8, cex=1.3,adj=0.1)
mtext("No whales", side=3, line=-3.4, cex=1,adj=0)
mtext("No birds", side=3, line=-4.9, cex=1,adj=0)

#Fall 2016
cruise<-filter(data, cruise=='4')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=2.5,add=TRUE)
SI.cgi(db.data, name="w_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="blue",lwd=2.5,add=TRUE)
SI.cgi(db.data, name="b_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="gold",lwd=2.5,add=TRUE)
mtext("G) Fall 2016", side=3, line=-1.8, cex=1.3,adj=0.1)
      
#Spring 2017
cruise<-filter(data, cruise=='9')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=2.5,add=TRUE)
mtext("H) Spring 2017", side=3, line=-1.8, cex=1.3,adj=0.1)
mtext("No whales", side=3, line=-3.4, cex=1,adj=0)
mtext("No birds", side=3, line=-4.9, cex=1,adj=0)

#Fall 2017
cruise<-filter(data, cruise=='5')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=2.5,add=TRUE)
SI.cgi(db.data, name="b_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="gold",lwd=2.5,add=TRUE)
mtext("I) Fall 2017", side=3, line=-1.8, cex=1.3,adj=0.1)
mtext("No whales", side=3, line=-3.4, cex=1,adj=0)

#Spring 2018 
cruise<-filter(data, cruise=='10')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=2.5,add=TRUE)
SI.cgi(db.data, name="w_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="blue",lwd=2.5,add=TRUE)
mtext("J) Spring 2018", side=3, line=-1.8, cex=1.3,adj=0.1)
mtext("No birds", side=3, line=-3.4, cex=1,adj=0)

#Summer 2018 
cruise<-filter(data, cruise=='11')
cruise<-dplyr::select(cruise, sl,w_total,b_total,lat,lon)
db.data = db.create(cruise)
db.data = db.locate(db.data,"lon","x",1)
db.data = db.locate(db.data,"lat","x",2)
db.data = db.locate(db.data,"sl","z",1)
db.data = db.locate(db.data,"w_total","z",2)
db.data = db.locate(db.data,"b_total","z",3)

plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
plot(sbnms, add=T)
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=3,add=TRUE)
SI.cgi(db.data, name="b_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="gold",lwd=2.5,add=TRUE)
SI.cgi(db.data, name="w_total",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="blue",lwd=2.5,add=TRUE)
mtext("K) Summer 2018", side=3, line=-1.8, cex=1.3,adj=0.1)

#legend
plot(sub_mask, col="white",axes=FALSE,box=FALSE, legend=FALSE) #plot a blank plot
#plot.new() #plot an empty plot
#legend(x=-70.65,y=42.74,bg="transparent",pt.cex=4,bty="n",y.intersp=4, fill = c("red","blue","gold"), legend = c("", "",""),xpd=TRUE)
legend(x=-70.6,y=42.74,bg="transparent",pt.cex=2,bty="n",y.intersp=4,pch=c(16,16,16),col=c("red","blue","gold"), legend = c("", "",""),xpd=TRUE)#legend with circles and not squares

sandeel <- readPNG("./data/sand eel.png") #read in animal pngs
mn <- readPNG("./data/Mn.png")
shear <- readPNG("./data/shearwater.png")

rasterImage(image=sandeel, xleft=-70.5,ybottom=42.58,xright=-70.1,ytop=42.71)
rasterImage(image=mn, xleft=-70.55,ybottom=42.34,xright=-70.08,ytop=42.58)
rasterImage(image=shear, xleft=-70.5,ybottom=42.18,xright=-70.25,ytop=42.38)
dev.off()

###Supplementary figures

### Figures S1 - S3 - animal counts by site by cruise, change symbol sizes as appropriate (for sand lance, whales, or birds)

##assign symbol sizes (size points according to counts)
#assign symbol sizes for sl, birds and whales based on bird numbers (largest spread and most variation) for plotting
#I just chose these symbol sizes manually by eye for now-probably good enough

#bird symbol sizes
data <- mutate(data, b_sym_size = ifelse(b_total==0, 0.7, 
                                         ifelse(b_total>0 & b_total < 3,1.1, 
                                                ifelse(b_total>2 & b_total < 6,1.4,
                                                       ifelse(b_total>5 & b_total < 11,1.8, 
                                                              ifelse(b_total>10 & b_total < 31,2.3,
                                                                     ifelse(b_total>30 & b_total < 51,2.9,
                                                                            ifelse(b_total>50 & b_total < 101,3.9, 5.1))))))))
#whale symbol sizes
data <- mutate(data, w_sym_size = ifelse(w_total==0, 0.7, 
                                         ifelse(w_total>0 & w_total < 3,1.1, 
                                                ifelse(w_total>2 & w_total < 6,1.4,
                                                       ifelse(w_total>5 & w_total < 11,1.8, 
                                                              ifelse(w_total>10 & w_total < 31,2.3,
                                                                     ifelse(w_total>30 & w_total < 51,2.9,
                                                                            ifelse(w_total>50 & w_total < 101,3.9, 5.1))))))))
#sand lance symbol sizes
data <- mutate(data, sl_sym_size = ifelse(sl==0, 0.7, 
                                          ifelse(sl>0 & sl < 3,1.1, 
                                                 ifelse(sl>2 & sl < 6,1.4,
                                                        ifelse(sl>5 & sl < 11,1.8, 
                                                               ifelse(sl>10 & sl < 31,2.3,
                                                                      ifelse(sl>30 & sl < 51,2.9,
                                                                             ifelse(sl>50 & sl < 101,3.9, 5.1))))))))

tiff( "./figures/Collocation_manuscript/revisions_conscipractice/b_totals.tif", width = 10, height = 8, units = 'in', res = 300)

par(mfrow=c(2,6))
par(mar=c(0,0,4,0),oma = c(0, 0, 0, 0))

#Fall 2013
cruise<-filter(data, cruise=='1')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black") #distinguish sites not visited, or no valid sample
mtext("A) Fall 2013", side=3, line=0,adj=0,padj=1.8)

#Spring 2014
cruise<-filter(data, cruise=='6')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("B) Spring 2014", side=3, line=0,adj=0,padj=1.8)

#Fall 2014
cruise<-filter(data, cruise=='2')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("C) Fall 2014", side=3, line=0,adj=0,padj=1.8)

#Spring 2015
cruise<-filter(data, cruise=='7')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("D) Spring 2015", side=3, line=0,adj=0,padj=1.8)

#Fall 2015
cruise<-filter(data, cruise=='3')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("E) Fall 2015", side=3, line=0,adj=0,padj=1.8)

#Spring 2016
cruise<-filter(data, cruise=='8')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("F) Spring 2016", side=3, line=0,adj=0,padj=1.8)

#Fall 2016
cruise<-filter(data, cruise=='4')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("G) Fall 2016", side=3, line=0,adj=0,padj=1.8)

#Spring 2017
cruise<-filter(data, cruise=='9')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("H) Spring 2017", side=3, line=0,adj=0,padj=1.8)

#Fall 2017
cruise<-filter(data, cruise=='5')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("I) Fall 2017", side=3, line=0,adj=0,padj=1.8)

#Spring 2018
cruise<-filter(data, cruise=='10')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("J) Spring 2018", side=3, line=0,adj=0,padj=1.8)

#Summer 2018
cruise<-filter(data, cruise=='11')
cruise$site<-as.factor(cruise$site)
plot(sub_mask,breaks=c(75,165,255),col=c("gray78","gray54","gray27"),axes=FALSE, box=FALSE,legend=FALSE)
plot(sbnms, add=T) 
points(cruise$lon, cruise$lat,cex=cruise$b_sym_size, pch=16,col=mycol)
not_samp_sites<-anti_join(sites, cruise, by="site")
points(not_samp_sites$lon, not_samp_sites$lat, cex=0.7, pch=4,lwd=1.2, col="black")
mtext("K) Summer 2018", side=3, line=0,adj=0,padj=1.8)
#legend stuff
pointsize <- c(0,"1-2","3-5","6-10","11-30","31-50","51-100","101-200")
#legend(-70.145, 42.78,legend=pointsize, pt.cex =sort(unique(data$w_sym_size)),y.intersp=1.3,  col= mycol,pch=16)
plot.new()
legend("topright",legend=pointsize, pt.cex =sort(unique(data$b_sym_size)),y.intersp=1.6,  col= mycol,pch=16,bty="n")
dev.off()


###Figure S4 - combined center of gravity and inertia for all cruises combined

tiff( "./figures/Collocation_manuscript/revisions_conscipractice/totalGIC.tif",width = 7, height = 10, units = 'in', res = 300) #create a tiff device, sized and 300 dpi res
#plot(depth_mask, col=cols,xlab="Longitude", ylab="Latitude",cex.lab=1.3,cex.axis=1.3, axis.args=list(cex.axis=1.3),legend.args=list(cex=1.3,text='Depth (m)', side=3,adj=0.15))
plot(sub_mask,breaks=c(16,75,165,255),col=c("gray100","gray80","gray55"),axes=FALSE,box=FALSE, legend=FALSE)
#legend("topright",bg="transparent",pt.cex=1.3,bty="n",y.intersp=3, fill = c("#4DAF4A","#377EB8","#FF7F00"), legend = c("", "", ""))
plot(sbnms, add=T) 
SI.cgi(db.data, name="sl",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="red",lwd=2,add=TRUE)
SI.cgi(db.data, name="w",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="blue",lwd=2,add=TRUE)
SI.cgi(db.data, name="b",flag.plot=T,flag.inertia=T,flag.ellipse=T,col="gold",lwd=2,add=TRUE)

rasterImage(image=sandeel, xleft=-70.18,ybottom=42.72,xright=-70.09,ytop=42.75)
rasterImage(image=mn, xleft=-70.19,ybottom=42.66,xright=-70.09,ytop=42.7)
rasterImage(image=shear, xleft=-70.18,ybottom=42.60,xright=-70.1,ytop=42.65)
dev.off() #end plot
