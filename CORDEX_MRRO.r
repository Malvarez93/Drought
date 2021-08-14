#Libraries
library(ncdf4) #CDF Files
library(xts) #Time Series
library(ggplot2) #Plotting
library(reshape2) #Manipulation
library(gridExtra) #Plot Arrangement
library(sp) #Geographical Distance
library(geosphere)

#Create an Array with the Files
directory <- "C:/Users/malva/Documents/1 - Flood Risk Management/Semester 3 - UPC/3 - Drought/Assignment/Data/Runoff"
setwd(directory)
files <- list.files(directory,pattern = ".nc")

#Read Initial File
ncin <- nc_open(files[1], write=FALSE, readunlim=TRUE, verbose=TRUE, auto_GMT=TRUE, suppress_dimvals=FALSE)

#Determine Data Extent
lat <- ncvar_get(ncin,"lat")
lon <- ncvar_get(ncin,"lon")

#Size of the Data Mesh
mesh_size <- dim(lat)

#Set my City (Lat/Lon) <- Belén, Heredia, Costa Rica
lat_city = 9.989828
lon_city = -84.184075

#Create Mesh Points Coordinates
coords_matrix <- cbind(as.vector(lon),as.vector(lat))
mesh_points <- SpatialPoints(coords_matrix, proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
city_point <- SpatialPoints(matrix(c(lon_city, lat_city), nrow=1, ncol=2), proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)

#Plot Data
#plot(mesh_points, axes=TRUE, las=1)
#points(city_point, col='red', pch=20, cex=3) #Plot City Location (within the Data grid)

#ID of the closest point to city in the Mesh
closest_point <- which.min(distGeo(mesh_points, city_point))

#Contruct a Matrix with the point's ID
point_ID_matrix <- matrix(1:length(lat),nrow = nrow(lat),ncol = ncol(lat))

# Get index i, j of the City (Point)
ij <- which( point_ID_matrix == closest_point, arr.ind=T )

#First Iteration
t <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
obsdatesfile <- as.Date(t, origin = '1949-12-01') #Check time units

#Sample file at City location
timeseriesfile <- ncvar_get(ncin, varid = 'mrro', start = c(ij[1],ij[2],1), count = c(1,1,-1))

#Create XTS Timeseries
timeseries_mrro <- xts(timeseriesfile,obsdatesfile)

#Close and Remove file
nc_close(ncin)
rm(ncin)

#LOOP 

for (i in 2:length(files)){ 
	#Read File +1
	ncin<-nc_open(files[i])	
	
	# Get Time +1
	t <- ncvar_get(ncin,"time")
	
	# Get Time Units +1
	tunits <- ncatt_get(ncin,"time","units")
	
	# Get Dates +1
	obsdatesfile <- as.Date(t, origin = '1949-12-01')
	# Sample File +1
	timeseriesfile <- ncvar_get(ncin, varid = 'mrro', start = c(ij[1],ij[2],1), count = c(1,1,-1))

	# Append the TS of the First Iteration with all the Files
	timeseries_mrrofile <- xts(timeseriesfile,obsdatesfile)
	timeseries_mrro <- rbind(timeseries_mrro,timeseries_mrrofile)
	nc_close(ncin)
	rm(ncin)
}

#Convert TS of Precipitation from kg/(m2s) to mm/day 
timeseries_mmd <- timeseries_pr*(24*60*60)
 
#Create Yearly Values
yearlyMRRO<-apply.yearly(timeseries_mmd,sum)

#Statistical Values of the Precipitation
summary(yearlyMRRO["1951-01-31/1979-12-31"])
summary(yearlyMRRO["1970-01-01/1999-12-31"])
summary(yearlyMRRO["1990-01-01/2019-12-31"])
summary(yearlyMRRO["2020-01-01/2049-12-31"])
summary(yearlyMRRO["2040-01-01/2069-12-31"])
summary(yearlyMRRO["2060-01-01/2099-12-31"])

#Convert MMRO to DF
yearlyMRRO_df<-fortify(yearlyMRRO)
colnames(yearlyMRRO_df)[2]="Runoff"
colnames(yearlyMRRO_df)[1]="Date"

#Create new column with Year
yearlyMRRO_df$Year<-as.numeric(format(yearlyMRRO_df$Date,"%Y"))

#Create Water Budget
yearlyP_df <- read.csv("C:/Users/malva/Documents/1 - Flood Risk Management/Semester 3 - UPC/3 - Drought/Assignment/Data/Data_PP.csv")
WaterBudget_DF <- data.frame(yearlyMRRO_df$Year, yearlyMRRO_df$Runoff, yearlyP_df$Precipitation)
colnames(WaterBudget_DF)[1]="Year"
colnames(WaterBudget_DF)[2]="Runoff"
colnames(WaterBudget_DF)[3]="Precipitation"

#Regression Analysis
#Plot yearly values for the periods required
hlay<-cbind(c(1,1,2,2),c(NA,3,3,NA))
hlay<-t(hlay)
#Historical Period
P1951<-ggplot(subset(yearlyP_df,Year<1980),aes(x=Year,y=Precipitation))+geom_bar(stat='identity',size=0.5,fill="steelblue3",color="grey85")+theme_bw()+labs(title="Annual Precipitation 1951-1980", y="Precipitation (mm)", x="Year")+geom_hline(yintercept = mean(yearlyP_df$Precipitation[yearlyP_df$Year<1980]),color="red")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,9000)+scale_x_continuous(breaks=seq(1950,2100,by=2))
P1970<-ggplot(subset(yearlyP_df,Year<2000&Year>=1970),aes(x=Year,y=Precipitation))+geom_bar(stat='identity',size=0.5,fill="steelblue3",color="grey85")+theme_bw()+labs(title="Annual Precipitation 1970-2000", y="Precipitation (mm)",, x="Year")+geom_hline(yintercept = mean(yearlyP_df$Precipitation[yearlyP_df$Year<1980]),color="red")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,9000)+scale_x_continuous(breaks=seq(1950,2100,by=2))
P1990<-ggplot(subset(yearlyP_df,Year<2020&Year>=1990),aes(x=Year,y=Precipitation))+geom_bar(stat='identity',size=0.5,fill="steelblue3",color="grey85")+theme_bw()+labs(title="Annual Precipitation 1990-2020", y="Precipitation (mm)", x="Year")+geom_hline(yintercept = mean(yearlyP_df$Precipitation[yearlyP_df$Year<1980]),color="red")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,9000)+scale_x_continuous(breaks=seq(1950,2100,by=2))
#Plot Historical
grid.arrange(P1951,P1970,P1990,layout_matrix=hlay)

#Future Period
P2020<-ggplot(subset(yearlyP_df,Year<2050&Year>=2020),aes(x=Year,y=Precipitation))+geom_bar(stat='identity',size=0.5,fill="steelblue4",color="grey85")+theme_bw()+labs(title="Annual Precipitation 2020-2050", y="Precipitation mm/year", x="Year")+geom_hline(yintercept = mean(yearlyP_df$Precipitation[yearlyP_df$Year<1980]),color="red")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,9000)+scale_x_continuous(breaks=seq(1950,2100,by=2))
P2040<-ggplot(subset(yearlyP_df,Year<2070&Year>=2040),aes(x=Year,y=Precipitation))+geom_bar(stat='identity',size=0.5,fill="steelblue4",color="grey85")+theme_bw()+labs(title="Annual Precipitation 2040-2070", y="Precipitation mm/year", x="Year")+geom_hline(yintercept = mean(yearlyP_df$Precipitation[yearlyP_df$Year<1980]),color="red")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,9000)+scale_x_continuous(breaks=seq(1950,2100,by=2))
P2060<-ggplot(subset(yearlyP_df,Year<2100&Year>=2060),aes(x=Year,y=Precipitation))+geom_bar(stat='identity',size=0.5,fill="steelblue4",color="grey85")+theme_bw()+labs(title="Annual Precipitation 2070-2100", y="Precipitation mm/year", x="Year")+geom_hline(yintercept = mean(yearlyP_df$Precipitation[yearlyP_df$Year<1980]),color="red")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,9000)+scale_x_continuous(breaks=seq(1950,2100,by=2))
#Plot Future
grid.arrange(P2020,P2040,P2060,layout_matrix=hlay)