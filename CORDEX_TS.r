#Libraries
library(ncdf4) #CDF Files
library(xts) #Time Series
library(ggplot2) #Plotting
library(reshape2) #Manipulation
library(gridExtra) #Plot Arrangement
library(sp) #Geographical Distance
library(geosphere)

#Create an Array with the Files
directory <- "C:/Users/malva/Documents/1 - Flood Risk Management/Semester 3 - UPC/3 - Drought/Assignment/Data/Temperature"
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
timeseriesfile <- ncvar_get(ncin, varid = 'ts', start = c(ij[1],ij[2],1), count = c(1,1,-1))

#Create XTS Timeseries
timeseries_ts <- xts(timeseriesfile,obsdatesfile)

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
	timeseriesfile <- ncvar_get(ncin, varid = 'ts', start = c(ij[1],ij[2],1), count = c(1,1,-1))

	# Append the TS of the First Iteration with all the Files
	timeseries_tsfile <- xts(timeseriesfile,obsdatesfile)
	timeseries_ts <- rbind(timeseries_ts,timeseries_tsfile)
	nc_close(ncin)
	rm(ncin)
}

#Convert TS from K to C 
timeseries_C <- timeseries_ts - 273.15

#Create Monthly Values
monthlyT <- apply.monthly(timeseries_C,mean)
 
#Create Yearly Values
yearlyT <- apply.yearly(timeseries_C,mean)

#Statistical Values for each Period
summary(yearlyT["1951-01-31/1979-12-31"])
summary(yearlyT["1970-01-01/1999-12-31"])
summary(yearlyT["1990-01-01/2019-12-31"])
summary(yearlyT["2020-01-01/2049-12-31"])
summary(yearlyT["2040-01-01/2069-12-31"])
summary(yearlyT["2060-01-01/2099-12-31"])

#Convert TS to Data Frame
yearlyT_df<-fortify(yearlyT)
colnames(yearlyT_df)[2]="Temperature"
colnames(yearlyT_df)[1]="Date"
monthlyT_df <-fortify(monthlyT)
colnames(monthlyT_df)[2]="Temperature"
colnames(monthlyT_df)[1]="Date"

#Create new column with the number of month for monthlyT_df
monthlyT_df$Month<-as.numeric(format(monthlyT_df$Date,'%m'))
#Create new column with the number of year for yearlyT_df
yearlyT_df$Year<-as.numeric(format(yearlyT_df$Date,"%Y"))

#Plot yearly values for the analyzed periods
hlay<-cbind(c(1,1,2,2),c(NA,3,3,NA))
hlay<-t(hlay)

#Historical Yearly Mean Plots
P1951 <- ggplot(subset(yearlyT_df,Year<1980),aes(x=Year,y=Temperature)) + geom_line(size=1.2,color="steelblue4") + theme_bw() + labs(title="Annual Temperature 1951-1980", y="Temperature °C", x="Year") + ylim(19,23) + geom_hline(yintercept = mean(yearlyT_df$Temperature[yearlyT_df$Year<1980]),color="red") + theme(plot.title = element_text(size=11,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9))

P1970 <- ggplot(subset(yearlyT_df,Year<2000&Year>=1970),aes(x=Year,y=Temperature)) + geom_line(size=1.2,color="steelblue4") + theme_bw() + labs(title="Annual Temperature 1970-2000", y="Temperature °C", x="Year") + ylim(19,23) + geom_hline(yintercept = mean(yearlyT_df$Temperature[yearlyT_df$Year<2000&yearlyT_df$Year>=1970]),color="red") + theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9))

P1990 <- ggplot(subset(yearlyT_df,Year<2020&Year>=1990),aes(x=Year,y=Temperature)) + geom_line(size=1.2,color="steelblue4") + theme_bw() + labs(title="Annual Temperature 1990-2020", y="Temperature °C", x="Year") + ylim(19,23) + geom_hline(yintercept = mean(yearlyT_df$Temperature[yearlyT_df$Year<2020&yearlyT_df$Year>=1990]),color="red") + theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9))
#Plot Historical Period
grid.arrange(p51,p70,p90,layout_matrix=hlay)

#Future Yearly Mean Plots
P2020 <- ggplot(subset(yearlyT_df,Year<2050&Year>=2020),aes(x=Year,y=Temperature)) + geom_line(size=1.2,color="tomato4") + theme_bw() + labs(title="Annual Temperature 2020-2050", y="Temperature °C", x="Year") + ylim(19,25) + geom_hline(yintercept = mean(yearlyT_df$Temperature[yearlyT_df$Year<2050&yearlyT_df$Year>=2020]),color="red") + theme(plot.title = element_text(size=12, hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9))

P2040 <- ggplot(subset(yearlyT_df,Year<2070&Year>=2040),aes(x=Year,y=Temperature)) + geom_line(size=1.2,color="tomato4") + theme_bw() + labs(title="Annual Temperature 2040-2070", y="Temperature °C", x="Year") + ylim(19,25) + geom_hline(yintercept = mean(yearlyT_df$Temperature[yearlyT_df$Year<2070&yearlyT_df$Year>=2040]),color="red") + theme(plot.title = element_text(size=12, hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9))

P2060 <- ggplot(subset(yearlyT_df,Year<2100&Year>=2060),aes(x=Year,y=Temperature)) + geom_line(size=1.2,color="tomato4") + theme_bw() + labs(title="Annual Temperature 2060-2100", y="Temperature °C", x="Year") + ylim(19,25) + geom_hline(yintercept = mean(yearlyT_df$Temperature[yearlyT_df$Year<2100&yearlyT_df$Year>=2060]),color="red") + theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9))
#Plot Future
grid.arrange(p2020,p2040,p2060,layout_matrix=hlay)

#Create an Array with mean temperature values for each month
T51_80 <- as.numeric(with(subset(monthlyT_df,Date<="1979-12-31"&Date>="1951-01-31"),by(Temperature,Month,mean)))
T70_00 <- as.numeric(with(subset(monthlyT_df,Date<="1999-12-31"&Date>="1970-01-01"),by(Temperature,Month,mean)))
T90_20 <- as.numeric(with(subset(monthlyT_df,Date<="2019-12-31"&Date>="1990-01-01"),by(Temperature,Month,mean)))
T20_50 <- as.numeric(with(subset(monthlyT_df,Date<="2049-12-31"&Date>="2020-01-01"),by(Temperature,Month,mean)))
T40_70 <- as.numeric(with(subset(monthlyT_df,Date<="2069-12-31"&Date>="2040-01-31"),by(Temperature,Month,mean)))
T60_100 <- as.numeric(with(subset(monthlyT_df,Date<="2099-12-31"&Date>="2060-01-31"),by(Temperature,Month,mean)))
Months <- c(1,2,3,4,5,6,7,8,9,10,11,12)

#Create data frame with monthly averages
MonthlyT_AVG_df <- data.frame(Months,T51_80,T70_00,T90_20,T20_50,T40_70,T60_100)

#Name Columns
colnames(MonthlyT_AVG_df)[2]="1951-1980"
colnames(MonthlyT_AVG_df)[3]="1970-2000"
colnames(MonthlyT_AVG_df)[4]="1990-2020"
colnames(MonthlyT_AVG_df)[5]="2020-2050"
colnames(MonthlyT_AVG_df)[6]="2040-2070"
colnames(MonthlyT_AVG_df)[7]="2070-2100"

#Prepare Data for ggplot
Month_AVG_Melt<-melt(MonthlyT_AVG_df,id.vars = 'Months')
colnames(Month_AVG_Melt)[1]="Month"
colnames(Month_AVG_Melt)[2]="Period"
colnames(Month_AVG_Melt)[3]="Temperature"

#Plot all monthly averages for periods of study
ggplot(Month_AVG_Melt,aes(x=Month,y=Temperature,colour=Period))+geom_line(size=1.2)+theme_bw()+labs(title="Monthly Temperature", y="Temperature °C", x="Month")+scale_x_continuous(name="Month",limits = c(1,12),breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))

#Regression Analysis
#Tendency for the Historical Period
lm(Temperature~Year,data = subset(yearlyT_df,Year<2006))
#Tendency for the Fistorical Period
lm(Temperature~Year,data = subset(yearlyT_df,Year>=2006))
#Plot trend line for both scenarios
ggplot(yearlyT_df,aes(x=Year,y=Temperature))+geom_line(color="indianred3",size=1.2)+geom_smooth(data = yearlyT_df[yearlyT_df$Year<2006,],method = lm,se=F)+geom_smooth(data = yearlyT_df[yearlyT_df$Year>=2006,],method = lm,se=F)+theme_bw()+labs(title="Annual Temperature", y="Temperature °C", x="Year")