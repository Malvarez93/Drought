#Libraries
library(ncdf4) #CDF Files
library(xts) #Time Series
library(ggplot2) #Plotting
library(reshape2) #Manipulation
library(gridExtra) #Plot Arrangement
library(sp) #Geographical Distance
library(geosphere)

#Create an Array with the Files
directory <- "C:/Users/malva/Documents/1 - Flood Risk Management/Semester 3 - UPC/3 - Drought/Assignment/Data/Precipitation"
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
timeseriesfile <- ncvar_get(ncin, varid = 'pr', start = c(ij[1],ij[2],1), count = c(1,1,-1))

#Create XTS Timeseries
timeseries_pr <- xts(timeseriesfile,obsdatesfile)

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
	timeseriesfile <- ncvar_get(ncin, varid = 'pr', start = c(ij[1],ij[2],1), count = c(1,1,-1))

	# Append the TS of the First Iteration with all the Files
	timeseries_prfile <- xts(timeseriesfile,obsdatesfile)
	timeseries_pr <- rbind(timeseries_pr,timeseries_prfile)
	nc_close(ncin)
	rm(ncin)
}

#Convert TS of Precipitation from kg/(m2s) to mm/day 
timeseries_mmd <- timeseries_pr*(24*60*60)

#Create Monthly Values
monthlyP<-apply.monthly(timeseries_mmd,sum)
 
#Create Yearly Values
yearlyP<-apply.yearly(timeseries_mmd,sum)

#Statistical Values of the Precipitation
summary(yearlyP["1951-01-31/1979-12-31"])
summary(yearlyP["1970-01-01/1999-12-31"])
summary(yearlyP["1990-01-01/2019-12-31"])
summary(yearlyP["2020-01-01/2049-12-31"])
summary(yearlyP["2040-01-01/2069-12-31"])
summary(yearlyP["2060-01-01/2099-12-31"])

#Conver PR to DF
dailyP_df<-fortify(timeseries_mmd)
colnames(dailyP_df)[2]="Precipitation"
colnames(dailyP_df)[1]="Date"
monthlyP_df <-fortify(monthlyP)
colnames(monthlyP_df)[2]="Precipitation"
colnames(monthlyP_df)[1]="Date"
yearlyP_df<-fortify(yearlyP)
colnames(yearlyP_df)[2]="Precipitation"
colnames(yearlyP_df)[1]="Date"

#Create new column with the number of Month
monthlyP_df$Month<-as.numeric(format(monthlyP_df$Date,'%m'))
#Create new column with Year
yearlyP_df$Year<-as.numeric(format(yearlyP_df$Date,"%Y"))
#Create new column with year to the dailyP_df
dailyP_df$Year<-as.numeric(format(dailyP_df$Date,"%Y"))

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

#Create values with average values for each month
P51_80<-as.numeric(with(subset(monthlyP_df,Date<="1979-12-31"&Date>="1951-01-31"),by(Precipitation,Month,mean)))
P70_00<-as.numeric(with(subset(monthlyP_df,Date<="1999-12-31"&Date>="1970-01-01"),by(Precipitation,Month,mean)))
P90_20<-as.numeric(with(subset(monthlyP_df,Date<="2019-12-31"&Date>="1990-01-01"),by(Precipitation,Month,mean)))
P20_50<-as.numeric(with(subset(monthlyP_df,Date<="2049-12-31"&Date>="2020-01-01"),by(Precipitation,Month,mean)))
P40_70<-as.numeric(with(subset(monthlyP_df,Date<="2069-12-31"&Date>="2040-01-31"),by(Precipitation,Month,mean)))
P60_100<-as.numeric(with(subset(monthlyP_df,Date<="2099-12-31"&Date>="2060-01-31"),by(Precipitation,Month,mean)))
Months<-c(1,2,3,4,5,6,7,8,9,10,11,12)

#Create data frame with monthly averages
MonthlyP_AVG_DF<-data.frame(Months,P51_80,P70_00,P90_20,P20_50,P40_70,P60_100)

#Name variables
colnames(MonthlyP_AVG_DF)[2]="1951-1980"
colnames(MonthlyP_AVG_DF)[3]="1970-2000"
colnames(MonthlyP_AVG_DF)[4]="1990-2020"
colnames(MonthlyP_AVG_DF)[5]="2020-2050"
colnames(MonthlyP_AVG_DF)[6]="2040-2070"
colnames(MonthlyP_AVG_DF)[7]="2070-2100"

#Prepare data for ggplot
MonthlyP_AVG_Melt<-melt(MonthlyP_AVG_DF,id.vars = 'Months')
colnames(MonthlyP_AVG_Melt)[1]="Month"
colnames(MonthlyP_AVG_Melt)[2]="Period"
colnames(MonthlyP_AVG_Melt)[3]="Precipitation"

#plot all monthly averages for periods of study
ggplot(MonthlyP_AVG_Melt,aes(x=Month,y=Precipitation,colour=Period))+geom_line(size=1.2)+theme_bw()+labs(title="Monthly Precipitation", y="Average Total Precipitation (mm)", x="Month")+scale_x_continuous(name="Month",limits = c(1,12),breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))+scale_color_manual(values=c("grey80","grey50","grey20","royalblue","royalblue1","royalblue2"))

#Create daily annual maximum for all periods
MAXDaily_AN<-data.frame(seq(1951,2100),as.numeric(with(subset(dailyP_df),by(Precipitation,Year,max))))
colnames(MAXDaily_AN)[1]="Year"
colnames(MAXDaily_AN)[2]="Precipitation"

#Extreme Analysis (Daily Maximyms)
#Split the DF into the Historical and Future Periods
MAXDaily_AN50_05<-subset(MAXDaily_AN,Year<2006)
MAXDaily_AN06_100<-subset(MAXDaily_AN,Year>=2006)
#Create GEV to obtain tendency
EXT51_05<-fevd(MAXDaily_AN50_05$Precipitation,type = "GEV")
summary(EXT51_05)
plot(EXT51_05)
return.level.fevd(EXT51_05,return.period = c(2,5,10,25,50,100,500))
EXT06_100<-fevd(MAXDaily_AN06_100$Precipitation,type = "GEV")
summary(EXT06_100)
plot(EXT06_100)
return.level.fevd(EXT06_100,return.period = c(2,5,10,25,50,100,500))

#Extreme Analysis (Droughts)
#Load mam.r Script
source("C:/Users/malva/Documents/1 - Flood Risk Management/Semester 3 - UPC/3 - Drought/PRACTICAL CASE Project/R_scripts/MAM.r")
#Droughts 5 Years
AM5_P <- amfunction(yearlyP_df,5)
AM5_P <- AM5_P[-c(1),] #Deletes Erroneous Value
AM5_P_DF<-as.data.frame(AM5_P)
colnames(AM5_P_DF)[1]="Year"
colnames(AM5_P_DF)[2]="MAM"
#Divide Yearly Precipitation into Historical and Future Data Sets
YearlyP_Hist<-subset(yearlyP_df,Date<"2005-12-18")
YearlyP_Fut<-subset(yearlyP_df,Date>="2005-12-18")
#Calculate MAM 5 Years
AM5_Hist<-amfunction(YearlyP_Hist,5)
AM5_Hist<-AM5_Hist[-c(1),]
MAM5_Hist<-mean(AM5_Hist[,2])
AM5_Fut<-amfunction(YearlyP_Fut,5)
AM5_Fut<-AM5_Fut[-c(1),]
MAM5_Fut<-mean(AM5_Fut[,2])

#Plot MAM Each Period (Medina) + Mean Longer Period
#Historical
PMAM51<-ggplot(subset(AM5_P_DF,Year<1980),aes(x=Year,y=MAM))+geom_line(size=1.2,color="dodgerblue3")+theme_bw()+labs(title="MAM (5 yrs) 1951-1980", y="Annual P. (mm/year)", x="Year")+geom_hline(yintercept =MAM5_Hist,color="red")+theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,6500)+scale_x_continuous(breaks=seq(1950,2100,by=2))
PMAM70<-ggplot(subset(AM5_P_DF,Year<2000&Year>=1970),aes(x=Year,y=MAM))+geom_line(size=1.2,color="dodgerblue3")+theme_bw()+labs(title="MAM (5 yrs) 1970-2000", y="Annual P. (mm/year)", x="Year")+geom_hline(yintercept =MAM5_Hist,color="red")+theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,6500)+scale_x_continuous(breaks=seq(1950,2100,by=2))
PMAM90<-ggplot(subset(AM5_P_DF,Year<2020&Year>=1990),aes(x=Year,y=MAM))+geom_line(size=1.2,color="dodgerblue3")+theme_bw()+labs(title="MAM (5 yrs) 1990-2200", y="Annual P. (mm/year)", x="Year")+geom_hline(yintercept =MAM5_Hist,color="red")+theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,6500)+scale_x_continuous(breaks=seq(1950,2100,by=2))
#Plot Historical
grid.arrange(PMAM51,PMAM70,PMAM90,layout_matrix=hlay)

#Future
PMAM20<-ggplot(subset(AM5_P_DF,Year<2050&Year>=2020),aes(x=Year,y=MAM))+geom_line(size=1.2,color="dodgerblue4")+theme_bw()+labs(title="MAM (5 yrs) 2020-2050", y="Annual P. (mm/year)", x="Year")+geom_hline(yintercept =MAM5_Fut,color="red")+theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,6500)+scale_x_continuous(breaks=seq(1950,2100,by=2))
PMAM40<-ggplot(subset(AM5_P_DF,Year<2070&Year>=2040),aes(x=Year,y=MAM))+geom_line(size=1.2,color="dodgerblue4")+theme_bw()+labs(title="MAM (5 yrs) 2050-2070", y="Annual P. (mm/year)", x="Year")+geom_hline(yintercept =MAM5_Fut,color="red")+theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,6500)+scale_x_continuous(breaks=seq(1950,2100,by=2))
PMAM60<-ggplot(subset(AM5_P_DF,Year<2100&Year>=2060),aes(x=Year,y=MAM))+geom_line(size=1.2,color="dodgerblue4")+theme_bw()+labs(title="MAM (5 yrs) 2070-2100", y="Annual P. (mm/year)", x="Year")+geom_hline(yintercept =MAM5_Fut,color="red")+theme(plot.title = element_text(size=12,hjust = 0.5),axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),axis.text.x = element_text(angle=90, hjust=1))+ylim(0,6500)+scale_x_continuous(breaks=seq(1950,2100,by=2))
#Plot Historical
grid.arrange(PMAM20,PMAM40,PMAM60,layout_matrix=hlay)


#Season Analysis
#Classify as dry and wet season (according to the month)
monthlyP_df$Season<-ifelse(monthlyP_df$Month<=4|monthlyP_df$Month==12,"D","W")
monthlyP_df$Year<-as.numeric(format(monthlyP_df$Date,"%Y"))
#data frame with seasons
DRY<-subset(monthlyP_df,Season=="D")
WET<-subset(monthlyP_df,Season=="W")
Yearly_Dry<-aggregate(Precipitation~Season+Year,data =  DRY,sum)
Yearly_Wet<-aggregate(Precipitation~Season+Year,data =  WET,sum)
#Summary to extract mean values of dry and wet
summary(subset(Yearly_Dry,Yearly_Dry$Year<1979))
summary(subset(Yearly_Dry,Yearly_Dry$Year>=1970&Yearly_Dry$Year<2000))
summary(subset(Yearly_Dry,Yearly_Dry$Year>=1990&Yearly_Dry$Year<2020))
summary(subset(Yearly_Dry,Yearly_Dry$Year>=2020&Yearly_Dry$Year<2050))
summary(subset(Yearly_Dry,Yearly_Dry$Year>=2050&Yearly_Dry$Year<2070))
summary(subset(Yearly_Dry,Yearly_Dry$Year>=2060&Yearly_Dry$Year<2100))
summary(subset(Yearly_Wet,Yearly_Wet$Year<1979))
summary(subset(Yearly_Wet,Yearly_Wet$Year>=1970&Yearly_Wet$Year<2000))
summary(subset(Yearly_Wet,Yearly_Wet$Year>=1990&Yearly_Wet$Year<2020))
summary(subset(Yearly_Wet,Yearly_Wet$Year>=2020&Yearly_Wet$Year<2050))
summary(subset(Yearly_Wet,Yearly_Wet$Year>=2050&Yearly_Wet$Year<2070))
summary(subset(Yearly_Wet,Yearly_Wet$Year>=2060&Yearly_Wet$Year<2100))