library(data.table)
library(ggplot2)
library(patchwork)
library(broom)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(sf)
library(openxlsx)
library(car)
library(lubridate)
library(stats)
library(gt)

root = "C:\\Users\\cmeri\\OneDrive - Dartmouth College\\Research\\Runoff_Scaling"

#### Functions ####

scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "10^", l)
  # return this as an expression
  parse(text=l)
}

load_GRDC = function(number,file_list){
  prefix = ifelse(grepl('Reference',file_list,fixed=T),'Reference/','Data/')
  position = which(text_files == as.character(file_list[number]))
  file_name = file_list[position]
  file = fread(file_name)
  string = strsplit(file_list[position],'_')[[1]][3]
  site = sub(prefix,'',string)
  file = file[,site_number:=site]
}

PowerFit = function (data_table,cluster_list){
  cluster = data_table[cluster==cluster_list]
  lmfit = lm(data=cluster,formula = log10(MeanQ) ~ log10(area))
  lmSummary = glance(lmfit)
  coefficients1 = setDT(tidy(lmfit))
  coefficients = coefficients1[,c('cluster','k','c','StError_c','StError_k','p_value'):=.(cluster_list,coefficients1[1,10^(estimate)],coefficients1[2,estimate],coefficients1[2,std.error],coefficients1[1,std.error],coefficients1[2,p.value])]
  coefficients = coefficients[,k:=10^k]
  Statistics = unique(coefficients[,c('cluster','k','c','StError_c','StError_k','p_value')],by='cluster')
  #Statistics = Statistics[,Name:=Group_Names[cluster_list]]
}

#### Gauging Stations ####

# Load All Stations
HydroRegions = unique(fread(paste0(root,'\\Data\\Runoff_HydroRegions.csv'), colClasses = c("site_number" = "character")),by='site_number')
GRDC_Stations = setDT(read.xlsx(paste0(root,'\\Data\\grdc_climatesensitive_stations.xlsx')))
GRDC_Reference_Stations = setDT(read.xlsx(paste0(root,'\\Data\\grdc_reference_stations.xlsx')))

# Adjust units
GRDC_Stations = GRDC_Stations[,area:=area*1000000]
GRDC_Reference_Stations = GRDC_Reference_Stations[,area:=area*1000000]

# Rename for Consistency across tables
setnames(HydroRegions,'LATITUDE','lat')
setnames(HydroRegions, 'LONGITUDE','long')
setnames(HydroRegions, 'GEEArea','area')
setnames(GRDC_Stations,'grdc_no','site_number')
setnames(GRDC_Reference_Stations,'grdc_no','site_number')

# Remove US and CA sites not in HCDN or WSC Unregulated
GRDC_Reference_Stations = GRDC_Reference_Stations[!(country%in%c('US','CA')) | (startsWith(nat_id,'15'))]
GRDC_Stations = GRDC_Stations[!(country%in%c('US','CA')) | (startsWith(nat_id,'15'))]

# Remove Reference Stations in Countries with Climate Sensitve Stations and Duplicate Stations
GRDC_Reference_Stations = GRDC_Reference_Stations[!(site_number%in%GRDC_Stations$site_number)]
GRDC_Reference_Stations = GRDC_Reference_Stations[!(country%in%GRDC_Stations$country)]

# Remove Web Scraped Dam Results
write.csv(unique(GRDC_Reference_Stations[,'river'],by='river'),paste0(root,'\\Data\\Dams\\ReferenceStations_list.csv'))

Dam_Stations = fread(paste0(root,'\\Data\\Dams\\DammedStations_list.csv')) # Mono River Removed due to web scraping issues
Dam_Stations = merge(Dam_Stations,GRDC_Reference_Stations[,c('river','site_number')], by='river')

# Toggle Dam Stations
#GRDC_Reference_Stations = GRDC_Reference_Stations[site_number%in%(Dam_Stations[Dam=='No',site_number])]

# Create list of all Stations
HydroRegions = HydroRegions[,country:=ifelse(nchar(site_number)<8,'CA','US')]
Global_Sites = rbind(HydroRegions[,c('site_number','area','country','lat','long')],
                     GRDC_Stations[,c('site_number','area','country','lat','long')],
                     GRDC_Reference_Stations[,c('site_number','area','country','lat','long')])
Global_Sites = Global_Sites[area>0]
Global_Sites[,site_number:=paste('ID',site_number,sep='-')]
write.csv(Global_Sites,paste0(root,'\\Data\\Global_Sites.csv'))
Global_Sites = Global_Sites[,site_number:=gsub('ID-','',site_number)]

#### Discharge Data ####
# US and Canada
USCA_Discharge = readRDS(paste0(root,'\\Data\\DischargeData.RDS'))
setnames(USCA_Discharge,'LATITUDE','lat')
setnames(USCA_Discharge, 'LONGITUDE','long')

# GRDC Reference
text_files = list.files("C:\\Users\\cmeri\\OneDrive - Dartmouth College\\Research\\Runoff_Scaling\\Data\\GRDC_Reference", pattern='txt', full.names=T)
RefDischarge_Data = lapply(1:length(text_files), load_GRDC, file_list = text_files)
RefDischarge_Data = rbindlist(RefDischarge_Data)
RefDischarge_Data = RefDischarge_Data[Calculated>0]
RefDischarge_Data = RefDischarge_Data[,c('Month','year'):=.(month(as.Date(`YYYY-MM-DD`)),ifelse(month(as.Date(`YYYY-MM-DD`)) > 9, as.numeric(year(as.Date(`YYYY-MM-DD`)) + 1), as.numeric(year(as.Date(`YYYY-MM-DD`)))))]
setnames(RefDischarge_Data,'Calculated','mean_Qcms')
RefDischarge_Data = RefDischarge_Data[site_number%in%GRDC_Reference_Stations$site_number]

# GRDC Climate
text_files <- list.files("C:\\Users\\cmeri\\OneDrive - Dartmouth College\\Research\\Runoff_Scaling\\Data\\GRDC_Data", pattern='txt', full.names=T)
ClimDischarge_Data = lapply(1:length(text_files), load_GRDC, file_list = text_files)
ClimDischarge_Data = rbindlist(ClimDischarge_Data)
ClimDischarge_Data = ClimDischarge_Data[Calculated>0]
ClimDischarge_Data = ClimDischarge_Data[,c('Month','year'):=.(month(as.Date(`YYYY-MM-DD`)),ifelse(month(as.Date(`YYYY-MM-DD`)) > 9, as.numeric(year(as.Date(`YYYY-MM-DD`)) + 1), as.numeric(year(as.Date(`YYYY-MM-DD`)))))]
setnames(ClimDischarge_Data,'Calculated','mean_Qcms')
ClimDischarge_Data = ClimDischarge_Data[site_number%in%GRDC_Stations$site_number]

# Merge Discharge Data
AllDischarge = rbind(USCA_Discharge[,c('site_number','year','mean_Qcms')],ClimDischarge_Data[,c('site_number','year','mean_Qcms')],RefDischarge_Data[,c('site_number','year','mean_Qcms')])

# Prepare Discharge Data for clustering
GRDC_Discharge = rbind(ClimDischarge_Data[,c('site_number','year','Month','mean_Qcms')],RefDischarge_Data[,c('site_number','year','Month','mean_Qcms')])

GRDC_Discharge = GRDC_Discharge[,YearPeak:=max(mean_Qcms),by=c('site_number','year')]
GRDC_Discharge = GRDC_Discharge[,MonthPeak:=max(mean_Qcms),by=c('site_number','Month','year')]
GRDC_Discharge = unique(GRDC_Discharge,by=c('site_number','year','Month'))

GRDC_Discharge = GRDC_Discharge[,NormalizedQ:=MonthPeak/YearPeak]
GRDC_Discharge = GRDC_Discharge[,MeanNormQ:=mean(NormalizedQ,na.rm=T),by=c('site_number','Month')]
GRDC_Discharge = unique(GRDC_Discharge,by=c('site_number','Month'))

WideData = reshape(GRDC_Discharge[,c('site_number','Month','MeanNormQ')],idvar = 'site_number',timevar = 'Month',direction = 'wide')

setnames(WideData,'MeanNormQ.1','jan')
setnames(WideData,'MeanNormQ.2','feb')
setnames(WideData,'MeanNormQ.3','mar')
setnames(WideData,'MeanNormQ.4','apr')
setnames(WideData,'MeanNormQ.5','may')
setnames(WideData,'MeanNormQ.6','jun')
setnames(WideData,'MeanNormQ.7','jul')
setnames(WideData,'MeanNormQ.8','aug')
setnames(WideData,'MeanNormQ.9','sep')
setnames(WideData,'MeanNormQ.10','oct')
setnames(WideData,'MeanNormQ.11','nov')
setnames(WideData,'MeanNormQ.12','dec')

# Add Latitude and Longitude
WideData = merge(WideData,Global_Sites[,c('site_number','lat','long')],by='site_number')

#### Google Earth Engine Elevation & Precipitation ####

# Elevation
Elevation = fread(paste0(root,'\\Data\\Elevation\\Global_Elevations.csv'), colClasses = c('site_number'='character'))
setnames(Elevation,'be75','elevation',skip_absent = T)
setnames(Elevation,'bedrock','elevation',skip_absent = T) # For variable source dataset
Elevation = Elevation[,site_number:=gsub('ID-','',site_number)]

WideData = merge(WideData,Elevation[,c('site_number','elevation')])

# Precipitation
Precipitation = fread(paste0(root,'\\Data\\Gauge_Precipitation.csv'), colClasses = c('shedID' = 'character'))
setnames(Precipitation,'shedID','site_number')
Precipitation = Precipitation[,site_number:=gsub('ID-','',site_number)]
Precipitation = Precipitation[,Precip:=mean(annPcpt),by='site_number']
Precipitation = unique(Precipitation,by='site_number')

WideData = merge(WideData,Precipitation[,c('site_number','Precip')])

#### Clustering ####
# Load US and Canada data 
ClusteringData =  readRDS(paste0(root,'\\Data\\HydroRegionClusteringData.RDS'))
ClusteringData = ClusteringData[,c(1:14,17)]
setnames(ClusteringData,'LATITUDE','lat')
setnames(ClusteringData,'LONGITUDE','long')

# Add Site Elevation and Precipitation
ClusteringData = merge(ClusteringData, Precipitation[,c('site_number','Precip')],by='site_number')
ClusteringData = merge(ClusteringData, Elevation[,c('site_number','elevation')],by='site_number')

# Merge US & CA dataset with GRDC dataset
ClusteringData = rbind(ClusteringData,WideData)
ClusteringData = setorder(ClusteringData,'site_number')
ClusteringData = na.omit(ClusteringData)

saveRDS(ClusteringData,paste0(root,'\\Data\\Global_ClusteringData.RDS'))

# Record Sites used in analysis
setorder(ClusteringData, site_number)
site_number = ClusteringData$site_number

# Scale Clustering Data
Scaled_Clustering = scale(ClusteringData[,-1])

# kmeans clustering
Scaled_Clustering = data.frame(Scaled_Clustering)

rownames(Scaled_Clustering) <- site_number

Scaled_Clustering = na.omit(Scaled_Clustering)

#Adjust Geographic Importance
#Scaled_Clustering$lat = Scaled_Clustering$lat*1.5
#Scaled_Clustering$long = Scaled_Clustering$long*1.5

set.seed(0)
Kmeans = kmeans(Scaled_Clustering,25,40)
Kmeans = data.table(site_number=names(Kmeans$cluster), cluster=Kmeans$cluster)

# Add clusters to unscaled data
ClusteringData = merge(ClusteringData,Kmeans,by='site_number')

# Visualize Clusters
AllClusters_sf = st_as_sf(ClusteringData,coords = c('long','lat'))
AllClusters_sf = st_set_crs(AllClusters_sf,4326)
AllClusters_sf = st_transform(AllClusters_sf,crs=("+proj=robin +lon_0=0w"))

world = ne_countries(type = 'countries',scale = 'medium', returnclass = 'sf')
world = st_transform(world, crs=("+proj=robin +lon_0=0w"))

AllCluster_map = ggplot(data = NULL) +  
  geom_sf(data=world, fill='white',color = "gray60")+
  geom_sf(data=AllClusters_sf, aes(color=as.character(cluster)),size=1)+
  scale_color_manual(values = c("#5f86b7", "#4ed31b", "#c257d3", "#569244", "#eb04dc", "#90d796", "#631e62", "#bce333", "#5010c8", "#f6932e", "#1d333f", "#f4cacb", "#670501", "#6ad5eb", "#f82f65", "#1cf1a3", "#ad5f78", "#02531d", "#ffa8ff", "#4e4809", "#7d76e4", "#e2c627", "#ff2a0d", "#169294", "#a0581c"))+
  labs(x='Longitude',y='Latitude',color='Cluster')

#### Scaling ####
# Prepare discharge data
AllDischarge = AllDischarge[,MeanQ:=mean(mean_Qcms),by='site_number']
AllDischarge = merge(AllDischarge,Global_Sites[,c('site_number','country','area')],by='site_number')
AllDischarge = merge(AllDischarge, ClusteringData[,c('site_number','cluster')],by='site_number')

ScalingData = unique(AllDischarge[,c(1,4:7)],by='site_number')

# Run Scaling Function
Statistics = lapply(1:max(ScalingData$cluster),PowerFit,data_table=ScalingData)
Statistics = rbindlist(Statistics)

Scaling_facet= ggplot(unique(AllDischarge,by='site_number'),aes(area,mean_Qcms,color=country))+
  facet_wrap(vars(cluster),scales='free')+
  geom_point()+
  scale_x_log10(labels=scientific)+
  scale_y_log10(labels=scientific)+
  theme_bw()

c = ggplot(Statistics,aes(cluster,c))+
  geom_hline(yintercept=1, linetype='dashed')+
  geom_errorbar(aes(ymin=(c - StError_c), ymax=c + StError_c, width=0.25))+
  geom_point()+
  scale_x_continuous(breaks=c(1:25))+
  #scale_y_continuous(limits = c(0.25,1.4))+
  labs(x='Cluster',y='Scaling Parameter, c')+
  theme_bw()

AvePrecip = unique(ClusteringData[,meanPrecip:=mean(Precip),by='cluster'],by='cluster')
AvePrecip = merge(Statistics,AvePrecip[,c('cluster','meanPrecip')],by='cluster')
ggplot(AvePrecip,aes(meanPrecip,c))+
  #geom_hline(yintercept=1, linetype='dashed')+
  geom_errorbar(aes(ymin=(c - StError_c), ymax=c + StError_c, width=0.05))+
  geom_point()+
  scale_x_continuous(breaks=c(1:25))+
  #scale_y_continuous(limits = c(0.25,1.4))+
  labs(x='Cluster Average Precipitation',y='Scaling Parameter, c')+
  theme_bw()

print(mean(Statistics$c))

#### Format export ####
# Convert data.table to data frame
Stats_df <- as.data.frame(Statistics)

# Create a gt table from the data frame
Formatted_Table <- gt(data = Stats_df)

# Perform any table formatting or styling using the gt package functions
Formatted_Table <- Formatted_Table %>%
  tab_header(title = "Scaling Results") %>%
  tab_style(style = cell_text(weight = "bold"),locations = cells_title()) %>%
  opt_horizontal_padding(scale=3) %>%
  cols_label(StError_c = 'St.E. c', StError_k = 'St.E. k', p_value = 'p-value') %>%
  fmt_number(columns = c('k','c','StError_c','StError_k','p_value'),n_sigfig = 3) %>%
  fmt_scientific(columns = c('k','StError_c','StError_k','p_value'),decimals = 2)

# Print the table
print(Formatted_Table)
as_latex(Formatted_Table)
