library(data.table)
library(ggplot2)
library(dataRetrieval)
library(tidyhydat)
library(smwrBase)
library(e1071)
library(broom)
library(gstat)
library(car)
library(ggpubr)
library(gstat)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)
library(ggspatial)
library(gt)

root = "C:\\Users\\cmeri\\OneDrive - Dartmouth College\\Research\\Runoff_Scaling"

#### Functions ####
Calculate_Q = function(site_list,Qi,table) {
  RI = (1/Qi)
  Data = table[site_number==site_list]
  mean = mean(Data[,log(mean_Qcms)],na.rm=T)
  median = median(Data[,log(mean_Qcms)],na.rm = T)
  sd = sd(Data[,log(mean_Qcms)],na.rm = T)
  print(site_list)
  skew = skewness(Data[,log(mean_Qcms)],2,na.rm = T)
  Flow = qlpearsonIII((1-RI),mean,sd,skew)
  Q = data.table(site_number=as.character(site_list), Q=as.numeric(Flow))
  setnames(Q,'Q',paste0('Q',Qi))
}

PowerFit = function (cluster_list,data_table,Q_column){
  cluster = data_table[cluster==cluster_list]
  lmfit = lm(data=cluster,formula = log10(get(Q_column)) ~ log10(area))
  lmSummary = glance(lmfit)
  coefficients1 = setDT(tidy(lmfit))
  coefficients = coefficients1[,c('cluster','k','c','StError_c','StError_k','p_value'):=.(cluster_list,coefficients1[1,10^(estimate)],coefficients1[2,estimate],coefficients1[2,std.error],coefficients1[1,std.error],coefficients1[2,p.value])]
  Statistics = coefficients[,c('cluster','k','c','StError_c','StError_k','p_value')]
}

Apply_Q = function(Q_column) {
  Q_ = Q_column
  Statistics = unique(rbindlist(lapply(1:max(Annual$cluster),PowerFit,data_table=Annual,Q_column=Q_)),by='cluster')
  Statistics = Statistics[,Qi:=Q_]
}

std.error = function(x) sd(x)/length(x)

Calculate_c = function(group,data,Q) {
  df = data[Group==group & Qi==Q]
  grouplm = lm(data=df, c ~ AvRE)
  intercept = (coef(grouplm))[[1]]
  slope = (coef(grouplm))[[2]]
  Modeled = copy(Annual[Group==group])
  Modeled_c = Modeled[Group==group,Calc_c:=(slope*RunoffRatio)+intercept]
  return(Modeled_c)
}

Apply_Qc = function(Q) {
  Q_ = Q
  ModeledValues = unique(rbindlist(lapply(c('A','B','C'), Calculate_c, data=RegressionData, Q=Q_)),by='site_number')
  ModeledValues = ModeledValues[,Qi:=Q]
}

krige_Q = function(Q,data) {
  df = data[Qi==Q]
  sf_data = st_as_sf(df,coords = c('LONGITUDE','LATITUDE'),crs=4326)
  sf_data = st_transform(sf_data,'ESRI:102008')
  
  # Fit Variogram
  lzn.vgm <- variogram(Calc_c ~ 1, data = sf_data)
  lzn.fit <- fit.variogram(lzn.vgm, vgm(c("Exp", "Sph")), fit.kappa = TRUE)
  
  plot(lzn.vgm, lzn.fit)
  
  # Krige
  kriged <- krige(sf_data$Calc_c ~ 1, as_Spatial(sf_data$geometry), newdata=North_America_grid, model=lzn.fit)
  kriged <- st_as_sf(kriged)
  
  Kriged_Map = st_intersection(kriged,North_America)
  print(paste('Kriging completed for',Q,sep=' '))
  return(Kriged_Map)
}

#### Data Preparation ####
Group_Names = c('1' = 'Northern Parallel',
                '2' = 'Central Canada',
                '3' = 'Rocky Mountains',
                '4' = 'Appalachians',
                '5' = 'Northern Pacific Coast',
                '6' = 'Southern Plains',
                '7' = 'Northern Plains',
                '8' = 'Central East',
                '9' = 'Northern Atlantic',
                '10' = 'Southwest',
                '11' = 'Southern Pacific Coast',
                '12' = 'Western Canada',
                '13' = 'Rocky Lowland',
                '14' = 'Pacific Northwest',
                '15' = 'Great Lakes',
                '16' = 'Southeast')

Annual = fread(paste0(root,'\\Data\\Runoff_HydroRegions.csv'), colClasses = c("site_number" = "character"))
HydroSites = unique(Annual,by='site_number')
HydroSites = HydroSites[site_number!='08LF033']

HCDN_sites = HydroSites[nchar(site_number)>=8]
WSC_sites = HydroSites[nchar(site_number)<8] # 08LF033 has insufficient observations

# Read data
ScalingData = fread(paste0(root,'\\Data\\Runoff_HydroRegions.csv'), colClasses = c("site_number" = "character"))
setnames(ScalingData,'GEEArea','area')

# Annual Data
Annual = Annual[,c('site_number','year','mean_Qcms','LATITUDE','LONGITUDE','GEEArea','RunoffRatio','cluster')]
Annual = Annual[,mean_Qcms:=mean(mean_Qcms),by='site_number']
Annual = Annual[,RunoffRatio:=mean(RunoffRatio),by='site_number']
Annual = unique(Annual,by='site_number')
setnames(Annual, 'GEEArea', 'area')

#### Qi Series ####
Q2 = rbindlist(lapply(HydroSites$site_number, Calculate_Q, Qi=2, table=ScalingData))
Q5 = rbindlist(lapply(HydroSites$site_number, Calculate_Q, Qi=5, table=ScalingData))
Q10 = rbindlist(lapply(HydroSites$site_number, Calculate_Q, Qi=10, table=ScalingData))
Q50 = rbindlist(lapply(HydroSites$site_number, Calculate_Q, Qi=50, table=ScalingData))
Q100 = rbindlist(lapply(HydroSites$site_number, Calculate_Q, Qi=100, table=ScalingData))

# Merge Q values with Annual Data
Annual = merge(Annual,Q2,by='site_number')
Annual = merge(Annual,Q5,by='site_number')
Annual = merge(Annual,Q10,by='site_number')
Annual = merge(Annual,Q50,by='site_number')
Annual = merge(Annual,Q100,by='site_number')

#### Scaling ####
Statistics = rbindlist(lapply(c('mean_Qcms','Q2','Q5','Q10','Q50','Q100'), Apply_Q))

# Group
RegressionData = copy(Annual)
RegressionData = RegressionData[,c('AvArea','AvRE'):=.(mean(area),mean(RunoffRatio)),by='cluster']
RegressionData = RegressionData[,-c(1:7,9:13)]
RegressionData = unique(RegressionData,by='cluster')

RegressionData = merge(RegressionData, Statistics[,c('cluster','c','Qi')],by='cluster')
RegressionData = RegressionData[,Name:=Group_Names[cluster]]

RegressionData = RegressionData[Name%in%c('Southern Plains','Northern Parallel','Western Canada'),Group:='A']
RegressionData = RegressionData[Name%in%c('Southwest','Northern Plains','Central Canada','Appalachians'),Group:='B']
RegressionData = RegressionData[Name%in%c('Southeast','Central East','Southern Pacific Coast','Rocky Lowland','Rocky Mountains','Great Lakes','Northern Atlantic','Northern Pacific Coast','Pacific Northwest'),Group:='C']

# Add Uncertainties
Annual = Annual[,c('RE_se'):=.(std.error(RunoffRatio)),by='cluster']
RegressionData = merge(RegressionData,unique(Annual[,c('cluster','RE_se')],by='cluster'),by='cluster')
RegressionData = merge(RegressionData,Statistics[,c('cluster','StError_c','Qi')],by=c('cluster','Qi'))

ggplot(RegressionData[Qi=='Q100'],aes(AvRE,c))+
  stat_smooth(data=RegressionData[Group%in%c('A','B')],method='lm',aes(color=Group),se=F)+
  stat_cor(data=RegressionData[Group%in%c('A','B')],aes(color=Group,label=..rr.label..),label.x.npc = 0.75,label.y.npc = 0.3)+
  stat_regline_equation(data=RegressionData[Group%in%c('A','B')],aes(color=Group),label.x.npc = 0.75,label.y.npc = 0.1)+
  geom_hline(yintercept = 1,linetype='dashed')+
  geom_errorbar(aes(ymin=c-StError_c, ymax=c+StError_c),alpha=0.3)+
  geom_errorbarh(aes(xmin=AvRE-RE_se,xmax=AvRE+RE_se),alpha=0.3)+
  geom_point()+
  geom_text(aes(label=Name),vjust=-1)+
  guides(color='none')+
  theme_bw()

# Calculate c for each site using groupings
Annual = merge(Annual,unique(RegressionData[,c('cluster','Group')],by='cluster'),by='cluster')

ModeledValues = rbindlist(lapply(c('mean_Qcms','Q2','Q5','Q10','Q50','Q100'), Apply_Qc))
ModeledValues = ModeledValues[Group=='C',Calc_c:=1]

#### Krige ####
# Map Data
US = ne_countries(country = 'United States of America',type = 'countries',scale = 'medium', returnclass = 'sf')
Canada = ne_countries(country = 'Canada',type = 'countries',scale = 'medium', returnclass = 'sf')
North_America = rbind(US,Canada)  
North_America = st_set_crs(North_America,4326)
North_America = st_transform(North_America, 'ESRI:102008')
North_America_grid = st_make_grid(North_America,cellsize = 50000)

# Interpolation Data
Kriged_Map = lapply(c('mean_Qcms','Q2','Q5','Q10','Q50','Q100'), krige_Q, data=ModeledValues)

Mean_krige = rbindlist(Kriged_Map[1])
Q2_krige = rbindlist(Kriged_Map[2])
Q5_krige = rbindlist(Kriged_Map[3])
Q10_krige = rbindlist(Kriged_Map[4])
Q50_krige = rbindlist(Kriged_Map[5])
Q100_krige = rbindlist(Kriged_Map[6])

#### Maps ####
MeanMap = ggplot()+
  geom_sf(data=st_as_sf(Mean_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(55, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill='Kriged Scaling Parameter (c)',title='Mean Annual Discharge')

Q2Map = ggplot()+
  geom_sf(data=st_as_sf(Q2_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(55, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill='Kriged Scaling Parameter (c)',title='Q2')

Q5Map = ggplot()+
  geom_sf(data=st_as_sf(Q5_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(55, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill='Kriged Scaling Parameter (c)',title='Q5')

Q10Map = ggplot()+
  geom_sf(data=st_as_sf(Q10_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(55, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill='Kriged Scaling Parameter (c)',title='Q10')

Q50Map = ggplot()+
  geom_sf(data=st_as_sf(Q50_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(55, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill='Kriged Scaling Parameter (c)',title='Q50')

Q100Map = ggplot()+
  geom_sf(data=st_as_sf(Q100_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(55, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill='Kriged Scaling Parameter (c)',title='Q100')
