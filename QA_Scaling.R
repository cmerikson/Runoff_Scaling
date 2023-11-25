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
library(openxlsx)
library(stats)
library(qqplotr)
library(patchwork)
library(scales)

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
  coefficients = coefficients1[,c('cluster','k','c','StError_c','StError_k','p_value','r_squared'):=.(cluster_list,coefficients1[1,10^(estimate)],coefficients1[2,estimate],coefficients1[2,std.error],coefficients1[1,std.error],coefficients1[2,p.value],summary(lmfit)$r.squared)]
  Statistics = coefficients[,c('cluster','k','c','StError_c','StError_k','p_value','r_squared')]
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
  r2 = summary(grouplm)$r.squared
  Modeled = copy(Annual[Group==group])
  Modeled = Modeled[Group==group,c('Slope','Intercept','R_squared'):=.(slope,intercept,r2)]
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
  kriged <- krige(sf_data$Calc_c ~ 1, as_Spatial(sf_data$geometry), newdata=North_America_grid, model=lzn.fit, maxdist=300000)
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

# Summary Statistics
SiteStatistics = copy(ScalingData)
SiteStatistics = SiteStatistics[,Length:=.N,by='site_number']
SiteStatistics = SiteStatistics[,AvClusterLength:=mean(Length),by='cluster']
SiteStatistics = SiteStatistics[,SDClusterLength:=sd(Length),by='cluster']
SiteStatistics = SiteStatistics[,SiteCount:=length(unique(site_number)),by='cluster']

SiteStatistics = SiteStatistics[,cluster:=Group_Names[cluster]]
SiteStats_df <- as.data.frame(unique(SiteStatistics,by='cluster'))

Site_Table <- gt(data = SiteStats_df[,c('cluster','AvClusterLength','SDClusterLength','SiteCount')])

Site_Table <- Site_Table %>%
  tab_style(style = cell_text(weight = "bold"),locations = cells_title()) %>%
  opt_horizontal_padding(scale=3) %>%
  cols_label(cluster = 'Cluster', AvClusterLength = 'Mean Length of Record', SDClusterLength = 'Standard Deviation of Record Length',SiteCount= 'Number of Gauging Sites') %>%
  cols_align(align='center') %>%
  fmt_number(columns = c('AvClusterLength','SDClusterLength'),n_sigfig = 2) #%>%
  #fmt_scientific(columns = c('k','p_value'),decimals = 2) %>%

print(Site_Table)
as_latex(Site_Table)
  
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

ggplot(RegressionData[Qi=='mean_Qcms'],aes(AvRE,c))+
  stat_smooth(data=RegressionData[Group%in%c('A','B') & Qi=='mean_Qcms'],method='lm',aes(color=Group),se=F)+
  stat_cor(data=RegressionData[Group%in%c('A','B') & Qi=='mean_Qcms'],aes(color=Group,label=..rr.label..),label.x.npc = 0.75,label.y.npc = 0.3,size=6)+
  stat_regline_equation(data=RegressionData[Group%in%c('A','B') & Qi=='mean_Qcms'],aes(color=Group),label.x.npc = 0.75,label.y.npc = 0.2,size=6)+
  geom_hline(yintercept = 1,linetype='dashed')+
  geom_errorbar(aes(ymin=c-StError_c, ymax=c+StError_c),alpha=0.3)+
  geom_errorbarh(aes(xmin=AvRE-RE_se,xmax=AvRE+RE_se),alpha=0.3)+
  geom_point()+
  geom_text(aes(label=Name),vjust=-1,size=5)+
  scale_color_manual(values=c('blue','red'))+
  #guides(color='none')+
  labs(x='Average Runoff Efficiency', y=expression(paste('Scaling Exponent, ', italic('(c)'),sep = ' ')))+
  theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=16),legend.text = element_text(size=16),legend.title = element_text(size=18))

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
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(2, 'cm'), width = unit(2,'cm'),
                         pad_x = unit(42, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill=expression(paste('Scaling Exponent, ', italic('(c)'),sep = ' ')))+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16),legend.text = element_text(size=16),legend.title = element_text(size=18))

Q2Map = ggplot()+
  geom_sf(data=st_as_sf(Q2_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(2, 'cm'), width = unit(2,'cm'),
                         pad_x = unit(42, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill=expression(paste('Scaling Exponent, ', italic('(c)'),sep = ' ')))+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16),legend.text = element_text(size=16),legend.title = element_text(size=18))

Q5Map = ggplot()+
  geom_sf(data=st_as_sf(Q5_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(2, 'cm'), width = unit(2,'cm'),
                         pad_x = unit(42, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill=expression(paste('Scaling Exponent, ', italic('(c)'),sep = ' ')))+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16),legend.text = element_text(size=16),legend.title = element_text(size=18))

Q10Map = ggplot()+
  geom_sf(data=st_as_sf(Q10_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(2, 'cm'), width = unit(2,'cm'),
                         pad_x = unit(42, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill=expression(paste('Scaling Exponent, ', italic('(c)'),sep = ' ')))+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16),legend.text = element_text(size=16),legend.title = element_text(size=18))

Q50Map = ggplot()+
  geom_sf(data=st_as_sf(Q50_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(2, 'cm'), width = unit(2,'cm'),
                         pad_x = unit(42, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill=expression(paste('Scaling Exponent, ', italic('(c)'),sep = ' ')))+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16),legend.text = element_text(size=16),legend.title = element_text(size=18))

Q100Map = ggplot()+
  geom_sf(data=st_as_sf(Q100_krige),aes(fill=var1.pred),color=NA)+
  scale_fill_viridis(option='turbo', direction = -1,limits=c(0.45,1),oob=scales::squish)+
  coord_sf(xlim = c(-3000000, 3000000), ylim = c(-2000000,3000000), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(2, 'cm'), width = unit(2,'cm'),
                         pad_x = unit(42, "pt"), pad_y = unit(25, "pt"),
                         style = north_arrow_nautical()) +
  labs(x='Longitude',y='Latitude',fill=expression(paste('Scaling Exponent, ', italic('(c)'),sep = ' ')))+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16),legend.text = element_text(size=16),legend.title = element_text(size=18))

#### Tables ####
Statistics = Statistics[,cluster:=Group_Names[cluster]]
Stats_df <- as.data.frame(Statistics[Qi=='Q100'])

# Create a gt table from the data frame
Formatted_Table <- gt(data = Stats_df)

# Perform any table formatting or styling using the gt package functions
Formatted_Table <- Formatted_Table %>%
  tab_style(style = cell_text(weight = "bold"),locations = cells_title()) %>%
  opt_horizontal_padding(scale=3) %>%
  cols_label(StError_c = 'St.E. c', StError_k = 'St.E. k', p_value = 'p-value',r_squared=paste0('R','\U00B2') ,Qi = 'Flow Magnitude') %>%
  cols_align(align='center') %>%
  fmt_number(columns = c('k','c','StError_c','StError_k','p_value','r_squared'),n_sigfig = 2) %>%
  fmt_scientific(columns = c('k','p_value'),decimals = 2) %>%
  text_case_match("mean_Qcms" ~ "Mean Annual")

# Print the table
print(Formatted_Table)
as_latex(Formatted_Table)
gtsave(Formatted_Table,paste0(root,'\\Output\\ScalingTable.png'))

#Regression table
ModeledValues = ModeledValues[,cluster:=Group_Names[cluster]]
Mod_df <- as.data.frame(unique(ModeledValues[Group%in%c('A','B'),c('Group','Slope','Intercept','Qi','R_squared')]),by=c('Qi','Group'))

# Create a gt table from the data frame
gt_mod <- gt(data = Mod_df)

# Perform any table formatting or styling using the gt package functions
gt_mod <- gt_mod %>%
  tab_header(title = 'Table 1') %>%
  tab_style(style = cell_text(weight = "bold"),locations = cells_title()) %>%
  #opt_horizontal_padding(scale=1) %>%
  cols_label(R_squared = paste0('R','\U00B2'), Qi = 'Flow Magnitude') %>%
  cols_align(align='center') %>%
  fmt_number(columns = c('R_squared','Slope','Intercept'),n_sigfig = 2) %>%
  text_case_match("mean_Qcms" ~ "Mean Annual") %>%
  tab_footnote(html(paste('<b>Regression Equations are of the form c = (Slope * R\U2091) + Intercept.<b>
                      <br>
                      <b>Group A includes Northern Parallel, Southern Plains, and Western Canada.<b>
                      <br>
                      <b>Group B includes Appalachians, Central Canada, Northern Plains and Southwest.<b>')))

# Print the table
print(gt_mod)
as_latex(gt_mod)
gtsave(gt_mod,paste0(root,'\\Output\\CoeffTable.png'))

#### Bankfull Geometry ####

#Trampush data
ChannelSize = setDT(read.xlsx(paste0(root,'\\Data\\ChannelSize_Data.xlsx')))
setnames(ChannelSize,'Wbf.[m]','Width.(m)')

# State Codes
State_Codes = na.omit(setDT(copy(stateCd)))
State_Codes = State_Codes[!(STUSAB%in%c('AK','HI','GU','PR','MP','AS'))]

USGS_sites = data.table()
for (i in unique(State_Codes$STUSAB)) {
  USGS_data = whatNWISdata(stateCd=State_Codes[STUSAB==i,STUSAB],parameterCd="00060")
  USGS_sites = rbind(USGS_sites,USGS_data)
}
USGS_sites = USGS_sites[,c(1:6)]

write.csv(USGS_sites,paste0(root,'\\Data\\USGS_site_information.csv'))
USGS_sites = read.csv(paste0(root,'\\Data\\USGS_site_information.csv'))
setnames(USGS_sites,'station_nm','Site_Name')
setnames(USGS_sites,'dec_lat_va','Latitude')
setnames(USGS_sites,'dec_long_va','Longitude')

# Phillips data
Phillips = setDT(read.xlsx(paste0(root,'\\Data\\Phillips_et_al_2022.xlsx')))
Phillips = Phillips[!(is.na(`Drainage.area.(km^2)`))]

# Assign Location and Cluster
Phillips = merge(Phillips,USGS_sites[,c('Site_Name','Latitude','Longitude')],by='Site_Name')

Phillips_sf = st_as_sf(Phillips,coords = c('Latitude','Longitude'))
HydroSites_sf = st_as_sf(HydroSites,coords = c('LATITUDE','LONGITUDE'))

Phillips = Phillips[,cluster:=HydroSites[st_nearest_feature(Phillips_sf,HydroSites_sf),cluster]]

Phillips = Phillips[,`Scale.Type`:=ifelse(cluster%in%c(1,2,4,6,7,10,12),'NonLinear','Linear')]
setnames(Phillips,'Drainage.area.(km^2)','DrainageArea')
Phillips = Phillips[,cluster_name:=Group_Names[cluster]]

# Combine Datasets
Width = rbind(Phillips[,c('Site_Name','DrainageArea','Width.(m)','cluster_name','Scale.Type')],ChannelSize[,c('Site_Name','DrainageArea','Width.(m)','cluster_name','Scale.Type')])
Width = unique(Width,by=c('Site_Name','DrainageArea'))
Width = Width[,Count:=.N,by='Site_Name']

# Check Normality
ggplot(ChannelSize[Scale.Type=='Linear'], mapping = aes(sample=log(Bankfull_Area))) +
  stat_qq_band() +      # plots uncertainty band 
  stat_qq_line() +      # plots line 
  stat_qq_point()

ggplot(ChannelSize[Scale.Type=='Linear'])+
  geom_histogram(aes(x=log(Bankfull_Area)))+
  theme_bw()

breaks = 10^(-10:10)
minor_braks = rep(1:9,21)*(10^rep(-10:10,each=9))
ggplot(ChannelSize,aes(DrainageArea,Bankfull_Area,color=Scale.Type))+
  geom_point()+
  stat_smooth(method='lm',se=F)+
  stat_cor(aes(label=..rr.label..),label.x.npc = 0.75,label.y.npc = 0.35)+
  #stat_regline_equation(label.x.npc = 0.75,label.y.npc = 0.2)+
  scale_x_continuous(trans = 'log10',breaks=breaks,minor_breaks=minor_braks)+
  scale_y_continuous(trans = 'log10',breaks=breaks,minor_breaks=minor_braks)+
  scale_color_manual(values=c('cornflowerblue','darkred'))+
  annotation_logticks(color = 'gray20')+
  annotate('text',label=(expression(y==0.75*x^(0.65))),x=525,y=0.8,color='cornflowerblue')+
  annotate('text',label=(expression(y==0.34*x^(0.8))),x=500,y=0.5,color='darkred')+
  labs(x=bquote("Drainage Area "~(km^2)),y=bquote("Bankfull Channel Area "~(m^2)),color='Scaling Group')+
  theme_bw()

# t-test
Linear = lm(Bankfull_Area ~ DrainageArea,ChannelSize[Scale.Type=='Linear'])
NonLinear = lm(Bankfull_Area ~ DrainageArea,ChannelSize[Scale.Type=='NonLinear'])

# Extract slope and standard error for linear
slope_linear <- coef(Linear)[2]
se_slope_linear <- summary(Linear)$coefficients[2, "Std. Error"]

# Extract slope and standard error for group 2
slope_nonlinear <- coef(NonLinear)[2]
se_slope_nonlinear <- summary(NonLinear)$coefficients[2, "Std. Error"]

# Calculate the standard error of the difference between slopes
se_difference <- sqrt(se_slope_linear^2 + se_slope_nonlinear^2)

# Calculate the t-statistic for comparing slopes
t_stat <- (slope_linear - slope_nonlinear) / se_difference

# Degrees of freedom for the t-test
df <- min(nrow(ChannelSize[Scale.Type=='Linear']), nrow(ChannelSize[Scale.Type=='NonLinear'])) - 2

# Calculate the p-value
p_value <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
 
#### Scaling Plot ####
scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

Annual = Annual[,Name:=Group_Names[cluster]]

ggplot(Annual[Name=='Northern Atlantic' & area>10^7],aes(x=area,y=mean_Qcms))+
  stat_function(fun = function(x) Statistics[cluster==9 & Qi=='Q2', k] * (x)^Statistics[cluster==9 & Qi=='mean_Qcms', c],linewidth=2)+
  #stat_function(fun = function(x) (3*(10^-6)) * (x)^0.5,linetype='dashed')+
  #stat_function(fun = function(x) (3*(10^-3)) * (x)^0.5,linetype='dashed')+
  stat_function(fun = function(x) Statistics[cluster==4 & Qi=='Q2', k] * (x)^Statistics[cluster==4 & Qi=='mean_Qcms', c],color='#a7e831')+
  stat_function(fun = function(x) Statistics[cluster==12 & Qi=='Q2', k] * (x)^Statistics[cluster==12 & Qi=='mean_Qcms', c],color='#2cf52b')+
  stat_function(fun = function(x) Statistics[cluster==7 & Qi=='Q2', k] * (x)^Statistics[cluster==7 & Qi=='mean_Qcms', c],color='#0f1f5f')+
  stat_function(fun = function(x) Statistics[cluster==6 & Qi=='Q2', k] * (x)^Statistics[cluster==6 & Qi=='mean_Qcms', c],color='#d0cc36')+
  stat_function(fun = function(x) Statistics[cluster==1 & Qi=='Q2', k] * (x)^Statistics[cluster==1 & Qi=='mean_Qcms', c],color='#c79ae2')+
  stat_function(fun = function(x) Statistics[cluster==10 & Qi=='Q2', k] * (x)^Statistics[cluster==10 & Qi=='mean_Qcms', c],color='#6ceac0')+
  stat_function(fun = function(x) Statistics[cluster==2 & Qi=='Q2', k] * (x)^Statistics[cluster==2 & Qi=='mean_Qcms', c],color='#378811')+
  scale_x_log10(labels=scientific,limits=c(10^7,10^10))+
  scale_y_log10(labels=scientific)+
  annotate('text',label='Southwest',x=10^9,y=0.1,angle=17)+
  annotate('text',label='Northern Plains',x=10^8,y=0.4,angle=18)+
  annotate('text',label='Southern Plains',x=10^9,y=4,angle=17)+
  annotate('text',label='Central Canada',x=10^9.3,y=15,angle=20)+
  annotate('text',label='Northern Parallel',x=10^7.5,y=2,angle=18)+
  annotate('text',label='Western Canada',x=10^7.2,y=0.7,angle=23)+
  annotate('text',label='Appalachians',x=10^9.5,y=50,angle=23)+
  labs(x=bquote("Drainage Area "~(m^2)),y='Mean Annual Discharge (cms)')+
  theme_bw()

ggplot(Annual[Name%in%c('Northern Atlantic','Southwest')],aes(x=area/1000000,y=mean_Qcms))+
  geom_point(aes(color=Name))+
  stat_function(fun = function(x) setDT(tidy(lm(data=Annual[cluster==9],formula = log10((mean_Qcms)) ~ log10(area/1000000))))[1,10^estimate] * x^(setDT(tidy(lm(data=Annual[cluster==9],formula = log10((mean_Qcms)) ~ log10(area/1000000))))[2,estimate]))+
  stat_function(fun = function(x) setDT(tidy(lm(data=Annual[cluster==10],formula = log10((mean_Qcms)) ~ log10(area/1000000))))[1,10^estimate] * x^(setDT(tidy(lm(data=Annual[cluster==10],formula = log10((mean_Qcms)) ~ log10(area/1000000))))[2,estimate]))+
  scale_x_log10()+
  scale_y_log10(labels=scales::comma)+
  scale_color_manual(values=c("#5d1800", "#1ed598"))+
  labs(x=bquote("Drainage Area "~(km^2)),y='Mean Annual Discharge (cms)',color='Hydro-Region')+
  theme_bw()

# Runoff Efficiency
RE_data = readRDS("C:\\Users\\Cmeri\\OneDrive - Dartmouth College\\Research\\Runoff_Ratio\\Code_Exports\\Clustered_AllData.Rds")
RE_data = unique(RE_data[,AvRE:=mean(RunoffRatio),by='site_number'],by='site_number')
RE_data = RE_data[,Scaling:=ifelse(cluster%in%c(1,2,4,6,7,10,12),'NonLinear','Linear')]

ggplot(RE_data,aes(GEEArea,AvRE,color=Scaling))+
  #facet_wrap(vars(cluster),labeller = as_labeller(Group_Names))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_x_log10()+
  scale_y_continuous(limits = c(0,1))+
  #scale_color_manual(values=c("#c79ae2", "#378811", "#cb1775", "#a7e831", "#c6dbae", "#d0cc36", "#0f1f5f", "#5648d3", "#5d1800", "#6ceac0", "#f24219", "#2cf52b", "#b00bd9", "#0b4512", "#fa7ee3", "#6a7d54"),labels = c('1' = 'Northern Parallel','2' = 'Central Canada','3' = 'Rocky Mountains','4' = 'Appalachians','5' = 'Northern Pacific Coast','6' = 'Southern Plains','7' = 'Northern Plains','8' = 'Central East','9' = 'Northern Atlantic','10' = 'Southwest','11' = 'Southern Pacific Coast','12' = 'Western Canada','13' = 'Rocky Lowland','14' = 'Pacific Northwest','15' = 'Great Lakes','16' = 'Southeast'), breaks=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'))+
  theme_bw()

# Recurrence Interval Plot
Statistics = Statistics[Qi!='mean_Qcms',RI:=as.numeric(gsub('Q','',Qi))]
Statistics = Statistics[,c_slope:=coef(lm(c~RI))[2],by='cluster']

Centroids = setDT(read.xlsx(paste0(root,'\\Data\\ClusterCentroids_shift.xlsx')))

GrobSample = ggplot(Statistics[Qi!='mean_Qcms' & cluster%in%c(1,10)],aes(RI,c))+
  geom_point()+
  stat_smooth(method='lm',se=F,aes(color=ifelse(c_slope>0, 'Positive','Negative')))+
  scale_color_manual(values = c('Positive'='red','Negative'='blue'))+
  labs(x='Recurrence Interval',y='Scaling Exponent (c)',color='Magnitude')+
  theme_bw()

c_slope_map <- function(cluster_number){
  long <- Centroids[cluster == cluster_number]$Cent_lon
  lat <- Centroids[cluster == cluster_number]$Cent_lat
  plot_sel <- Statistics[cluster == cluster_number]
  return(annotation_custom(grob = ggplotGrob(
    ggplot(plot_sel,aes(RI,c))+
      geom_point()+
      stat_smooth(method='lm',se=F,aes(color=ifelse(c_slope>0, 'red','blue')))+
      scale_color_manual(values = c('red'='red','blue'='blue'))+
      labs(x='Recurrence Interval',y='Scaling Exponent (c)',title = Group_Names[cluster_number])+
      theme_minimal()+
      theme(legend.position = 'none',plot.title = element_text(size=8,hjust = 0.5,vjust = 0.5),
            axis.text.x = element_text(size=4), axis.text.y = element_blank(),
            axis.title = element_text(size=5))
  ),
  xmin = long - 4.5, xmax = long + 4.5,
  ymin = lat - 3.5, ymax = lat + 3.5))
}

extract_legend <- function(GrobObject) {
  step1 <- ggplot_gtable(ggplot_build(GrobObject))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

Map_c_slope = function(){
  bar_sel = lapply(c(1:16),c_slope_map)
  
  Map = ggplot(data = NULL) +  
    geom_sf(data=ne_countries(country = 'United States of America',type = 'countries',scale = 'medium', returnclass = 'sf'), fill='white',color = "gray60") +
    geom_sf(data=ne_countries(country = 'Canada',type = 'countries',scale = 'medium', returnclass = 'sf'), fill='white',color = "gray60") + 
    coord_sf(xlim = c(-142, -50), ylim = c(24.5, 60), expand = FALSE)+
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(55, "pt"), pad_y = unit(25, "pt"),
                           style = north_arrow_nautical()) +
    labs(x='Longitude',y='Latitude')
  
  Combined = Map + bar_sel + inset_element(extract_legend(GrobSample),left = 0.9,bottom = 0.1,right = 0.98,top = 0.4)
  
  return(Combined)
}

Map_c_slope()
