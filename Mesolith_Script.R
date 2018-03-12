# Start R from within a GRASS session

setwd("F:/R-Code_Projects/Mesolithic-PPM")

## Load Required Packages ##

library(rgrass7)


# Importing raster files (covariates) and region from GRASS

cadore<-readRAST(c("elevation","erosion_class","lithology_reclass","distpaths","slope",
                   "landuse_class","insolation","aspect_cat"),cat=c(F,F,F,F,F,F,F,T),plugin=F)
str(cadore)

# Importing sites from GRASS 

prj_area<-readVECT("region_def",plugin=F)

sites_vect<-readVECT("siti",plugin=F)
head(data.frame(sites_vect))

# Save GRASS files loaded into R

save(cadore,prj_area,sites_vect,file="Meso.RData")


## Load required packages

library(spatstat)
library(MASS)
library(maptools)
library(MuMIn)


# Assessing collinearity between covariates at site locations

sites_var<-data.frame(sites_vect)[,5:12]
corr_sites<-cor(sites_var)
# No significant collinearity between the selected covariates


# Setting the working region

region<-as(as(prj_area,"SpatialPolygons"),"owin")


# Convert covariates to objects of class image and create list

covar<-list(El=as.im(cadore["elevation"]),Er=as.im(cadore["erosion_class"]),Li=as.im(cadore["lithology_reclass"]),
            Pa=as.im(cadore["distpaths"]),Sl=as.im(cadore["slope"]),Lu=as.im(cadore["landuse_class"]),
            In=as.im(cadore["insolation"]),As=as.im(cadore["aspect_cat"]))


# Converting site locations to point pattern process

sites<-data.frame(x=coordinates(sites_vect)[,1],y=coordinates(sites_vect)[,2])
sites_ppp<-ppp(sites[,1],sites[,2],window=region)


# Create null point pattern model

model_null<-ppm(sites_ppp,~1)


# Create model 1: spatial pattern depends on archaeological biases

model_1<-ppm(sites_ppp,~Er+Lu+Li,covariates=covar)


# Create model 2: spatial pattern depends on Mesolithic settlement strategy 

model_2 <- ppm(sites_ppp, ~El+In+Sl+As,covariates=covar)


# Create model 3 (all the covariates) and model 4 (all the covariates + dist. from paths)

model_3<-ppm(sites_ppp, ~Er+Lu+Li+As+El+In+Sl,covariates=covar)
model_4<-ppm(sites_ppp, ~Er+Lu+Li+Pa+As+El+In+Sl,covariates=covar)


# Create model list

models<-list(model0=model_null,model1=model_1,model2=model_2,model3=model_3,model4=model_4)


# Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC)
# stepwise variable selection

AIC_models<-lapply(models,stepAIC)
# Model 1 keeps Er & Lu
# Model 2 keeps El, In & Sl
# Model 3 keeps Er, Lu, El & In & Sl
# Model 4 keeps Er, Lu, Li, Pa, El & Sl

n_sites<-length(unique(sites_vect$ID))
BIC_models<-lapply(models,stepAIC,k=log(n_sites))
# Model 1 keeps Er & Lu
# Model 2 keeps El, In & Sl
# Model 3 keeps Er, Lu & Sl
# Model 4 keeps Lu, Pa, El & Sl


# Compare AIC scores for the different models

AIC_scores<-lapply(AIC_models,AICc)
BIC_scores<-lapply(BIC_models,BIC)


# AIC weight for diferent models

AIC.BIC.weight<-function(x){
  for(i in 1:length(x)){
    x.vect<-as.numeric(c(x[1:i]))}
  delta<-x.vect-min(x.vect)
  L<-exp(-0.5*delta)
  result<-round(L/sum(L),digits=7)
  return(result)
  }

write.table(AIC.BIC.weight(AIC_scores),file="AIC_Models_weights.txt")
write.table(AIC.BIC.weight(BIC_scores),file="BIC_Models_weights.txt")


# Export model summary

write.table(capture.output(print(AIC_models)),file="AIC_Models.txt")
write.table(capture.output(print(BIC_models)),file="BIC_Models.txt")


# Plot images

heat.colors(3)

tiff("Model1.tif",width=1440,height=480)
par(mar=c(1.5,1.5,1.5,1.5))
par(mfrow=c(1,3))
par(cex=1.5)
plot.im(covar$Er,main="erosion",col=c("#FFFF00FF","#FF0000FF"),ribsep=0.05)
plot.im(covar$Lu,main="landuse",col=c("#FFFF00FF","#FF8000FF","#FF0000FF"),ribsep=0.05)
plot.im(covar$Li,main="lithology",col=c("#FFFF00FF","#FF8000FF","#FF0000FF"),ribsep=0.05)
dev.off()

tiff("Model2.tif",width=960,height=960)
par(mar=c(1.5,1.5,1.5,1.5))
par(mfrow=c(2,2))
par(cex=1.5)
plot.im(covar$El,main="elevation",col=rev(heat.colors(250)),ribsep=0.05)
plot(covar$Sl,main="slope",col=rev(heat.colors(90)),ribsep=0.05)
plot(covar$In,main="insolation time",col=rev(heat.colors(15)),ribsep=0.05)
plot(covar$As,main="aspect",col=heat.colors(4),ribsep=0.05)
dev.off()

tiff("Additional_Covar.tif",width=480,height=480)
par(mar=c(1.5,1.5,1.5,1.5))
par(mfrow=c(1,1))
par(cex=1.5)
plot.im(covar$Pa,main="distance from paths (e-03)",col=heat.colors(250),ribscale=0.001,ribsep=0.05)
dev.off()


# Additional figures

library(GISTools)

tiff("Map.tif",width=960,height=960)
par(mar=c(0,0,0,0))
image(cadore,"elevation",col=terrain.colors(10))
plot(sites_vect,pch=20,col="black",cex=1.5,add=T)
map.scale(1718725,5138140,len=10000,units="Kilometer",ndivs=4,subdiv=2.5)
north.arrow(1718725,5140150,len=1000,lab="North",col="grey")
dev.off()