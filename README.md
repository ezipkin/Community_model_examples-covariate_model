# Community model examples
# Part II - Covariate model
Based on work from: *Zipkin E.F., Royle J.A., Dawson D.K., and Bates S. 2010. Multi-species occurrence models to evaluate the effects of conservation and management actions. Biological Conservation. 143: 479-484.*

## **Description:**
The hierarchical community model is a multi-species approach to obtain community information, such as species or assemblage richness, by estimating individual species occurrence probabilities. The fundamental idea behind the community approach is that collective information on all observed species can inform probabilities of detection and occurrence for both observed and unobserved species, even those that are rare or elusive. This results in an improved composite analysis of the community and increased precision in species specific estimates of occurrence. The hierarchical model can be specified to incorporate habitat and sampling effects that influence occurrence and detection. Thus the community approach can provide the best possible estimates of species richness and other metrics of interest across a heterogeneous landscape, while accounting for variation in occurrence and detection among species..

This repo provides R and WinBUGS code to implement a model with covariates using this modeling framework. The code here estimates species-specific detection and occupancy assuming covariate effects on both processes. It also estimates species richness using data augmentation.

#It is designed to estimate static species-specific occupancy and detection 
#with site specific habitat and sampling covariates using the community model.
#The occurence data are in the file "occ data.csv". 
#The covariate data are in the files "habitat.csv" (occurence) and "date.csv" (detection).
#Species are grouped into one of three categories in the file "groups.csv"

See: https://www.mbr-pwrc.usgs.gov/site/communitymodeling/home/ for more information.

## **Data:**
occ data.csv - Bird species detection nondetection data. Each column is an observation containing information on 1) the site where the observation occured (Point), 2) the time (Time) and date (Date) of the detection, 3) the species that was detected (Species), and 4) which survey replicate (Rep) the observation data was collected (as determined by the unique site/date combination).

habitat.csv - Habitat covariate data used to estimate species-specific occupancy.

date.csv - Sampling covariate data used to estimate species-specific detection.

group.csv - Species grouping data used to estimate assemblage richness.

## **Code:**
covariate model code.R - R code to run the multi-species occupancy model with covariates on occupancy (habitat type - CATO,FCW, understory foliage cover - ufc, tree basal area - ba) and detection (habitat type - CATO,FCW, and survey date - date). Contains code to import and reshape the data, create the BUG model file, and run the model file in WinBUGS. There is also code to process the results and make some figures.
