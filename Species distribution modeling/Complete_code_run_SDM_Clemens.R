library(sp)
library(terra)
library(maps)
library(sdm)
library(raster)
library(ncdf4)

#######1st step: prepare species occurrence data

# Load IUCN Data
example_amphibian_data <- terra::vect('C:/Users/Claus/MasterThesis/Data/AMPHIBIANS/AMPHIBIANS.shp')

# Choose a specific species
specific_species <- example_amphibian_data[example_amphibian_data$sci_name == 'Caecilia tentaculata', ]

# Calculate the extent of the specific species
species_extent <- ext(specific_species)

# Add relative extent for plotting
extent_buffer <- 0.3  

# Calculate xlim and ylim for plotting species occurrence
xlim <- c(species_extent$xmin - (species_extent$xmax - species_extent$xmin) * extent_buffer,
          species_extent$xmax + (species_extent$xmax - species_extent$xmin) * extent_buffer)
ylim <- c(species_extent$ymin - (species_extent$ymax - species_extent$ymin) * extent_buffer,
          species_extent$ymax + (species_extent$ymax - species_extent$ymin) * extent_buffer)

# Define resolution of the raster (z.B. 0.5°)
resolution <- 0.5  # in Grad

# Create a raster in the occurrence region of the species
r_custom_resolution <- terra::rast(ext = species_extent, res = resolution)

# Rasterize the specific species on the predefined raster
specific_species_custom_resolution <- terra::rasterize(specific_species, r_custom_resolution)

# transform raster in SpatialPointsDataFrame
raster_points <- terra::as.points(specific_species_custom_resolution)
raster_points_spdf <- as(raster_points, "Spatial")
names(raster_points_spdf) <- "species"


#Plot information and maps of the species occurrence data:
#print(raster_points_spdf)
#maps::map('world', xlim = xlim, ylim = ylim)
#plot(specific_species_custom_resolution, add = TRUE, alpha = 0.6, legend = FALSE)

#######2nd step: prepare environmental data

# Load Data (prepared from Claus)
nc_file <- "C:/Users/Claus/MasterThesis/Data/ClimateVariables/bioclimatic_variables_old.nc"
climate_data <- terra::rast(nc_file)

# Close the NetCDF file
nc_close(nc)

# Get the name of the layers
layer_names <- names(climate_data)

# Filter the layers for the year 2022 data
selected_layers <- grep("2022", layer_names, value = TRUE)

# Create a new SpatRaster object with the filtered layers
climate_data_2022 <- climate_data[[selected_layers]]

# Change the names of the columns (to not run into problems later)
names(climate_data_2022)[names(climate_data_2022) == "MAT_year=2022"] <- "MAT"
names(climate_data_2022)[names(climate_data_2022) == "MTWM_year=2022"] <- "MTWM"
names(climate_data_2022)[names(climate_data_2022) == "MTCM_year=2022"] <- "MTCM"
names(climate_data_2022)[names(climate_data_2022) == "AP_year=2022"] <- "AP"
names(climate_data_2022)[names(climate_data_2022) == "PDQ_year=2022"] <- "PDQ"

# > climate_data_2022
# class       : SpatRaster
# dimensions  : 721, 1440, 5  (nrow, ncol, nlyr)
# resolution  : 0.25, 0.25  (x, y)
# extent      : -180.125, 179.875, -90.125, 90.125  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs
# sources     : bioclimatic_variables_old.nc:MAT
#               bioclimatic_variables_old.nc:MTWM
#               bioclimatic_variables_old.nc:MTCM
#               ... and 2 more source(s)
# varnames    : MAT
#               MTWM
#               MTCM
#               ...
# names       : MAT, MTWM, MTCM, AP, PDQ
# unit        :    ,     ,     ,  m,   m


# Transform the Climate data to RasterStack data format 
climate_data_2022_raster <- as(climate_data_2022, "Raster")
climate_data_2022_stack <- stack(climate_data_2022_raster)
# now make climate_data_2022_stack into Raster* object


# Plot the climate data

#plot all climate variables
#plot(climate_data_2022_stack)

###### 3rd step: Combine the two data sets using the sdm package

# Create an sdmData object with 3 times the number of background points then occurrence points
num_features <- nrow(raster_points_spdf) 
sdm_data <- sdmData(formular= species~., train = raster_points_spdf, predictors = climate_data_2022_stack, bg = list(n=num_features * 3, method ='gRandom'))

# Build and train model with sdm Data

# List of possible methods:
#'glm','gam','maxent','rf','tree','fda','mars','svm','brt'
methods <- c('bioclim','gam')

# Define test percentage for testing data (splitting data set in test and train data)
test_precent <- 30
# Runs of subsampling
n <- 2

# Run the model
sdm_model <- sdm(species~.,data=sdm_data,methods=methods, replication='sub', test.precent=test_precent, n=n)

# Consider performance measures ROC for model run
roc(sdm_model)

###### 4th step: Prediction of the habitat suitability into the whole study area
# add slot of name "ccp" to climate_data_2022_stack
pred_model <- predict(sdm_model,newdata=climate_data_2022_stack)
#Plot 
plot(pred_model)

### Ensemble approach:
e1 <- ensemble(sdm_model,newdata=climate_data_2022_stack,setting=list(method='weighted',stat='AUC'))
plot(e1)

### Plot niche space

niche_plot <- niche(climate_data_2022_stack,e1,n=c('MTCM','MAT'))

###### Some plots
#Plot species
maps::map('world', xlim = xlim, ylim = ylim)
plot(specific_species, col = 'red', add = TRUE)
plot(specific_species_custom_resolution,col="grey", add = TRUE, alpha = 0.6, legend = FALSE)
plot(climate_data_2022_stack)

#Plot prediction 
maps::map('world')
plot(e1, add = T, alpha = 0.6)



