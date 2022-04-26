# PACKAGES ####
library(sp)
library(raster)
library(ctmm)
library(GISTools)

# ANALYSIS ####
# Read job number from command line
args = commandArgs(trailingOnly=TRUE)
ind_file <- args[1]
print(ind_file)

# Load data and create telemetry object
df <- fread(ind_file)
aid <- df$individual.local.identifier[1]

l <- as.telemetry(df)

# Load habitat raster (trees) and crop it to save some time and RAM
r1 <- raster("data/treecover2010.tif")
e <- extent(min(l$longitude) - 0.2, max(l$longitude) + 0.2, min(l$latitude) - 0.2, max(l$latitude) + 0.2)
r2 <- crop(r1, e)

# Record start time to monitor how long replicates take to compute
sTime <- Sys.time()
print(sTime)

# Create training and test sets ###
tl <- length(l$timestamp)
fh <- tl/2
sh <- tl/2+1
train <- l[1:fh,]
test <- l[sh:tl,]

# ANALYSIS ####
# Fit the movement model to the test data
svf <- variogram(train)
guess <- ctmm.guess(train, variogram=svf, interactive=FALSE)
fit <- ctmm.select(train, guess, trace=2)
summary(fit) 
print("Fitted movement model")
  
# Calculate the habitat-naive AKDE UD ###
ud <- akde(train, fit, weights=TRUE)
ess <- summary(ud)$DOF["area"]
ud_area <- summary(ud)$CI[1,2]
print("UD created")

# Fit the RSF ###
rsf <- ctmm:::rsf.fit(train, UD=ud, R=list(), debias=TRUE, error=0.01)
summary(rsf)
print("Fitted RSF")

# Calculate the RSF-informed AKDE
ud_rsf <- akde(train, rsf, R=list(trees=r2))
ud_rsf_area <- summary(ud_rsf)$CI[1,2]

# Cross-validate the UDs ###

# Export the UDs as sp objects
ud95 <- SpatialPolygonsDataFrame.UD(ud, level.UD = 0.95)
ud50 <- SpatialPolygonsDataFrame.UD(ud, level.UD = 0.50)

ud95_rsf <- SpatialPolygonsDataFrame.UD(ud_rsf, level.UD = 0.95)
ud50_rsf <- SpatialPolygonsDataFrame.UD(ud_rsf, level.UD = 0.50)

# Export the test data as an sp object
test_sp <- SpatialPoints.telemetry(test)

# Determine which points are inside the 95% OD
pct95 <- round(unname(GISTools:::poly.counts(test_sp, ud95)[2])/nrow(test)*100, 2)
pct50 <- round(unname(GISTools:::poly.counts(test_sp, ud50)[2])/nrow(test)*100, 2)

pct95_rsf <- round(unname(GISTools:::poly.counts(test_sp, ud95_rsf)[2])/nrow(test)*100, 2)
pct50_rsf <- round(unname(GISTools:::poly.counts(test_sp, ud50_rsf)[2])/nrow(test)*100, 2)

eTime <- Sys.time()
  
# Create vector of results
results <- data.frame(aid, ind_file, ess, ud_area, pct95, pct50, ud_rsf_area, pct95_rsf, pct50_rsf)
  
# Store results in data.frame
write.table(results, 'results/cv_summary_tapir.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
  
# Print indicators of progress
print(aid)
print(eTime)
print(eTime - sTime)
