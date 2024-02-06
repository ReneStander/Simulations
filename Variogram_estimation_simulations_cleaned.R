# Simulations 
# Domain 1

library(sf)
library(dplyr)
library(ggplot2)
library(gstat)
library(raster)

####################################
# Domain
# Police precincts in Gauteng

wards <- read_sf("../../../../Data/Crime/police_bounds_crime.shp")
wards <- wards |> filter(Province == "Gauteng")
wards <- wards |> st_transform(22234)
wards <- wards |> st_geometry()
wards <- wards |> st_as_sf()
union <- st_union(wards)
centroids <- wards |> st_centroid()
####################################
# GRF simulation grid

domain <- wards |> as("Spatial")

grid <- raster(domain, resolution = c(2000,2000), crs = proj4string(domain))
grid.poly <- rasterToPolygons(grid)
grid.points <- rasterToPoints(grid, spatial = T)
grid.in <- intersect(grid.poly, domain)

field.sf <- st_as_sf(grid.points)
field.sf <- field.sf |> st_set_crs(22234)

contain <- st_contains(union, field.sf)
field.sf <- field.sf[contain[[1]],]

field <- field.sf |> as("Spatial") |> as.data.frame()
colnames(field) <- c("x", "y")

poi.in <- st_contains(wards, field.sf)
####################################
# Parameters

# Simulation
range <- c(30000, 50000, 70000, 90000)
nugget <- c(10, 15, 20)
psill <- c(200, 250, 350)
Beta <- c(200, 500, 1000)
mod <- c("Sph", "Exp", "Gau")

simulation.df <- expand.grid(range = range, 
                             nugget = nugget,
                             psill = psill,
                             beta = Beta,
                             model = mod,
                             rep = 1:3)

# Variogram estimation
num_points <- c(300, 400, 500, 600)
adap <- c(3,4,5,6,7)
iter <- c(250,500,1000)

estimation.df <- expand.grid(numpoints = num_points,
                             iter = iter)


####################################
# Simulation

results.sim <- tibble(
  range = numeric(0),
  nugget = numeric(0),
  psill = numeric(0),
  beta = numeric(0),
  model = character(0),
  rep = numeric(0),
  variance = numeric(0),
  num_points = numeric(0),
  adap = numeric(0),
  iter = numeric(0),
  nug.sph = numeric(0),
  sill.sph = numeric(0),
  ran.sph = numeric(0),
  nug.exp = numeric(0),
  sill.exp = numeric(0),
  ran.exp = numeric(0),
  nug.gau = numeric(0),
  sill.gau = numeric(0),
  ran.gau = numeric(0),
  n.sph = numeric(0),
  n.exp = numeric(0),
  n.gau = numeric(0)
)

for(ii in 1:nrow(simulation.df)){
  dsim <- simulation.df[ii,]
  
  theo <- vgm(psill = dsim$psill, range = dsim$range, 
              nugget = dsim$nugget, model = as.character(dsim$model))
  
  RDT_modelling <- gstat(formula=z~1, ## We assume that there is a constant trend in the data
                      locations=~x+y,
                      dummy=T,    ## Logical value to set to True for unconditional simulation
                      beta=dsim$beta,  ## Necessity to set the average value over the field
                      model=theo, 
                      nmax=100)
  grf <- predict(RDT_modelling, newdata=field, nsim=1)
  
  rgf <- as.data.frame(grf)
  coordinates(rgf) <- ~x+y
  
  sim.sf <- st_as_sf(rgf)
  sim.sf <- sim.sf |> st_set_crs(22234)
  
  wards$simm <- NA
  for(i in 1:length(poi.in)){
    din <- poi.in[[i]]
    wards$simm[i] <- mean(sim.sf$sim1[din])
  }
  
  for(jj in 1:nrow(estimation.df)){
    print(c(ii,jj))
    
    dd <- estimation.df[jj,]

    variog <- LatticeVario(wards, "simm",
                           cent = centroids,
                           nran = dd$numpoints,
                           nreg = dd$numpoints,
                           niter = dd$iter)

    resu <- variog[[2]]
    
    to_concat <- cbind(dsim, var(wards$simm), dd, resu)
    results.sim <- rbind(results.sim,to_concat)
    
  }
  
}



















