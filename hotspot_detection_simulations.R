# Simulation study
# RENE STANDER
# Cleaned

library(png)
library(dplyr)
library(Matrix)
library(sf)
library(rlist)
library(spatstat)
library(igraph)
library(tidygraph)

# Import shape files for the region
# Unfortunately these cannot be shared
boundaries <- read_sf("../Shapefiles/MN_SA_2011.shp")
boundaries_dc <- read_sf("../Shapefiles/DC_SA_2011.shp")
boundaries_pr <- read_sf("../Shapefiles//PR_SA_2011.shp")
fs <- boundaries_dc %>% filter(PR_NAME == "Free State") %>% select(DC_NAME) %>% unique()
boundaries_fs <- boundaries %>% filter(DC_NAME %in% fs$DC_NAME)


num_reg <- c(boundaries_fs, boundaries_dc)
hotspot <- c(0.05, 0.3, 0.6)
num_hot <- c("one", "two", "two_sep")

points_list <- list()

c <- 0


for(i in 1:2){
  
  if(i == 1){
    data <- boundaries_fs
  } else{
    data <- boundaries_dc
  }
  
  neigh_sparse <- create_neigh_mat_lat(data)
  neigh_sparse_2 <- neigh_sparse
  for(p in 1:nrow(data)){
    neigh_sparse_2[p,] <- neigh_sparse_2[p,]/sum(neigh_sparse_2[p,])
  }
  
  for(h in hotspot){
    for(k in num_hot){
      for(m in 1:500){
        
        print(c)
        c <- c + 1
        
        results$num[c] <- m
        results$num_reg[c] <- nrow(data)
        results$hotspot[c] <- h
        results$num_hot[c] <- k
        
        hom <- st_as_sf(st_sample(data$geometry, nrow(data)*50))
        ts <- sample.int(nrow(data), 1)
        clust <- st_as_sf(st_sample(data$geometry[ts], size = nrow(data)*50*h))
        
        results$hot1_true[c] <- ts
        
        if(k == "two"){
          ts2 <- sample(which(neigh_sparse[ts,] == 1),1)
          results$hot2_true[c] <- ts2
          clust1<- st_as_sf(st_sample(data$geometry[ts2], size = nrow(data)*50*h))
          pat <- rbind(hom, clust, clust1)
          t <- c(ts, ts2)
        } else if(k == "two_sep"){
          ts2 <- sample(which(neigh_sparse[ts,] != 1),1)
          results$hot2_true[c] <- ts2
          clust1<- st_as_sf(st_sample(data$geometry[ts2], size = nrow(data)*50*h))
          pat <- rbind(hom, clust, clust1)
          t <- c(ts, ts2)
        } else{
          pat <- rbind(hom, clust)
          t <- ts
        }
        
        
        
        pat <- st_transform(pat, st_crs(data$geometry))
        data$points <- rowSums(as.matrix(st_contains(data$geometry, pat)))
        data$points_by_area <- round(data$points/data$AREA_SQKM,5)
        
        points_list <- list.append(points_list, data$points_by_area)
        
      }
    }
  }
}


points_list_df <- as.data.frame(points_list[[1]])
colnames(points_list_df) <- "points"
points_list_df$num <- 1

for(i in 2:length(points_list)){
  print(i)
  d <- as.data.frame(points_list[[i]])
  colnames(d) <- "points"
  d$num <- i
  
  points_list_df <- rbind(points_list_df, d)
}

