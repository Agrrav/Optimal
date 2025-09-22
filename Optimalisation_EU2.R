##### libraries #####
library(readxl)
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(caret)
library(imputeTS)
library(data.table)
library(clusterCrit)
library(dtwclust)
library(TSrepr)
library(TSclust)
library(NbClust)
library(LPWC)
library(clValid)
library(mclust)
library(patchwork)

#testing
library(dbscan)
library(fpc)

# loading datasets
setwd("C:/Users/Karolina/Documents/Doktorat/Time_series_DT")
House_EU <- as.data.frame(read_excel("Data/Housing_prices_EU_c.xlsx"))
ID <- House_EU[,1]
rownames(House_EU) <- ID
data <- House_EU[,-1]

################################################################################
##### visualisation ############################################################

graph_func <- function(data){
  view <- as.data.frame(cbind(seq(1,ncol(data)),t(data)))
  colnames(view) <- c("year", ID)
  
  # Convert sample_data from wide form to long form
  view_data <- reshape2::melt(view, id.vars = "year")
  
  # Plot the final data
  ggplot(view_data,                            
         aes(x = year,
             y = value,
             col = variable)) + geom_line()
}

################################################################################
##### Clustering ###############################################################
################################################################################

### Distance matrixes
##### distance matrices and data
process <- preProcess(t(data), method=c("center", "scale"))
data_norm <- as.data.frame(t(predict(process, t(data))))
process <- preProcess(t(data), method=c("range"))
data_range <- as.data.frame(t(predict(process, t(data))))

dist_list <- list(dtw_distmat_raw <- as.matrix(proxy::dist(data, method = "dtw", 
                                                           upper = TRUE, diag = TRUE)),
                  dtw_distmat_norm <- as.matrix(proxy::dist(data_norm, method = "dtw", 
                                                            upper = TRUE, diag = TRUE)),
                  sbd_distmat_norm <- as.matrix(proxy::dist(data_norm, method = "sbd", 
                                                            upper = TRUE, diag = TRUE)))

graph_func(data)
graph_func(data_norm)

################################################################################
##### MDS

MDS_func <- function(distmat, id = ID){
  
  fit1 <- cmdscale(distmat, eig = TRUE, k=2)
  
  mds_data1 <- data.frame(ID = id,
                          x = fit1$points[,1],
                          y = fit1$points[,2])
  
  ggplot(mds_data1,                            
         aes(x = x,
             y = y,
             label = ID)) + 
    geom_point(size = 5) + 
    labs(title = "Metric MDS",
         x = "Coordinate 1", y = "Coordinate 2") + 
    geom_text(hjust = 1, vjust = 2)
}

MDS_func(dist_list[[1]])
MDS_func(dist_list[[2]])
MDS_func(dist_list[[3]])

################################################################################
##### Cluster number ###########################################################
################################################################################

CH_values <- list()
av_CH <- list()
max_CH <- list()
median_CH <- list()
k <- nrow(data)

clusterings_tmp1 <- lapply(c(2:(k-1)), function(x)
  tsclust(data_norm, type = "p", k = x,
          distance = "dtw_basic", centroid = "mean",
          seed = 42,
          trace = TRUE,
          control = partitional_control(nrep = 100L)))

clusterings_tmp2 <- lapply(c(2:(k-1)), function(x)
  tsclust(data_norm, type = "p", k = x,
          distance = "dtw_basic", centroid = "pam",
          seed = 42,
          trace = TRUE,
          control = partitional_control(nrep = 100L)))

clusterings_tmp3 <- lapply(c(2:(k-1)), function(x)
  tsclust(data_norm, type = "p", k = x,
          distance = "dtw_basic", centroid = "dba",
          seed = 42,
          trace = TRUE,
          control = partitional_control(nrep = 100L)))

clust_CH <- list(mean = clusterings_tmp1,
                 pam = clusterings_tmp2,
                 dba = clusterings_tmp3)
rm(clusterings_tmp1)
rm(clusterings_tmp2)
rm(clusterings_tmp3)

temp_avg_CH <- c()
temp_max_CH <- c()
temp_median_CH <- c()

for(j in 1:3){
  CH_values[[j]] <- list()
  for(i in 1:(k-2)){
    CH_values[[j]][[i]] <- list()
    CH_values[[j]][[i]] <- sapply(seq_along(clust_CH[[j]][[i]]), function(x) 
      intCriteria(as.matrix(data_norm), as.integer(clust_CH[[j]][[i]][[x]]@cluster), c("Calinski_Harabasz")))
    
    temp_avg_CH[i] <- mean(unlist(CH_values[[j]][[i]]))
    temp_max_CH[i] <- max(unlist(CH_values[[j]][[i]]))
    temp_median_CH[i] <- median(unlist(CH_values[[j]][[i]]))
  }
  av_CH[[j]] <- temp_avg_CH
  max_CH[[j]] <- temp_max_CH
  median_CH[[j]] <- temp_median_CH
}

CH_plots <- list()
for(j in 1:3){
  CH_plots[[j]] <- list()
  CH_plots[[j]][[1]] <- ggplot(data.table(Clusters = 2:(k-1), CHindex = av_CH[[j]]),
                               aes(Clusters, CHindex)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    theme_bw() + 
    ggtitle("Average")
    
  CH_plots[[j]][[2]] <- ggplot(data.table(Clusters = 2:(k-1), CHindex = max_CH[[j]]),
                               aes(Clusters, CHindex)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    theme_bw() + 
    ggtitle("Maximum")
  
  CH_plots[[j]][[3]] <- ggplot(data.table(Clusters = 2:(k-1), CHindex = median_CH[[j]]),
                               aes(Clusters, CHindex)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    theme_bw() + 
    ggtitle("Median")
}

CH_plots[[1]][[1]] + CH_plots[[1]][[2]] + CH_plots[[1]][[3]]
CH_plots[[2]][[1]] / CH_plots[[2]][[2]] / CH_plots[[2]][[3]]
CH_plots[[3]][[1]] / CH_plots[[3]][[2]] / CH_plots[[3]][[3]]

## BOXPLOT

for(j in 1:3){
  CH_values[[j]] <- list()
  for(i in 1:(k-2)){
    CH_values[[j]][[i]] <- list()
    CH_values[[j]][[i]] <- sapply(seq_along(clust_CH[[j]][[i]]), function(x) 
      intCriteria(as.matrix(data_norm), as.integer(clust_CH[[j]][[i]][[x]]@cluster), c("Calinski_Harabasz")))
    
    temp_avg_CH[i] <- mean(unlist(CH_values[[j]][[i]]))
    temp_max_CH[i] <- max(unlist(CH_values[[j]][[i]]))
    temp_median_CH[i] <- median(unlist(CH_values[[j]][[i]]))
  }
  av_CH[[j]] <- temp_avg_CH
  max_CH[[j]] <- temp_max_CH
  median_CH[[j]] <- temp_median_CH
}

for(j in 1:3){
  CH_values[[j]] <- list()
  for(i in 1:(k-2)){
    CH_values[[j]][[i]] <- list()
    CH_values[[j]][[i]] <- sapply(seq_along(clust_CH[[j]][[i]]), function(x) 
      intCriteria(as.matrix(data_norm), as.integer(clust_CH[[j]][[i]][[x]]@cluster), c("Calinski_Harabasz")))
    
    temp_avg_CH[i] <- mean(unlist(CH_values[[j]][[i]]))
    temp_max_CH[i] <- max(unlist(CH_values[[j]][[i]]))
    temp_median_CH[i] <- median(unlist(CH_values[[j]][[i]]))
  }
  av_CH[[j]] <- temp_avg_CH
  max_CH[[j]] <- temp_max_CH
  median_CH[[j]] <- temp_median_CH
}


my_list <- CH_values[[1]]
names(my_list) <- c(LETTERS)

df <- bind_rows(
  lapply(names(my_list), function(group) {
    data.frame(
      value = unlist(my_list[[group]]),
      group = group
    )
  })
)

ggplot(df, aes(x = group, y = value)) +
    geom_boxplot(coef = Inf) +
    theme_bw() + 
    ggtitle("Calinski_Harabasz")

###############################################################################
pairs.matrix <- matrix(0, k, k)
colnames(pairs.matrix) <- ID
rownames(pairs.matrix) <- ID

for (s in 1:100){
  cluster <- clust_CH[[1]][[2]][[s]]@cluster
  
  for (i in 1:(length(cluster)-1))
    for (j in (i+1):length(cluster))
      if (cluster[i] == cluster[j]) pairs.matrix[i,j] <- pairs.matrix[i,j] + 1
}

pairs.matrix/100

hc1 <- hclust(as.dist(t(1-pairs.matrix/100)), method = "complete")
plot(hc1)

fit1 <- cmdscale(dist_list[[2]], eig = TRUE, k=2)
mds_data1 <- data.frame(ID = ID,
                        x = fit1$points[,1],
                        y = fit1$points[,2],
                        group = as.factor(cutree(hc1, k = 3)))

mds_plot1 <- ggplot(mds_data1,                            
                   aes(x = x,
                       y = y,
                       label = ID,
                       col = group)) + 
              geom_point(size = 5) + 
              labs(title = "Initial result",
                   x = "Coordinate 1", y = "Coordinate 2") + 
              geom_text(hjust = 1, vjust = -1)

write.csv(pairs.matrix,"pairs.matrix.csv", row.names = TRUE)

###############################################################################
#### Bootstrap ################################################################
###############################################################################
pairs.matrix.boot <- list()
pairs.matrix.boot2 <- matrix(0, k, k)
colnames(pairs.matrix.boot2) <- ID
rownames(pairs.matrix.boot2) <- ID

for (s in 1:k){
  pairs.matrix.temp <- matrix(0, k-1, k-1)
  colnames(pairs.matrix.temp) <- ID[-s]
  rownames(pairs.matrix.temp) <- ID[-s]
  data_temp <- data_norm[-s,]
  
  for (l in 1:100){
    clust_result_temp <- tsclust(data_temp, type = "p", k = 3,
                                 distance = "dtw", centroid = "mean",
                                 trace = TRUE,
                                 control = partitional_control(nrep = 1L))
    
    cluster_temp <- clust_result_temp@cluster
    
    for (i in 1:(length(cluster_temp)-1))
      for (j in (i+1):length(cluster_temp))
        if (cluster_temp[i] == cluster_temp[j]) pairs.matrix.temp[i,j] <- pairs.matrix.temp[i,j] + 1
    
    if(s == 1){
      cluster <- c(0,clust_result_temp@cluster)
    } else if(s == k){
      cluster <- c(clust_result_temp@cluster,0)
    }else {
      take_out <- seq(1,s,1)
      cluster <- c(clust_result_temp@cluster[take_out],0,clust_result_temp@cluster[-take_out])
    }
    for (i in 1:(length(cluster)-1))
      for (j in (i+1):length(cluster))
        if (cluster[i] == cluster[j]) pairs.matrix.boot2[i,j] <- pairs.matrix.boot2[i,j] + 1
  }
  
  pairs.matrix.boot[[s]] <- pairs.matrix.temp
}

hc2 <- hclust(as.dist(t(1-pairs.matrix.boot2/(k*100))), method = "complete")
plot(hc2)


f1 <- cutree(hc2, k = 3)
f1[which(f1 == 1)] <- 4
f1[which(f1 == 2)] <- 1
f1[which(f1 == 4)] <- 2

fit1 <- cmdscale(dist_list[[2]], eig = TRUE, k=2)
mds_data1 <- data.frame(ID = ID,
                        x = fit1$points[,1],
                        y = fit1$points[,2],
                        group = as.factor(f1))

mds_plot2 <- ggplot(mds_data1,                            
                   aes(x = x,
                       y = y,
                       label = ID,
                       col = group)) + 
              geom_point(size = 5) + 
              labs(title = "Modified result",
                   x = "Coordinate 1", y = "Coordinate 2") + 
              geom_text(hjust = 1, vjust = -1)

par(mfrow=c(1,2))
plot(hc1)
plot(hc2)

sapply(seq_along(pairs.matrix.boot), function(x)
  write.csv(pairs.matrix.boot[[x]], paste0("~/Doktorat/Time_series_DT/pairs.matrix.boot/pairs.matrix.boot",x,".csv"), row.names = TRUE))
write.csv(pairs.matrix.boot2, "~/Doktorat/Time_series_DT/pairs.matrix.boot/pairs.matrix.boot.clean.csv", row.names = TRUE)

###############################################################################
###############################################################################

mds_plot1 + mds_plot2

factor1 <- cutree(hc1, k = 3)
factor2 <- cutree(hc2, k = 3)

plot_data <- cbind(data, factor1, factor2, factor1 == factor2); plot_data

g1_1 <- data[factor1 == 1,]
g1_2 <- data[factor1 == 2,]
g1_3 <- data[factor1 == 3,]

g2_1 <- data[factor2 == 1,]
g2_2 <- data[factor2 == 2,]
g2_3 <- data[factor2 == 3,]


graph_func2 <- function(data, f){
  view <- as.data.frame(cbind(seq(1,ncol(data)),t(data)))
  colnames(view) <- c("year", ID[f])
  
  # Convert sample_data from wide form to long form
  view_data <- reshape2::melt(view, id.vars = "year")
  
  # Plot the final data
  ggplot(view_data,                            
         aes(x = year,
             y = value,
             col = variable)) + geom_line()
}

graph_func2(g2_1, f = which(factor2 == 1)) +
graph_func2(g2_2, f = which(factor2 == 2)) +
graph_func2(g2_3, f = which(factor2 == 3))
