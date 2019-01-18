library(plyr)
library(randomForest)
setwd("~/GIT/darts-heterogeneous/data_img/data/")

# machines <- c("ccsl", "f3", "supermicro", "debian")
machine <- 1
versions <- c( "gpu_streamNot")
data_gpu_streamNot <- data.frame()    
data_gpu_streams <- data.frame()

data <- data.frame(read.csv(paste("./nvprof_postpose_version2.csv", sep="")))
df_exec_time <- data.frame(read.csv(paste("./postpose_weak_execution3.csv", sep="")))

flag_streams = FALSE
flag_streamNot = FALSE

for (version in versions){
    if (version == "gpu_streams"){
        dataTemp <- data[data['IsStream'] == 1,]
        data_exec_time_temp <- df_exec_time[df_exec_time['IsStream'] == 1,]
    } else{
        dataTemp <- data[data['IsStream'] != 1,]
        data_exec_time <- df_exec_time[df_exec_time['IsStream'] != 1,]
    }

    dataTemp$size <- paste(as.character(dataTemp$size.x.), 
                           as.character(dataTemp$size.y.),
                           as.character(dataTemp$size.z.), sep = '-')
    # dataTemp2 <- dataTemp[dataTemp$metric %in% c("total_CUDA_functions", "total_running_time"),]
    sizes <- count(dataTemp, c("size"))[1]
    
    for (size in c(1:dim(sizes)[1])){
        if (!(as.character(sizes[size,]) ==  "100-100-50" | as.character(sizes[size,]) ==  "1000-1000-50" |
              as.character(sizes[size,]) ==  "200-200-50" | as.character(sizes[size,]) ==  "400-400-50" |
              as.character(sizes[size,]) ==  "800-800-50")){ 
            size_dim    <- strsplit(x = sizes[size,],split = '-')
            
            iterations <- max(1, dataTemp$nTile[dataTemp['size'] == sizes[size,]])
            dataTemp$nTile[dataTemp['size'] == sizes[size,]]
    
            temp_exec_time <- data_exec_time$totalExeTime[data_exec_time$size.x. == size_dim[[1]][1] & 
                                data_exec_time$size.y. == size_dim[[1]][2] &
                                data_exec_time$size.z. == size_dim[[1]][3]]
            
            temp_nTile <- data_exec_time$nTile[data_exec_time$size.x. == size_dim[[1]][1] & 
                                                              data_exec_time$size.y. == size_dim[[1]][2] &
                                                              data_exec_time$size.z. == size_dim[[1]][3] ]
            if (version == "gpu_streams"){
                png(paste(sizes[size,],'.png', sep=""))
                plot(temp_nTile, temp_exec_time, main = paste(sizes[size,], sep = ""))
                dev.off()
            }
            
                
            if(version == "gpu_streams"){
                for (iter in c(0:iterations)){
                    dataTemporal <- as.array(c(machine, as.numeric(size_dim[[1]][1]),as.numeric(size_dim[[1]][2]),as.numeric(size_dim[[1]][3]), iter, 
                                               dataTemp[dataTemp$size == sizes[size,] & dataTemp$nTile == iter,'minExe']))
                    for (column in c('avgExe', 'maxExe', 'totalExe', 'callNum')){
                        dataTemporal <- c(dataTemporal, dataTemp[dataTemp$size == sizes[size,] & dataTemp$nTile == iter,column])
                    }
                    
                    metric_names <- as.character(dataTemp$metric[dataTemp$size == sizes[size,] & dataTemp$nTile == iter])
                    names_dataTemporal <- c('machine', 'dim_X', 'dim_Y', 'dim_Z', 'nTile', 
                                            paste(metric_names, '_minExe', sep=""), 
                                            paste(metric_names, '_avgExe', sep=""),
                                            paste(metric_names, '_maxExe', sep=""),
                                            paste(metric_names, '_totalExe', sep=""),
                                            paste(metric_names, '_callNum', sep=""))
                    names(dataTemporal) <- names_dataTemporal
                    data_gpu_streams <- rbind(data_gpu_streams, dataTemporal)
                    if(flag_streams != TRUE){
                        names(data_gpu_streams) <- names_dataTemporal
                        flag_streams = TRUE
                    }
                }
            } else {
                    dataTemporal <- as.array(c(machine, as.numeric(size_dim[[1]][1]),as.numeric(size_dim[[1]][2]),as.numeric(size_dim[[1]][3]), 
                                                   dataTemp[dataTemp$size == sizes[size,],'minExe']))
                        for (column in c('avgExe', 'maxExe', 'totalExe', 'callNum')){
                            dataTemporal <- c(dataTemporal, dataTemp[dataTemp$size == sizes[size,],column])
                        }
                        
                        metric_names <- as.character(dataTemp$metric[dataTemp$size == sizes[size,]])
                        names_dataTemporal <- c('machine', 'dim_X', 'dim_Y', 'dim_Z', 
                                                paste(metric_names, '_minExe', sep=""), 
                                                paste(metric_names, '_avgExe', sep=""),
                                                paste(metric_names, '_maxExe', sep=""),
                                                paste(metric_names, '_totalExe', sep=""),
                                                paste(metric_names, '_callNum', sep=""))
                        names(dataTemporal) <- names_dataTemporal
                        data_gpu_streams <- rbind(data_gpu_streams, dataTemporal)
                        if(flag_streams != TRUE){
                            names(data_gpu_streams) <- names_dataTemporal
                            flag_streams = TRUE
                        }
                    
                data_gpu_streamNot <- rbind(data_gpu_streamNot, dataTemporal)
                    if(flag_streamNot != TRUE){
                        names(data_gpu_streamNot) <- names_dataTemporal
                        flag_streamNot = TRUE
                    }
                }
            }
        }
    }


