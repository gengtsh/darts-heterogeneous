library(plyr)
library(randomForest)
setwd("~/GIT/darts-heterogeneous/data_img/data/stencil3d7pt/profile/")

# machines <- c("ccsl", "f3", "supermicro", "debian")
machines <- c("ccsl", "debian","f3", "supermicro")
versions <- c("gpu_streams", "gpu_streamNot")
data_gpu_streamNot <- data.frame()    
data_gpu_streams <- data.frame()
for (version in versions){
    flag_streams = FALSE
    flag_streamNot = FALSE
    for (machine in c(1:length(machines))){
    
        dataTemp <- data.frame(read.csv(paste("./ScriptAndOutput/GPUExecutionTimeSAndNS/nvprofTraces_", machines[machine], "_",version, ".csv", sep="")))
        # dataTemp <- dataTemp[-c(5,6,7,8)]
        dataTemp$sum <- as.numeric(dataTemp$sum)
        dataTemp$size <- as.character(dataTemp$size)
        # dataTemp2 <- dataTemp[dataTemp$metric %in% c("total_CUDA_functions", "total_running_time"),]
        sizes <- count(dataTemp, c("size"))[1]
        iterations <- count(dataTemp, c("iteration"))[1]
        
        for (size in c(1:dim(sizes)[1])){
            size_dim <- strsplit(x = sizes[size,],split = '_')
            for (iter in c(1:dim(iterations)[1])){
                dataTemporal <- as.array(c(machine, as.numeric(size_dim[[1]][1]),as.numeric(size_dim[[1]][2]),as.numeric(size_dim[[1]][3]), iter, 
                                           dataTemp[dataTemp$size == sizes[size,] & dataTemp$iteration == iter,'sum']))
                for (column in c('min', 'max', 'mean', 'invocations')){
                    dataTemporal <- c(dataTemporal, dataTemp[dataTemp$size == sizes[size,] & dataTemp$iteration == iter,column])
                }
                
                metric_names <- dataTemp[dataTemp$size == sizes[size,] & dataTemp$iteration == iter,'metric']
                names_dataTemporal <- c('machine', 'dim_X', 'dim_Y', 'dim_Z', 'iteration', paste(metric_names, '_sum', sep=""), paste(metric_names, '_min', sep=""),
                                        paste(metric_names, '_max', sep=""),paste(metric_names, '_mean', sep=""),
                                        paste(metric_names, '_invocations', sep=""))
                names(dataTemporal) <- names_dataTemporal
                
                if(version == "gpu_streams"){
                    data_gpu_streams <- rbind(data_gpu_streams, dataTemporal)
                    if(flag_streams != TRUE){
                        names(data_gpu_streams) <- names_dataTemporal
                        flag_streams = TRUE
                    }
                }
            
                if(version == "gpu_streamNot"){
                    data_gpu_streamNot <- rbind(data_gpu_streamNot, dataTemporal)
                    if(flag_streamNot != TRUE){
                        names(data_gpu_streamNot) <- names_dataTemporal
                        flag_streamNot = TRUE
                    }
                }
            }
        }
    }
}


# if(version == "gpu_streamNot")
#     write.csv(data_gpu_streamNot, file = "./ScriptAndOutput/GPUExecutionTimeSAndNS/data_gpu_streamNot.csv", row.names=FALSE)
# if(version == "gpu_streams")
#     write.csv(data_gpu_streams, file = "./ScriptAndOutput/GPUExecutionTimeSAndNS/data_gpu_streams.csv", row.names=FALSE)

data_gpu_streams['exec_time'] <- 0
data_gpu_streamNot['exec_time'] <- 0

for (version in versions){
    for (machine in c(1:length(machines))){
        dataTemp <- data.frame(read.csv(paste("./ScriptAndOutput/GPUExecutionTimeSAndNS/weak_execution_", machines[machine], "_",version, ".csv", sep="")))
        names(dataTemp) <- c('size', 'volumn', c(1:10))
        
        size_dim <- strsplit(x = as.character(dataTemp$size),split = "*", fixed=TRUE)
        for (size in c(1:length(size_dim))){
            for (iter in c(1:9)){
                
                if(version == "gpu_streams"){
                    data_gpu_streams$exec_time[data_gpu_streams$dim_X == size_dim[[size]][2] &
                                                   data_gpu_streams$dim_Y == size_dim[[size]][3] &
                                                   data_gpu_streams$dim_Z == size_dim[[size]][1] &
                                                   data_gpu_streams$iteration == iter &
                                                   data_gpu_streams$machine == machine
                                                   ] <-  dataTemp[[size,iter+2]]
                    
                }
                if(version == "gpu_streamNot"){
                    data_gpu_streamNot$exec_time[data_gpu_streamNot$dim_X == size_dim[[size]][2] &
                                                     data_gpu_streamNot$dim_Y == size_dim[[size]][3] &
                                                     data_gpu_streamNot$dim_Z == size_dim[[size]][1] &
                                                     data_gpu_streamNot$iteration == iter &
                                                     data_gpu_streamNot$machine == machine
                                               ] <-  dataTemp[[size,iter+2]]
                }
                
            }
        }
        
    }
}

if(version == "gpu_streamNot")
    write.csv(data_gpu_streamNot, file = "./ScriptAndOutput/GPUExecutionTimeSAndNS/data_gpu_streamNot.csv", row.names=FALSE)
if(version == "gpu_streams")
    write.csv(data_gpu_streams, file = "./ScriptAndOutput/GPUExecutionTimeSAndNS/data_gpu_streams.csv", row.names=FALSE)
