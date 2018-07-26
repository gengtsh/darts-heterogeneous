library(plyr)
setwd("~/GIT/darts-heterogeneous/data_img/data/stencil3d7pt/profile/")

machines <- c("ccsl", "f3", "supermicro", "debian")
versions <- c( "gpu_streams", "gpu_streamNot")

data_gpu_streams <- data.frame()
data_gpu_streamNot <- data.frame()    
for(version in c(versions)){
    for (machine in c(machines)){
        print(machine)
        dataTemp <- data.frame(read.csv(paste("./ScriptAndOutput/GPUExecutionTimeSAndNS/nvprofTraces_", machine, "_",version, ".csv", sep="")))
        dataTemp <- dataTemp[-c(5,6,7,8)]
        dataTemp$sum <- as.numeric(dataTemp$sum)
        # dataTemp2 <- dataTemp[dataTemp$metric %in% c("total_CUDA_functions", "total_running_time"),]
        sizes <- count(dataTemp, c("size"))[1]
        iterations <- count(dataTemp, c("iteration"))[1]
        
        for (size in 1:dim(sizes)[1]){
            for(iter in 1:9){
                dataTemporal <- as.array(t(dataTemp[dataTemp$size == as.character(sizes[size,]) & dataTemp$iteration == iter,"sum"]))
                
                if(version == "gpu_streams")
                    data_gpu_streams <- as.matrix(rbind(data_gpu_streams, c(Machines=machine, Sizes=as.character(sizes[size,]), Iteration=iter, dataTemporal)))
                
                if(version == "gpu_streamNot")
                    data_gpu_streamNot <- as.matrix(rbind(data_gpu_streams, c(Machines=machine, Sizes=as.character(sizes[size,]), Iteration=iter, dataTemporal)))
                
                
            }
        }
    }
    if(version == "gpu_streams")
    colnames(data_gpu_streams) <- c( "Machines","Sizes", "Iteration" , 
                     t(as.character(dataTemp[dataTemp$size == "400_400_400" & dataTemp$iteration == 1,"metric"])))
    
    if(version == "gpu_streamNot")
        colnames(data_gpu_streamNot) <- c( "Machines","Sizes", "Iteration" , 
                                         t(as.character(dataTemp[dataTemp$size == "400_400_400" & dataTemp$iteration == 1,"metric"])))
}


dataTS <- data.frame(read.csv(paste("./ScriptAndOutput/GPUExecutionTimeSAndNS/nvprofTraces_overlap.csv", sep="")))
