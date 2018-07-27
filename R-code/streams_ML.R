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


PP_data_gpu <- data.frame()
for(metric in c("cudadriver","gpu_running","kernel","overlap","transDtoH","transHtoD")){
    dataTemp <- data.frame(Metrics=metric, read.csv(paste("./ScriptAndOutput/GPUExecutionTimeSAndNS/nvprof_postpose/nvprofTraces_", metric, ".csv", sep="")))
    
    for(machine in machines){
            dataNew <- dataTemp[dataTemp$sever == machine, 1:5]
            dataIter <- dataTemp[dataTemp$sever == machine, 6:15]
            PP_data_gpu_temporal <- data.frame()
            for (i in 1:10){
                PP_data_gpu_temporal <- rbind(PP_data_gpu_temporal, cbind(dataNew, iteration=i, time=dataIter[,paste("its.", i, sep="")]))
            }
            PP_data_gpu <- rbind(PP_data_gpu, PP_data_gpu_temporal)
    }
}


PP_data_gpu <- PP_data_gpu[PP_data_gpu$iteration != 10,]
PP_data_gpu$time <- PP_data_gpu$time/max(PP_data_gpu$time)

par("mar")
par(mar=c(1,1,1,1), mfrow=c(6,2))
# png("./plot_postpone_data.png", width = 600, height = 800)

for(metric in c("cudadriver","gpu_running","kernel","overlap","transDtoH","transHtoD")){
    for(version in versions){
        
        trainingSet <- PP_data_gpu[PP_data_gpu$stream.Y.NO == version & PP_data_gpu$Metrics == metric &
                                       PP_data_gpu$iteration <= 5 ,]
        
        testSet <- PP_data_gpu[PP_data_gpu$stream.Y.NO == version & PP_data_gpu$Metrics == metric &
                                       PP_data_gpu$iteration >= 7 & PP_data_gpu$iteration <= 9 ,]
        
        trainingSet$Metrics <- NULL
        trainingSet$stream.Y.NO <- NULL
        trainingSet$size <- NULL
        trainingSet$sever <- NULL
        trainingSet[apply(trainingSet, 2, is.infinite)] <- 0
        trainingSet[apply(trainingSet, 2, is.na)] <- 0
        trainingSet <- trainingSet[, apply(trainingSet, 2, function(v) var(v, na.rm = TRUE) != 0)]
        
        testSet$Metrics <- NULL
        testSet$stream.Y.NO <- NULL
        testSet$size <- NULL
        testSet$sever <- NULL
        testSet[apply(testSet, 2, is.infinite)] <- 0
        testSet[apply(testSet, 2, is.na)] <- 0
        testSet <- testSet[, apply(testSet, 2, function(v) var(v, na.rm = TRUE) != 0)]
        durationReal <- testSet$time
        testSet$time = NULL
        
        fit.model <- lm((trainingSet$time)~., data=trainingSet)
        predictTime <- predict(fit.model, testSet)

        accuracy <- predictTime/durationReal
        
        boxplot(accuracy, main=paste(metric, " in version ",version, sep="" ))
    }
}
# dev.off()

data_Overlap <- data.frame(read.csv(paste("./ScriptAndOutput/GPUExecutionTimeSAndNS/nvprofTraces_overlap.csv", sep="")))


data_gpu_streams <- as.data.frame(data_gpu_streams)
data_gpu_streams$total_CUDA_functions <- NULL


# fit <- lm(data_gpu_streams$total_running_time~., data = data_gpu_streams)