library(plyr)
library(randomForest)
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
    colnames(data_gpu_streams) <- c( "Machines","size", "Iteration" , 
                     t(as.character(dataTemp[dataTemp$size == "400_400_400" & dataTemp$iteration == 1,"metric"])))
    
    if(version == "gpu_streamNot")
        colnames(data_gpu_streamNot) <- c( "Machines","size", "Iteration" , t(as.character(dataTemp[dataTemp$size == "400_400_400" & dataTemp$iteration == 1,"metric"])))
}

data_gpu_streamNot <- data.frame(data_gpu_streamNot)
data_gpu_streams <- data.frame(data_gpu_streams)


PP_data_gpu <- data.frame()
for(metric in c("cudadriver","gpu_running","kernel","overlap","transDtoH","transHtoD")){
    dataTemp <- data.frame(Metrics=metric, read.csv(paste("./ScriptAndOutput/GPUExecutionTimeSAndNS/nvprof_postpose/nvprofTraces_", metric, ".csv", sep="")))
    
    for(machine in machines){
            dataNew <- dataTemp[dataTemp$sever == machine, 1:5]
            dataIter <- dataTemp[dataTemp$sever == machine, 6:15]
            PP_data_gpu_temporal <- data.frame()
            for (i in 1:9){
                PP_data_gpu_temporal <- rbind(PP_data_gpu_temporal, cbind(dataNew, iteration=i, time=dataIter[,paste("its.", i, sep="")]))
            }
            PP_data_gpu <- rbind(PP_data_gpu, PP_data_gpu_temporal)
    }
}
colnames(PP_data_gpu)[3] <- "Machines"

# data_gpu_streams <- data_gpu_streams[order(data_gpu_streams$size),]
# data_gpu_streamNot<- data_gpu_streamNot[order(data_gpu_streamNot$size),]

# data_gpu_streamNot <- data_gpu_streamNot[data_gpu_streamNot$Sizes != "1000_1000_1000",]

# PP_data_gpu<-PP_data_gpu[order(PP_data_gpu$size),]
# 
# for(size in count(PP_data_gpu$size)[1]){
#     if (!is.na(size.)){
#         XYZ <- strsplit(as.character(size), "*")
#         PP_data_gpu$size[PP_data_gpu$size == as.character(size)] <- paste(XYZ[1], "_", XYZ[2], "_", XYZ[3], sep = "")
# }

for( version in versions){
    
    for(metric in c("cudadriver","gpu_running","kernel","overlap","transDtoH","transHtoD")){
        for(machine in machines){
            for(iter in 1:9){
                temporal.PP_data_gpu <- PP_data_gpu[PP_data_gpu$stream.Y.NO == version & PP_data_gpu$Metrics == metric
                                   & PP_data_gpu$Machines == machine & PP_data_gpu$iteration == iter,]
                
                temporal.data_gpu_streamNot <- data_gpu_streamNot[data_gpu_streamNot$Machines == machine & data_gpu_streamNot$Iteration == iter,]
                
                combineAll <- matrix()
                for( size_iter in as.character(temporal.data_gpu_streamNot$Sizes)){
                    
                    
                    
                    if (!is.null(as.character(temporal.data_gpu_streamNot$Sizes[temporal.data_gpu_streamNot$Sizes == size_iter])) &
                        !is.null(as.character(temporal.PP_data_gpu$size[temporal.PP_data_gpu$size == size_iter]))){
                            print("Yes!")
                        print(as.character(temporal.data_gpu_streamNot$Sizes[temporal.data_gpu_streamNot$Sizes == size_iter]))
                        print(as.character(temporal.PP_data_gpu$size[temporal.PP_data_gpu$size == size_iter]))
                    } 
                    # combineAll <- temporal.data_gpu_streamNot[temporal.data_gpu_streamNot$Sizes == size_iter]
                }
            }
        }
    }
    
}


# PP_data_gpu$time <- PP_data_gpu$time/max(PP_data_gpu$time)


for (ML in c("lm", "glm", "RF", "svm" )){
png(paste("./plot_postpone_data-", ML, ".png", sep=""), width = 600, height = 800)
par("mar")
par(mar=c(1,1,1,1), mfrow=c(5,2))


for(metric in c("cudadriver","gpu_running","kernel","transDtoH","transHtoD")){
    for(version in versions){

        trainingSet <- PP_data_gpu[PP_data_gpu$stream.Y.NO == version & PP_data_gpu$Metrics == metric &
                                       PP_data_gpu$iteration <= 3 ,]
        
        testSet <- PP_data_gpu[PP_data_gpu$stream.Y.NO == version & PP_data_gpu$Metrics == metric &
                                       PP_data_gpu$iteration >= 8 & PP_data_gpu$iteration <= 9 ,]
        
        trainingSet$Metrics <- NULL
        trainingSet$stream.Y.NO <- NULL
        trainingSet$size <- NULL
        trainingSet$Machines <- NULL
        trainingSet[apply(trainingSet, 2, is.infinite)] <- 0
        trainingSet[apply(trainingSet, 2, is.na)] <- 0
        trainingSet <- trainingSet[, apply(trainingSet, 2, function(v) var(v, na.rm = TRUE) != 0)]
        
        testSet$Metrics <- NULL
        testSet$stream.Y.NO <- NULL
        testSet$size <- NULL
        testSet$Machines <- NULL
        testSet[apply(testSet, 2, is.infinite)] <- 0
        testSet[apply(testSet, 2, is.na)] <- 0
        testSet <- testSet[, apply(testSet, 2, function(v) var(v, na.rm = TRUE) != 0)]
        durationReal <- testSet$time
        testSet$time = NULL
        
        if (ML == "lm")
            fit.model <- lm(log(trainingSet$time, 2)~., data=log(trainingSet + 0.0000000001, 2))
        
        if (ML == "glm")
            fit.model <- glm(log(trainingSet$time, 2)~., data=log(trainingSet + 0.0000000001, 2))
        
        if (ML == "RF")
            fit.model <- randomForest(log(trainingSet$time, 2)~., data=log(trainingSet + 0.0000000001, 2),  ntree=50)
        
        if (ML == "svm")
            fit.model <- svm(log(trainingSet$time, 2)~., data=log(trainingSet + 0.0000000001, 2), kernel="linear")
        
        predictTime <- predict(fit.model, log(testSet+ 0.0000000001, 2))
        predictTime <- 2^predictTime
        accuracy <- predictTime/durationReal
        boxplot(accuracy, main=paste(metric, " in version ",version, sep="" ))
    }
}
dev.off()
}

# data_Overlap <- data.frame(read.csv(paste("./ScriptAndOutput/GPUExecutionTimeSAndNS/nvprofTraces_overlap.csv", sep="")))