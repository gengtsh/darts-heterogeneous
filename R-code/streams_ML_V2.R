library("randomForest")
library("plyr")
library("corrplot")
library("e1071")
library("ggplot2")
library("data.table")
library("cluster")
library("dendextend")

setwd("~/GIT/darts-heterogeneous/data_img/data/")

data_gpu_streams <- read.csv(file = "data_gpu_streams.csv")
data_gpu_streamNot <- read.csv(file = "data_gpu_streamNot.csv")

# machines <- c("ccsl", "f3", "supermicro", "debian")


# plot(data_gpu_streams$nTile, log(data_gpu_streams$exec_time))

nums <- sapply(data_gpu_streams, is.numeric)
data_gpu_streams <- data_gpu_streams[, nums]
data_gpu_streams[apply(data_gpu_streams, 2, is.infinite)] <- 0
data_gpu_streams[apply(data_gpu_streams, 2, is.na)] <- 0

tempFeatures <- data_gpu_streams[, apply(data_gpu_streams, 2, function(v) var(v, na.rm = TRUE) != 0)]



corFeatures <- cor(normalizeLogMax(getElement(tempFeatures, "exec_time")), apply(tempFeatures, 2, normalizeLogMax),
                   method = "spearman", use = "complete.obs")

tempFeatures$duration <- NULL
corFeatures <- corFeatures[, colnames(corFeatures) != "exec_time"]

    for (threshCorr in c(0.25, .5, .75)){
    
        tempData <- data.frame()
        tempData <- subset(tempFeatures[which(abs(corFeatures) >= threshCorr)])
        # varImp(tempData)
        
        col <- colorRampPalette(c("blue", "yellow", "red"))(20)
        png(filename = paste("./heatMap", "-Thresh_", threshCorr, ".png", sep=""), width = 1600, height = 800)
        heatmap(x = cor(apply(tempData, 2, normalizeLogMax), method = "spearman", use = "complete.obs"), col = col, symm = TRUE)
        dev.off()
        
        # png(filename = paste("./images/phase2/correlation/corClustring_All_App_GPUs", "-Thresh=", threshCorr, ".png", sep=""), width = 1600, height = 800)
        # corrplot(cor(apply(tempData, 2, normalizeLogMax),
        #              method = "spearman", use = "complete.obs"), type = "upper", order = "hclust", hclust.method="average")
        # dev.off()
        
        if(length(tempData) > 10){
        hcFeatures <- hclust(as.dist(1-abs(cor(apply(tempData, 2, normalizeLogMax),
                                       method = "spearman", use = "complete.obs"))), method = "average")
        
        # plot(hcFeatures)
        
        # roc_imp <- filterVarImp(x = tempData, y = tempDuration)
        
        for(numberFeatures in c(5, 10)){
            cutedTree <- cutree(hcFeatures, k = numberFeatures)
            png(filename = paste("./Cluster-Thresh_", threshCorr, "-NParam_", numberFeatures, ".png", sep=""), 
            width = 1600, height = 800)
            
            dend <- as.dendrogram(hcFeatures)
            dend %>% color_branches(k=numberFeatures) %>% plot(horiz=TRUE, 
                                                           main = paste( " Thresh=", 
                                                                         threshCorr, " NParam=", numberFeatures, sep=""))
            
            # add horiz rect
            dend %>% rect.dendrogram(k=numberFeatures,horiz=TRUE)
            # add horiz (well, vertical) line:
            abline(v = heights_per_k.dendrogram(dend)[paste(numberFeatures, sep = "")], 
               lwd = 2, lty = 2, col = "blue")
            # text(50, 50, table(cutedTree))
            dev.off()
            
            parNameTemp <- vector()
            
            for(numberCluster in 1:numberFeatures){
                Tempvariance <-  apply(apply(tempData[cutedTree == numberCluster],2, normalizeLogMax), 2,var)
                parNameTemp[numberCluster] <- names(sort(Tempvariance)[length(Tempvariance)])
            }
            Data <- tempData[parNameTemp]
        }
    }
}