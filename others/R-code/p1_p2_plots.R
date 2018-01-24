
setwd("/home/marcos/Doctorate/GIT/darts-heterogeneous/others/data/")

for(i in c("p1", "p2")){
    for(j in c("ccss", "debian", "f4", "hive", "supermicro")){
        if(i == "p1"){
            png(filename = paste("./images/",i, "-", j, ".png", sep = ""), 
                width = 1200, height = 800)
            par(mfrow = c(4, 2), oma=c(0,0,2,0))
            for(k in c(1000, 3000, 5000, 25000)){
                data <- read.csv(paste("./p1_p2/", j,"_p1/StencilCudaCpu_", k, ".dat", sep=""), header = F, sep = ",")
                dim(data)
                
                plot(rowMeans(data[4:93])~data[,1],xlab = "No. Threads", ylab="nanosecond", main=paste("Size ", k, sep = ""),
                     cex.lab = 2, cex.main = 2, cex = 1.5, cex.axis = 2)
                plot((rowMeans(data)/rowMeans(data[1,]))~data[,1],xlab = "No. Threads", ylab="speedup", main=paste("Size ", k, sep = ""), 
                     cex.lab = 2, cex.main = 2, cex = 2.5, cex.axis = 2)
            }
            title(paste(i, " ", j), outer=TRUE,cex.main=2)
            dev.off()
        }
        if(i == "p2"){
            if(j == "ccss") threads <- 7
            if(j == "debian") threads <- 11
            if(j == "f4") threads <- 31
            if(j == "hive") threads <- 11
            if(j == "supermicro") threads <- 39              
            png(filename = paste("./images/",i, "-", j, ".png", sep = ""), 
                width = 600, height = 1800)
            par(mfrow = c(8, 2), oma=c(0,0,2,0))
            for(k in c(17000, 25000)){
                for(l in c(1000, 2000, 4000, 8000)){
                    data <- read.csv(paste("./p1_p2/", j,"_p2/StencilCudaHybrid4_",k, "_", l, "_", threads, "_1.dat", sep=""), header = F, sep = ",")
                    dim(data)
                    
                    plot(rowMeans(data[4:93])~data[,5],xlab = "No. Samples", ylab="nanosecond", main=paste("Size ", k, "GPU rows=", l, sep = ""),cex.lab = 2,
                         cex.main = 2,
                         cex = 1.5,
                         cex.axis = 2)
                    plot((rowMeans(data)/rowMeans(data[1,]))~data[,5],xlab = "No. Samples", ylab="speedup", main=paste("Size ", k, "GPU rows=", l, sep = ""), cex.lab = 2,
                         cex.main = 2,
                         cex = 1.5,
                         cex.axis = 2)
                }
            }
            title(paste(i, " ", j), outer=TRUE, cex.main=2)
            dev.off()
        }
    }
}

