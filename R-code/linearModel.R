library(caret)
df <- data.frame(read.csv("~/GIT/darts-heterogeneous/others/data/p1_p2/p2_more_data/weak_execution.csv"))

df <- df[df$Server_Name %in% c("supermicro", "fatnode", "debian", "hive", "ccsl"),]

# df$Server_Name <- NULL
df$CPU_Type <- NULL
df$GPU_Type <- NULL
df$SP.SM <- NULL
df$Available_GPU_Mem.GB. <- NULL
df$CPU_Mem.GB. <- NULL
# df$GPU_Clock.GHz. <- NULL
# df$PCIe  <- NULL
# df$NumOfCE  <- NULL
# df$CC  <- NULL
# df$GPU_L2_Cache.MB.  <- NULL
# df$GPU_Mem.GB.  <- NULL
# df$NumOfSM  <- NULL
# df$CPU_Threads  <- NULL
# 
# df$CPU.GPU_Workload_Ratio  <- NULL
# df$Initial_GPU_Workload  <- NULL
# df$Total_Workload  <- NULL

# df$Total_Workload <- NULL
# df$Initial_GPU_Workload <- NULL
# df$CPU.GPU_Workload_Ratio <- NULL


x <- df$Exe_Time
h<-hist(x, breaks=10, col="gray", xlab="Execution Time", 
        main="Histogram with Normal Curve of the Execution Times") 
xfit<-seq(min(x),max(x),length=40) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)


df$Exe_Time <- log(df$Exe_Time)


for(host in c("supermicro", "fatnode", "debian", "hive", "ccsl")){
    trainingSet <- df[df$Server_Name != host,]
    testSet <- df[df$Server_Name == host,]
    
    trainingSet$Server_Name <- NULL
    testSet$Server_Name <- NULL
    
    dataReal <- testSet$Exe_Time
    testSet$Exe_Time <- NULL
    
    
    fit.lm <- lm(trainingSet$Exe_Time ~. , data = trainingSet )
    dataPredict <- predict(fit.lm, testSet)
    # plot(fit.lm)
    # alias(fit.lm)
    # print(host)
    # print(summary(fit.lm))
    # print(varImp(fit.lm, scale = TRUE,sort = TRUE,type = 1))
    
    
    qqnorm(rstandard(fit.lm)) #, which=2, caption = "Madre", main = "Quantile-Quantile Plot of the\n execution times", cex.id = 1,qqline = TRUE)
    qqline(rstandard(fit.lm))
    
    qqnorm(rstandard(fit.lm),col="blue", main = paste("Nomal Q-Q of the respose variable", sep="" ))
    qqline(rstandard(fit.lm), col="red")
    
    qqnorm(residuals(fit.lm),col="blue", main = paste("Nomal Q-Q of the respose variable", sep="" ))
    qqline(residuals(fit.lm), col="red")
    
    boxplot(fit.lm$fitted.values/trainingSet$Exe_Time, main = paste("Learning model without ",host, sep="" ))
    plot(fit.lm$fitted.values~trainingSet$Exe_Time, type="l", main = paste("Learning model without ",host, sep="" ))
    

    
    plot(dataPredict, cex.lab = 1.5, cex.axis=2, xlab="Mesaured", ylab="Predited", cex=1.5, type="l", col="blue",
         main = paste("Predicted(blue) VS Measured(red) - ",host, sep="" ), 
         ylim = c(min(dataPredict, dataReal), max(dataPredict, dataReal)))
    lines(dataReal , col="red", cex=2.5)
    
    boxplot(dataPredict/dataReal, main = paste("Learning model - ",host, sep="" ))
    
    print(mape <- mean(abs((as.matrix(dataReal)  - dataPredict)/dataPredict))*100)

}


fit.rf <- randomForest(Exe_Time ~ ., data=df, mtry=5, ntree=100,
                          keep.forest=TRUE, importance=TRUE)

varImpPlot(fit.rf, scale = TRUE,sort = TRUE,type = 1)







fit.svm <- svm(Exe_Time ~., data = df, kernel="linear",scale = TRUE, cost=0.0001)
fit.svm
summary(fit.svm)
qqnorm(residuals(fit.svm))
qqline(residuals(fit.svm))


######################  
createCounter <- function(value) {function(i) {value <<- value + i}}
counter <- createCounter(0)
rss <- matrix(nrow = 15, ncol = 3)

trainFeatures <- df
trainFeatures$Server_Name <- NULL

for (resultY in c("", "log2", "poly")) {
    for (FeaturesX in c("", "log2",  "exp", "sqrt", "poly")) {
        count <- counter(1)
        
        fit <- lm(as.formula(
            paste(resultY, "(", "Exe_Time",  ifelse(resultY == "log2", "+ 0.0000000001", ""),") ~",
                  paste(FeaturesX, "(",ifelse(FeaturesX == "exp", "normalizeMax(",""),  
                        colnames(trainFeatures[!names(trainFeatures) == "Exe_Time"]), 
                        ifelse(FeaturesX == "poly", ",2,raw = TRUE",""),
                        ifelse(FeaturesX == "log2", "+ 0.0000000001", ""),
                        ifelse(FeaturesX == "exp", ")", "") , ")",
                        collapse = " + " ))), data =  trainFeatures)
        
        print(count)
        print(varImp(fit, scale = TRUE))
        
        rss[count, 1] <- summary(fit)$r.squared
        rss[count, 2] <- resultY
        rss[count, 3] <- FeaturesX
        
        qqnorm(
            residuals(fit),
            ylab = "Studentized Residual",
            xlab = "t Quantiles",
            main = paste(" Model: ", resultY, "(y) = ", FeaturesX,"X_i", sep = ""),
            cex.lab = 2,
            cex.main = 2,
            cex = 1.5,
            cex.axis = 2
        )
        qqline(residuals(fit), col = 2, lwd = 5)
    }
}

rss


summary(fit.lm)
summary(fit.lm)$r.squared

qqnorm(residuals(fit.lm))
qqline(residuals(fit.lm))

plot(importance(mtcars.rf, scale = TRUE))
varImpPlot(mtcars.rf, scale = FALSE,sort = TRUE)




as.dist(1 - abs(cor(df, method = "spearman", use = "complete.obs")))

plot(hclust(as.dist(1 - abs(cor(df, method = "spearman", use = "complete.obs"))), method = "average"))



nSocket = 2;
CPU_clock = 3.4;
CPU_cores =32;
CPU_L3_Cache=25;


for( wl in seq(17000, 49000, 2000)){
    pred = predict_best_wl(wl,nSocket,CPU_clock,CPU_cores,CPU_L3_Cache);
    if (wl == 17000){
        m2 = pred;
    }else{
        m2= cat(" ",m2,pred);
    }
}

# predict_best_wl(17000,nSocket,CPU_clock,CPU_cores,CPU_L3_Cache)

predict_best_wl <- function(wl,nSocket,CPU_clock,CPU_cores,CPU_L3_Cache){
    min_exe_time = Inf;
    best_g_wl = Inf;
    best_cg_ratio=Inf;
    for(g_wl in seq(1000, 8000, 1000)){
        for(cg_ratio in seq(0.25, 2, 0.25)){
            execution_time = 2^((-4.280e+00)*nSocket  - (1.093e+00) *CPU_clock + (7.307e-01)*CPU_cores -(5.239e-01)*CPU_L3_Cache + (8.458e-05)*wl + (3.619e-05)*g_wl -(9.041e-03)*cg_ratio + (3.094e+01))
            if(execution_time < min_exe_time){
                min_exe_time = execution_time;
                best_g_wl = g_wl;
                best_cg_ratio = cg_ratio;
            }
        }
    }
    pred = c(min_exe_time, best_g_wl, best_cg_ratio)
}



