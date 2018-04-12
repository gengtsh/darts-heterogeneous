df <- data.frame(read.csv("~/GIT/darts-heterogeneous/others/data/p1_p2/p2_more_data/weak_execution.csv"))

df <- df[df$Server_Name %in% c("supermicro", "fatnode", "debian", "hive", "ccsl"),]



df$Server_Name <- NULL
df$CPU_Type <- NULL
df$GPU_Type <- NULL
df$SP.SM <- NULL
df$Available_GPU_Mem.GB. <- NULL
df$CPU_Mem.GB. <- NULL
df$GPU_Clock.GHz. <- NULL
df$PCIe  <- NULL
df$NumOfCE  <- NULL
df$CC  <- NULL
df$GPU_L2_Cache.MB.  <- NULL
df$GPU_Mem.GB.  <- NULL
df$NumOfSM  <- NULL
df$CPU_Threads  <- NULL

df$CPU.GPU_Workload_Ratio  <- NULL
df$Initial_GPU_Workload  <- NULL
df$Total_Workload  <- NULL

# df$Total_Workload <- NULL
# df$Initial_GPU_Workload <- NULL
# df$CPU.GPU_Workload_Ratio <- NULL


trainingSet <- df[df$CPU_Threads != 32,]
testSet <- df[df$CPU_Threads == 32,]

dataReal <- testSet$Exe_Time
testSet$Exe_Time <- NULL


fit.lm <- lm(Exe_Time ~., data = trainingSet )
summary(fit.lm)
qqnorm(residuals(fit.lm))
qqline(residuals(fit.lm))

boxplot(fit.lm$fitted.values/trainingSet$Exe_Time)
plot(fit.lm$fitted.values~trainingSet$Exe_Time, type="l")


dataPredict <- predict(fit.lm, testSet)
plot(dataPredict, cex.lab = 1.5, cex.axis=2, xlab="Mesaured", ylab="Predited", cex=1.5, type="l", col="blue")
lines(dataReal , col="red", cex=2.5)


varImp(fit.lm, scale = TRUE,sort = TRUE,type = 1)



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

