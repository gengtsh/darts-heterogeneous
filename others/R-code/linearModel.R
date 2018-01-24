df <- data.frame(read.csv("~/GIT/darts-heterogeneous/others/data/p1_p2/p2/weak_execution.csv"))

df$Server_Name <- NULL
df$CPU_Type <- NULL
df$GPU_Type <- NULL


trainingSet <- df[df$CPU_Threads != 32,]
testSet <- df[df$CPU_Threads == 32,]

dataReal <- testSet$Initial_GPU_Workload
testSet$Initial_GPU_Workload <-  NULL

fit <- lm(Initial_GPU_Workload~., data = df )

dataPredict <- predict(fit, testSet)

plot(dataReal ~ dataPredict, cex.lab = 1.5, cex.axis=2, xlab="Mesaured", ylab="Predited", cex=1.5)




str(summary(fit)) 
summary(fit)$r.squared


qqnorm(residuals(fit))
qqline(residuals(fit))




hclust(df)
