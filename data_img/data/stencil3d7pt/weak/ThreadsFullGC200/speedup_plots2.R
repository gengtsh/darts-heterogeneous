library(ggplot2)
library(stats)
#setwd("/home/marcos/GIT/darts-heterogeneous/others/data/weak/")
currPath = getwd()
setwd(currPath)

# cbbPalette <- c("#000000", "#F0E442", "#56B4E9", "#009E73", "#E69F00", "#0072B2", "#D55E00", "#CC79A7")
# cbbPalette <- gray(1:4/ 12)#c("red", "blue", "darkgray", "orange","black","brown", "lightblue","violet")

df <- data.frame()
for(j in c("ccsl","debian","f4", "supermicro")){
#for(j in c("ccsl", "debian", "f3", "supermicro")){
    if(j == "ccsl") threads <- 7
    if(j == "debian") threads <- 11
    if(j == "f4") threads <- 31
    if(j == "hive") threads <- 11
    if(j == "supermicro") threads <- 39       
    
    
    
    #data <- read.csv(paste("./", j, "_", threads, "_1_weak_speedup.dat", sep=""), header = T, sep = ",")
    data <- read.csv(paste("./", j,"_weak_speedup.dat", sep=""), header = T, sep = ",")
    names(data) <- c("size", "CPU-Sequential", "GPU-only", "DARTS-CPU", "DARTS-GPU", "DARTS-DAWL")
    
    
    if(j == "f4") j <- "Fatnode"
    dataFrame <- data.frame(size=data$size, machine=j, apps=c(rep(names(data)[2], dim(data)[1]),
                                                                rep(names(data)[3], dim(data)[1]),
                                                                rep(names(data)[4], dim(data)[1]),
                                                                rep(names(data)[5], dim(data)[1]),
                                                                rep(names(data)[6], dim(data)[1])),
                            speedup=c(data[,2], data[,3], data[,4], data[,5], data[,6]))
    
    df <- rbind(df, dataFrame)
}


df <- df[df$apps != "CPU-Sequential",]

df$size <- factor(df$size, levels = c("50*200*200", "100*200*200", "200*200*200", "200*800*800", 
                                      "800*400*400", "400*800*800",    "800*800*800", "800*1000*1000", "1000*1000*1000"))

df$apps <- factor(df$apps, levels = c("CPU-Sequential", "DARTS-DAWL", "DARTS-CPU", "DARTS-GPU", "GPU-only"))

Graph <- ggplot(data=df, aes(x=size, y=speedup, group=apps, col=apps, pch=apps, linetype = apps)) + 
  geom_line(size=5)+
  geom_point(cex=7) +
  scale_linetype_manual(values=c("dotted", "longdash","solid",  "twodash" ))+
  xlab("\n Size of the Problem") + 
  theme_bw() +
#    scale_colour_grey() +
  ylab("Speedup(baseline = CPU-Sequential)" ) +
  theme(plot.title = element_text(family = "Times", face="bold", size=50)) +
  theme(axis.title = element_text(family = "Times", face="bold", size=50)) +
  theme(axis.text  = element_text(family = "Times", face="bold", size=50, colour = "Black")) +
#  scale_x_continuous(breaks=seq(0,50000,5000)) +
  theme(axis.text.x= element_text(family = "Times", face="bold", size=30, colour = "Black", angle=15, hjust=0.85)) +
  theme(legend.title  = element_text(family = "Times", face="bold", size=0)) +
  theme(legend.text  = element_text(family = "Times", face="bold", size=50)) +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.key=element_rect(size=50),
        legend.key.size = unit(10, "lines")) +
  guides(col = guide_legend(nrow = 10)) +
  # facet_grid(.~machine, scales="free") +
  facet_wrap(~machine, ncol=2, scales="free_x", nrow=2) +
  theme(strip.text = element_text(size=50))
ggsave(paste("./speedUp-Stencil3D2.pdf",sep=""), Graph, device = pdf, height=30, width=40)



#df <-  df[df$size >= 17000,]
#Graph <- ggplot(data=df, aes(x=size, y=speedup, group=apps, col=apps, pch=apps, linetype = apps)) + 
#  geom_line(size=1.5)+
#  geom_point(cex=3.5) +
#  xlab("Size of the Problem") + 
#  theme_bw() +
#    scale_colour_grey() +
#  ylab("Speedup(baseline = CPU-Sequence)" ) +
#  theme(plot.title = element_text(family = "Times", face="bold", size=40)) +
#  theme(axis.title = element_text(family = "Times", face="bold", size=30)) +
#  theme(axis.text  = element_text(family = "Times", face="bold", size=15, colour = "Black")) +
#  scale_x_continuous(breaks=seq(0,50000,5000)) +
#  theme(axis.text.x= element_text(family = "Times", face="bold", size=15, colour = "Black", angle=0, hjust=1)) +
#  theme(legend.title  = element_text(family = "Times", face="bold", size=0)) +
#  theme(legend.text  = element_text(family = "Times", face="bold", size=16)) +
#  theme(legend.direction = "horizontal", 
#        legend.position = "bottom",
#        legend.key=element_rect(size=5),
#        legend.key.size = unit(3, "lines")) +
#  guides(col = guide_legend(nrow = 1)) +
#  #facet_grid(.~machine, scales="free") +
#  #facet_wrap(~machine, ncol=1, scales="free_x") +
#  facet_wrap(~machine, ncol=1, scales="free") +
#  theme(strip.text = element_text(size=20))
#ggsave(paste("./speedUpZoom.pdf",sep=""), Graph, device = pdf, height=18, width=9)#