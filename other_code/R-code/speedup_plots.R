library(ggplot2)
setwd("/home/marcos/GIT/darts-heterogeneous/others/data/weak/")

cbbPalette <- c("blue", "#E69F00", "red", "#CC79A7", "", "#0072B2", "#D55E00", "#CC79A7")

df <- data.frame()
for(j in c("ccsl", "debian", "f4", "hive", "supermicro")){
    if(j == "ccsl") threads <- 7
    if(j == "debian") threads <- 11
    if(j == "f4") threads <- 31
    if(j == "hive") threads <- 11
    if(j == "supermicro") threads <- 39       
    
    
    
    data <- read.csv(paste("./", j, "_", threads, "_1_weak_speedup.dat", sep=""), header = T, sep = ",")
    
    names(data) <- c("size", "Sequential", "CUDA", "DartsCpu", "DartsGpu", "DartsHybrid3")
    
    
    if(j == "f4") j <- "Fatnode"
    dataFrame <- data.frame(size=data$size, machine=j, apps=c(rep(names(data)[2], dim(data)[1]),
                                                                rep(names(data)[3], dim(data)[1]),
                                                                rep(names(data)[4], dim(data)[1]),
                                                                rep(names(data)[5], dim(data)[1]),
                                                                rep(names(data)[6], dim(data)[1])),
                            speedup=c(data[,2], data[,3], data[,4], data[,5], data[,6]))
    
    df <- rbind(df, dataFrame)
}

df <- df[df$apps != "Sequential",]

Graph <- ggplot(data=df, aes(x=size, y=speedup, group=apps, col=apps, pch=apps, linetype = apps)) + 
    geom_line(size=1.5) +
    geom_point(cex=3.5) +
    xlab("Size of the Problem") + 
    theme_bw() +
    scale_colour_manual(values=cbbPalette) +
    ylab("Speedup" ) +
    theme(plot.title = element_text(family = "Times", face="bold", size=40)) +
    theme(axis.title = element_text(family = "Times", face="bold", size=30)) +
    theme(axis.text  = element_text(family = "Times", face="bold", size=15, colour = "Black")) +
    scale_x_continuous(breaks=seq(0,50000,5000)) +
    theme(axis.text.x= element_text(family = "Times", face="bold", size=15, colour = "Black", angle=0, hjust=1)) +
    theme(legend.title  = element_text(family = "Times", face="bold", size=0)) +
    theme(legend.text  = element_text(family = "Times", face="bold", size=18)) +
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.key=element_rect(size=5),
          legend.key.size = unit(4, "lines")) +
    guides(col = guide_legend(nrow = 1)) +
     # facet_grid(.~machine, scales="free") +
    facet_wrap(~machine, ncol=1, scales="free_x") +
    theme(strip.text = element_text(size=20))
ggsave(paste("./speedUp.pdf",sep=""), Graph, device = pdf, height=18, width=9)



df <-  df[df$size >= 17000,]
Graph <- ggplot(data=df, aes(x=size, y=speedup, group=apps, col=apps, pch=apps, linetype = apps)) + 
    geom_line(size=1.5)+
    geom_point(cex=3.5) +
    xlab("Size of the Problem") + 
    theme_bw() +
    scale_colour_manual(values=cbbPalette) +
    ylab("Speedup" ) +
    theme(plot.title = element_text(family = "Times", face="bold", size=40)) +
    theme(axis.title = element_text(family = "Times", face="bold", size=30)) +
    theme(axis.text  = element_text(family = "Times", face="bold", size=15, colour = "Black")) +
    scale_x_continuous(breaks=seq(0,50000,5000)) +
    theme(axis.text.x= element_text(family = "Times", face="bold", size=15, colour = "Black", angle=0, hjust=1)) +
    theme(legend.title  = element_text(family = "Times", face="bold", size=0)) +
    theme(legend.text  = element_text(family = "Times", face="bold", size=18)) +
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.key=element_rect(size=5),
          legend.key.size = unit(4, "lines")) +
    guides(col = guide_legend(nrow = 1)) +
    # facet_grid(.~machine, scales="free") +
    facet_wrap(~machine, ncol=1, scales="free_x") +
    theme(strip.text = element_text(size=20))
ggsave(paste("./speedUpZoom.pdf",sep=""), Graph, device = pdf, height=18, width=9)

