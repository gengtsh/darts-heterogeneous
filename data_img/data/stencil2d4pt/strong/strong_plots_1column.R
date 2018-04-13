library(ggplot2)
library(stats)

#setwd("/home/marcos/GIT/darts-heterogeneous/others/data/strong/")
cbbPalette <- c("#000000", "#F0E442", "#56B4E9", "#009E73", "#E69F00", "#0072B2", "#D55E00", "#CC79A7")

df <- data.frame()
for(machine in c("ccsl", "debian", "f4", "hive", "supermicro")){
    for (sched in c("scatter", "compact")){
        for(app in c("StencilCUDA", "StencilCudaCpu", "StencilCudaHybrid3")){
            data <- read.csv(paste("./", machine, "_strong_G1C1_G0.5M_allsize_", sched, "/", app, "_1_strong_execution.dat", sep=""), header = T, sep = ",",stringsAsFactors = FALSE)
            
            temp <- data.frame(cbind(rep(4:max(as.numeric(data[,1])),1), app,  ifelse(machine == "f4", "Fatnode", machine), sched,                          
                                     rep(seq(from=1000, to=((dim(data)[2]-1)*2000)-1000, by = 2000), each=max(as.numeric(data[,1]))-3), 
                                     unlist(data[,2:dim(data)[2]],use.names = FALSE)))
            
            names(temp) <- c("Threads", "App", "machine", "sched", "size", "time")
            df <- rbind(df, temp)
        }
    }
}

names(df)

df$Threads <- as.numeric(as.character(df$Threads))
df$size <- as.character(df$size)
df$time <- as.numeric(as.character(df$time))
df$sched <- as.character(df$sched)
df$App <- as.character(df$App)

df <- df[df$size %in% c( "11000", "17000", "25000", "31000", "35000"),]
df <- df[df$App %in% c("StencilCudaHybrid3", "StencilCudaCpu"),]

index <- df$App == "StencilCudaHybrid3"
df$App[index] <- "EDRT-DAWL"

index <- df$App == "StencilCudaCpu"
df$App[index] <- "EDRT-CPU"


dfTemp <- df[df$machine %in% c("hive"),]
Graph <- ggplot(data=dfTemp, aes(x=Threads, y=time, group=size, col=size, pch=size, linetype = size)) + 
    geom_line(size=1.5) +
    geom_point(cex=3) + 
    xlab("Number of Threads") + 
    theme_bw() +
    ggtitle("Machine Hive") +
    theme(plot.title=element_text(hjust=0.5)) +
    # scale_colour_manual(values=cbbPalette) +
    scale_colour_grey() +
    ylab("Speedup(baseline: CPU-Sequential)" ) +
    theme(plot.title = element_text(family = "Times", face="bold", size=20)) +
    theme(axis.title = element_text(family = "Times", face="bold", size=15)) +
    theme(axis.text  = element_text(family = "Times", face="bold", size=15, colour = "Black")) +
    scale_x_continuous(breaks=seq(4,40,4)) +
    theme(axis.text.x= element_text(family = "Times", face="bold", size=15, colour = "Black", angle=0, hjust=1)) +
    theme(legend.title  = element_text(family = "Times", face="bold", size=0)) +
    theme(legend.text  = element_text(family = "Times", face="bold", size=15)) +
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.key=element_rect(size=0.1),
          legend.key.size = unit(4, "lines")) +
    guides(col = guide_legend(nrow = 1)) +
    facet_grid(sched~App, scales="free") +
    # facet_wrap(sched~App, scales="free_x",nrow = 2) +
    theme(strip.text = element_text(size=20))
ggsave(paste("./speedUp-strong-Hive.pdf",sep=""), Graph, device = pdf, height=6, width=9)
    

dfTemp <- df[df$machine %in% c("Fatnode"),]
Graph <- ggplot(data=dfTemp, aes(x=Threads, y=time, group=size, col=size, pch=size, linetype = size)) + 
    geom_line(size=1.5) +
    geom_point(cex=1.5) + 
    xlab("Number of Threads") + 
    theme_bw() +
    ggtitle("Machine Fatnode") +
    theme(plot.title=element_text(hjust=0.5)) +
    # scale_colour_manual(values=cbbPalette) +
    scale_colour_grey() +
    ylab("Speedup(baseline: CPU-Sequential)" ) +
    theme(plot.title = element_text(family = "Times", face="bold", size=20)) +
    theme(axis.title = element_text(family = "Times", face="bold", size=15)) +
    theme(axis.text  = element_text(family = "Times", face="bold", size=15, colour = "Black")) +
    scale_x_continuous(breaks=seq(4,40,4)) +
    theme(axis.text.x= element_text(family = "Times", face="bold", size=15, colour = "Black", angle=0, hjust=1)) +
    theme(legend.title  = element_text(family = "Times", face="bold", size=0)) +
    theme(legend.text  = element_text(family = "Times", face="bold", size=15)) +
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.key=element_rect(size=1),
          legend.key.size = unit(5, "lines")) +
    guides(col = guide_legend(nrow = 1)) +
    facet_grid(sched~App, scales="free") +
    # facet_wrap(sched~App, scales="free_x",nrow = 2) +
    theme(strip.text = element_text(size=15))
ggsave(paste("./speedUp-strong-Fatnode.pdf",sep=""), Graph, device = pdf, height=6, width=9)
