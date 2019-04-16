library(ggplot2)
#setwd("/home/marcos/GIT/darts-heterogeneous/data_img/data/stencil2d4pt/weak/threads_num_equals_to_machines_threads_G2000C2000_2/")

#cbbPalette <- gray.colors(5, start = 0.1, end = 0.9)
cbbPalette <- gray.colors(6, start = 0.1, end = 0.7)

#setwd("~/GIT/darts-heterogeneous/data_img/data/stencil2d4pt/weak/threads_number_equals_to_machines_threads_w_static/f4_weak_2GB_c")

temp <- read.csv("31_1_weak_speedup.csv",header = TRUE)

#speedup <- c(temp$StencilCudaHybrid3S12, temp$StencilCudaHybrid3S21,
#            temp$StencilCudaHybrid3S22,temp$StencilCudaHybrid3S24, temp$StencilCudaHybrid3S44, temp$StencilCudaHybrid3_IDAWL)
speedup <- c(
             temp$StencilCudaHybrid3S22,temp$StencilCudaHybrid3S24, temp$StencilCudaHybrid3S42,temp$StencilCudaHybrid3S44, temp$StencilCudaHybrid3_IDAWL)

size <- c(temp$X.size, temp$X.size, temp$X.size,
          temp$X.size, temp$X.size) 

apps <- c(
          rep("DARTS_Static_22", 25),rep("DARTS_Static_24", 25),rep("DARTS_Static_42", 25),rep("DARTS_Static_44", 25), rep("DARTS_IDAWL",25))

df <- data.frame(cbind(size, apps, speedup))

df <- df[df$apps != "StencilSEQ",]
df$speedup <- as.numeric(as.character(df$speedup))
df$size <- as.numeric(as.character(df$size))
Graph <- ggplot(data=df, aes(x=size, y=speedup, group=apps, col=apps, pch=apps, linetype = apps)) + 
    geom_line(size=1.5) +
    geom_point(cex=3.5) +
    xlab("Size of the Problem") + 
    theme_bw() +
    scale_colour_manual(values=cbbPalette) +
    ylab("Speedup (baseline = CPU-Sequential)" ) +
    theme(plot.title = element_text(family = "Times", face="bold", size=40)) +
    theme(axis.title = element_text(family = "Times", face="bold", size=20)) +
    theme(axis.text  = element_text(family = "Times", face="bold", size=15, colour = "Black")) +
    theme(axis.text.x= element_text(family = "Times", face="bold", size=15, colour = "Black", angle=0, hjust=1)) +
    theme(legend.title  = element_text(family = "Times", face="bold", size=0)) +
    theme(legend.text  = element_text(family = "Times", face="bold", size=16)) +
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.key=element_rect(size=5),
          legend.key.size = unit(3, "lines")) +
    guides(col = guide_legend(nrow = 2))

#ggsave(paste("./speedUp_dynamic.pdf",sep=""), Graph, device = pdf, height=7, width=10)
ggsave("speedUp_dynamic.pdf", height=7, width=10)
