library("ggpubr")
library("caret")
library("ggplot2")
library("plyr")

df_supermicro <- data.frame(read.csv("~/GIT/darts-heterogeneous/data_img/data/stencil3d7pt/profile/supermicro/nvprofTraces.csv"))
df_ccsl <- data.frame(read.csv("~/GIT/darts-heterogeneous/data_img/data/stencil3d7pt/profile/ccsl/nvprofTraces.csv"))
df_debian <-     data.frame(read.csv("~/GIT/darts-heterogeneous/data_img/data/stencil3d7pt/profile/debian/nvprofTraces.csv"))
df_f3 <-     data.frame(read.csv("~/GIT/darts-heterogeneous/data_img/data/stencil3d7pt/profile/f3/nvprofTraces.csv"))

df <- rbind(df_supermicro, df_ccsl, df_debian,df_f3)


Traces <- c("cuDeviceGetAttribute", "cudaSetupArgument", "cudaMallocHost", "cudaGetDeviceCount", "cudaGetDeviceProperties", "cudaMemGetInfo", 
            "cudaLaunch", "cudaMalloc", "cudaMemcpy", "CUDA memcpy HtoD", "CUDA memcpy DtoH", "cudaFree", "\"gpu_stencil37_hack1_cp_slices", 
            "cudaLaunch (gpu_stencil37_hack1_cp_rows", "gpu_stencil37_hack1_cp_rows", "cudaLaunch (gpu_stencil37_hack1_cp_cols", "gpu_stencil37_hack1_cp_cols",
            "cudaLaunch (gpu_stencil37_hack2", "gpu_stencil37_hack2")

ZXY <- c("200_200_50", "200_200_100", "200_200_200", "400_400_800", "800_800_200", "800_800_400", "800_800_800")

df$size <- as.character(df$size)
df$size <- factor(df$size, levels = c(" 200_200_50", " 200_200_100", " 200_200_200", " 400_400_800", " 800_800_200", " 800_800_400", " 800_800_800"))

df$metric <- as.character(df$metric)
df$mean <- df$mean/10^9

df_temp <- df[df$metric %in% c(" CUDA memcpy HtoD", " CUDA memcpy DtoH", " cudaLaunch", " cudaMemcpy"),]

Graph1 <- ggplot(data=df_temp, aes(x=size, y=mean, group=machine, col=machine)) + 
    geom_line(size=1.25 , aes(linetype = machine)) + 
    geom_point(size=2, aes(shape = machine)) +  
    theme_bw() +
    ylab("Mean (sec)") + 
    theme(legend.position="none") +
    theme(legend.direction = "horizontal",
    legend.position = "bottom",
    legend.key=element_rect(size=0),
    legend.key.size = unit(3, "lines")) +
    facet_wrap(~metric, ncol=1,  scales="free_y") +
    theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
    theme(legend.title=element_blank())+
    theme(strip.text = element_text(size=14, family = "Times", face="bold")) +
    scale_colour_grey()
ggsave(paste("~/GIT/darts-heterogeneous/data_img/Img/tracesPlots-Memcpy", ".png" ,sep=""), Graph1, height=12, width=6)


df_temp <- df[df$metric %in% c(" gpu_stencil37_hack1_cp_slices", " gpu_stencil37_hack1_cp_rows", " gpu_stencil37_hack1_cp_cols", " gpu_stencil37_hack2"),]

Graph2 <- ggplot(data=df_temp, aes(x=size, y=mean, group=machine, col=machine)) + 
    geom_line(size=1.25,aes(linetype = machine)) + 
    geom_point(size=2, aes(shape = machine)) + 
    theme_bw() +
    theme(legend.position="none") +
    ylab("Mean (sec)") + 
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key=element_rect(size=0),
          legend.key.size = unit(3, "lines")) +
    facet_wrap(~metric, ncol=1, scales="fixed") +
    theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
    scale_y_log10()  +
    theme(legend.title=element_blank()) +
    theme(strip.text = element_text(size=14, family = "Times", face="bold")) +
    scale_colour_grey()
ggsave(paste("~/GIT/darts-heterogeneous/data_img/Img/tracesPlots-kernels", ".png" ,sep=""), Graph2, height=12, width=6)

pdf("~/GIT/darts-heterogeneous/data_img/Img/plot-Memcpy-kernels.pdf", width = 10, height = 14)
ggarrange(Graph1, Graph2, nrow=1, ncol=2, labels=c(LETTERS[1:2]), legend="bottom")
dev.off()

png("~/GIT/darts-heterogeneous/data_img/Img/plot-Memcpy-kernels.png", width = 600, height = 1000, units="px")
ggarrange(Graph1, Graph2, nrow=1, ncol=2, labels=c(LETTERS[1:2]), legend="bottom")
dev.off()
