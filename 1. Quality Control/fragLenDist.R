library(ggplot2)

## Please set your project directory
PROJECT_DIR = "/%PATH%/mouse_olfactory_bulb_atac"

df = read.table("fragLenStat.txt")
colnames(df) = c("length", "count")

ggplot(df, aes(x=length, y=count))+
    geom_bar(stat = 'identity', fill = '#a3ceec')+
    theme(
        panel.background = element_blank(),
    )

ggsave("fragLenDist.pdf", width = 12, height = 8)
