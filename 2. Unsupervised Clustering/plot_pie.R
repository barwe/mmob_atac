df = data.frame(
    ct = c("IN", "EX", "OEC", "End", "AST", "OSN", "Purk", "OPC", "OLI", "Gran", "MG"),
    vl = c(61,10,5,3,4,2,2,1,2,7,3)/100,
    stringsAsFactors = F
)
rownames(df) = df$ct
df = df[c('IN','Purk','EX','Gran','OEC','OPC','OLI','AST','MG','End'),]
df$ct = factor(df$ct, levels = c('IN','Purk','EX','Gran','OEC','OPC','OLI','AST','MG','End'))

library(ggplot2)

ggplot(df, aes(x="", y = vl, fill = ct))+
    geom_bar(stat = "identity", width = 1)+
    coord_polar(theta = "y")+
    labs(x = "", y = "", title = "") + 
    scale_fill_discrete(breaks=df$ct, labels=paste0(df$ct,' (',df$vl*100,'%)'))+
    theme(
        panel.background = element_blank(),
        axis.ticks = element_blank())
ggsave("各个cluster细胞比例的饼图.pdf", width=6, height=4)
