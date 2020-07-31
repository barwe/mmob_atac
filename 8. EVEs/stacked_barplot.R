width.df = element.types
local({
    width.df = width.df[order(rowSums(width.df)),]
    ggdf = data.frame(
        XN = rep(nrow(width.df):1, each = ncol(width.df)),
        XS = rep(rev(rownames(width.df)), each = ncol(width.df)),
        L = factor(rep(rev(colnames(width.df)), nrow(width.df)), levels = rev(colnames(width.df))),
        Y = rev(c(t(width.df))), stringsAsFactors = F
    )
    ggdf$Y = ifelse(ggdf$Y==0,NA,ggdf$Y)
    x_labels = unique(ggdf$XS)[order(unique(ggdf$XN))]
    ggplot(ggdf, aes(x = XN, y = Y, fill = L)) + 
        geom_bar(stat = 'identity', position = 'stack') +
        coord_flip() + 
        labs(x = '', y = '', title = '') +
        #geom_text(aes(label = LABEL),size = 3,position = position_stack()) +
        scale_x_continuous(breaks = 1:21, labels = x_labels) +
        scale_y_continuous(breaks = 1, labels = '') +
        guides(fill = guide_legend(reverse = TRUE)) +
        theme(
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            legend.position = 'right',
            legend.title = element_blank()
        )
})
ggsave("E:/project/MMOB/[07]EVE/stacked_barplot/stacked_barplot.pdf", width = 10, height = 8)
