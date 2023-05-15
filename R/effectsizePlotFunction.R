plot_DataEFFECT_vs_RefEFFECT<-function(input.data, effPlotPath, plot.title.text)
{


  eff.ggplot <- ggplot(data = input.data, aes(x = EFFECT.x, y = EFFECT.y,colour = as.factor(HQ))) +
    geom_point( size = .8) +
    #geom_rug() +
    geom_smooth(method=lm, aes(fill= as.factor(HQ)),fullrange=TRUE)+
    labs(title="Effect-Size comparison plot",
         x="reference effect-size",
         y="reported effect-size",
         subtitle= plot.title.text) +
    theme_classic() +
    #  scale_x_continuous(limits=c(-1,1)) +
    # scale_y_continuous(limits=c(-1,1)) +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size=8, face = "bold")
          ,strip.text.y = element_text(size=8, face = "bold")
          ,legend.position = 'none'
          ,plot.title=element_text(size=10, face="bold",hjust = 0.5)
          ,plot.subtitle = element_text(size=8, hjust=0.5, face="italic", color="darkblue")
    ) +
    annotate("text",
             x = min(input.data$EFFECT.x) * 0.9,
             y = max(input.data$EFFECT.y) * 0.95,
             # x=Inf,
             # y=-Inf,
             # hjust=0,
             # vjust=0,
             label=sprintf('italic(r) == %.3f' ,  .QC$thisStudy$effect.rho_4),
             parse= TRUE)+
    scale_color_manual(values=c('red','blue'))



  ggsave(plot=eff.ggplot,
         device = .QC$graphic.device,
         filename = effPlotPath,
         units = c('mm'),
         width = 180,
         height =120,
         dpi = 300)

  print_and_log(sprintf('Effect-Size comparison plot saved! %s variants.',nrow(input.data)),
                'info')

  rm(eff.ggplot)
  invisible(gc())
}
