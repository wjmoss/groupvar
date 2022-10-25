library(ggplot2)

my.theme <- theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=11),
        axis.title = element_text(size=15),
        legend.position="bottom",
        legend.title = element_text(face = "bold", size=15),
        legend.text = element_text(size=13))

