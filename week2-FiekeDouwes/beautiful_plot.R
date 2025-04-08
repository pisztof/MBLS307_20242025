library("tidyverse")

mpg_df <- mpg

png("output/beautiful_plot.png")

ggplot(mpg_df) + 
  geom_boxplot(aes(y=displ, x = class), outlier.shape = NA) +
  geom_jitter(aes(y=displ, x= class, color = manufacturer))+
  labs(x = "Class of the car", y = "Engine size, L")+
  lims(y = c(0,10))+
  theme_classic()
dev.off()
