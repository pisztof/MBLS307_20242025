library(ggplot2)

ggplot(mpg, aes(x=class,y=displ)) + 
  geom_boxplot(outliers=FALSE) +
  geom_jitter(aes(color=manufacturer)) +
  ylim(c(0,10)) +
  xlab("Class of the car") +
  ylab("Engine size, L") +
  theme_classic()

ggsave("Plots/ggplot_exercise.png")