install.packages("ggplot2")
library(ggplot2)
install.packages("tidyverse")
library(tidyverse)

mpg_df <- mpg
epic_figure <- ggplot(aes(y=displ, x=class),data = mpg_df)+
  geom_boxplot(outliers=FALSE)+
  geom_jitter(aes(color=manufacturer), width=0.3)+
  xlab("Class of car") +
  ylab("Engine size (L)") +
  theme_classic()+
  ylim(0,10)

epic_figure


ggsave(filename = "exercise_plot_making.pdf", plot = epic_figure, device = "pdf", 
       path = "data", width = 250, height = 150, units = "mm", 
       dpi = 500, limitsize = TRUE)

mpg_df %>% select(displ) %>% pull() %>% mean()
