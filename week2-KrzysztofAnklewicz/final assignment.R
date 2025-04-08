library(tidyverse)
library(ggplot2)
library(readxl)

# please note that I've taken the liberty to replace all NA values with blanks in Excel. This could have been handled in the code, had I had the patience.
data_tumour <- read_xlsx("proteomicsdata.xlsx", sheet="Tumor_Proteome_SISII_Log2Abunda")
data_tcell <- read_xlsx("proteomicsdata.xlsx", sheet="Tcell_Proteome_SISII_Log2Abund")
data_newsynth <- read_xlsx("proteomicsdata.xlsx", sheet="NewSynth_Proteome_SISII_Log2Abu")

head(data_tumour) # inspect the dataframe

# the following code graphs the plot for the tumour dataframe. Graphing all three at once can be done by wrapping all this in a for loop.
long_df <- data_tumour %>% # pivot to a longer dataframe to compare intensities per sample per time point
  pivot_longer(
    cols = !c(Protein.ID, Protein.name, Gene.name, SILAC.label), 
    names_to = "sample", 
    values_to = "log2_intensity"
  ) %>%
  separate(sample, into = c("cell_line", "timepoint", "replicate"), sep = "_", remove = FALSE) # separate sample names into categories

log2f_df <- long_df %>%    # take the intensities from each replicate
  group_by(across(all_of("Protein.ID")), cell_line, timepoint) %>%
  filter(timepoint %in% c("T0", "T6h")) %>%
  pivot_wider(
    id_cols = c(all_of("Protein.ID"), cell_line, replicate),
    names_from = timepoint,
    values_from = log2_intensity) %>% 
  mutate(log2f = T6h - T0) # from log laws, log(a)-log(b) = log(a/b)
  
head(log2f_df) # the table has 4 columns: Protein.ID, cell_line, T0, T6h

log2f_df <- log2f_df %>%  # filter the dataframe to only include proteins with at least 2 data points
  filter(!is.na(log2f)) %>%
  group_by(across(all_of("Protein.ID"))) %>%
  filter(n() >= 2) %>%
  ungroup()

high_fc <- log2f_df %>% # takes the average of all data points for any protein from any cell line
  group_by(across(all_of("Protein.ID"))) %>%
  summarize(mean_fc = mean(log2f, na.rm = TRUE), .groups = "drop") # average of log values is valid, gives log2 of the geometric mean

top_n <- 10  # choose how many top hits you want to show per plot

upregulated <- high_fc %>%   # get top_n proteins with the highest log2f  
  arrange(desc(mean_fc)) %>%
  slice_head(n = top_n) %>%
  pull(!!sym("Protein.ID"))

downregulated <- high_fc %>% # get top_n proteins with the lowest log2f
  arrange(mean_fc) %>%
  slice_head(n = top_n) %>%
  pull(!!sym("Protein.ID"))

up_data <- log2f_df %>% # retrieves data points from the main dataframe
  filter(!!sym("Protein.ID") %in% upregulated) %>%
  mutate(direction = "Upregulated")

down_data <- log2f_df %>%
  filter(!!sym("Protein.ID") %in% downregulated) %>%
  mutate(direction = "Downregulated")

# plotting the data
plot_data <- bind_rows(up_data, down_data)

ggplot(plot_data, aes(x = !!sym("Protein.ID"), y = log2f, fill = direction)) +
  #geom_violin(trim = FALSE, alpha = 0.7) + # tempting as it was, it's too ugly
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.5) +
  facet_wrap(~direction, scales = "free_x") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Log2 Fold Change (T6h/T0)",
    subtitle = "Tumour Cells",
    x = "Protein ID",
    y = "Log2 Fold Change (T6h/T0)",
    fill = "Direction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # rotate names to avoid overcrowding
    strip.text = element_text(size = 14, face = "bold")
  )
