library("tidyverse")

the_data <- read_delim(file="summary_sealed_bid.csv",delim=";")

the_data_l <- pivot_longer(data=the_data
        ,cols=c("mean_belligerence","mean_bravery")
        ,names_to="trait"
        ,values_to="trait_values")

ggplot(data=the_data_l
        ,mapping=aes(x=d, y=trait_values)) +
    geom_point(mapping=aes(colour=trait)) +
    facet_grid(~h)

ggsave(filename="plot_bel_brav.pdf")
