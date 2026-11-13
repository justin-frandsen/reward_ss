Subj100Run02 <- read.csv("~/MATLAB/repos/reward_ss/data/bx_data/Subj100Run02.csv", header=TRUE)

library(dplyr)

Subj100Run02 <- Subj100Run02 %>% 
  mutate(position_type = ifelse(target_position == 1 | target_position == 2, "wall", 
                                ifelse(target_position == 3 | target_position == 4, "counter", 
                                       ifelse(target_position == 5 | target_position == 6, "floor", NA))))

Subj100Run02 %>%
  count(reward_amount)

Subj100Run02 %>%
  group_by(target_shape_idx, position_type) %>% 
  count(reward_amount) %>%
  mutate(prop = n / sum(n))
