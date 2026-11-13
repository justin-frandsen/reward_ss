#-------------------------------------------------------------------------
# Script: graphs.R
# Description: Graphs for reward_ss (same style as region_suppression
#              graphs.R): violin per condition with the mean and
#              within-subject 95% CI on top. Saved as .png in ../figures.
#-------------------------------------------------------------------------
source("reward_ss_analysis.R")

library(gghalves)

palatte_num <- 7

nicelimits <- function(x) {
  range(scales::extended_breaks(only.loose = TRUE)(x))
}

# data: one row per subject x RewardRegion; dv: column name (string)
condition_violin <- function(data, dv, ylab, percent = FALSE) {
  data <- data %>% filter(!is.na(.data[[dv]])) %>% droplevels()

  # Compute within-subjects summary stats
  summary_within <- Rmisc::summarySEwithin(data = as.data.frame(data),
                                           measurevar = dv,
                                           withinvars = "RewardRegion",
                                           idvar = "sub_num",
                                           na.rm = TRUE)

  p <- ggplot(data, aes(x = RewardRegion, y = .data[[dv]])) +
    geom_half_violin(aes(fill = RewardRegion), side = "r", alpha = 0.6,
                     position = position_nudge(x = 0.07)) +
    geom_line(aes(group = sub_num), color = "grey70", alpha = 0.5) +
    geom_point(data = summary_within, aes(color = RewardRegion), size = 3) +
    geom_errorbar(data = summary_within,
                  aes(color = RewardRegion,
                      ymin = .data[[dv]] - ci,
                      ymax = .data[[dv]] + ci),
                  width = 0.1) +
    labs(x = "Reward of the target's region", y = ylab) +
    scale_fill_brewer(type = "qual", palette = palatte_num) +
    scale_color_brewer(type = "qual", palette = palatte_num) +
    theme_classic() +
    theme(legend.position = "none")

  if (percent) {
    p + scale_y_continuous(labels = scales::percent, limits = nicelimits)
  } else {
    p + scale_y_continuous(limits = nicelimits)
  }
}

figure_folder <- "../figures"
dir.create(figure_folder, showWarnings = FALSE)
save_fig <- function(plot, name) {
  ggsave(file.path(figure_folder, name), plot, width = 4, height = 4, dpi = 300)
  plot
}

response_time_violin <- condition_violin(bx_summary, "meanRT", "Mean Response Time (ms)")
save_fig(response_time_violin, "rt_by_reward_region.png")

accuracy_violin <- condition_violin(bx_summary, "meanACC", "Accuracy", percent = TRUE)
save_fig(accuracy_violin, "accuracy_by_reward_region.png")

if (exists("fixation_summary")) {
  target_first_fix_violin <- condition_violin(fixation_summary, "mean_target_first_fix",
                                              "Percentage of first fixations on target", percent = TRUE)
  print(save_fig(target_first_fix_violin, "target_first_fix_by_reward_region.png"))

  time_to_target_violin <- condition_violin(fixation_summary, "mean_time_to_target",
                                            "Time to first target fixation (ms)")
  print(save_fig(time_to_target_violin, "time_to_target_by_reward_region.png"))

  ordinal_fix_violin <- condition_violin(fixation_summary, "mean_Target_fixation_number",
                                         "Fixation number of first target fixation")
  print(save_fig(ordinal_fix_violin, "target_fixation_number_by_reward_region.png"))
}
