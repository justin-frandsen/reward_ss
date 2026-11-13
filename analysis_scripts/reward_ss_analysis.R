#-------------------------------------------------------------------------
# Script: reward_ss_analysis.R
# Description: Behavioral + eye movement analysis for reward_ss.
#              Based on region_suppression/analysis_scripts/
#              region_suppression_analysis.R. The question here is whether
#              search is faster when the target appears in its high-reward
#              region than in a low-reward region, so the design is a paired
#              comparison (High vs Low) within subject.
#
# Usage:
# - Working directory = analysis_scripts/
# - Run eye_preprocessing.R first if you want the eye movement section.
#
# Design (reward_ss.m + setup/randomizor_reward.m on the main branch):
#   Each subject has 4 target shapes. Every shape is tied to one region of the
#   scene (target_shape_association: 1 = wall, 2 = counter, 3 = floor; the 4th
#   shape repeats one of the three regions).
#   condition 0     = target is in its associated region -> HIGH reward if correct
#   condition 1, 2  = target is in one of the other two regions -> LOW reward
#   Each condition is 1/3 of the trials, so every shape is equally likely to
#   be in each region: only the reward differs, not the location probability.
#   Region of the target on a trial = ceiling(target_position / 2).
#   Run 1 is practice (no reward), runs 2-7 are the main runs.
#-------------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(lmerTest)    # for p-values in lmer
library(emmeans)
library(effectsize)

# ======= SETTINGS =====================
data_root     <- Sys.getenv("RWD_DATA_ROOT", "../data")
bx_folder     <- file.path(data_root, "bx_data")
fixation_file <- file.path(data_root, "eye_data", "fixations_search_period.csv")
ia_file       <- file.path(data_root, "eye_data", "interest_areas.csv")

required_main_runs <- 6   # main runs a subject needs to be included
min_accuracy     <- 0.80   # overall accuracy cutoff (main runs)
min_rt           <- 200    # ms, hard lower bound
sd_cutoff        <- 3      # RTs beyond this many SDs are removed
roi_pad          <- 0      # pixels added around each shape rect when scoring fixations
excluded_subs    <- c()    # e.g. c(120)

# ======= IMPORT BEHAVIORAL DATA =====================
bx_files <- list.files(bx_folder, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
if (length(bx_files) == 0) stop("No CSV files found in folder: ", bx_folder)

# everything is read as text first so empty columns do not get guessed as different types
raw_imported_bx_files <- map_dfr(bx_files, ~ read_csv(.x, col_types = cols(.default = col_character()))) %>%
  type_convert(col_types = cols(phase = col_character(), scene_file = col_character(),
                                response_key = col_character(), timestamp = col_character()),
               na = c("", "NA", "NaN"))

region_names <- c("Wall", "Counter", "Floor")

all_imported_bx_files <- raw_imported_bx_files %>%
  filter(run_num != 1) %>%   # run 1 is practice
  mutate(RewardRegion = case_when(condition == 0 ~ "High",
                                  condition %in% c(1, 2) ~ "Low"),
         RewardRegion = factor(RewardRegion, levels = c("High", "Low")),
         target_region = factor(region_names[ceiling(target_position / 2)], levels = region_names),
         correct = as.integer(accuracy == 1))   # accuracy: 1 correct, 0 incorrect, -1 no response

# sanity check: condition 0 must be the trials where the target is in its associated region
stopifnot(with(all_imported_bx_files,
               all((condition == 0) == (ceiling(target_position / 2) == target_shape_association))))

# ======= RUN SUMMARY & ACCURACY =====================
unique_run_summary <- raw_imported_bx_files %>%
  group_by(sub_num) %>%
  summarise(unique_runs = n_distinct(run_num[run_num != 1]),   # main runs only
            overall_accuracy = mean(accuracy[run_num != 1] == 1, na.rm = TRUE),
            .groups = "drop")
print(unique_run_summary, n = Inf)

# ======= BEHAVIORAL DATA CLEANUP & FILTERING =======
all_bx_files <- all_imported_bx_files %>%
  left_join(unique_run_summary, by = "sub_num") %>%
  # --- Basic inclusion filters ---
  filter(unique_runs >= required_main_runs,
         overall_accuracy > min_accuracy,
         !sub_num %in% excluded_subs) %>%
  # --- Remove RT outliers within participant and condition ---
  group_by(sub_num, RewardRegion) %>%
  mutate(rt = ifelse(correct == 0, NA, rt),
         rt = ifelse(rt <= min_rt, NA, rt),   # Hard lower bound
         mean_rt = mean(rt, na.rm = TRUE),
         sd_rt   = sd(rt, na.rm = TRUE),
         rt = ifelse(rt > mean_rt + sd_cutoff * sd_rt | rt < mean_rt - sd_cutoff * sd_rt, NA, rt)) %>%
  ungroup() %>%
  select(-mean_rt, -sd_rt) %>%
  mutate(log_rt = log(rt),
         run_c   = run_num - mean(unique(run_num)),   # centered run, for the learning model
         sub_num = as.factor(sub_num),
         run_num = as.factor(run_num))

removed_trials <- all_bx_files %>%
  group_by(sub_num) %>%
  summarise(total_trials = n(),
            kept_trials = sum(!is.na(rt)),
            removed_trials = sum(is.na(rt)))
print(removed_trials, n = Inf)

# ======= SUMMARIES (one row per subject x condition) =======
bx_summary <- all_bx_files %>%
  group_by(sub_num, RewardRegion) %>%
  summarise(meanRT = mean(rt, na.rm = TRUE),
            medianRT = median(rt, na.rm = TRUE),
            meanACC = mean(correct, na.rm = TRUE),
            n_trials = n(),
            .groups = "drop")

print(table(bx_summary$RewardRegion)) # checks to make sure that we have balanced data

# ---- helper: High vs Low paired t-test + effect size ----
# (with two levels the repeated measures ANOVA is the same test: F = t^2)
run_paired_test <- function(data, dv) {
  wide <- data %>% filter(!is.na(.data[[dv]])) %>%
    select(sub_num, RewardRegion, all_of(dv)) %>%
    pivot_wider(names_from = RewardRegion, values_from = all_of(dv)) %>%
    filter(!is.na(High), !is.na(Low))   # keep only subjects who have both levels

  cat("\n==================", dv, ": High vs Low ==================\n")
  cat(sprintf("N = %d | High M = %.3f (SD %.3f) | Low M = %.3f (SD %.3f) | High - Low = %.3f | %d of %d subjects High < Low\n",
              nrow(wide), mean(wide$High), sd(wide$High), mean(wide$Low), sd(wide$Low),
              mean(wide$High - wide$Low), sum(wide$High < wide$Low), nrow(wide)))
  test <- t.test(wide$High, wide$Low, paired = TRUE)
  print(test)
  print(effectsize::cohens_d(wide$High, wide$Low, paired = TRUE))   # d_z
  invisible(test)
}

# ======= RT & ACCURACY =======
t_RT       <- run_paired_test(bx_summary, "meanRT")
t_medianRT <- run_paired_test(bx_summary, "medianRT")
t_ACC      <- run_paired_test(bx_summary, "meanACC")

# ======= SINGLE-TRIAL LOG RT =======
# target_region is a covariate: each subject's 4th shape repeats one region, so
# High trials are not perfectly balanced over wall/counter/floor.
lmer_data <- all_bx_files %>% filter(!is.na(log_rt))
lmer_RT <- lmer(log_rt ~ RewardRegion + target_region + (1 + RewardRegion | sub_num) + (1 | scene_idx),
                data = lmer_data)
print(summary(lmer_RT))
print(emmeans(lmer_RT, pairwise ~ RewardRegion, type = "response"))

# ======= DOES THE REWARD EFFECT GROW WITH LEARNING? =======
run_summary <- all_bx_files %>%
  group_by(sub_num, run_num, RewardRegion) %>%
  summarise(meanRT = mean(rt, na.rm = TRUE), .groups = "drop")

cat("\n--- Mean RT by run (average of subject means) ---\n")
run_summary %>%
  group_by(run_num, RewardRegion) %>%
  summarise(meanRT = mean(meanRT), .groups = "drop") %>%
  pivot_wider(names_from = RewardRegion, values_from = meanRT) %>%
  mutate(High_minus_Low = High - Low) %>%
  print()

lmer_RT_run <- lmer(log_rt ~ RewardRegion * run_c + target_region + (1 + RewardRegion | sub_num) + (1 | scene_idx),
                    data = lmer_data)
print(summary(lmer_RT_run))

#======================= EYETRACKING ANALYSIS =========================
if (!file.exists(fixation_file) || !file.exists(ia_file)) {
  message("[INFO] No fixation file found (", fixation_file, "). Run eye_preprocessing.R to get the eye movement analysis.")
} else {
  fixations <- read_csv(fixation_file, show_col_types = FALSE) %>%
    mutate(sub_num = as.factor(sub_num), run_num = as.factor(run_num))

  # ---- distractor rects ([left top right bottom]): the bx files only hold the
  # target rect, so the three distractor rects come from the eyetracker messages ----
  distractor_rects <- read_csv(ia_file, show_col_types = FALSE) %>%
    filter(str_detect(ia_label, "^NonCritDistBox")) %>%
    mutate(sub_num = as.factor(sub_num), run_num = as.factor(run_num),
           d = str_extract(ia_label, "\\d+$")) %>%
    select(sub_num, run_num, trial_num, d, left, top, right, bottom) %>%
    pivot_wider(names_from = d, values_from = c(left, top, right, bottom), names_glue = "dist{d}_{.value}")

  in_rect <- function(x, y, left, top, right, bottom) {
    !is.na(left) & x >= left - roi_pad & x <= right + roi_pad & y >= top - roi_pad & y <= bottom + roi_pad
  }

  fixations_with_ROIs <- fixations %>%
    inner_join(all_bx_files %>%
                 select(sub_num, run_num, trial_num, RewardRegion, condition, rt,
                        matches("_rect")),
               by = c("sub_num", "run_num", "trial_num")) %>%
    left_join(distractor_rects, by = c("sub_num", "run_num", "trial_num")) %>%
    filter(!is.na(rt)) %>%   # same trials as the RT analysis (correct, not an outlier)
    mutate(in_TargetBox = in_rect(fix_x, fix_y, target_rect_1, target_rect_2, target_rect_3, target_rect_4),
           in_DistBox =
             in_rect(fix_x, fix_y, dist1_left, dist1_top, dist1_right, dist1_bottom) |
             in_rect(fix_x, fix_y, dist2_left, dist2_top, dist2_right, dist2_bottom) |
             in_rect(fix_x, fix_y, dist3_left, dist3_top, dist3_right, dist3_bottom),
           current_roi = case_when(in_TargetBox ~ "TargetBox",
                                   in_DistBox ~ "DistBox",
                                   TRUE ~ NA_character_))

  # check that the eye files line up with the behavioral files:
  # the search period in the EDF should be as long as the RT
  timing_check <- fixations_with_ROIs %>% distinct(sub_num, run_num, trial_num, rt, search_duration)
  cat(sprintf("\n[CHECK] search period vs RT: median |difference| = %.1f ms, max = %.1f ms, %d trials\n",
              median(abs(timing_check$search_duration - timing_check$rt)),
              max(abs(timing_check$search_duration - timing_check$rt)), nrow(timing_check)))

  # ---- Summarize fixations per trial ----
  trial_fixation_summary <- fixations_with_ROIs %>%
    group_by(sub_num, run_num, trial_num, RewardRegion, rt) %>%
    arrange(fix_index, .by_group = TRUE) %>%
    summarise(
      total_fixations        = n(),
      fixations_on_Target    = sum(in_TargetBox),
      fixations_on_Dist      = sum(in_DistBox),
      # ---- First shape fixated (in time order) ----
      first_fixation_ROI     = current_roi[which(!is.na(current_roi))[1]],
      # ---- Fixation number when target was first fixated ----
      Target_fixation_number = fix_index[which(in_TargetBox)[1]],
      # ---- Time from search onset to the start of the first target fixation (ms) ----
      time_to_target         = pmax(fix_start[which(in_TargetBox)[1]], 0),
      .groups = "drop") %>%
    mutate(target_first_fix  = as.integer(!is.na(first_fixation_ROI) & first_fixation_ROI == "TargetBox"),
           # time from first landing on the target to the response
           verification_time = rt - time_to_target)

  cat(sprintf("[CHECK] target was fixated on %.1f%% of analyzed trials\n",
              100 * mean(!is.na(trial_fixation_summary$Target_fixation_number))))

  fixation_summary <- trial_fixation_summary %>%
    group_by(sub_num, RewardRegion) %>%
    summarise(mean_target_first_fix       = mean(target_first_fix, na.rm = TRUE),
              mean_Target_fixation_number = mean(Target_fixation_number, na.rm = TRUE),
              mean_time_to_target         = mean(time_to_target, na.rm = TRUE),
              mean_verification_time      = mean(verification_time, na.rm = TRUE),
              mean_total_fixations        = mean(total_fixations, na.rm = TRUE),
              n_trials                    = n(),
              .groups = "drop")

  # ---- Target measures: High vs Low ----
  t_First_fix    <- run_paired_test(fixation_summary, "mean_target_first_fix")
  t_ordinal_fix  <- run_paired_test(fixation_summary, "mean_Target_fixation_number")
  t_time_target  <- run_paired_test(fixation_summary, "mean_time_to_target")
  t_verification <- run_paired_test(fixation_summary, "mean_verification_time")
  t_total_fix    <- run_paired_test(fixation_summary, "mean_total_fixations")
}
