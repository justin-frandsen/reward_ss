#-------------------------------------------------------------------------
# Script: eye_preprocessing.R  (reward_ss)
# Description: Replaces the Data Viewer fixation report. Converts the EDF
#              files to ASC (with SR Research's edf2asc utility), then pulls
#              every fixation that falls in the SEARCH_PERIOD of each trial
#              and saves one long .csv (one row per fixation).
#
# Usage:
# - Run with the working directory set to analysis_scripts/ (like E1.R/E2.R).
# - edf2asc comes with the free EyeLink Developers Kit (the same kit the
#   Eyelink toolbox in Psychtoolbox needs, so the tracking PCs have it).
#   EDF is a closed binary format, so this one step cannot be done in pure R.
#   If edf2asc is not installed on this computer, convert on a computer that
#   has it and copy the .asc files into data/asc_data/.
# - EDF files that are missing the start of the file are repaired on the fly
#   (see STEP 0); the files in data/edf_data are never modified.
#
# NOTE: the experiment script only sends END_TIME SEARCH_PERIOD when a
#   response is made, so timed-out trials (no response) have no search period
#   and are dropped here. They are excluded from the RT analysis anyway.
#
# Output 1: data/eye_data/fixations_search_period.csv
#   sub_num, run_num, trial_num : match the behavioral files
#   fix_index            : order of the fixation within the search period
#                          (1 = fixation in progress at search onset, same as
#                          CURRENT_FIX_INDEX in a Data Viewer interest period)
#   fix_x, fix_y         : average gaze position (pixels)
#   fix_start, fix_end   : ms relative to search display onset
#   fix_dur              : full fixation duration (ms)
#   fix_dur_in_period    : duration clipped to the search period (ms)
#   started_before_onset : TRUE if the fixation began before the search display
#
# Output 2: data/eye_data/interest_areas.csv (one row per shape per trial)
#   The bx files only store the target rect, so the distractor rects are
#   taken from the "!V IAREA RECTANGLE" messages the experiment sends.
#   sub_num, run_num, trial_num, ia_id, left, top, right, bottom, ia_label
#   (ia_label: TargetBox, NonCritDistBox1..3)
#-------------------------------------------------------------------------
library(tidyverse)

data_root  <- Sys.getenv("RWD_DATA_ROOT", "../data")
edf_folder <- file.path(data_root, "edf_data")
asc_folder <- file.path(data_root, "asc_data")
eye_folder <- file.path(data_root, "eye_data")
out_file   <- file.path(eye_folder, "fixations_search_period.csv")
ia_file    <- file.path(eye_folder, "interest_areas.csv")

period_start_msg <- "START_TIME SEARCH_PERIOD"
period_end_msg   <- "END_TIME SEARCH_PERIOD"

dir.create(asc_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(eye_folder, showWarnings = FALSE, recursive = TRUE)

# ======= STEP 0: repair EDF files that lost their first chunk =====================
# About half of the EDF files in this experiment arrived without the start of
# the file (the text header plus the first few records), and edf2asc refuses
# them ("End of file Exception"). Everything after the missing chunk is intact,
# so a readable copy is rebuilt as:
#   text header of any intact EDF  +  this file from its first clean record boundary
# - If the file starts before the first recording block (in the calibration
#   messages) nothing else is needed and every trial is recovered.
# - If it starts inside trial 1's recording, the recording-start block of the
#   file's next trial is put in front so the samples can be decoded. Trial 1
#   has then lost its TRIALID message and is dropped in step 2.
edf_has_header <- function(edf) {
  identical(readBin(edf, "raw", 11), charToRaw("SR_RESEARCH"))
}

# position of the first occurrence of a byte pattern in a raw vector (NA if none)
find_bytes <- function(bytes, pattern, from = 1) {
  hits <- which(bytes == pattern[1])
  hits <- hits[hits >= from & hits <= length(bytes) - length(pattern) + 1]
  for (k in seq_along(pattern)[-1]) hits <- hits[bytes[hits + k - 1] == pattern[k]]
  if (length(hits) == 0) NA_integer_ else hits[1]
}

repair_headerless_edf <- function(donor_edf, bad_edf, out_edf) {
  msg_tag    <- as.raw(c(0x18, 0x00, 0x21, 0x00))   # message event
  sample_tag <- as.raw(c(0x6f, 0xc1, 0x00))         # full (non delta-coded) sample
  start_tag  <- as.raw(c(0x0f, 0x00, 0x21, 0x00))   # recording-start block

  donor  <- readBin(donor_edf, "raw", 4096)
  header <- donor[1:(find_bytes(donor, charToRaw("ENDP:\n")) + 5)]
  bad    <- readBin(bad_edf, "raw", file.size(bad_edf))

  i_msg    <- find_bytes(bad, msg_tag)
  i_sample <- find_bytes(bad, sample_tag)
  i_start  <- find_bytes(bad, start_tag)
  if (is.na(i_msg) || is.na(i_sample) || is.na(i_start)) return(FALSE)
  body <- bad[min(i_msg, i_sample):length(bad)]   # first clean record boundary
  if (i_sample < i_start) {
    # the file begins inside trial 1's recording (its start block is gone):
    # borrow the next trial's start block (it ends where the next message begins)
    e <- find_bytes(bad, msg_tag, from = i_start)
    if (is.na(e)) return(FALSE)
    body <- c(bad[i_start:(e - 1)], body)
  }
  writeBin(c(header, body), out_edf)
  TRUE
}

# ======= STEP 1: EDF -> ASC =====================
convert_edf_to_asc <- function(edf_folder, asc_folder) {
  # recursive: in case a transfer saved the file inside a folder of the same name
  edf_files <- list.files(edf_folder, pattern = "\\.edf$", full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  asc_names <- paste0(tools::file_path_sans_ext(basename(edf_files)), ".asc")
  todo      <- edf_files[!file.exists(file.path(asc_folder, asc_names))]
  if (length(todo) == 0) {
    message("[INFO] No new EDF files to convert.")
    return(invisible())
  }

  edf2asc <- Sys.which("edf2asc")
  if (edf2asc == "") {
    warning(length(todo), " EDF file(s) have no .asc yet and edf2asc was not found on this computer. ",
            "Install the EyeLink Developers Kit, or convert elsewhere and copy the .asc files into ",
            asc_folder, call. = FALSE)
    return(invisible())
  }

  # donor for the text header of headerless files (see repair_headerless_edf)
  has_header <- map_lgl(edf_files, edf_has_header)
  donor      <- edf_files[has_header][1]

  for (edf in todo) {
    src <- edf
    if (!edf_has_header(edf)) {
      # headerless transfer: convert a repaired copy, the original EDF is never changed
      src <- file.path(asc_folder, basename(edf))
      ok  <- !is.na(donor) && tryCatch(repair_headerless_edf(donor, edf, src), error = function(e) FALSE)
      if (!ok) {
        warning("Could not repair headerless EDF ", basename(edf), call. = FALSE)
        next
      }
    }
    # -e = events only (fixations, saccades, messages; no samples) -> small files
    # -y = overwrite an existing .asc
    system2(edf2asc, c("-e", "-y", shQuote(src)), stdout = FALSE, stderr = FALSE)
    made <- paste0(tools::file_path_sans_ext(src), ".asc")  # edf2asc writes next to the EDF
    if (file.exists(made)) {
      if (src != edf) file.remove(src) else file.rename(made, file.path(asc_folder, basename(made)))
      message("[INFO] Converted ", basename(edf), if (src != edf) " (header repaired)" else "")
    } else {
      if (src != edf) file.remove(src)
      warning("edf2asc did not produce an .asc for ", basename(edf),
              " (file may be a partial transfer).", call. = FALSE)
    }
  }
}

# ======= STEP 2: ASC -> fixations in the search period =====================
parse_asc_file <- function(asc_file) {
  ids <- str_match(basename(asc_file), "^RS(\\d+)R(\\d+)")
  if (is.na(ids[1, 1])) {
    warning("Skipping ", basename(asc_file), ": name does not look like RS<sub>R<run>.asc", call. = FALSE)
    return(NULL)
  }
  lines <- readLines(asc_file, warn = FALSE)

  # ---- messages: "MSG <time> <text>" ----
  msg <- str_match(lines[startsWith(lines, "MSG")], "^MSG\\s+(\\d+)\\s+(.*)$")
  msg <- tibble(time = as.numeric(msg[, 2]), text = str_trim(msg[, 3])) %>% filter(!is.na(time))

  trial_starts <- msg %>%
    filter(str_detect(text, "^TRIALID\\s+\\d+")) %>%
    transmute(trial_num = as.integer(str_match(text, "^TRIALID\\s+(\\d+)")[, 2]), trial_time = time)
  if (nrow(trial_starts) == 0) {
    warning("No TRIALID messages in ", basename(asc_file), call. = FALSE)
    return(NULL)
  }

  # first period-start / period-end message after each TRIALID (and before the next one)
  next_trial_time <- c(trial_starts$trial_time[-1], Inf)
  first_msg_between <- function(pattern, from, to) {
    t <- msg$time[msg$text == pattern & msg$time >= from & msg$time < to]
    if (length(t) == 0) NA_real_ else t[1]
  }
  trials <- trial_starts %>%
    mutate(search_onset = map2_dbl(trial_time, next_trial_time, ~ first_msg_between(period_start_msg, .x, .y)),
           search_end   = map2_dbl(trial_time, next_trial_time, ~ first_msg_between(period_end_msg, .x, .y))) %>%
    filter(!is.na(search_onset), !is.na(search_end))   # drops aborted trials

  # ---- fixations: "EFIX <eye> <start> <end> <dur> <x> <y> <pupil>" ----
  fix <- str_split_fixed(str_squish(lines[startsWith(lines, "EFIX")]), " ", 8)
  fix <- tibble(eye = fix[, 2],
                st  = suppressWarnings(as.numeric(fix[, 3])),
                et  = suppressWarnings(as.numeric(fix[, 4])),
                dur = suppressWarnings(as.numeric(fix[, 5])),
                x   = suppressWarnings(as.numeric(fix[, 6])),
                y   = suppressWarnings(as.numeric(fix[, 7]))) %>%
    filter(!is.na(st), !is.na(x), !is.na(y))

  # binocular recording: keep one eye (right, same as the experiment script)
  if (n_distinct(fix$eye) > 1) fix <- fix %>% filter(eye == "R")
  if (nrow(fix) == 0 || nrow(trials) == 0) return(NULL)

  # ---- interest areas: "!V IAREA RECTANGLE <id> <left> <top> <right> <bottom> <label>" ----
  ia <- str_match(msg$text, "^!V IAREA RECTANGLE\\s+(\\d+)\\s+(-?\\d+)\\s+(-?\\d+)\\s+(-?\\d+)\\s+(-?\\d+)\\s+(\\S+)")
  ia <- tibble(time = msg$time, ia_id = as.integer(ia[, 2]),
               left = as.numeric(ia[, 3]), top = as.numeric(ia[, 4]),
               right = as.numeric(ia[, 5]), bottom = as.numeric(ia[, 6]), ia_label = ia[, 7]) %>%
    filter(!is.na(ia_id)) %>%
    mutate(trial_idx = findInterval(time, trial_starts$trial_time)) %>%
    filter(trial_idx > 0) %>%   # sent before the first TRIALID in the file (repaired files)
    mutate(trial_num = trial_starts$trial_num[trial_idx],
           sub_num = as.integer(ids[1, 2]), run_num = as.integer(ids[1, 3])) %>%
    select(sub_num, run_num, trial_num, ia_id, left, top, right, bottom, ia_label)

  # ---- keep fixations that overlap each trial's search period ----
  fixations <- trials %>%
    mutate(fixations = map2(search_onset, search_end, ~ fix %>% filter(et > .x, st < .y))) %>%
    unnest(fixations) %>%
    group_by(trial_num) %>%
    arrange(st, .by_group = TRUE) %>%
    mutate(fix_index = row_number()) %>%
    ungroup() %>%
    transmute(sub_num   = as.integer(ids[1, 2]),
              run_num   = as.integer(ids[1, 3]),
              trial_num,
              fix_index,
              eye,
              fix_x     = x,
              fix_y     = y,
              fix_start = st - search_onset,
              fix_end   = et - search_onset,
              fix_dur   = dur,
              fix_dur_in_period    = pmin(et, search_end) - pmax(st, search_onset),
              started_before_onset = st < search_onset,
              search_duration      = search_end - search_onset)

  list(fixations = fixations, interest_areas = ia)
}

# ======= RUN =====================
convert_edf_to_asc(edf_folder, asc_folder)

asc_files <- list.files(asc_folder, pattern = "\\.asc$", full.names = TRUE, ignore.case = TRUE)
if (length(asc_files) == 0) stop("No .asc files found in ", asc_folder)

parsed         <- map(asc_files, parse_asc_file) %>% compact()
fixations      <- map_dfr(parsed, "fixations")
interest_areas <- map_dfr(parsed, "interest_areas")

write_csv(fixations, out_file)
write_csv(interest_areas, ia_file)
message("[INFO] Saved ", nrow(fixations), " fixations from ", length(asc_files), " file(s): ", out_file)
message("[INFO] Saved ", nrow(interest_areas), " interest areas: ", ia_file)

# quick check: every subject/run should have its full set of trials
fixations %>%
  group_by(sub_num, run_num) %>%
  summarise(trials_with_fixations = n_distinct(trial_num), fixations = n(), .groups = "drop") %>%
  print(n = Inf)
