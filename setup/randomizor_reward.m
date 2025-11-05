%-------------------------------------------------------------------------
% Script: randomizor_drinking.m
% Author: Justin Frandsen
% Date: 2025/08/05 yyyy/mm/dd
% Description: Prerandomizor for the drinking_ss experiment.
%
% For each subject:
% - Assigns scenes to runs (balanced across 6 runs)
% - Randomizes scene orders, distractors, targets, and conditions
% - Generates trial-specific parameters for 7 experimental runs
% - Outputs a struct per subject with trial information
%-------------------------------------------------------------------------

%% CONFIGURATION
SAVE_OUTPUT = true;  % Set to true to save output .mat file

% Column constants for scene matrix
SCENE_ID   = 1;
REP        = 2;
RUN        = 3;
TARGET     = 4;
CONDITION  = 5;

% Experiment parameters
total_subs              = 1000;   % Total subjects to generate
total_runs              = 7;     % Number of experimental runs per subject
number_trials_per_run   = 72;    % Trials per run
total_trials            = (total_runs-1) * number_trials_per_run; % Total trials
total_scenes            = 108;   % Number of unique scenes in each type
total_reps_per_scene    = 4;     % Number of times each scene is shown

% Condition configurations
distractor_inds         = [1 2 3];                     % Indices of distractor shapes
condition_inds          = [0 1 2 0 1 2];               % conditions

% Load and filter shapes from stimulus folder
all_shapes = dir('../stimuli/shapes/transparent_black/*');
all_shapes = all_shapes(~ismember({all_shapes.name}, {'.','..','.DS_Store'}));
shape_inds = 1:length(all_shapes);  % Shape indices from filtered directory

% Safety checks
assert(mod(total_scenes * total_reps_per_scene*2, total_runs-1) == 0, ...
    'Scene repetitions must divide evenly into runs.');

% Generate all permutations of [0 0 1 1] (used for target directions)
t_directions = unique(perms([0 0 1 1]), 'rows');

%% Initialize output
randomizor_matrix = struct();
rng('shuffle');  % Seed RNG for per-subject randomness

fprintf('[INFO] Beginning randomization for %d subjects...\n', total_subs);

%% Subject loop
for sub_num = 1:total_subs
    sub_struct_name = sprintf('subj%d', sub_num);  % e.g., 'subj1'
    subject_struct = struct();

    % Initialize scene matrix for this subject
    % Columns: [scene_id, rep, run, distractor, target, condition]
    scene_randomizor = zeros(total_scenes * total_reps_per_scene, 5);

    % Populate scene ID and repetition
    row = 1;
    for scene = 1:total_scenes
        for rep = 1:total_reps_per_scene
            scene_randomizor(row, SCENE_ID:REP) = [scene, rep];
            row = row + 1;
        end
    end

    % Assign randomized run numbers (1–6) in balanced blocks
    scene_randomizor = assign_balanced_runs(scene_randomizor, total_runs, RUN);

    %% Define stimuli for this subject
    all_targets = randsample(shape_inds, 4);  % 3 targets
    all_distractors = setdiff(shape_inds, all_targets);         % Remaining distractors

    target_associations = [1 2 3 randsample(3, 1)];

    %% Create distractor triplets for runs 5–6
    triplets = nchoosek(all_distractors, 3);
    all_triplets = [];

    % Generate all permutations for each triplet
    for i = 1:size(triplets, 1)
        all_triplets = [all_triplets; perms(triplets(i, :))];
    end

    trial_distractor_inds = shuffle_matrix(all_triplets, [1 2 3], [3 3 3]);

    % Assign run numbers to triplets (72 per run × 6 runs: 2 through 7)
    row = 1;
    for run = 2:7
        trial_distractor_inds(row:row+71, 4) = run;
        row = row + 72;
    end

    %% Trial run loop
    current_target = 1;  % Start with target index 1


    for run = 1:total_runs
        run_struct = struct();
        run_name = sprintf('run%d', run); % this is because the first run is practice

        % Extract this run’s rows
        if run == 1
            %practice run
            this_run_scene = scene_randomizor(scene_randomizor(:, RUN) == 1, :);
            % Practice run: Randomize scenes, no distractors/conditions
            target_choice =[1 2 3 4 4 3 2 1];
            practice_run_matrix = zeros(8, 5);
            counter = 1;
            for i = 1:8
                practice_run_matrix(i, SCENE_ID) = counter;
                practice_run_matrix(i, REP) = 0; %this doesn't matter for practice
                practice_run_matrix(i, TARGET) = target_choice(i);
                practice_run_matrix(i, CONDITION) = 0; %always valid for practice
                practice_run_matrix(i, RUN) = 1; %practice run
                counter = counter + 1;
                if counter > 4
                    counter = 1;
                end
            end

            practice_run_matrix = shuffle_matrix(practice_run_matrix, ...
                [SCENE_ID TARGET], [1 1], 10000);

            numRows = size(trial_distractor_inds, 1);
            randIdx = randperm(numRows, 8);  % random permutation of row indices
            run_distractors = trial_distractor_inds(randIdx, :);

            all_possible_locations = zeros(8, 4);
            loc1 = [1 2];
            loc2 = [3 4];
            loc3 = [5 6];
            for trial = 1:8
                rand1 = loc1(randperm(2));
                rand2 = loc2(randperm(2));
                rand3 = loc3(randperm(2));
                
                possible_locations = [rand1(1), rand2(1), rand3(1), 0];  % first 3 are random from each group
                
                remaining_distractors = [...
                                        setdiff(loc1, rand1(1)), ...
                                        setdiff(loc2, rand2(1)), ...
                                        setdiff(loc3, rand3(1))...
                                        ];
                
                possible_locations(4) = remaining_distractors(randi(length(remaining_distractors)));  % random from remaining
                all_possible_locations(trial, :) = possible_locations;
            end
        elseif run > 1
            % Experimental runs 1–6
            this_run_scene = scene_randomizor(scene_randomizor(:, RUN) == run, :);

            this_run_scene = this_run_scene(randperm(size(this_run_scene, 1)), :);

            % Assign targets and conditions in blocks of 8
            con_len = length(condition_inds);
            for i = 1:number_trials_per_run/con_len
                idx = (i-1)*con_len + 1;
                this_run_scene(idx:idx+(con_len-1), TARGET) = current_target;
                this_run_scene(idx:idx+(con_len-1), CONDITION) = condition_inds(randperm(con_len));
                current_target = mod(current_target, 4) + 1;
            end

            % Final shuffle
            this_run_scene = shuffle_matrix(this_run_scene, [SCENE_ID TARGET], [2 3], 10000);

            % Use second-half distractors
            run_distractors = trial_distractor_inds(trial_distractor_inds(:,4) == run, :);

            all_possible_locations = zeros(72, 4);
            loc1 = [1 2];
            loc2 = [3 4];
            loc3 = [5 6];
            for trial = 1:number_trials_per_run
                rand1 = loc1(randperm(2));
                rand2 = loc2(randperm(2));
                rand3 = loc3(randperm(2));
            
                possible_locations = [rand1(1), rand2(1), rand3(1), 0];  % first 3 are random from each group
                
                remaining_distractors = [...
                                        setdiff(loc1, rand1(1)), ...
                                        setdiff(loc2, rand2(1)), ...
                                        setdiff(loc3, rand3(1))...
                                        ];
                
                possible_locations(4) = remaining_distractors(randi(length(remaining_distractors)));  % random from remaining
                all_possible_locations(trial, :) = possible_locations;
            end
        end

        % Generate target direction matrix
        rep_count = length(this_run_scene) / size(t_directions, 1);
        if run == 1
            rep_count = 2;  % Practice run uses each direction once
        end
        full_directions = repmat(t_directions, rep_count, 1);
        full_directions = full_directions(randperm(size(full_directions, 1)), :);

        % Store all run-specific info
        run_struct.targets                = all_targets;
        run_struct.distractors            = all_distractors;
        if run == 1
            run_struct.scene_randomizor = practice_run_matrix;
        else
            run_struct.scene_randomizor = this_run_scene;
        end
        run_struct.t_directions                      = full_directions;
        run_struct.target_associations               = target_associations;
        run_struct.this_run_distractors              = run_distractors;
        run_struct.all_possible_locations            = all_possible_locations;

        subject_struct.(run_name) = run_struct;
    end

    % Save full subject struct
    randomizor_matrix.(sub_struct_name) = subject_struct;

    % Print progress every 50 subjects
    if mod(sub_num, 50) == 0
        fprintf('[INFO] Completed subject %d/%d\n', sub_num, total_subs);
    end
end

fprintf('[INFO] Randomization complete.\n');

% Save output if configured
if SAVE_OUTPUT
    save('../trial_structure_files/randomizor.mat', 'randomizor_matrix');
    fprintf('[INFO] Output saved to ../trial_structure_files/randomizor.mat\n');
end

%% Local function for assigning run numbers in a balanced way
function scene_randomizor = assign_balanced_runs(scene_randomizor, total_runs, RUN, practice_run)
    if nargin < 4
        practice_run = true;
    end

    if practice_run
        % Exclude run 1 (practice) from balancing
        possible_runs = 2:total_runs;
        number_main_runs = total_runs - 1;
    else
        possible_runs = 1:total_runs;
        number_main_runs = total_runs;
    end

    num_rows = size(scene_randomizor, 1);

    % Check that the number of rows is divisible by the number of runs
    assert(mod(num_rows, number_main_runs) == 0, ...
        'Total number of rows (%d) must be divisible by total_runs (%d).', ...
        num_rows, number_main_runs);

    % Assign a random permutation of run numbers to each group of total_runs rows
    for i = 1:(num_rows / number_main_runs)
        idx = (i - 1) * number_main_runs + 1;
        scene_randomizor(idx : idx + number_main_runs - 1, RUN) = possible_runs(randperm(number_main_runs));
    end
end
