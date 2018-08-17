% Prediction Error Study parameters and folder structure

subjects = [2, 6:8];

%%% Filename and folder structure informations. folders will be created automatically!
study_folder = 'P:\Project_Sezen\data\';

% everything from here is according to the general pipeline, changes only recommended if you know
% the whole structure
raw_data_folder = '0_raw_data\';
mobi_data_folder = '1_mobi_data\';
raw_EEGLAB_data_folder = '2_raw_EEGLAB\';
spatial_filters_folder = '3_spatial_filters\';
spatial_filters_folder_AMICA = '3-1_AMICA\';
single_subject_analysis_folder = '4_single_subject_analysis\';
single_subject_analysis_folder_ERSPs = 'ERSPs\';
single_subject_analysis_folder_ERPs = 'ERPs\';
study_level = '5_study_level\';

single_subject_analysis_folder_epochs_1 = '';
merged_filename = 'merged_locs_EEG.set';
interpolated_filename = 'resampled_filtered_interpolated.set';
segments_filename = 'segments.set';
FH_cleaning_output_filename = 'filtered_clean.set';
amica_filename_input = 'resampled_filtered_interpolated_avRef.set';
amica_filename_output = 'postICA_1.set';
warped_dipfitted_filename = 'warped_dipfitted.set';
copy_weights_interpolate_avRef_filename = 'interpolated_avRef_ICA_weights.set';

epochs_filename = 'epochs.set';
study_1_filename = strcat('predError', single_subject_analysis_folder_epochs_1, '.study');

% mocap data naming
mocap_data_fname = 'merged_locs_EEG_mocap.set';

% filter frequencies
filter_lowCutoffFreqAMICA = 1;
filter_highCutoffFreqAMICA = [];
lowCutoffFreqERSP_preprocessing = [];
highCutoffFreqERSP_preprocessing = [];
lowCutoffFreqERP_preprocessing = 0.2;
highCutoffFreqERP_preprocessing = 35;

channel_locations_filename = [];
resample_freq = [];

%%% AMICA
% what is amica: https://sccn.ucsd.edu/~jason/amica_a.pdf and in general to
% better understand ICA: http://pressrelease.brainproducts.com/independent-component-analysis-demystified/

% on some PCs AMICA may crash before the first iteration if the number of
% threads and the amount the data does not suit the algorithm. Jason Palmer
% has been informed, but no fix so far. just roll with it. if you see the
% first iteration working there won't be any further crashes. in this case
% just press "close program" or the like and the bemobil_spatial_filter
% algorithm will AUTOMATICALLY reduce the number of threads and start AMICA
% again. this way you will always have the maximum number
% of threads that should be used for AMICA. check in the
% task manager how many threads you have theoretically available and think
% how much computing power you want to devote for AMICA. on the bpn-s1
% server, 12 is half of the capacity and can be used. be sure to check with
% either Ole or your supervisor and also check the CPU usage in the task
% manager before!
% 4 threads are most effective for single subject speed, more threads don't
% really shorten the calculation time much. best efficiency is using just 1
% thread and have as many matlab instances open as possible (limited by the
% CPU usage). Remember your RAM limit in this case.
max_threads = 8;
num_models = 1;

% warp electrodemontage and run dipfit
RV_threshold = 15;
remove_outside_head = 'on';
number_of_dipoles = 1;

% epoching
epochs_1_boundaries = [-1  2];
% TODO, this has to be specified and events have to be renamed after event type parsing
epochs_1_event = {'box:touched'}; 

% study
STUDY_1_components_to_use = [];

% ERSPs_1
n_times_1 = 1000;
trial_normalization_1 = 'off';
baseline_start_end_1 = [-1 0];

% fft options
fft_cycles = [3 0.5];
fft_freqrange = [3 100];
fft_padratio = 2;
fft_freqscale = 'log';
fft_alpha = NaN;
fft_powbase = NaN;
fft_c_type   = 'ersp'; % 'itc' 'both'
n_freqs = 98;

% .icaersp file options
savetrials_icaersp = 'off';
parameters_icaersp = { 'cycles', fft_cycles, 'padratio', fft_padratio, 'alpha', fft_alpha, 'freqscale', fft_freqscale};
prefix_icaersp = 'comp';
experiment_conditions_to_test_icaersp = [];
design_name_icaersp = 'design1';

% precluster
STUDY_1_clustering_weights = struct('dipoles', 6, 'scalp_topographies', 0, 'spectra', 1, 'ERSPs', 3);
STUDY_1_clustering_freqrange = [3 100];

% repeated clustering, TODO set Talairach of peak interest
outlier_sigma = 3;
STUDY_1_n_clust = 50;
n_iterations = 10000;
STUDY_1_cluster_ROI_talairach = struct('x', 0, 'y', -45, 'z', 10);
STUDY_1_quality_measure_weights = [2,-3,-1,-1,-3,-1];
do_clustering = true;
do_multivariate_data = true;
STUDY_1_filepath_clustering_solutions = '\clustering_solutions\box_touch\';
filename_clustering_solutions = 'solutions';
filepath_multivariate_data = '';
filename_multivariate_data = 'multivariate_data';