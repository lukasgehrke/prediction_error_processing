% Prediction Error Study parameters and folder structure
subjects = [2:11];

%%% Filename and folder structure informations. folders will be created automatically!
study_folder = 'P:\Project_Sezen\data\';
eeglab_path = 'M:\Toolboxes_Skripts_and_Coding_examples\eeglab-by-marius\eeglab14_1_0b';

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
single_subject_analysis_folder_ERP_images = 'ERP_images\';
study_level = '5_study_level\';

merged_filename = 'merged_EEG_incl_EMS.set';
interpolated_filename = 'locs_resampled_filtered_interpolated_incl_EMS.set';
FH_cleaning_output_filename = 'filtered_clean_incl_EMS.set';
amica_filename_output = 'postICA_1_incl_EMS.set';
warped_dipfitted_filename = 'warped_dipfitted_incl_EMS.set';
copy_weights_interpolate_avRef_filename = 'interpolated_avRef_ICA_weights_incl_EMS.set';

% filter frequencies
filter_lowCutoffFreqAMICA = 1;
filter_highCutoffFreqAMICA = [];

channel_locations_filename = [];
resample_freq = 250;

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
max_threads = 12;
num_models = 1;

% warp electrodemontage and run dipfit
RV_threshold = 15;
remove_outside_head = 'on';
number_of_dipoles = 1;
