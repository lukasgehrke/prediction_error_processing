% set Study parameters for ERPs of event: "box:touched"
addpath 'P:\Project_Sezen\data_processing';
study_params_PE_incl_EMS;

% filtering for ERPs
lowCutoffFreqERP_preprocessing = 0.2;
highCutoffFreqERP_preprocessing = 35;

% epoch duration/boundaries
epochs_boundaries = [-0.3  0.7];
% specified epoching event, must be present in EEG.event.type
epochs_event = {'box:touched'}; 

% filenames and folders
epochs_filename = 'epochs_incl_EMS.set';
single_subject_analysis_folder_epochs = 'ERP_box_touched';
study_filename = strcat('predError_', single_subject_analysis_folder_epochs, '_incl_EMS', '.study');

% study parameters
STUDY_components_to_use = [];

% precluster
STUDY_clustering_weights = struct('dipoles', 6, 'scalp_topographies', 0, 'spectra', 1, 'ERPs', 3);
STUDY_clustering_freqrange = [3 100];

% repeated clustering, TODO set Talairach of peak interest
outlier_sigma = 3;
STUDY_n_clust = 50;
n_iterations = 10000;
STUDY_cluster_ROI_talairach = struct('x', 0, 'y', -45, 'z', 10);
STUDY_quality_measure_weights = [2,-3,-1,-1,-3,-1];
do_clustering = true;
do_multivariate_data = true;
STUDY_filepath_clustering_solutions = '\clustering_solutions\box_touch\';
filename_clustering_solutions = 'solutions';
filepath_multivariate_data = '';
filename_multivariate_data = 'multivariate_data';