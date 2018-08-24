% set Study parameters for ERPs of event: "box:touched"
study_params_PE_vis_vibro;

% filtering for ERPs
% todo: ask @marius or @laurens for parameter selection
lowCutoffFreqERP_preprocessing = 0.1;
highCutoffFreqERP_preprocessing = 40;

% epoch duration/boundaries
epochs_boundaries = [-1  2];
% specified epoching event, must be present in EEG.event.type
epochs_event = {'box:touched'}; 

% filenames and folders
epochs_filename = 'epochs.set';
single_subject_analysis_folder_epochs = 'ERP_box_touched';
study_filename = strcat('predError_', single_subject_analysis_folder_epochs, '.study');

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

%% (TODO) STEP L: remove undesired source projections out of data and filter for ERPs
% question to answer: can we filter after copying ICA weights or should it be done
% before?

%% STEP M: Epoch Loop
% makes epochs of the data set specified. Currently this is changed each time, that's not good. Also
% experiment-specific event fields are entered in a script where you need to specify some stuff in
% again. BE AWARE WHAT'S THE CURRENT STATE OF THE OTHER SCRIPT BEFORE STARTING THIS!

input_path = [study_folder single_subject_analysis_folder];
output_path = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERPs single_subject_analysis_folder_epochs];

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path '\' num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', copy_weights_interpolate_avRef_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    [EEG, created_epochs_indices] = pop_epoch( EEG, epochs_event, epochs_boundaries, 'newname',...
        'epochs', 'epochinfo', 'yes');
    
    mkdir(output_filepath);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew',[output_filepath '\' epochs_filename],'gui','off');
     
%     % event fields are inserted here
%     [epoch_info, EEG.event] = insert_event_fields_visual_maze(EEG.event, subject);
%     
%     % the epoching sometimes doesnt create all epochs due to data set time limits. only the used epochs
%     % should be in the epoch info struct.
%     all_epochs_info = epoch_info(created_epochs_indices);
%     
%     % save all info about each epoch in a mat file
%     save([output_filepath '\' epochs_info_filename], 'all_epochs_info')
    
end

%% STEP N: Epoch cleaning
% Status 04/07/2018: 

input_path = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERPs single_subject_analysis_folder_epochs];
output_path = input_path;

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    
    disp(['Subject #' num2str(subject)]);
    input_filepath = [input_path '\' num2str(subject)];
    output_filepath = [output_path '\' num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    oriEEG = pop_loadset('filename', epochs_filename, 'filepath', input_filepath);
    oriEEG = eeg_checkset( oriEEG );
    
    % start FH epoch cleaning
    %%% define folder names; folders for saving will be automatically created, if not existing already
    datapath_specifications.datapath_original_files = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERPs single_subject_analysis_folder_epochs '\' num2str(subject) '\'];
    datapath_specifications.datapath_save_new_files=datapath_specifications.datapath_original_files;         %%% keep last \; path for saving in single subject indices of bad epochs (no EEG data)
    datapath_specifications.datapath_save_figure_badEpochs=datapath_specifications.datapath_original_files;  %%% keep last \; path for saving in single subject figure of epoch rejection
    datapath_specifications.datapath_load_selectedChannels=[];  %%% indicate folder name for selected channels/components; e.g., if you want to analyze different IC per subject; leave empty for using all channels

    %%% file names
    filename_specifications.filename_EEG_raw='epochs.set';  %%% name of EEG.set to load
    filename_specifications.filename_saveBadEpochIndices='epochs.set'; %%% file name for saving in single subject the info/indices of rejected epochs
    filename_specifications.filename_selectedChannels=[];    %%% indicate file name for selected channels/components of a single subject; name of the .mat file; variable in the file should be named identical;  1D vector for each subject, no string; digits indicate selected channels/components

    %%% settings for automatic epoch cleaning
    automatic_cleaning_settings.cleanTimeWarpedEpochs=false;  %%% =true for epochs of different lengths, defined by timewarp; =false for epochs of same length by EEGLAB epoching
    automatic_cleaning_settings.baselineLatencyTimeWarp=[]; %%% [] for no time warp
    automatic_cleaning_settings.timeWarpEventSequence=[];  %%% warp with respect to marker names; [] for no time warp

    automatic_cleaning_settings.clean_epochs_ICA_use_components=false;  %%% =false for cleaning epochs sensor level ("channels" are EEG channels); =true for cleaning IC epochs ("channels" are components)

    automatic_cleaning_settings.deselect_channels_for_epoch_cleaning=true; %%% =true for using all available sensors/components; "false" for selecting specified channels/components per subject; cf. "datapath_load_selectedChannels" and "filename_selectedChannels"
    if ~automatic_cleaning_settings.deselect_channels_for_epoch_cleaning
        automatic_cleaning_settings.selected_channels_for_epoch_cleaning_all=[];  %%% empty [] for loading INDIVIDUAL channel/component indices per subject; else: specify [digits] to indicate fixed channels/components for all subjects
    end
    automatic_cleaning_settings.crit_keep_epochs_percent=0.9;  %%% keep 90 % of epochs, i.e. reject 10 % of the "worst" epochs
    automatic_cleaning_settings.crit_percent_sample_epoch=0.1;  %%% [] for nothing; else remove epochs if they have more than x % of samples with zero or NaN ("flat line")
    automatic_cleaning_settings.weighting_factor_epoch_cleaning=[1 2 2]; %%% method I mean of epochs, method II = channel heterogeneity --> SD across channels of mean_epochs; method III = channel heterogeneity --> Mahal. distance of mean_epochs across channels; recommended: put method I not at zero, because Mahal. might not yield results if data set is extremely short 

    % 2-4) load the EPOCHED data (EEG.set) that should be cleaned
    %%%     clean, save figures and indices of bad epochs (indices saved only! No changing or saving of EEG.data)

    [auto_epoch_cleaning]=wrapper_automatic_epoch_cleaning_EEG(datapath_specifications,filename_specifications,automatic_cleaning_settings);
    
    % copy cleaning results and save dataset
    oriEEG.etc.auto_epoch_cleaning = auto_epoch_cleaning;
    pop_saveset( oriEEG, 'filename', epochs_filename, 'filepath', output_filepath);
    close(gcf);
end

%% STEP O: create EEGLAB study structure (first rejects bad epochs)
% this is group level analyses specific and comes at the end of the data
% collection phase
% Mike X. Cohen: Part VI, Chapter 34, 35

input_path = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERPs single_subject_analysis_folder_epochs];
output_path = [study_folder study_level];

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, 'option_saveversion6', 0,...
    'option_single', 1, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1,...
    'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
    'option_checkversion', 1, 'option_chat', 1);

STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

% load all subjects into EEGLAB
for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    
    EEG = pop_loadset('filename', epochs_filename, 'filepath', input_filepath);
    % reject bad epochs
    trialrej = zeros(1,length(EEG.epoch));
    trialrej(EEG.etc.auto_epoch_cleaning.indices_bad_epochs) = 1;
    EEG = pop_rejepoch( EEG, trialrej, 0);
    
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
end

eeglab redraw
disp('All study files loaded. Creating STUDY...')

% create command which subjects and ICs to load into the study
command = cell(0,0);
for set = 1:length(subjects)
    command{end+1} = {'index' set 'subject' num2str(subjects(set)) };
end
if ~isempty(STUDY_components_to_use)
    for set = 1:length(subjects)
        command{end+1} = {'index' set 'comps' STUDY_components_to_use };
    end
end

% create study
[STUDY ALLEEG] = std_editset( STUDY, ALLEEG, 'name',study_filename,'commands',command,...
    'updatedat','on','savedat','off','rmclust','on' );
[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];

% save study
disp('Saving STUDY...')
mkdir(output_path)
[STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename',study_filename,'filepath',output_path);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
eeglab redraw
disp('...done')

%% STEP P: precompute topographies and spectra

input_path = [study_folder study_level];

if ~exist('ALLEEG','var')
    eeglab;
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, 'option_saveversion6', 0,...
        'option_single', 1, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1,...
        'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
        'option_checkversion', 1, 'option_chat', 1);
end
if isempty(STUDY)
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    [STUDY ALLEEG] = pop_loadstudy('filename', study_filename, 'filepath', input_path);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    
    eeglab redraw
end

% TODO add correct function call for this study: ERPs 
% precompute component measures except ERSPs
disp('Precomputing topographies and spectra.')
[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, 'components','recompute','on',...
    'erp', 'on', 'scalp','on','spec','on','specparams',{'specmode' 'fft' 'logtrials' 'off'},...
    'ersp', 'on', 'itc', 'on', 'erspparams',{ 'cycles' [ 3 0.5 ], 'alpha', 0.01, 'padratio' 1 });

% save study
disp('Saving STUDY...')
mkdir(output_path)
[STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename',study_filename,'filepath',output_path);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
eeglab redraw
disp('...done')

%% STEP Q: Precluster EEGLAB study

input_path = [study_folder study_level];

% what does this do?
% input_path_latencies = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERSPs single_subject_analysis_folder_epochs];
% load([input_path_latencies timewarp_name_1 '_latencyMeans'],'latencyMeans')

if ~exist('ALLEEG','var')
    eeglab;
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, 'option_saveversion6', 0,...
        'option_single', 1, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1,...
        'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
        'option_checkversion', 1, 'option_chat', 1);
end
if isempty(STUDY)
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    [STUDY ALLEEG] = pop_loadstudy('filename', study_filename, 'filepath', input_path);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    
    eeglab redraw
end

% create preclustering array that is used for clustering and save study
[STUDY, ALLEEG, EEG] = bemobil_precluster(STUDY, ALLEEG, EEG, STUDY_clustering_weights, STUDY_clustering_freqrange,...
    [-1000 2000], study_filename, input_path);

%% STEP R: Repeated clustering EEGLAB study

input_path = [study_folder study_level];
input_path_latencies = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERSPs single_subject_analysis_folder_epochs];

full_path_clustering_solutions = [STUDY_filepath_clustering_solutions num2str(STUDY_n_clust) '-cluster_' num2str(outlier_sigma)...
    '-sigma_' num2str(STUDY_clustering_weights.dipoles) '-dipoles_' num2str(STUDY_clustering_weights.spectra) '-spec_'...
    num2str(STUDY_clustering_weights.scalp_topographies) '-scalp_' num2str(STUDY_clustering_weights.ERSPs) '-ersp_' num2str(n_iterations) '-iterations'];

% load([input_path_latencies timewarp_name_1 '_latencyMeans'],'latencyMeans')

if ~exist('ALLEEG','var')
    eeglab;
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, 'option_saveversion6', 0,...
        'option_single', 1, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1,...
        'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
        'option_checkversion', 1, 'option_chat', 1);
end
if isempty(STUDY)
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    [STUDY ALLEEG] = pop_loadstudy('filename', study_filename, 'filepath', input_path);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    
    eeglab redraw
end

% cluster the components repeatedly and use a region of interest and
% quality measures to find the best fitting solution
[STUDY, ALLEEG, EEG] = bemobil_repeated_clustering_and_evaluation(STUDY, ALLEEG, EEG, outlier_sigma,...
    STUDY_n_clust, n_iterations, STUDY_cluster_ROI_talairach, STUDY_quality_measure_weights, false,...
    false, input_path, study_filename, [input_path full_path_clustering_solutions],...
    filename_clustering_solutions, [input_path full_path_clustering_solutions], filename_multivariate_data);