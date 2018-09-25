% AFTER ICA continue here
% set Study parameters
study_params_PE_vis_vibro;

%% STEP J (OPTIONAL): Warping of locations and dipole fitting

% ToDo: to be used in a later step of the analysis, section will be relocated

% see Mike X. Cohen: Part IV, 24 "Basics of Single Dipole..."
% renames the specified channels, warps the chanlocs on a standard head model and fits dipoles for
% each IC below the threshold of residual variance

% input_path = [study_folder spatial_filters_folder spatial_filters_folder_AMICA];
% output_path = input_path;

% if ~exist('ALLEEG','var'); eeglab; end
% pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

% for subject = subjects
%     disp(['Subject #' num2str(subject)]);
    
%     input_filepath = [input_path num2str(subject)];
%     output_filepath = [output_path num2str(subject)];
%     
%     STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
%     EEG = pop_loadset('filename', amica_filename_output, 'filepath', input_filepath);
%     EEG = eeg_checkset( EEG );
%     [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
%     warping_channel_names = [];
    
    % do the warp and dipfit
%     bemobil_dipfit( EEG , ALLEEG, CURRENTSET, warping_channel_names, eeglab_path, RV_threshold,...
%         remove_outside_head, number_of_dipoles, warped_dipfitted_filename, output_filepath);
    
% end

%% STEP K: Create 1st Level data sets
% by copying the spatial filter information (including dipfit info) into a
% fresh data set which is just channel-domain interpolated and average
% referenced, but NOT cleaned in the time domain. No filter is applied yet,
% do this as needed

input_path = [study_folder raw_EEGLAB_data_folder];
spatial_filtering_weights_path = [study_folder spatial_filters_folder spatial_filters_folder_AMICA];
output_path = [study_folder single_subject_analysis_folder];

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    spatial_filter_filepath = [spatial_filtering_weights_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    % load AMICA file (first file in EEGLAB) including dipole fitting
     % AMICA_EEG = pop_loadset('filename', warped_dipfitted_filename, 'filepath', spatial_filter_filepath);
     % AMICA_EEG = eeg_checkset( AMICA_EEG );
     % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, AMICA_EEG, 0 );
    
    % load AMICA file (first file in EEGLAB) NOT including dipole fitting
    AMICA_EEG = pop_loadset('filename', amica_filename_output, 'filepath', spatial_filter_filepath);
    AMICA_EEG = eeg_checkset( AMICA_EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, AMICA_EEG, 0 );
    
    % switch case here for IC selection and rejection
    switch subject
        case 2
            eyeICs = [1 2 26]; %26 has some muscle artefacts as well.
        case 3
            eyeICs = []; %the same error message: surface error: not a valid SurfaceCData
        case 4
            eyeICs = [1 2];
        case 5
            eyeICs = [2 5 8 9 12]; %after new ICA, unfortunately the same ICs are generated.
        case 6
            eyeICs = [1 3 4]; %3 and 4 are almost identical
        case 7
            eyeICs = [1 2 3];
        case 8
            eyeICs = [1 2 3];
        case 9
            eyeICs = [1];
        case 10
            eyeICs = [2 5 6]; %5 and 6 are almost identical
        case 11
            eyeICs = [1 2]; %30 ICs only
        case 12
            eyeICs = [1 4 6 8]; % idx of eye components to be removed
    end
    
    % load preprocessed, interpolated, average referenced file
    EEG = pop_loadset('filename', interpolated_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    EEG.etc.eye_ICs = eyeICs;
    
    % parse events
    EEG = parse_events_PredError(EEG);    
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % copy AMICA weights and NO dipfit info
    bemobil_copy_spatial_filter(EEG, ALLEEG, CURRENTSET, AMICA_EEG, true, false, copy_weights_interpolate_avRef_filename , output_filepath);
    
end
