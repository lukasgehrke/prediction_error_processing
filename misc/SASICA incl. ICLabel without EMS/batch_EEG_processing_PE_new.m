%% STEP K: Create 1st Level data sets with SASICA & ICLabel ix of bad comps
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
    
    % load AMICA file (first file in EEGLAB) NOT including dipole fitting
    AMICA_EEG = pop_loadset('filename', amica_filename_output, 'filepath', spatial_filter_filepath);
    AMICA_EEG = eeg_checkset( AMICA_EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, AMICA_EEG, 0 );
    
    % Find eye component indices using SASICA 
    [AMICA_EEG,settings_SASICA] = eeg_SASICA(AMICA_EEG,SASICA_settings);
    
    % Classify components using IClabel
    AMICA_EEG = pop_iclabel(AMICA_EEG); % pop_viewprops( EEG, 0, [1:144], {'freqrange', [2 80]}, {}, 1, 'ICLabel' )
    
    % load preprocessed, interpolated, average referenced file
    EEG = pop_loadset('filename', interpolated_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    
    % add eye comps indeces determined on the AMICA set by SASICA to the
    % EEG dataset
    EEG.etc.sasica.components_rejected = find(AMICA_EEG.reject.gcompreject);
    EEG.etc.sasica.SASICA_settings = SASICA_settings;
    EEG.etc.iclabel = AMICA_EEG.etc.ic_classification;
    
    % parse events
    EEG = parse_events_PredError(EEG);    
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % copy AMICA weights and NO dipfit info
    bemobil_copy_spatial_filter(EEG, ALLEEG, CURRENTSET, AMICA_EEG, true, false, copy_weights_interpolate_avRef_filename , output_filepath);
    
end