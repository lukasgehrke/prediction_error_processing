% set Study parameters
study_params_PE_incl_EMS;

%% STEP E: Merge all datasets ALL conditions! 
% This merges all EEG data files into one big file

input_path = [study_folder raw_EEGLAB_data_folder];
output_path = input_path;

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    fnames = dir([input_path '\' num2str(subject)]);
    for i=1:length(fnames)
        if ~contains(fnames(i).name, '_EEG.set') || contains(fnames(i).name, 'Training') || contains(fnames(i).name, 'merged')
            fnames(i).name = [];
        end
    end
    fnames = {fnames.name};
    fnames = fnames(~cellfun('isempty',fnames));
            
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', fnames, 'filepath', input_filepath);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
    
    % merges all files currently loaded in EEGLAB into one file and stores
    % the original filenames in EEG.etc.appended_files
    [ALLEEG EEG CURRENTSET] = bemobil_merge(ALLEEG,EEG,CURRENTSET,[1:length(ALLEEG)], merged_filename, output_filepath);
    
end

%% STEP F: Locations, Resampling, high-pass filter, Channel rejection, interpolation of channels, Re-referencing
% Done manually in this case, but you can replace this loop with a loop
% that uses automatic cleaning!
% @Sezen: Mike X. Cohen, Part II, 7: Preprocessing

input_path = [study_folder raw_EEGLAB_data_folder '\'];
output_path = input_path;

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

%%% executing the cleaning for all subjects
for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', merged_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % resampling code
    EEG = pop_resample(EEG, resample_freq);    
    
    % 1Hz Highpass (effectivelz 1-124 Hz bandpass with 250Hz sampling rate)
    [ALLEEG, EEG, CURRENTSET] = bemobil_filter(ALLEEG, EEG, CURRENTSET, filter_lowCutoffFreqAMICA, filter_highCutoffFreqAMICA);
        
    % change electrode labels and load default electrode positions
    for i = 1:length(EEG.chanlocs)
        EEG.chanlocs(i).labels = erase(EEG.chanlocs(i).labels, 'brainvision_rda_bpn-c012_');
    end
    % insert FCz channel
    EEG=pop_chanedit(EEG, 'insert',64,'changefield',{64 'labels' 'FCz'});
    EEG=pop_chanedit(EEG, 'lookup','P:\\Lukas_Gehrke\\toolboxes\\eeglab\\plugins\\dipfit2.4\\standard_BESA\\standard-10-5-cap385.elp','rplurchanloc');
    EEG = eeg_checkset( EEG );
    
    % Todo: do manual channel rejection, add subjects from vis vibro script
    switch subject
        case 2
            removed_chans = [4 16];
        case 3
            removed_chans = [9 10 55 60];
        case 4
            removed_chans = [41];
        case 5
            removed_chans = [1 33 41 42];
        case 6
            removed_chans = [9 16 43 46 10 14];
            %FC4 only bad after 2600ms, not yet removed
        case 7
            removed_chans = [17 32 49];
        case 8
            removed_chans = [41 42 62 63 9 17 55];
        case 9
            removed_chans = [12 41 46];
        case 10
            removed_chans = [42 45 41 33 17];
        case 11
            removed_chans = [22];
        case 12
            removed_chans = [2 22 31 64]; % idx of channels to be removed
        case 13
            removed_chans = [7, 16, 22, 30, 46, 48];
        case 14
            removed_chans = [2, 29, 30, 3, 12];
        case 15
            removed_chans = [8, 12, 32, 40, 46, 5, 50];
        case 16
            removed_chans = [45, 60, 57, 41, 42, 32, 28, 16];
    end
    
    % if interested try the below methods for automatic channel rejection
%     % automatic bad channel rejection
%     reject_chans_EEG = clean_channels(rerefEEG, 0.5);
%     % find removed chans for interpolation
%     try
%         removed_chans = find(~reject_chans_EEG.etc.clean_channel_mask);
%     catch
%         disp('No channels need to be interpolated');
%     end
%     % different approach
%     % [newEEG removed_chans] = pop_rejchan(EEG, 'elec',[1:64] ,'threshold',5,'norm','on','measure','kurt');

    if ~isempty(removed_chans)
        [ALLEEG, EEG, CURRENTSET] = bemobil_interp( EEG , ALLEEG, CURRENTSET, removed_chans);
        % -> interpolated bad channels
        % the function stores the interpolated channels in EEG.etc and saves a
        % data set
    end
    
    % Compute average reference, after adding additional channel for averef
    %EEG = fullRankAveRef( EEG );
    EEG = pop_reref( EEG, [1:64] ,'refloc',struct('labels',{'FCz'},'type',{'EEG'},'theta',{0},'radius',{0.12662},'X',{32.9279},'Y',{0},'Z',{78.363},'sph_theta',{0},'sph_phi',{67.208},'sph_radius',{85},'urchan',{65},'ref',{''},'datachan',{0}),'keepref','on');
    
    % save dataset
    mkdir(output_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(output_filepath);
    if ismember(interpolated_filename, {dir_files.name})
        warning([interpolated_filename ' file already exists in: ' output_filepath '. File will be overwritten...']);
    end
    
    pop_saveset( EEG, 'filename', interpolated_filename, 'filepath', output_filepath);
end

disp('Channel based cleaning and interpolation done!')

%% STEP G: Automatic time domain cleaning on the channel level
% Possible: Segmentation of data set to throw out pre- and post-experiment 
% times as well as breaks for AMICA
% you can delete events before this if they hinder your vision, this data
% set is only for AMICA:

input_path = [study_folder raw_EEGLAB_data_folder];
output_path = [study_folder spatial_filters_folder spatial_filters_folder_AMICA];

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects

    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    oriEEG = pop_loadset('filename', interpolated_filename, 'filepath', input_filepath);
    oriEEG = eeg_checkset( oriEEG );
    
    %%% specify folder names
    datapath_specifications.datapath_original_files=input_filepath; %%% keep last \; single subject folder
    datapath_specifications.datapath_save_files=output_filepath;     %%% keep last \; path for saving updated EEG
    datapath_specifications.datapath_save_figures=output_filepath;   %%% keep last \; path for saving figures of cleaning

    %%% specify file names
    filename_specifications.file_name_original_EEG=interpolated_filename;   %%% loads "fresh" EEG (raw, unfiltered)
    filename_specifications.filename_saveBadEpochIndices='';  

    automatic_cleaning_settings.cleaned_data_type='sensor data'; %%% ICA not implemented yet; usually bad segments found on sensor level are also fine for IC later on

    %%% select channels that should be considered for cleaning
    automatic_cleaning_settings.selected_sensor_channels_for_cleaning=[]; %%% [] use all available channels for cleaning, else specify [1 2 ...]; currently same channels for all subjects alike
    %%% select channel(s) for cleaning before vs. after
    automatic_cleaning_settings.chan_select_plot_before_after=[5];  %%% [] use all available channels for cleaning

    %%% define frequency band (band-pass filter) only for cleaning 
    automatic_cleaning_settings.band_artifact_cleaning=[1 40];  %%% in Hz; [] for no filter
    automatic_cleaning_settings.band_stop_artifact_cleaning=[]; %%% [] for nothing; else e.g. [48 52] for removal of line artifacts (notch filter)
    automatic_cleaning_settings.band_filtorder=2;               %%% for IIR Butterworth filter; since filtfilt is used, the effective order is double (here 4)

    %%% further settings for the cleaning algorithm
    automatic_cleaning_settings.analyzed_files_N=1;  %%% so far: analyze only 1 file! (no appended different conditions)
    automatic_cleaning_settings.crit_all=[0.85]; %%% e.g., 0.9=90% keep amount of epochs; 1 value if only 1 file (no appended recordings); else indicate, e.g., 4 values for 4 appended files 
    automatic_cleaning_settings.wind_ms=1000;    %%% in ms, epochs for finding large artifacts
    automatic_cleaning_settings.crit_keep_sec=4; %%% in seconds; for NOT removing any additional data, set: automatic_cleaning_settings.crit_keep_sec=automatic_cleaning_settings.wind_ms/1000; else: value should be multiple of "wind_ms" for additionally removing data, i.e., keep uninterrupted "good" data segments not shorter than this value
    automatic_cleaning_settings.crit_percent_sample_epoch=0.1;  %%% [] for nothing; e.g., 0.1 = 10%; remove epochs if they have more than x % of samples with zero or NaN ("flat line")
    automatic_cleaning_settings.weighting_factor_epoch_cleaning_methods=[1 1 2];     %%% method I mean of epochs, method II = channel heterogeneity --> SD across channels of mean_epochs; method III = channel heterogeneity --> Mahal. distance of mean_epochs across channels; recommended: put method I not at zero, because Mahal. might not yield results if data set is extremely short 
    automatic_cleaning_settings.visual_inspection_mode=false;  %%% =false if visual threshold rejection after automatic cleaning should not be applied; in this case, bad segments from previous automatic artifact rejection are taken
    if ~automatic_cleaning_settings.visual_inspection_mode
        automatic_cleaning_settings.threshold_visual_reject=zeros(1,automatic_cleaning_settings.analyzed_files_N);
    end
    
    % 2-5) loading, cleaning, and saving (results, figures) 
    %%% only indices are saved in the struct "auto_continuous_cleaning"
    %%% so far: before/after comparisons are saved as .jpg only, because
    %%% otherwise usually files are too large with continuous data. However, if
    %%% you open the "wrapper_automatic_cleaning_continuous_EEG" and run each
    %%% section manually, you can save the same figure as .fig (it's just disabled by
    %%% now). Please do not enable it in the wrapper itself! 
    addpath 'P:\Project_Friederike\2017_spot_rotation\1_analysis\analysis_Matlab_diaries_scripts_releases\3_release_internal_use_only\scripts\'
    [auto_continuous_cleaning]=wrapper_automatic_cleaning_continuous_EEG(datapath_specifications,filename_specifications,automatic_cleaning_settings);
    
    % copy cleaning results and save dataset
    oriEEG.etc.auto_continuous_cleaning = auto_continuous_cleaning;
    pop_saveset( oriEEG, 'filename', FH_cleaning_output_filename, 'filepath', output_filepath);
    close(gcf);
end

%% STEP H: ICA loop 1
% what is amica: https://sccn.ucsd.edu/~jason/amica_a.pdf and in general to
% better understand ICA: http://pressrelease.brainproducts.com/independent-component-analysis-demystified/

% BEST TO USE 1 THREAD! So far, no crashes ocurred with 1 thread, and it's
% most efficient. Start several MATLAB instances and run AMICA in parallel
% for different subjects. Remember your CPU and RAM limit!

% ToDo before ICA:

% 1. preprocess data
% ('preprocessed_poss_wrong_chanlocs.set')

% 2. interpolate bad channel locations (only if necessary)
% 3. interchange channels back to correct (only if necessary)
% ('preprocessed.set')

% 4. interpolate bad channels
% ('interpolated.set')

% 5. average reference
% ('interpolated_avRef.set')

% 6. take only experiment segments
% ('segments.set')

% 7. manually and/or automatically clean for artifacts
% ('pre_ICA_1.set')

% -> Then off you go :)

%1: 4     5     6     7
%2: 12    13    14    15    16
%3: 18    20    21    22    23
%4: 25    26    31    32

input_path = [study_folder spatial_filters_folder spatial_filters_folder_AMICA];
output_path = input_path;
if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', FH_cleaning_output_filename , 'filepath', input_filepath);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
    
    % reject data from previously computed FH channel data rejection
    invalid_segments_index = EEG.etc.auto_continuous_cleaning.invalid_segments_final_start_stop_sample;
    EEG = eeg_eegrej(EEG, invalid_segments_index);
    
    % rank is reduced by the number of interpolated channels
    rankEEG = rank(EEG.data');
    channelSubset = loc_subsets(EEG.chanlocs, rankEEG);
    EEG = pop_select( EEG,'channel', channelSubset{1});
    EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
    
    % running signal decomposition with values specified above
    [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, ...
        CURRENTSET, true, num_models, max_threads, rankEEG, [], ...
        amica_filename_output, [output_filepath]);
    
end

%% STEP I (OPTIONAL): Warping of locations and dipole fitting

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

%% STEP J: Create 1st Level data sets with SASICA & ICLabel ix of bad comps
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
