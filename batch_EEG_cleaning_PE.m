% set Study parameters
study_params_PE_vis_vibro;

%% STEP G: Remove segments of non-experiment data
% now, reject data not in experiment condition, i.e. between blocks,
% before and after starting the experiment

input_path = [study_folder raw_EEGLAB_data_folder];
output_path = input_path;

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

%%% executing the cleaning for all subjects
for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', interpolated_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    EEG = parse_events_PredError(EEG);
    
    % extract data between first spawn and last touch
    trials_spawn = [1 101 201];
    if subject == 3
        trials_spawn = [1 102 201]; % data collection started after the first trial was completed
    end
    ind_spawn = find(strcmp('spawned', {EEG.event.box}));
    trials_touch = [100 200 300];
    ind_touch = find(strcmp('touched', {EEG.event.box}));
    spawns = [];
    touches = [];

    for i = 1:length(trials_spawn)
        ind_tr = find(strcmp(num2str(trials_spawn(i)) ,{EEG.event.trial_nr}));
        spawns = [spawns intersect(ind_spawn, ind_tr)]; 

        ind_tr = find(strcmp(num2str(trials_touch(i)) ,{EEG.event.trial_nr}));
        touches = [touches intersect(ind_touch, ind_tr)];
    end
    out_dat_tmp = sort([spawns, touches]);
    for i = 1:length(out_dat_tmp)
        out_dat(i) = EEG.event(out_dat_tmp(i)).latency;
    end
    out_dat = [1, out_dat, EEG.pnts];
    out_dat = reshape(out_dat, [2,length(spawns)+1])';

    EEG = eeg_eegrej(EEG, out_dat);
    pop_saveset( EEG, 'filename', segments_filename, 'filepath', output_filepath);
    clear out_dat
end

disp('Remove segments cleaning done!')

%% STEP H: Automatic time domain cleaning on the channel level
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
    
    oriEEG = pop_loadset('filename', segments_filename, 'filepath', input_filepath);
    oriEEG = eeg_checkset( oriEEG );
    
    %%% specify folder names
    datapath_specifications.datapath_original_files=input_filepath; %%% keep last \; single subject folder
    datapath_specifications.datapath_save_files=output_filepath;     %%% keep last \; path for saving updated EEG
    datapath_specifications.datapath_save_figures=output_filepath;   %%% keep last \; path for saving figures of cleaning

    %%% specify file names
    filename_specifications.file_name_original_EEG=segments_filename;   %%% loads "fresh" EEG (raw, unfiltered)
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
    automatic_cleaning_settings.crit_all=[0.90]; %%% e.g., 0.9=90% keep amount of epochs; 1 value if only 1 file (no appended recordings); else indicate, e.g., 4 values for 4 appended files 
    automatic_cleaning_settings.wind_ms=1000;    %%% in ms, epochs for finding large artifacts
    automatic_cleaning_settings.crit_keep_sec=6; %%% in seconds; for NOT removing any additional data, set: automatic_cleaning_settings.crit_keep_sec=automatic_cleaning_settings.wind_ms/1000; else: value should be multiple of "wind_ms" for additionally removing data, i.e., keep uninterrupted "good" data segments not shorter than this value
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

%% STEP I: ICA loop 1
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
