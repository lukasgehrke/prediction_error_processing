% Specific batch script for Prediction Error Experiment
% using EEGLAB & MoBILAB & Marius Klug's MoBILAB extensions (https://github.com/MariusKlug/mobilab)
% and BeMoBIL pipeline functions and experiment specific scripts &
% functions

% Versions:
% MATLAB Version: 9.2.0.538062 (R2017a)
% MATLAB License Number: FreeForAll
% Operating System: Microsoft Windows 10 Pro Version 10.0 (Build 17134)

% EEGLAB

% MoBILAB

%% START: set Study parameters
subjects = 2; 

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
preprocessed_filename = 'preprocessed.set';
interpolated_filename = 'interpolated.set';
interpolated_avRef_filename = 'interpolated_avRef.set';
filtered_filename = 'filtered.set';
segments_filename = 'segments.set';
copy_weights_interpolate_avRef_filename = 'interpolated_avRef_ICA_weights.set';
amica_filename_input = 'filtered_clean.set';
amica_filename_output = 'postICA_1.set';
warped_dipfitted_filename = 'warped_dipfitted.set';
merged_amica_filename = 'merged_ICA.set';

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
max_threads = 4;
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

%% STEP A: MoBILAB import loop
% what is mobilab? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3942646/
% here the .xdf files of the experiment are loaded into mobilab and the
% according folders are created.

input_path = [study_folder raw_data_folder];
output_path = [study_folder mobi_data_folder];

if ~exist('mobilab','var'); eeglab; runmobilab; end

for subject = subjects
    
    disp(['Subject #' num2str(subject)]);
    fnames = dir([input_path '\' num2str(subject)]);
    
    for file=1:length(fnames)
        if contains(fnames(file).name, '.xdf')
            disp(['Importing file: ' fnames(file).name '...']);

            input_filepath = [input_path '\' num2str(subject) '\' fnames(file).name];
            output_filepath = [output_path '\' num2str(subject) '\' fnames(file).name(1:end-4) '_MoBI'];

            mobilab.allStreams = dataSourceXDF(input_filepath,output_filepath);
            disp('...done.');
        end
    end
end

%% STEP B: MoBILAB processing loop
% Movement data preprocessing
% the _MoBI data folders have to have been cleaned to before the first processing

input_path = [study_folder mobi_data_folder];

if ~exist('mobilab','var'); eeglab; runmobilab; end

for subject = subjects
    
    disp(['Subject #' num2str(subject)]);
    fnames = dir([input_path '\' num2str(subject)]);
    
    for file = 1:length(fnames)
        if contains(fnames(file).name, 'MoBI')

            input_filepath = [input_path num2str(subject) '\' fnames(file).name];
            mobilab.allStreams = dataSourceMoBI(input_filepath);

            % find hand and head RBs                
            for i=1:length(mobilab.allStreams.item)
                if strcmp(mobilab.allStreams.item{i}.name, 'rigid_handr_BPN-C043')
                    hand = i;
                elseif strcmp(mobilab.allStreams.item{i}.name, 'rigid_head_BPN-C043')
                    head = i;
                end
            end
            
            % quaternion values sometimes flip their sign for mathematical reasons
            % which is bad for filtering. but the values stay the same if they are
            % flipped back, so this is done to allow filtering.
            unflipHand = mobilab.allStreams.item{hand}.unflipSigns();
            unflipHead = mobilab.allStreams.item{head}.unflipSigns();

            % lowpass filtering with the specified cutoff frequency
            lowpassHand = unflipHand.lowpass(6);
            lowpassHead = unflipHead.lowpass(6);

            % quaternion orientation values are transformed to euler angle to be
            % interpretable for humans
            eulerHand = lowpassHand.quaternionsToEuler();
            eulerHead = lowpassHead.quaternionsToEuler();

            % 3 time derivatives are calculated (velocity, acceleration, and jerk)
            eulerHand.timeDerivative(3);
            eulerHead.timeDerivative(3);
                       
            % ugly, throwout 2 channels repeatedly for 4 channel vel data
            % because
            % mobilab.allStreams.item{hand+9}.throwOutChannels(2:3); not
            % working
            % TODO redo throwOutChannels with all channels but 1 which is
            % then overwritten with 3D magnitude
            mobilab.allStreams.item{hand+9}.throwOutChannels(2);
            mobilab.allStreams.item{hand+9}.children{1}.throwOutChannels(2);
            
            % create event markers based on 3D hand velocity thresholding
            % get velocity data of hand RB
            vel = mobilab.allStreams.item{hand+9};
            vel_xyz_newChannel = sqrt(vel.data(:,1).^2 + vel.data(:,2).^2 + vel.data(:,3).^2);
            
            % add channel to the matching length RB stream
            % not working: vel.addChannels(1, vel_xyz_newChannel);
            mobilab.allStreams.item{hand+9}.children{1}.children{1}.data(:,1) = vel_xyz_newChannel;
            
            % explore velocity data for later parameter selection in
            % onset/offset detection
%             figure;histogram(vel_xyz_newChannel);
%             figure;plot(sort(vel_xyz_newChannel));
%             prctile(vel_xyz_newChannel, [90:1:100]);
%             figure;histogram(mobilab.allStreams.item{hand+9}.data(:,1), 20);
            
            mobilab.allStreams.item{hand+9}.children{1}.children{1}.createEventsFromMagnitude(1,...
                'movements','hand_movement:start hand_movement:end',0.7,0,90,3,200);
            
            % for head; just exploration, no hypotheses/interest as of yet
%             head_euler = mobilab.allStreams.item{head+11}.data(:,4);
%             figure;histogram(abs(head_euler));
%             figure;plot(sort(abs(head_euler)));
%             prctile(abs(head_euler), [1:1:100])
%             figure;histogram(mobilab.allStreams.item{head+9}.data(:,4), 20);
            
            %mobilab.gui
        end
    end
end

disp('mobilab processing done!')

%% STEP C: MoBILAB export loop
% the whole mobilab data has to be exported to EEGLAB in order to be
% processed further. EEG data, mocap data, experimental markers, and mocap
% markers are aligned and exported as .set data file

input_path = [study_folder mobi_data_folder];
output_path = [study_folder raw_EEGLAB_data_folder];

if ~exist('ALLEEG','var'); eeglab; end
if ~exist('mobilab','var'); runmobilab; end

pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    
    disp(['Subject #' num2str(subject)]);
    output_filepath = [output_path num2str(subject)];
    mkdir(output_filepath); % create folder
    
    fnames = dir([input_path '\' num2str(subject)]);
    
    for file = 1:length(fnames)    
        if contains(fnames(file).name, 'MoBI')
        
            % get filenames in folder to save the files with postfix _MOBI
            input_filepath = [input_path num2str(subject) '\' fnames(file).name];
            mobilab.allStreams = dataSourceMoBI(input_filepath);

            disp(['Exporting Subject #' num2str(subject) ': ' fnames(file).name '...']);

            % the first vector are the streams that should be included as data streams, the second the
            % streams, from which the markers should be taken. 1 is the EEG
            % data, 2 the experimental event marker stream and the others are
            % mocap data and marker. this has to first be checken in mobilab.
            % mocap marker from head vel and acc, and from PS vel if present.

            % find streams to export, see mobilab export loop
            for i=1:length(mobilab.allStreams.item)
                switch mobilab.allStreams.item{i}.name
                    case 'rigid_handr_BPN-C043'
                        hand = i;
                    case 'rigid_head_BPN-C043'
                        head = i;
                    case 'brainvision_rda_bpn-c012'
                        eeg = i;
                    case 'unity_markers_prederror_BPN-C043'
                        exp_marker = i;
                    case 'quat2eul_filt_unflip_rigid_handr_BPN-C043'
                        hand_dat = i;
                    case 'vel_quat2eul_filt_unflip_rigid_handr_BPN-C043'
                        hand_dat_deriv = i;
                    case 'throwOut_throwOut_vel_quat2eul_filt_unflip_rigid_handr_BPN-C043'
                        hand_mark = i;
                    case 'quat2eul_filt_unflip_rigid_head_BPN-C043'
                        head_dat = i;
                    case 'vel_quat2eul_filt_unflip_rigid_head_BPN-C043'
                        head_dat_deriv = i;                        
                end
            end
            
            EEG = mobilab.allStreams.export2eeglab( [eeg hand_dat hand_dat_deriv head_dat head_dat_deriv],...
                [exp_marker hand_mark]);

            % the exported data set has the suffix _MoBI, since it contains both
            % brain and body imaging data
            EEG = pop_saveset( EEG, 'filename', fnames(file).name, 'filepath', output_filepath);
            disp('...done');
        end
    end
    
end

%% STEP D: Separating EEG and Mocap loop
% since most of the EEGLAB functions are only useful for EEG data, the
% _MoBI data set has to be split up into EEG and mocap data

input_path = [study_folder raw_EEGLAB_data_folder];
output_path = input_path;
if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject) '\'];
    fnames = dir([input_path '\' num2str(subject)]);
    
    for file = 1:length(fnames)
        if contains(fnames(file).name, '.set')
            % clean EEGLAB before each iteration
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

            % load the _MoBI set which has EEG and mocap data combined
            EEG = pop_loadset('filename', fnames(file).name, 'filepath', input_filepath);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
            EEG = eeg_checkset( EEG );

            % split it and keep just the EEG channels, save it with _EEG suffix

            % check correct channels to kick out.
            % TODO add standard channel names of 64 electrode layout
            EEG = pop_select( EEG,'channel',{'brainvision_rda_bpn-c012_Fp1' 'brainvision_rda_bpn-c012_Fp2' 'brainvision_rda_bpn-c012_F7' 'brainvision_rda_bpn-c012_F3' 'brainvision_rda_bpn-c012_Fz' 'brainvision_rda_bpn-c012_F4' 'brainvision_rda_bpn-c012_F8' 'brainvision_rda_bpn-c012_FC5' 'brainvision_rda_bpn-c012_FC1' 'brainvision_rda_bpn-c012_FC2' 'brainvision_rda_bpn-c012_FC6' 'brainvision_rda_bpn-c012_T7' 'brainvision_rda_bpn-c012_C3' 'brainvision_rda_bpn-c012_Cz' 'brainvision_rda_bpn-c012_C4' 'brainvision_rda_bpn-c012_T8' 'brainvision_rda_bpn-c012_TP9' 'brainvision_rda_bpn-c012_CP5' 'brainvision_rda_bpn-c012_CP1' 'brainvision_rda_bpn-c012_CP2' 'brainvision_rda_bpn-c012_CP6' 'brainvision_rda_bpn-c012_TP10' 'brainvision_rda_bpn-c012_P7' 'brainvision_rda_bpn-c012_P3' 'brainvision_rda_bpn-c012_Pz' 'brainvision_rda_bpn-c012_P4' 'brainvision_rda_bpn-c012_P8' 'brainvision_rda_bpn-c012_PO9' 'brainvision_rda_bpn-c012_O1' 'brainvision_rda_bpn-c012_Oz' 'brainvision_rda_bpn-c012_O2' 'brainvision_rda_bpn-c012_PO10' 'brainvision_rda_bpn-c012_AF7' 'brainvision_rda_bpn-c012_AF3' 'brainvision_rda_bpn-c012_AF4' 'brainvision_rda_bpn-c012_AF8' 'brainvision_rda_bpn-c012_F5' 'brainvision_rda_bpn-c012_F1' 'brainvision_rda_bpn-c012_F2' 'brainvision_rda_bpn-c012_F6' 'brainvision_rda_bpn-c012_FT9' 'brainvision_rda_bpn-c012_FT7' 'brainvision_rda_bpn-c012_FC3' 'brainvision_rda_bpn-c012_FC4' 'brainvision_rda_bpn-c012_FT8' 'brainvision_rda_bpn-c012_FT10' 'brainvision_rda_bpn-c012_C5' 'brainvision_rda_bpn-c012_C1' 'brainvision_rda_bpn-c012_C2' 'brainvision_rda_bpn-c012_C6' 'brainvision_rda_bpn-c012_TP7' 'brainvision_rda_bpn-c012_CP3' 'brainvision_rda_bpn-c012_CPz' 'brainvision_rda_bpn-c012_CP4' 'brainvision_rda_bpn-c012_TP8' 'brainvision_rda_bpn-c012_P5' 'brainvision_rda_bpn-c012_P1' 'brainvision_rda_bpn-c012_P2' 'brainvision_rda_bpn-c012_P6' 'brainvision_rda_bpn-c012_PO7' 'brainvision_rda_bpn-c012_PO3' 'brainvision_rda_bpn-c012_POz' 'brainvision_rda_bpn-c012_PO4' 'brainvision_rda_bpn-c012_PO8'});

            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[fnames(file).name(1:end-4) '_EEG'],'savenew',[output_filepath fnames(file).name(1:end-4) '_EEG.set'],'gui','off');

            % go back to the first data set
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',1,'study',0);
            EEG = eeg_checkset( EEG );

            % split it and kick all EEG data channels out, save it with _mocap suffix
            % TODO add standard channel names of 64 electrode layout
            EEG = pop_select( EEG,'nochannel',{'brainvision_rda_bpn-c012_Fp1' 'brainvision_rda_bpn-c012_Fp2' 'brainvision_rda_bpn-c012_F7' 'brainvision_rda_bpn-c012_F3' 'brainvision_rda_bpn-c012_Fz' 'brainvision_rda_bpn-c012_F4' 'brainvision_rda_bpn-c012_F8' 'brainvision_rda_bpn-c012_FC5' 'brainvision_rda_bpn-c012_FC1' 'brainvision_rda_bpn-c012_FC2' 'brainvision_rda_bpn-c012_FC6' 'brainvision_rda_bpn-c012_T7' 'brainvision_rda_bpn-c012_C3' 'brainvision_rda_bpn-c012_Cz' 'brainvision_rda_bpn-c012_C4' 'brainvision_rda_bpn-c012_T8' 'brainvision_rda_bpn-c012_TP9' 'brainvision_rda_bpn-c012_CP5' 'brainvision_rda_bpn-c012_CP1' 'brainvision_rda_bpn-c012_CP2' 'brainvision_rda_bpn-c012_CP6' 'brainvision_rda_bpn-c012_TP10' 'brainvision_rda_bpn-c012_P7' 'brainvision_rda_bpn-c012_P3' 'brainvision_rda_bpn-c012_Pz' 'brainvision_rda_bpn-c012_P4' 'brainvision_rda_bpn-c012_P8' 'brainvision_rda_bpn-c012_PO9' 'brainvision_rda_bpn-c012_O1' 'brainvision_rda_bpn-c012_Oz' 'brainvision_rda_bpn-c012_O2' 'brainvision_rda_bpn-c012_PO10' 'brainvision_rda_bpn-c012_AF7' 'brainvision_rda_bpn-c012_AF3' 'brainvision_rda_bpn-c012_AF4' 'brainvision_rda_bpn-c012_AF8' 'brainvision_rda_bpn-c012_F5' 'brainvision_rda_bpn-c012_F1' 'brainvision_rda_bpn-c012_F2' 'brainvision_rda_bpn-c012_F6' 'brainvision_rda_bpn-c012_FT9' 'brainvision_rda_bpn-c012_FT7' 'brainvision_rda_bpn-c012_FC3' 'brainvision_rda_bpn-c012_FC4' 'brainvision_rda_bpn-c012_FT8' 'brainvision_rda_bpn-c012_FT10' 'brainvision_rda_bpn-c012_C5' 'brainvision_rda_bpn-c012_C1' 'brainvision_rda_bpn-c012_C2' 'brainvision_rda_bpn-c012_C6' 'brainvision_rda_bpn-c012_TP7' 'brainvision_rda_bpn-c012_CP3' 'brainvision_rda_bpn-c012_CPz' 'brainvision_rda_bpn-c012_CP4' 'brainvision_rda_bpn-c012_TP8' 'brainvision_rda_bpn-c012_P5' 'brainvision_rda_bpn-c012_P1' 'brainvision_rda_bpn-c012_P2' 'brainvision_rda_bpn-c012_P6' 'brainvision_rda_bpn-c012_PO7' 'brainvision_rda_bpn-c012_PO3' 'brainvision_rda_bpn-c012_POz' 'brainvision_rda_bpn-c012_PO4' 'brainvision_rda_bpn-c012_PO8'});

            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname', [fnames(file).name(1:end-4) '_mocap'],'savenew',[output_filepath fnames(file).name(1:end-4) '_mocap.set'],'gui','off');
        end
    end
end

%% STEP D2 Merge all datasets
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
        if ~contains(fnames(i).name, '_EEG.set')
            fnames(i).name = [];
        end
    end
    fnames = {fnames.name};
    fnames = fnames(~cellfun('isempty',fnames));
            
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', fnames, 'filepath', input_filepath);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
    
    % change electrode labels and load default electrode positions
    for i = 1:length(EEG.chanlocs)
        EEG.chanlocs(i).labels = erase(EEG.chanlocs(i).labels, 'brainvision_rda_bpn-c012_');
    end
    EEG=pop_chanedit(EEG, 'lookup','P:\\Lukas_Gehrke\\toolboxes\\eeglab\\plugins\\dipfit2.4\\standard_BESA\\standard-10-5-cap385.elp','rplurchanloc');
        
    % merges all files currently loaded in EEGLAB into one file and stores
    % the original filenames in EEG.etc.appended_files
    [ALLEEG EEG CURRENTSET] = bemobil_merge(ALLEEG,EEG,CURRENTSET,[1:length(ALLEEG)], merged_filename, output_filepath);
    
end

%% STEP E: Resampling, Channel rejection and interpolation of channels
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
    EEG = pop_resample(EEG, 250);    
    
    % Todo: do manual channel rejection
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
        [ALLEEG, EEG, CURRENTSET] = bemobil_interp( EEG , ALLEEG, CURRENTSET, removed_chans, interpolated_filename, output_filepath);
        % -> interpolated bad channels
        % the function stores the interpolated channels in EEG.etc and saves a
        % data set
    end
end

disp('Channel based cleaning and interpolation done!')

%% STEP F: Automatic time domain cleaning on the channel level
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
    
    oriEEG = pop_loadset('filename', interpolated_avRef_filename, 'filepath', input_filepath);
    oriEEG = eeg_checkset( oriEEG );
    
    %%% specify folder names
    datapath_specifications.datapath_original_files=input_filepath; %%% keep last \; single subject folder
    datapath_specifications.datapath_save_files=output_filepath;     %%% keep last \; path for saving updated EEG
    datapath_specifications.datapath_save_figures=output_filepath;   %%% keep last \; path for saving figures of cleaning

    %%% specify file names
    filename_specifications.file_name_original_EEG='interpolated_avRef.set';   %%% loads "fresh" EEG (raw, unfiltered)
    filename_specifications.filename_saveBadEpochIndices='';  

    automatic_cleaning_settings.cleaned_data_type='sensor data'; %%% ICA not implemented yet; usually bad segments found on sensor level are also fine for IC later on

    %%% select channels that should be considered for cleaning
    automatic_cleaning_settings.selected_sensor_channels_for_cleaning=[]; %%% [] use all available channels for cleaning, else specify [1 2 ...]; currently same channels for all subjects alike
    %%% select channel(s) for cleaning before vs. after
    automatic_cleaning_settings.chan_select_plot_before_after=[10];  %%% [] use all available channels for cleaning

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
    [auto_continuous_cleaning]=wrapper_automatic_cleaning_continuous_EEG(datapath_specifications,filename_specifications,automatic_cleaning_settings);
    
    % copy cleaning results and save dataset
    oriEEG.etc.auto_continuous_cleaning = auto_continuous_cleaning;
    pop_saveset( oriEEG, 'filename', amica_filename_input, 'filepath', output_filepath);
    close(gcf);
end

%% STEP G: Rereference loop
% Use average reference in this case, but can be something else, if
% desired!

input_path = [study_folder raw_EEGLAB_data_folder '\'];
output_path = input_path;

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', interpolated_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % Compute average reference
    EEG = pop_reref( EEG, []);
    
    % new data set in EEGLAB
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
    EEG = eeg_checkset( EEG );
    
    % save on disk
    EEG = pop_saveset( EEG, 'filename',interpolated_avRef_filename,'filepath', output_filepath);
    disp('...done');
    
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
end

disp('Rereferencing to average reference done!');

%% STEP H: Filter the data for AMICA
% uses 1Hz highpass

input_path = [study_folder raw_EEGLAB_data_folder];
output_path = [study_folder spatial_filters_folder spatial_filters_folder_AMICA];

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', interpolated_avRef_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    [ALLEEG, EEG, CURRENTSET] = bemobil_filter(ALLEEG, EEG, CURRENTSET, filter_lowCutoffFreqAMICA, filter_highCutoffFreqAMICA, filtered_filename, output_filepath);
     
end

disp('Filtering for AMICA done!');

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
    
    EEG = pop_loadset('filename', amica_filename_input , 'filepath', input_filepath);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
    
    % reject data from previously computed FH channel data rejection
    invalid_segments_index = EEG.etc.auto_continuous_cleaning.invalid_segments_final_start_stop_sample;
    EEG = eeg_eegrej(EEG, invalid_segments_index);
    
    % rank of the data set is reduced by 1, because of average reference,
    % and then reduced by the number of interpolated channels
    rank = EEG.nbchan - 1 - length(EEG.etc.interpolated_channels);
    
    % running signal decomposition with values specified above
    [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, ...
        CURRENTSET, true, num_models, max_threads, rank, [], ...
        amica_filename_output, [output_filepath]);
    
end

%% STEP J: Warping of locations and dipole fitting
% see Mike X. Cohen: Part IV, 24 "Basics of Single Dipole..."
% renames the specified channels, warps the chanlocs on a standard head model and fits dipoles for
% each IC below the threshold of residual variance

input_path = [study_folder spatial_filters_folder spatial_filters_folder_AMICA];
output_path = input_path;

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', amica_filename_output, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % do the warp and dipfit
    bemobil_dipfit( EEG , ALLEEG, CURRENTSET, warping_channel_names, eeglab_path, RV_threshold,...
        remove_outside_head, number_of_dipoles, warped_dipfitted_filename, output_filepath);
    
end

%% STEP K: Create data set for single subject analysis
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
    AMICA_EEG = pop_loadset('filename', warped_dipfitted_filename, 'filepath', spatial_filter_filepath);
    AMICA_EEG = eeg_checkset( AMICA_EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, AMICA_EEG, 0 );
    
    % load preprocessed, interpolated, average referenced file
    EEG = pop_loadset('filename', interpolated_avRef_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % copy AMICA weights and dipfit info
    bemobil_copy_spatial_filter(EEG, ALLEEG, CURRENTSET, AMICA_EEG, true, false, copy_weights_interpolate_avRef_filename , output_filepath);
    
end

%% STEP L: Epoch Loop
% makes epochs of the data set specified. Currently this is changed each time, that's not good. Also
% experiment-specific event fields are entered in a script where you need to specify some stuff in
% again. BE AWARE WHAT'S THE CURRENT STATE OF THE OTHER SCRIPT BEFORE STARTING THIS!

input_path = [study_folder single_subject_analysis_folder];
output_path = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERSPs single_subject_analysis_folder_epochs_1];

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', copy_weights_interpolate_avRef_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % add filter
        
    % fill eventfields by parsing EEG.event.type string
    for file = 1:length(EEG.event)
        if ~strcmp(EEG.event(file).type, 'boundary')
            current_event = cellstr(strsplit(EEG.event(file).type, ';'));
            for j = 1:length(current_event)
                key_val = cellstr(strsplit(current_event{j}, ':'));
                EEG.event(file).(key_val{1}) = key_val{2};
            end
        end
    end

    % change eventfield 'type' to "box:touched" for all touch events
    for file = 1:length(EEG.event)
        if contains(EEG.event(file).type, 'box:spawned')
            EEG.event(file+1).condition = EEG.event(file).condition;
            EEG.event(file+1).normal_or_conflict = EEG.event(file).normal_or_conflict;
            EEG.event(file+1).type = 'box:touched';
            
            % TODO BEWARE the below is only for pilot
            %EEG.event(i).normal_or_conflict = EEG.event(i-1).normal_or_conflict;
%         else
%             EEG.event(i) = [];
        end
    end
    
    [EEG, created_epochs_indices] = pop_epoch( EEG, epochs_1_event, epochs_1_boundaries, 'newname',...
        'epochs', 'epochinfo', 'yes');
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


%% export Reaction for a single loaded EEG set in eeglab GUI

for file = 1:length(EEG.epoch)
    rt(file) = EEG.epoch(file).eventreaction_time(1);
end
rt = str2double(rt)';
xlsxwrite('rt_out', rt);

%% STEP M: epoch cleaning
% Status 04/07/2018: 

input_path = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERSPs single_subject_analysis_folder_epochs_1];
output_path = input_path;

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    
    disp(['Subject #' num2str(subject)]);
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    oriEEG = pop_loadset('filename', epochs_filename, 'filepath', input_filepath);
    oriEEG = eeg_checkset( oriEEG );
    
    % start FH epoch cleaning
    %%% define folder names; folders for saving will be automatically created, if not existing already
    datapath_specifications.datapath_original_files = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERSPs single_subject_analysis_folder_epochs_1 num2str(subject) '\'];
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

    automatic_cleaning_settings.clean_epochs_ICA_use_components=true;  %%% =false for cleaning epochs sensor level ("channels" are EEG channels); =true for cleaning IC epochs ("channels" are components)

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

%% STEP N: Safety checks data quality

input_path = [study_folder single_subject_analysis_folder];
output_path = [study_folder single_subject_analysis_folder];

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    output_filepath = [output_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', copy_weights_interpolate_avRef_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % check all single subject decompositions (scalp map topographies) and save
    % the resulting .fig
    pop_topoplot(EEG,0, [1:size(EEG.icaweights,1)] ,'scalp maps',[8 8] ,0,'electrodes','off');
    savefig(gcf, [output_filepath 'scalp_map' num2str(subject)]);
    close(gcf);
end

%% epoch rejection

% reject bad epochs
trialrej = zeros(1,length(EEG.epoch));
trialrej(EEG.etc.auto_epoch_cleaning.indices_bad_epochs) = 1;
EEG = pop_rejepoch( EEG, trialrej, 0);

%% STEP O: create EEGLAB study structure
% this is group level analyses specific and comes at the end of the data
% collection phase
% Mike X. Cohen: Part VI, Chapter 34, 35

input_path = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERSPs single_subject_analysis_folder_epochs_1];
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
if ~isempty(STUDY_1_components_to_use)
    for set = 1:length(subjects)
        command{end+1} = {'index' set 'comps' STUDY_1_components_to_use };
    end
end

% create study
[STUDY ALLEEG] = std_editset( STUDY, ALLEEG, 'name',study_1_filename,'commands',command,...
    'updatedat','on','savedat','off','rmclust','on' );
[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];

% save study
disp('Saving STUDY...')
mkdir(output_path)
[STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename',study_1_filename,'filepath',output_path);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
eeglab redraw
disp('...done')

%% precompute topographies and spectra

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
    [STUDY ALLEEG] = pop_loadstudy('filename', study_1_filename, 'filepath', input_path);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    
    eeglab redraw
end

% precompute component measures except ERSPs
disp('Precomputing topographies and spectra.')
[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, 'components','recompute','on',...
    'erp', 'on', 'scalp','on','spec','on','specparams',{'specmode' 'fft' 'logtrials' 'off'},...
    'ersp', 'on', 'itc', 'on', 'erspparams',{ 'cycles' [ 3 0.5 ], 'alpha', 0.01, 'padratio' 1 });

% save study
disp('Saving STUDY...')
mkdir(output_path)
[STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename',study_1_filename,'filepath',output_path);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
eeglab redraw
disp('...done')

%% Precluster EEGLAB study

input_path = [study_folder study_level];

% what does this do?
% input_path_latencies = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERSPs single_subject_analysis_folder_epochs_1];
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
    [STUDY ALLEEG] = pop_loadstudy('filename', study_1_filename, 'filepath', input_path);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    
    eeglab redraw
end

% create preclustering array that is used for clustering and save study
[STUDY, ALLEEG, EEG] = bemobil_precluster(STUDY, ALLEEG, EEG, STUDY_1_clustering_weights, STUDY_1_clustering_freqrange,...
    [-1000 2000], study_1_filename, input_path);

%% Repeated clustering EEGLAB study

input_path = [study_folder study_level];
input_path_latencies = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERSPs single_subject_analysis_folder_epochs_1];

full_path_clustering_solutions = [STUDY_1_filepath_clustering_solutions num2str(STUDY_1_n_clust) '-cluster_' num2str(outlier_sigma)...
    '-sigma_' num2str(STUDY_1_clustering_weights.dipoles) '-dipoles_' num2str(STUDY_1_clustering_weights.spectra) '-spec_'...
    num2str(STUDY_1_clustering_weights.scalp_topographies) '-scalp_' num2str(STUDY_1_clustering_weights.ERSPs) '-ersp_' num2str(n_iterations) '-iterations'];

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
    [STUDY ALLEEG] = pop_loadstudy('filename', study_1_filename, 'filepath', input_path);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    
    eeglab redraw
end

% cluster the components repeatedly and use a region of interest and
% quality measures to find the best fitting solution
[STUDY, ALLEEG, EEG] = bemobil_repeated_clustering_and_evaluation(STUDY, ALLEEG, EEG, outlier_sigma,...
    STUDY_1_n_clust, n_iterations, STUDY_1_cluster_ROI_talairach, STUDY_1_quality_measure_weights, false,...
    false, input_path, study_1_filename, [input_path full_path_clustering_solutions],...
    filename_clustering_solutions, [input_path full_path_clustering_solutions], filename_multivariate_data);
