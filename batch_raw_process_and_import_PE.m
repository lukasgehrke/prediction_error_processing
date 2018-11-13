%% Specific batch script for Prediction Error Experiment
% using EEGLAB & MoBILAB & Marius Klug's MoBILAB extensions (https://github.com/MariusKlug/mobilab)
% and BeMoBIL pipeline functions and experiment specific scripts &
% functions

% Versions:
% MATLAB Version: 9.2.0.538062 (R2017a)
% MATLAB License Number: FreeForAll
% Operating System: Microsoft Windows 10 Pro Version 10.0 (Build 17134)

% EEGLAB
addpath 'M:\Toolboxes_Skripts_and_Coding_examples\eeglab-by-marius\eeglab14_1_0b'
eeglab

% MoBILAB

% START: set Study parameters
addpath 'P:\Project_Sezen\data_processing'
study_params_PE_vis_vibro;

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
            output_filepath = [output_path num2str(subject) '\' fnames(file).name(1:end-4) '_MoBI'];

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
            
            for i = 1:100
                try
                    if subject==7
                        mobilab.allStreams.deleteItem(6);
                    else
                        mobilab.allStreams.deleteItem(5);
                    end                    
                catch
                    break
                end
            end

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
            
            eulerHand.children{1}.createEventsFromMagnitude(3,...
                'movements','hand_movement:start hand_movement:end',2,0,65,5,200);
                       
%             % ugly, throwout 2 channels repeatedly for 4 channel vel data
%             % because
%             % mobilab.allStreams.item{hand+9}.throwOutChannels(2:3); not
%             % working
%             % TODO redo throwOutChannels with all channels but 1 which is
%             % then overwritten with 3D magnitude
%             mobilab.allStreams.item{hand+9}.throwOutChannels(2);
%             mobilab.allStreams.item{hand+9}.children{1}.throwOutChannels(2);
%             
%             % create event markers based on 3D hand velocity thresholding
%             % get velocity data of hand RB
%             vel = mobilab.allStreams.item{hand+9};
%             vel_xyz_newChannel = sqrt(vel.data(:,1).^2 + vel.data(:,2).^2 + vel.data(:,3).^2);
%             
%             % add channel to the matching length RB stream
%             % not working: vel.addChannels(1, vel_xyz_newChannel);
%             mobilab.allStreams.item{hand+9}.children{1}.children{1}.data(:,1) = vel_xyz_newChannel;
%             
%             % explore velocity data for later parameter selection in
%             % onset/offset detection
% %             figure;histogram(vel_xyz_newChannel);
% %             figure;plot(sort(vel_xyz_newChannel));
% %             prctile(vel_xyz_newChannel, [90:1:100]);
% %             figure;histogram(mobilab.allStreams.item{hand+9}.data(:,1), 20);
%             
%             mobilab.allStreams.item{hand+9}.children{1}.children{1}.createEventsFromMagnitude(1,...
%                 'movements','hand_movement:start -',0.7,0,90,3,200);
%             
%             % for head; just exploration, no hypotheses/interest as of yet
% %             head_euler = mobilab.allStreams.item{head+11}.data(:,4);
% %             figure;histogram(abs(head_euler));
% %             figure;plot(sort(abs(head_euler)));
% %             prctile(abs(head_euler), [1:1:100])
% %             figure;histogram(mobilab.allStreams.item{head+9}.data(:,4), 20);
            
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

%% STEP D: Separating EEG and Mocap loop & correct EEG age_of_sample
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
            EEG = pop_select( EEG,'channel',{'brainvision_rda_bpn-c012_Fp1' 'brainvision_rda_bpn-c012_Fp2' 'brainvision_rda_bpn-c012_F7' 'brainvision_rda_bpn-c012_F3' 'brainvision_rda_bpn-c012_Fz' 'brainvision_rda_bpn-c012_F4' 'brainvision_rda_bpn-c012_F8' 'brainvision_rda_bpn-c012_FC5' 'brainvision_rda_bpn-c012_FC1' 'brainvision_rda_bpn-c012_FC2' 'brainvision_rda_bpn-c012_FC6' 'brainvision_rda_bpn-c012_T7' 'brainvision_rda_bpn-c012_C3' 'brainvision_rda_bpn-c012_Cz' 'brainvision_rda_bpn-c012_C4' 'brainvision_rda_bpn-c012_T8' 'brainvision_rda_bpn-c012_TP9' 'brainvision_rda_bpn-c012_CP5' 'brainvision_rda_bpn-c012_CP1' 'brainvision_rda_bpn-c012_CP2' 'brainvision_rda_bpn-c012_CP6' 'brainvision_rda_bpn-c012_TP10' 'brainvision_rda_bpn-c012_P7' 'brainvision_rda_bpn-c012_P3' 'brainvision_rda_bpn-c012_Pz' 'brainvision_rda_bpn-c012_P4' 'brainvision_rda_bpn-c012_P8' 'brainvision_rda_bpn-c012_PO9' 'brainvision_rda_bpn-c012_O1' 'brainvision_rda_bpn-c012_Oz' 'brainvision_rda_bpn-c012_O2' 'brainvision_rda_bpn-c012_PO10' 'brainvision_rda_bpn-c012_AF7' 'brainvision_rda_bpn-c012_AF3' 'brainvision_rda_bpn-c012_AF4' 'brainvision_rda_bpn-c012_AF8' 'brainvision_rda_bpn-c012_F5' 'brainvision_rda_bpn-c012_F1' 'brainvision_rda_bpn-c012_F2' 'brainvision_rda_bpn-c012_F6' 'brainvision_rda_bpn-c012_FT9' 'brainvision_rda_bpn-c012_FT7' 'brainvision_rda_bpn-c012_FC3' 'brainvision_rda_bpn-c012_FC4' 'brainvision_rda_bpn-c012_FT8' 'brainvision_rda_bpn-c012_FT10' 'brainvision_rda_bpn-c012_C5' 'brainvision_rda_bpn-c012_C1' 'brainvision_rda_bpn-c012_C2' 'brainvision_rda_bpn-c012_C6' 'brainvision_rda_bpn-c012_TP7' 'brainvision_rda_bpn-c012_CP3' 'brainvision_rda_bpn-c012_CPz' 'brainvision_rda_bpn-c012_CP4' 'brainvision_rda_bpn-c012_TP8' 'brainvision_rda_bpn-c012_P5' 'brainvision_rda_bpn-c012_P1' 'brainvision_rda_bpn-c012_P2' 'brainvision_rda_bpn-c012_P6' 'brainvision_rda_bpn-c012_PO7' 'brainvision_rda_bpn-c012_PO3' 'brainvision_rda_bpn-c012_POz' 'brainvision_rda_bpn-c012_PO4' 'brainvision_rda_bpn-c012_PO8'});
            
            % apply alignment shift compensating for age_of_sample of EEG data
            shift = round(data_shift * EEG.srate);
            EEG.data = [EEG.data zeros(size(EEG.data, 1), shift)];
            EEG.data = EEG.data(:,shift+1:end);
            
            % reject shift buffers and obtain original data set length
            EEG = eeg_eegrej(EEG, [length(EEG.data)-shift length(EEG.data)]);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[fnames(file).name(1:end-4) '_EEG'],'savenew',[output_filepath fnames(file).name(1:end-4) '_EEG.set'],'gui','off');

            % go back to the first data set
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',1,'study',0);
            EEG = eeg_checkset( EEG );

            % split it and kick all EEG data channels out, save it with _mocap suffix
            EEG = pop_select( EEG,'nochannel',{'brainvision_rda_bpn-c012_Fp1' 'brainvision_rda_bpn-c012_Fp2' 'brainvision_rda_bpn-c012_F7' 'brainvision_rda_bpn-c012_F3' 'brainvision_rda_bpn-c012_Fz' 'brainvision_rda_bpn-c012_F4' 'brainvision_rda_bpn-c012_F8' 'brainvision_rda_bpn-c012_FC5' 'brainvision_rda_bpn-c012_FC1' 'brainvision_rda_bpn-c012_FC2' 'brainvision_rda_bpn-c012_FC6' 'brainvision_rda_bpn-c012_T7' 'brainvision_rda_bpn-c012_C3' 'brainvision_rda_bpn-c012_Cz' 'brainvision_rda_bpn-c012_C4' 'brainvision_rda_bpn-c012_T8' 'brainvision_rda_bpn-c012_TP9' 'brainvision_rda_bpn-c012_CP5' 'brainvision_rda_bpn-c012_CP1' 'brainvision_rda_bpn-c012_CP2' 'brainvision_rda_bpn-c012_CP6' 'brainvision_rda_bpn-c012_TP10' 'brainvision_rda_bpn-c012_P7' 'brainvision_rda_bpn-c012_P3' 'brainvision_rda_bpn-c012_Pz' 'brainvision_rda_bpn-c012_P4' 'brainvision_rda_bpn-c012_P8' 'brainvision_rda_bpn-c012_PO9' 'brainvision_rda_bpn-c012_O1' 'brainvision_rda_bpn-c012_Oz' 'brainvision_rda_bpn-c012_O2' 'brainvision_rda_bpn-c012_PO10' 'brainvision_rda_bpn-c012_AF7' 'brainvision_rda_bpn-c012_AF3' 'brainvision_rda_bpn-c012_AF4' 'brainvision_rda_bpn-c012_AF8' 'brainvision_rda_bpn-c012_F5' 'brainvision_rda_bpn-c012_F1' 'brainvision_rda_bpn-c012_F2' 'brainvision_rda_bpn-c012_F6' 'brainvision_rda_bpn-c012_FT9' 'brainvision_rda_bpn-c012_FT7' 'brainvision_rda_bpn-c012_FC3' 'brainvision_rda_bpn-c012_FC4' 'brainvision_rda_bpn-c012_FT8' 'brainvision_rda_bpn-c012_FT10' 'brainvision_rda_bpn-c012_C5' 'brainvision_rda_bpn-c012_C1' 'brainvision_rda_bpn-c012_C2' 'brainvision_rda_bpn-c012_C6' 'brainvision_rda_bpn-c012_TP7' 'brainvision_rda_bpn-c012_CP3' 'brainvision_rda_bpn-c012_CPz' 'brainvision_rda_bpn-c012_CP4' 'brainvision_rda_bpn-c012_TP8' 'brainvision_rda_bpn-c012_P5' 'brainvision_rda_bpn-c012_P1' 'brainvision_rda_bpn-c012_P2' 'brainvision_rda_bpn-c012_P6' 'brainvision_rda_bpn-c012_PO7' 'brainvision_rda_bpn-c012_PO3' 'brainvision_rda_bpn-c012_POz' 'brainvision_rda_bpn-c012_PO4' 'brainvision_rda_bpn-c012_PO8'});

            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname', [fnames(file).name(1:end-4) '_mocap'],'savenew',[output_filepath fnames(file).name(1:end-4) '_mocap.set'],'gui','off');
        end
    end
end

%% STEP E: Merge all datasets condition Visual and VisualVibro only! 
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
        if ~contains(fnames(i).name, '_EEG.set') || contains(fnames(i).name, 'Training') || contains(fnames(i).name, 'EMS')
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
        
    % Todo: do manual channel rejection
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
            % also has two frontal electrodes misplaced that need to be switched
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
            removed_chans = [7 16 40 46 48];
        case 14
            removed_chans = [2 3 7 16 28];
        case 15
            removed_chans = [5 6 12 33 34 46];
        case 16
            removed_chans = [28 29 41 45 60];
        case 17
            removed_chans = [1 2 3 22 28 36];
        case 18
            removed_chans = [15 17 26 30 45];
        case 19
            removed_chans = [15 22 26 46 55 59 60];        
        case 20
            removed_chans = [2 8 11 36 62];
    end
    
    removed_chans = sort(removed_chans);
    
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
    EEG = fullRankAveRef( EEG );
    
    % save dataset
    mkdir(output_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(output_filepath);
    if ismember(interpolated_filename, {dir_files.name})
        warning([interpolated_filename ' file already exists in: ' output_filepath '. File will be overwritten...']);
    end
    
    pop_saveset( EEG, 'filename', interpolated_filename, 'filepath', output_filepath);
end

disp('Channel based cleaning and interpolation done!')
