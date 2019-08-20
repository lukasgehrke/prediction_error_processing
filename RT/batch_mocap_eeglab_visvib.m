%% set Study parameters for ERPs of event: "box:touched"
params_mocap;

%% merge mocap files

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
        if ~contains(fnames(i).name, '_mocap.set') || contains(fnames(i).name, 'Training') || contains(fnames(i).name, 'merged')
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
    [ALLEEG EEG CURRENTSET] = bemobil_merge(ALLEEG,EEG,CURRENTSET,[1:length(ALLEEG)], merged_mocap_filename, output_filepath);
    
end

%% extract RTs

addpath 'P:\Lukas_Gehrke\toolboxes\utils_LG'

input_path = [study_folder raw_EEGLAB_data_folder];
output_path = 'P:\Project_Sezen\data_processing\RT';

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', merged_mocap_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
    % parse events
    EEG = parse_events_PredError(EEG);

    touches = find(strcmp('touched', {EEG.event.box}));
    
    RT{1,1} = 'subject';
    RT{2,1} = 'condition';
    RT{3,1} = 'congruency';
    RT{4,1} = 'RT';
        
    for i = 1:length(touches)
        try
            for j = 1:10
                if ~isempty(EEG.event(touches(i)+j).hand_movement)
                    
                    switch EEG.event(touches(i)+j).hand_movement
                        case 'end'
                            RT{4,i+1} = num2str((EEG.event(touches(i)+j).latency - EEG.event(touches(i)).latency) / EEG.srate);
                        case 'start'
                            RT{4,i+1} = '0';    
                    end
                    RT{2,i+1} = EEG.event(touches(i)).condition;
                    RT{3,i+1} = EEG.event(touches(i)).normal_or_conflict;
                    
                    break;
                end
            end
        end
    end
    RT(1,2:end) = {num2str(subject)};
    
    if subject==2
        RT_all = RT;
    else
        RT_all = [RT_all, RT(:,2:end)];
    end
    
    % save matrix to csv format per subject
    cell2csv([output_path '\RT_results_' num2str(subject) '.csv'], RT', ';');
    clear RT
end


% save matrix to csv format for all subject
cell2csv([output_path '\RT_results_ALL.csv'], RT_all', ';');
