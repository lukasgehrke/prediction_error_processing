%% set Study parameters
study_params_PE_incl_EMS;

%% extract IPQ scores and save as .csv
input_path = [study_folder raw_EEGLAB_data_folder];
output_filepath = 'P:\Project_Sezen\data_processing\IPQ';

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    disp(['Subject #' num2str(subject)]);
    
    input_filepath = [input_path num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', interpolated_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );
    EEG = parse_events_PredError(EEG);    
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % get indices of questions and conditions
    ipq1_ix = find(~cellfun(@isempty,{EEG.event.ipq_question_nr_1_answer}));
    ipq2_ix = find(~cellfun(@isempty,{EEG.event.ipq_question_nr_2_answer}));
    ipq3_ix = find(~cellfun(@isempty,{EEG.event.ipq_question_nr_3_answer}));
    ipq4_ix = find(~cellfun(@isempty,{EEG.event.ipq_question_nr_4_answer}));
    conds_ix = find(~cellfun(@isempty,{EEG.event.condition}));
    
    % extract values to a new matrix
    for i = 1:length(ipq1_ix)
        [m, ix] = min(abs(conds_ix-ipq1_ix(i)));
        answers{1,1} = 'subject';
        answers{i+1,1} = EEG.event(conds_ix(ix)).condition;
        answers{i+1,2} = EEG.event(ipq1_ix(i)).ipq_question_nr_1_answer;
        answers{i+1,3} = EEG.event(ipq2_ix(i)).ipq_question_nr_2_answer;
        answers{i+1,4} = EEG.event(ipq3_ix(i)).ipq_question_nr_3_answer;
        answers{i+1,5} = EEG.event(ipq4_ix(i)).ipq_question_nr_4_answer;
        
    end
    answers(1,2:end) = {num2str(subject)};
    
    % save matrix to csv format per subject
    cell2csv([output_filepath '\ipq_results_' num2str(subject) '.csv'], answers', ';');
    clear answers
end




