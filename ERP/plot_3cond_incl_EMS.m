% set Study parameters for ERPs of event: "box:touched"
study_params_PE_incl_EMS;

input_path = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERPs single_subject_analysis_folder_epochs];
output_path = input_path;

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

for subject = subjects
    
    disp(['Subject #' num2str(subject)]);
    input_filepath = [input_path '\' num2str(subject)];
    output_filepath = [output_path '\' num2str(subject)];
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', epochs_filename, 'filepath', input_filepath);
    EEG = eeg_checkset( EEG );

    % select epochs
    bad = EEG.etc.auto_epoch_cleaning.indices_bad_epochs;
    ind_ems = find(cellfun(@(x) any(strcmp(x, 'ems')), {EEG.epoch.eventcondition}));
    ind_ems_rej = setdiff(ind_ems, bad);
    ind_vis = find(cellfun(@(x) any(strcmp(x, 'visual')), {EEG.epoch.eventcondition}));
    ind_vis_rej = setdiff(ind_vis, bad);
    ind_vibro = find(cellfun(@(x) any(strcmp(x, 'vibro')), {EEG.epoch.eventcondition}));
    ind_vibro_rej = setdiff(ind_vibro, bad);

    % find smallest array
    good_epochs = min([length(ind_ems_rej) length(ind_vis_rej) length(ind_vibro_rej)]);

    % get N = good_epochs random samples from the good set of samples
    rand_ind_ems = sort(randsample(ind_ems_rej, good_epochs));
    rand_ind_vis = sort(randsample(ind_vis_rej, good_epochs));
    rand_ind_vibro = sort(randsample(ind_vibro_rej, good_epochs));
    
    % extract each condition data with good epochs from complete dataset
    EEG_ems = pop_selectevent(EEG, 'epoch', rand_ind_ems);
    EEG_vis = pop_selectevent(EEG, 'epoch', rand_ind_vis);
    EEG_vibro = pop_selectevent(EEG, 'epoch', rand_ind_vibro);
    
    % extract condition data from each condition dataset with good epochs
    EEG_ems_normal = pop_selectevent(EEG_ems, 'normal_or_conflict', 'normal');
    EEG_ems_conflict = pop_selectevent(EEG_ems, 'normal_or_conflict', 'conflict');
    EEG_vis_normal = pop_selectevent(EEG_vis, 'normal_or_conflict', 'normal');
    EEG_vis_conflict = pop_selectevent(EEG_vis, 'normal_or_conflict', 'conflict');
    EEG_vibro_normal = pop_selectevent(EEG_vibro, 'normal_or_conflict', 'normal');
    EEG_vibro_conflict = pop_selectevent(EEG_vibro, 'normal_or_conflict', 'conflict');
    
    % make plots single subject
    plot_erp({{EEG_vis_normal}, {EEG_vibro_normal}, {EEG_ems_normal}},...
        'Cz', 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'});
    savefig(gcf, [output_filepath '\normal_3cond']);
    
    plot_erp({{EEG_vis_conflict}, {EEG_vibro_conflict}, {EEG_ems_conflict}},...
        'Cz', 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'});
    savefig(gcf, [output_filepath '\conflict_3cond']);
end


