%% set Study parameters for ERPs of event: "box:touched"
study_params_ERP_box_touched_incl_EMS;

%% load EEG sets
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
    
    % extract vis and vis+vibro condition data, all trials
    all_trials_vis_vibro{subject}{1} = pop_selectevent(EEG, 'epoch', ind_vis);
    all_trials_vis_vibro{subject}{2} = pop_selectevent(EEG, 'epoch', ind_vibro);
    
    all_trials_vis_vibro_normal{subject}{1} = pop_selectevent(all_trials_vis_vibro{subject}{1}, 'normal_or_conflict', 'normal');
    all_trials_vis_vibro_conflict{subject}{1} = pop_selectevent(all_trials_vis_vibro{subject}{1}, 'normal_or_conflict', 'conflict');
    all_trials_vis_vibro_normal{subject}{2} = pop_selectevent(all_trials_vis_vibro{subject}{2}, 'normal_or_conflict', 'normal');
    all_trials_vis_vibro_conflict{subject}{2} = pop_selectevent(all_trials_vis_vibro{subject}{2}, 'normal_or_conflict', 'conflict');
    
    % extract each condition data with good epochs from complete dataset
    % index 1 = ems, 2 = vis, 3 = vibro
    reduced_all_conds{subject}{1} = pop_selectevent(EEG, 'epoch', rand_ind_ems);
    reduced_all_conds{subject}{2} = pop_selectevent(EEG, 'epoch', rand_ind_vis);
    reduced_all_conds{subject}{3} = pop_selectevent(EEG, 'epoch', rand_ind_vibro);
    
    % extract condition data from each condition dataset with good epochs
    normal{subject}{1} = pop_selectevent(reduced_all_conds{subject}{1}, 'normal_or_conflict', 'normal');
    conflict{subject}{1} = pop_selectevent(reduced_all_conds{subject}{1}, 'normal_or_conflict', 'conflict');
    normal{subject}{2} = pop_selectevent(reduced_all_conds{subject}{2}, 'normal_or_conflict', 'normal');
    conflict{subject}{2} = pop_selectevent(reduced_all_conds{subject}{2}, 'normal_or_conflict', 'conflict');
    normal{subject}{3} = pop_selectevent(reduced_all_conds{subject}{3}, 'normal_or_conflict', 'normal');
    conflict{subject}{3} = pop_selectevent(reduced_all_conds{subject}{3}, 'normal_or_conflict', 'conflict');
      
end

%% plotting

chans = {'Fz' 'Cz' 'Pz' 'Oz'};

% 'standard' ERP plotting
for i = 1:2%length(chans)
    % Difference ERP vis+vibro+ems matched number normal trials
    plot_erp({{normal{2}{1}, normal{3}{1}, normal{6}{1}, normal{7}{1}, normal{8}{1}},...
        {conflict{2}{1}, conflict{3}{1}, conflict{6}{1}, conflict{7}{1},conflict{8}{1}}},...
        chans{i}, 'labels', {'Visual+Vibro+EMS Normal', 'Visual+Vibro+EMS Conflict'}, ...
        'plotstd', 'fill', 'plotdiff', 1);
    
    % Vis, vis+vibro, vis+vibro+ems matched number normal trials
    plot_erp({{normal{2}{2}, normal{3}{2}, normal{6}{2}, normal{7}{2}, normal{8}{2}},...
        {normal{2}{3}, normal{3}{3}, normal{6}{3}, normal{7}{3}, normal{8}{3}},...
        {normal{2}{1}, normal{3}{1}, normal{6}{1}, normal{7}{1}, normal{8}{1}}},...
        chans{i}, 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
        'plotstd', 'fill');

    % Vis, vis+vibro, vis+vibro+ems matched number conflict trials
    plot_erp({{conflict{2}{2}, conflict{3}{2}, conflict{6}{2}, conflict{7}{2}, conflict{8}{2}},...
        {conflict{2}{3}, conflict{3}{3}, conflict{6}{3}, conflict{7}{3}, conflict{8}{3}},...
        {conflict{2}{1}, conflict{3}{1}, conflict{6}{1}, conflict{7}{1},conflict{8}{1}}},...
        chans{i}, 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
        'plotstd', 'fill');
    
    % plot 2 main condition visual and visual+vibro with both conflict and
    % normal trials
    plot_erp({{all_trials_vis_vibro_normal{2}{1}, all_trials_vis_vibro_normal{3}{1}, all_trials_vis_vibro_normal{6}{1}, all_trials_vis_vibro_normal{7}{1}, all_trials_vis_vibro_normal{8}{1}},...
        {all_trials_vis_vibro_normal{2}{2}, all_trials_vis_vibro_normal{3}{2}, all_trials_vis_vibro_normal{6}{2}, all_trials_vis_vibro_normal{7}{2}, all_trials_vis_vibro_normal{8}{2}},...
        {all_trials_vis_vibro_conflict{2}{1}, all_trials_vis_vibro_conflict{3}{1}, all_trials_vis_vibro_conflict{6}{1}, all_trials_vis_vibro_conflict{7}{1}, all_trials_vis_vibro_conflict{8}{1}},...
        {all_trials_vis_vibro_conflict{2}{2}, all_trials_vis_vibro_conflict{3}{2}, all_trials_vis_vibro_conflict{6}{2}, all_trials_vis_vibro_conflict{7}{2}, all_trials_vis_vibro_conflict{8}{2}}},...
        chans{i}, 'labels', {'Visual Normal', 'Visual+Vibro Normal', 'Visual Conflict', 'Visual+Vibro Conflict'}, ...
        'plotstd', 'fill');
    
    % contrast 2 main condition visual and visual+vibro
    plot_erp({{all_trials_vis_vibro_conflict{2}{1}, all_trials_vis_vibro_conflict{3}{1}, all_trials_vis_vibro_conflict{6}{1}, all_trials_vis_vibro_conflict{7}{1}, all_trials_vis_vibro_conflict{8}{1}},...
        {all_trials_vis_vibro_conflict{2}{2}, all_trials_vis_vibro_conflict{3}{2}, all_trials_vis_vibro_conflict{6}{2}, all_trials_vis_vibro_conflict{7}{2}, all_trials_vis_vibro_conflict{8}{2}}},...
        chans{i}, 'labels', {'Visual Conflict', 'Visual+Vibro Conflict'}, ...
        'plotstd', 'fill', 'permute', 100);
end

% first calculate difference ERP of individually baseline corrected
% condition ERPs, then plot across conditions

for subject = subjects
    conflict{subject}{1}.data(14,:,:) = conflict{subject}{1}.data(find(strcmp({EEG.chanlocs.labels}, chans{i})),:,:)...
        - mean(normal{subject}{1}.data(find(strcmp({EEG.chanlocs.labels}, chans{i})),:,:),3);
    conflict{subject}{3}.data(14,:,:) = conflict{subject}{3}.data(find(strcmp({EEG.chanlocs.labels}, chans{i})),:,:)...
        - mean(normal{subject}{3}.data(find(strcmp({EEG.chanlocs.labels}, chans{i})),:,:),3);
    conflict{subject}{2}.data(14,:,:) = conflict{subject}{2}.data(find(strcmp({EEG.chanlocs.labels}, chans{i})),:,:)...
        - mean(normal{subject}{2}.data(find(strcmp({EEG.chanlocs.labels}, chans{i})),:,:),3);
end
% Vis, vis+vibro, vis+vibro+ems first difference within feedback condition 
plot_erp({{conflict{2}{2}, conflict{3}{2}, conflict{6}{2}, conflict{7}{2}, conflict{8}{2}, conflict{10}{2}, conflict{11}{2}},...
    {conflict{2}{3}, conflict{3}{3}, conflict{6}{3}, conflict{7}{3}, conflict{8}{3}, conflict{10}{3}, conflict{11}{3}},...
    {conflict{2}{1}, conflict{3}{1}, conflict{6}{1}, conflict{7}{1}, conflict{8}{1}, conflict{10}{1}, conflict{11}{1}}},...
    'Cz', 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
    'plotstd', 'fill');


% Vis, vis+vibro, vis+vibro+ems first difference within feedback condition 
plot_erp({{normal{2}{2}, normal{3}{2}, normal{6}{2}, normal{7}{2}, normal{8}{2}, normal{10}{2}, normal{11}{2}},...
    {normal{2}{3}, normal{3}{3}, normal{6}{3}, normal{7}{3}, normal{8}{3}, normal{10}{3}, normal{11}{3}},...
    {normal{2}{1}, normal{3}{1}, normal{6}{1}, normal{7}{1}, normal{8}{1}, normal{10}{1}, normal{11}{1}}},...
    'Cz', 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
    'plotstd', 'fill');

% TODO make function that computes difference channel between 
% add data channel containing difference ERP
% first, add data channel to not lose original data
EEG.data(end+1,:) = 0;
EEG.nbchan = size(EEG.data,1);
if ~isempty(EEG.chanlocs)
    EEG.chanlocs(end+1).label = '';
end;

% make plot single subject
single_fig = 1;
chans = {'Cz' 'Fz' 'Pz' 'Oz'};
for i = 1:length(chans)
    if ~single_fig
        plot_erp({{EEG_vis}, {EEG_vibro}, {EEG_ems}},...
            chans{i}, 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
            'plotstd', 'fill');
        savefig(gcf, [output_filepath '\all_3cond_' chans{i}]);
    else
        f = figure;
        p = uipanel('Parent',f,'BorderType','none'); 
        p.Title = 'ERPs of box:touched events'; 
        p.TitlePosition = 'centertop'; 
        p.FontSize = 12;
        p.FontWeight = 'bold';
        subplot(3,4,i,'Parent',p);
        plot_erp({{EEG_vis}, {EEG_vibro}, {EEG_ems}},...
            chans{i}, 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
            'plotstd', 'fill', 'newfig', 0);
    end

    % normal trials
    if ~single_fig
        plot_erp({{EEG_vis_normal}, {EEG_vibro_normal}, {EEG_ems_normal}},...
            chans{i}, 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
            'plotstd', 'fill');
        savefig(gcf, [output_filepath '\normal_3cond_' chans{i}]);
    else
        subplot(3,4,i+4,'Parent',p);
        plot_erp({{EEG_vis_normal}, {EEG_vibro_normal}, {EEG_ems_normal}},...
            chans{i}, 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
            'plotstd', 'fill', 'newfig', 0);
    end

    % conflict trials
    if ~single_fig
        plot_erp({{EEG_vis_conflict}, {EEG_vibro_conflict}, {EEG_ems_conflict}},...
            chans{i}, 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
            'plotstd', 'fill');
        savefig(gcf, [output_filepath '\conflict_3cond_' chans{i}]);
    else
        subplot(3,4,i+8,'Parent',p);
        plot_erp({{EEG_vis_normal}, {EEG_vibro_normal}, {EEG_ems_normal}},...
            chans{i}, 'labels', {'Visual Only', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
            'plotstd', 'fill', 'newfig', 0);
    end
end


