%% set Study parameters for ERPs of event: "box:touched"
addpath 'P:\Project_Sezen\data_processing\ERP';
study_params_ERP_box_touched_incl_EMS;

% baseline
baseline = [-200 0];
data_shift = 0.063;

% filtering for peak detection
low = 0.2;
high = 10;

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
    
    % remove baseline
    EEG = pop_rmbase(EEG, baseline);
    
    %EEG = pop_eegfiltnew(EEG, low, high);

    % select epochs
    bad = EEG.etc.auto_epoch_cleaning.indices_bad_epochs;
    ind_ems = find(cellfun(@(x) any(strcmp(x, 'ems')), {EEG.epoch.eventcondition}));
    ind_ems_rej = setdiff(ind_ems, bad);
    ind_vis = find(cellfun(@(x) any(strcmp(x, 'visual')), {EEG.epoch.eventcondition}));
    ind_vis_rej = setdiff(ind_vis, bad);
    ind_vibro = find(cellfun(@(x) any(strcmp(x, 'vibro')), {EEG.epoch.eventcondition}));
    ind_vibro_rej = setdiff(ind_vibro, bad);
 
    % extract each condition data, ALL epochs from complete dataset
    % index 1 = ems, 2 = vis, 3 = vibro
    ems = pop_selectevent(EEG, 'epoch', ind_ems_rej);
    vis = pop_selectevent(EEG, 'epoch', ind_vis_rej);
    vibro = pop_selectevent(EEG, 'epoch', ind_vibro_rej);
    
    % extract condition data from each condition dataset with good epochs
    vis_norm = pop_selectevent(vis, 'normal_or_conflict', 'normal');
    vis_conf = pop_selectevent(vis, 'normal_or_conflict', 'conflict');
    vibro_norm = pop_selectevent(vibro, 'normal_or_conflict', 'normal');
    vibro_conf = pop_selectevent(vibro, 'normal_or_conflict', 'conflict');

%     % load 'good' data sets
%     vibro_norm = load([input_filepath '\CHI_vibro_norm_trials.mat']);
%     vibro_norm = vibro_norm.vibro_norm;
%     vibro_conf = load([input_filepath '\CHI_vibro_conf_trials.mat']);
%     vibro_conf = vibro_conf.vibro_conf;
%     vibro_norm = load([input_filepath '\CHI_vibro_norm_trials.mat']);
%     vibro_norm = vibro_norm.vibro_norm;
%     vibro_norm = load([input_filepath '\CHI_vibro_norm_trials.mat']);
%     vibro_norm = vibro_norm.vibro_norm;
    
    ems_norm = pop_selectevent(ems, 'normal_or_conflict', 'normal');
    ems_conf = pop_selectevent(ems, 'normal_or_conflict', 'conflict');
    
    % match trial number of vis and vibro condition to EMS
    conf_N = length(ems_conf.epoch);
    norm_N = length(ems_norm.epoch);
    
    % randomly select N conf and norm trials from Vis and Vibro conditions
    %rand_norm = sort(randsample(1:length(vis_norm.epoch), norm_N));
    rand_norm = load([output_filepath '\CHI_rand_norm_trials.mat']);
    rand_norm = rand_norm.rand_norm;
    %rand_conf = sort(randsample(1:length(vis_conf.epoch), conf_N));
    rand_conf = load([output_filepath '\CHI_rand_conf']);
    rand_conf = rand_conf.rand_conf;
   
    vis_norm_N = pop_selectevent(vis_norm, 'epoch', rand_norm);
    vis_conf_N = pop_selectevent(vis_conf, 'epoch', rand_conf);
    vibro_norm_N = pop_selectevent(vibro_norm, 'epoch', rand_norm);
    vibro_conf_N = pop_selectevent(vibro_conf, 'epoch', rand_conf);
    
    % save trial inds for CHI study
    %save([output_filepath '\CHI_rand_norm_trials.mat'], 'rand_norm');
%     save([output_filepath '\CHI_rand_conf.mat'], 'rand_conf');
    
    % extract condition data from each condition dataset with good epochs
    normal{subject}{1} = vis_norm_N;
    conflict{subject}{1} = vis_conf_N;
    normal{subject}{2} = vibro_norm_N;
    conflict{subject}{2} = vibro_conf_N;
    normal{subject}{3} = ems_norm;
    conflict{subject}{3} = ems_conf;
    
    clear EEG;
end
disp('All files loaded!');

%% plotting of single condition ERPs

input_path = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERPs single_subject_analysis_folder_epochs];
chans = {'Fp1' 'FCz' 'Fz' 'Cz' 'Pz' 'Oz'};

for j = 1:length(chans)
    plot_erp_LG({{conflict{2}{1}, conflict{3}{1}, conflict{6}{1}, conflict{7}{1},...
        conflict{8}{1}, conflict{11}{1}, conflict{12}{1}, conflict{13}{1},...
        conflict{14}{1}, conflict{15}{1}, conflict{16}{1}},...
        
        {conflict{2}{2}, conflict{3}{2}, conflict{6}{2}, conflict{7}{2},...
        conflict{8}{2}, conflict{11}{2}, conflict{12}{2}, conflict{13}{2},...
        conflict{14}{2}, conflict{15}{2}, conflict{16}{2}},...
        
        {conflict{2}{3}, conflict{3}{3}, conflict{6}{3}, conflict{7}{3},...
        conflict{8}{3}, conflict{11}{3}, conflict{12}{3}, conflict{13}{3},...
        conflict{14}{3}, conflict{15}{3}, conflict{16}{3}}},...
        chans{j}, 'labels', {'Visual', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
        'plotstd', 'fill', 'yticksize', 4);
    
    savefig(gcf, [input_path '\11subjects_incl_EMS_conflict_' num2str(baseline(1)) '-' num2str(baseline(2)) '_' chans{j}]);
    close(gcf);

    plot_erp_LG({{normal{2}{1}, normal{3}{1}, normal{6}{1}, normal{7}{1},...
        normal{8}{1}, normal{11}{1}, normal{12}{1}, normal{13}{1},...
        normal{14}{1}, normal{15}{1}, normal{16}{1}},...
        
        {normal{2}{2}, normal{3}{2}, normal{6}{2}, normal{7}{2},...
        normal{8}{2}, normal{11}{2}, normal{12}{2}, normal{13}{2},...
        normal{14}{2}, normal{15}{2}, normal{16}{2}},...
        
        {normal{2}{3}, normal{3}{3}, normal{6}{3}, normal{7}{3},...
        normal{8}{3}, normal{11}{3}, normal{12}{3}, normal{13}{3},...
        normal{14}{3}, normal{15}{3}, normal{16}{3}}},...
        chans{j}, 'labels', {'Visual', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
        'plotstd', 'fill', 'yticksize', 4);
    
    savefig(gcf, [input_path '\11subjects_incl_EMS_normal_' num2str(baseline(1)) '-' num2str(baseline(2)) '_' chans{j}]);
    close(gcf);
end

% % plot to validate difference plot
% plot_erp({{conflict{2}{1}, conflict{3}{1}, conflict{6}{1}, conflict{7}{1},...
%         conflict{8}{1}, conflict{11}{1}, conflict{12}{1}, conflict{13}{1},...
%         conflict{14}{1}, conflict{15}{1}, conflict{16}{1}},...
%         
%         {normal{2}{1}, normal{3}{1}, normal{6}{1}, normal{7}{1},...
%         normal{8}{1}, normal{11}{1}, normal{12}{1}, normal{13}{1},...
%         normal{14}{1}, normal{15}{1}, normal{16}{1}}},...        
%         chans{1}, 'labels', {'Visual Only Conflict', 'Visual Only Normal'},...
%         'plotstd', 'fill', 'plotdiff', 1, 'delaycorrection', 'yticksize', 4);

%% plotting of difference ERPs

input_path = [study_folder single_subject_analysis_folder single_subject_analysis_folder_ERPs single_subject_analysis_folder_epochs];
chans = {'Fp1' 'FCz' 'Fz' 'Cz' 'Pz' 'Oz'};

for j = 1:length(chans)

    c = find(strcmp({normal{2}{1}.chanlocs.labels}, chans{j}));

    for subject = subjects
        conflict{subject}{1}.data(c,:,1) = mean(conflict{subject}{1}.data(c,:,:),3)...
            - mean(normal{subject}{1}.data(c,:,:),3);
        conflict{subject}{2}.data(c,:,1) = mean(conflict{subject}{2}.data(c,:,:),3)...
            - mean(normal{subject}{2}.data(c,:,:),3);
        conflict{subject}{3}.data(c,:,1) = mean(conflict{subject}{3}.data(c,:,:),3)...
            - mean(normal{subject}{3}.data(c,:,:),3);

        diff{subject}{1} = pop_selectevent(conflict{subject}{1}, 'epoch', 1);
        diff{subject}{2} = pop_selectevent(conflict{subject}{2}, 'epoch', 1);
        diff{subject}{3} = pop_selectevent(conflict{subject}{3}, 'epoch', 1);
    end

    % Vis, vis+vibro, vis+vibro+ems first difference within feedback condition 
    plot_erp_LG({{diff{2}{1}, diff{3}{1}, diff{6}{1}, diff{7}{1},...
        diff{8}{1}, diff{11}{1}, diff{12}{1}, diff{13}{1},...
        diff{14}{1}, diff{15}{1}, diff{16}{1}},...

        {diff{2}{2}, diff{3}{2}, diff{6}{2}, diff{7}{2},...
        diff{8}{2}, diff{11}{2}, diff{12}{2}, diff{13}{2},...
        diff{14}{2}, diff{15}{2}, diff{16}{2}},...

        {diff{2}{3}, diff{3}{3}, diff{6}{3}, diff{7}{3},...
        diff{8}{3}, diff{11}{3}, diff{12}{3}, diff{13}{3},...
        diff{14}{3}, diff{15}{3}, diff{16}{3}}},...
        chans{j}, 'labels', {'Visual', 'Visual+Vibro', 'Visual+Vibro+EMS'}, ...
        'plotstd', 'fill', 'yticksize', 4);
    
    savefig(gcf, [input_path '\11subjects_incl_EMS_difference_ytick_4_' num2str(baseline(1)) '-' num2str(baseline(2)) '_' chans{j}]);
    close(gcf);
end

%% extract peaks and lats for stats analysis

chans = {'FCz' 'Fz' 'Cz' 'Pz' 'Oz'};
count = 1;
j = 3;

c = find(strcmp({diff{2}{1}.chanlocs.labels}, chans{j}));

for subject = subjects

    diff_filt{subject}{1} = pop_eegfiltnew(diff{subject}{1}, low, high);
    diff_filt{subject}{2} = pop_eegfiltnew(diff{subject}{2}, low, high);
    diff_filt{subject}{3} = pop_eegfiltnew(diff{subject}{3}, low, high);
    
    diff_filt_vis_erps(count,:) = diff_filt{subject}{1}.data(c,:);
    diff_filt_vibro_erps(count,:) = diff_filt{subject}{2}.data(c,:);
    diff_filt_ems_erps(count,:) = diff_filt{subject}{3}.data(c,:);
    
    diff_vis_erps(count,:) = diff{subject}{1}.data(c,:);
    diff_vibro_erps(count,:) = diff{subject}{2}.data(c,:);
    diff_ems_erps(count,:) = diff{subject}{3}.data(c,:);
    
    count = count +1;
end

% find peaks and locs of peaks
clear peaks_locs
%searchwindow_min = [95:115]; % 80 to 160 ms post stimulus
%searchwindow_min = [100:120]; % 100 to 180 ms post stimulus
%searchwindow_min = [105:125]; % 120 to 200 ms post stimulus
searchwindow_min = [100:150]; % 120 to 200 ms post stimulus
searchwindow_max = [100:200];

for j = 1:3
    
    peaks_locs{1,1} = 'Participant';
    peaks_locs{1,2} = 'Condition';
    peaks_locs{1,3} = 'Min_Peak_Amplitude';
    peaks_locs{1,4} = 'Min_Peak_Latency';
    peaks_locs{1,5} = 'Max_Peak_Amplitude';
    peaks_locs{1,6} = 'Max_Peak_Latency';
    
    switch j
        case 1
            cond = 'visual';
            lat_find = diff_filt_vis_erps;
            erps = diff_vis_erps;
        case 2
            cond = 'vibro';
            lat_find = diff_filt_vibro_erps;
            erps = diff_vibro_erps;
        case 3
            cond = 'ems';
            lat_find = diff_filt_ems_erps;
            erps = diff_ems_erps;
    end
            
    for k=1:size(diff_vis_erps, 1)

        peak_min = min(lat_find(k,searchwindow_min));
        peak_max = max(lat_find(k,searchwindow_max));
        
        loc_min = find(lat_find(k,:)==peak_min);
        loc_max = find(lat_find(k,:)==peak_max);
        
        peak_min = mean(erps(k,loc_min-2:loc_min+2));%erps(k,loc_min);%
        peak_max = mean(erps(k,loc_max-2:loc_max+2));%erps(k,loc_max)
        
        row = k+length(subjects)*(j-1)+1;
        
        peaks_locs{row,1} = subjects(k);
        peaks_locs{row,2} = cond;
        peaks_locs{row,3} = peak_min;
        peaks_locs{row,4} = ((loc_min / 250) - .3)*1000;
        peaks_locs{row,5} = peak_max;
        peaks_locs{row,6} = ((loc_max / 250) - .3)*1000;

    end
end

j = 3;
cell2csv([ 'P:\Project_Sezen\data_processing\ERP\peaks_locs_ERP_incl_EMS_' chans{j} '_' num2str(searchwindow_min(1)) '.csv'], peaks_locs, ';');

