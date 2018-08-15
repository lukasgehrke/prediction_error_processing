function [ EEG ] = parseEventsPredError( EEG )
%parseEventsPredError Parses long EEG.event.type string into separate
%EEG.event fields and their value. Long EEG.event.type string means string
%of the kind "ke1:value1;key2:value2;key3:value3...". This function catches
%some specifics in the event generation of PredError Experiment Ver1.0 (08/2018)

for i = 1:length(EEG.event)
    
    current_event = cellstr(strsplit(EEG.event(i).type, ';'));
    
    % clean up some mistakes made in writing markers to be automatically
    % parsed. Catch duplicate box:touched events in conditions vibro and ems
    % and add to previous touch event
    if i>1 && strcmp(EEG.event(i-1).type, 'box:touched')
        idx = i-1;
    else
        idx = i;
    end
    
    for j=1:length(current_event)
        key_val = cellstr(strsplit(current_event{j}, ':'));
        EEG.event(idx).(key_val{1}) = key_val{2};
        if j==1
            EEG.event(idx).type = strcat(key_val{1}, ':', key_val{2});
        end
    end
    
    if idx == i-1
        EEG.event(i).type = 'duplicate_event';
    end
    
    % clean up some mistakes made in writing markers to be automatically
    % parsed
    if isfield(EEG.event(idx), 'emsFeedback')
        EEG.event(idx).condition = 'ems';
    end
end

end

