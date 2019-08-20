function [ EEG ] = parseEventsPredError( EEG )
%parseEventsPredError Parses long EEG.event.type string into separate
%EEG.event fields and their value. Long EEG.event.type string means string
%of the kind "ke1:value1;key2:value2;key3:value3...". This function catches
%some specifics in the event generation of PredError Experiment Ver1.0 (08/2018)

if isempty(EEG.event)
    EEG.event = EEG.urevent;
end

% delete duplicate events in EMS condition
for i = 1:length(EEG.event)
   if strfind(EEG.event(i).type, 'box:touched')
       first_touch = i;
       second_touch = [];
       spawn = [];
       
       % find second touch event occuring before spawn of next trial
       try
           for j = 1:100
               if strfind(EEG.event(i+j).type, 'box:touched')
                   second_touch = i+j;
                   break
               end
           end
           
           for j = 1:100
               if strfind(EEG.event(i+j).type, 'box:spawned')
                   spawn = i+j;
                   break
               end
           end

           % if second touch found before next trial, delete first touch
           % event, use second touch event for parsing, contains EMS info
           if spawn > second_touch
               EEG.event(first_touch).type = 'duplicate_event';
           end
       end
   end
end

for i = 1:length(EEG.event)
    
    current_event = cellstr(strsplit(EEG.event(i).type, ';'));

    % parse all events
    try
        for j=1:length(current_event)
            key_val = cellstr(strsplit(current_event{j}, ':'));
            EEG.event(i).(key_val{1}) = key_val{2};
            if j==1
                EEG.event(i).type = strcat(key_val{1}, ':', key_val{2});
            end
        end
    end
    
    % maintain last 'box:spawned' event to overwrite condition param
    if strcmp(EEG.event(i).type, 'box:spawned')
        last_spawned = i;
    end
    
    if isfield(EEG.event(i), 'emsFeedback') && strcmp(EEG.event(i).emsFeedback, 'on')
        EEG.event(i).condition = 'ems';
        EEG.event(last_spawned).condition = 'ems';
    end
end