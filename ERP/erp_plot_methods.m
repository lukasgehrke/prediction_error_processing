% load an EEG (eeglab) set containing epochs
% so size(EEG.data) must be 3 dimensional

dat = EEG.data;
size(dat)

% select channel index to plot
channel_to_plot = channel_to_plot;

% at what latency (in samples) is the zero event?
% plots a vertical line there
zero_ix = 75;

figure;
ylim([-15 15]);
title('ERP')
xlabel('samples')
ylabel('amplitude')

fname = 'C:\Users\Lukas\Desktop\epoch_smooth_chi_2_cond.mp4';

vidfile = VideoWriter(fname,'MPEG-4');
open(vidfile);

hold on
for i = 1:size(dat,3)
    
%     if i < 300
%         p = plot(dat(channel_to_plot,:,i), '-g');
%     else
%         p = plot(dat(channel_to_plot,:,i), '-r');
%     end
    i

    if i <= round((size(dat,3)/2))
        p = plot(dat(channel_to_plot,:,i), '-b');
        p1 = plot(mean(dat(channel_to_plot,:,1:i),3), '-b');
    else
        p = plot(dat(channel_to_plot,:,i), '-r');
        p1 = plot(mean(dat(channel_to_plot,:,1:i),3), '-r');    
    end
    p.Color(4) = 0.02;
    
    p1.LineWidth = i/size(dat,3);
    p1.Color(4) = i/size(dat,3)/2;
   
    vline(zero_ix, '-r');
    
    frame = getframe(gcf);
    writeVideo(vidfile,frame);
    clear p p1 frame
end

hold off
close(vidfile)


%%

% raw
dat = EEG.data(channel_to_plot,400000:402000);
figure;
plot(dat);

% filtered
dat_filt = movmean(dat,5);
figure;
plot(dat_filt);

% time-locked
vline(250, '-r')
vline(750, '-r')
vline(1250, '-r')
vline(1750, '-r')

% epoched
wins = [150:400, channel_to_plot0:900, 1150:1400, 1channel_to_plot0:1901];
for i = 1:size(dat_filt,2)
    if intersect(i, wins)
        remove(i) = 0;
    else 
        remove(i) = 1;
    end
end
dat_wins = dat_filt;
dat_wins(logical(remove)) = NaN;
figure;
plot(dat_wins);
vline(250, '-r')
vline(750, '-r')
vline(1250, '-r')
vline(1750, '-r')

% mean
wins = [150:399, channel_to_plot0:899, 1150:1399, 1channel_to_plot0:1899];
epochs = dat_filt(wins);
epochs = reshape(epochs, 4, 250);
epochs = movmean(mean(epochs, 1),10);
figure;
plot(epochs);

