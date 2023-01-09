% generate plots of the matched filterbank & system model

function [fig, ax1, ax2] = plot_FB(bFB, uFB, fs)

    fig = figure;
    
    % matched filterbank (bipolar signal template, scaled to norm)
    ax1 = subplot(1,2,1);
    stoptime = size(bFB,2) / fs; % convert sample # to seconds
    img1 = imagesc(bFB, 'xdata', [0, stoptime]); % make x-axis in time instead of sample #
    title('Matched Filterbank (bFB, bipolar & scaled)');
    xlabel('time (sec)');
    ylabel('template #');
    colormap(ax1, 'jet');
    c1 = colorbar('southoutside');
    c1.Label.String = 'signal value';
    
    % system model (unipolar signal template)
    ax2 = subplot(1,2,2);
    stoptime = size(uFB,2) / fs; % convert sample # to seconds
    img2 = imagesc(uFB, 'xdata', [0, stoptime]); % make x-axis in time instead of sample #
    title('System Model (uFB, unipolar)');
    xlabel('time (sec)');
    ylabel('template #');
    colormap(ax2, 'summer');
    c2 = colorbar('southoutside');
    c2.Label.String = 'signal value';

end