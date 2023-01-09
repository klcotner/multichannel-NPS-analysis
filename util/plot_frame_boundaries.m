% * `timedata` = vector of timepoints for data [any units]
% * `data` = vector of resistance data [any units]
% * `Nstarts` = vector of ints
%       > start indices for all frames to plot
% * `Nstops` = vector of ints (same lengths as Nstarts)
%       > end indces for all frames to plot
% * `savepath` = char (optional): full path to where to save the .fig file
%   > default = ''
%   > if empty, won't save the figure to file

function fighandle = plot_frame_boundaries(timedata, data, Nstarts, Nstops, savepath)

    if nargin<5
        savepath = '';
    end
    
    % full-screen figure (normalized units)
    figpos_full = [0,0,1,1];
    
    fontsize = 18;

    %% plot frames based on Nstarts & Nstops
    disp('plotting frame boundaries...');
    
    fighandle = figure('units','normalized', 'position', figpos_full);
    
    plot(timedata, data, 'k');
    hold on;

    xlim([timedata(1), timedata(end)]);
    xlabel('time');
    ylabel('data');
    set(gca,'fontsize',fontsize);
    title('final frame boundaries');

    % plot the start of all frames
    for ii=1:length(Nstarts)
        startix = Nstarts(ii);

        xline(timedata(startix), 'r', 'linewidth',2);
    end
    % plot the end of the last frame
    xline(timedata(Nstops(end)), 'r', 'linewidth',2);
    
    % zoom in on the selected frames
    xlim( timedata([Nstarts(1),Nstops(end)]) );

    disp('finished plotting');
    
    %% save figure, if desired
    if ~isempty(savepath)
        disp('saving breakpoint plot...');
        exportgraphics(fighandle, savepath);
        disp('done saving breakpoint plot');
    end
    
end