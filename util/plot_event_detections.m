% plot detected events as rectangles overlaid on the data
% 
% if `model_baseline` is not provided or is empty, will plot the data without the baseline subtracted
%   > the event rectangles will still be flat, based on the baseline at the event center
% `plotnums` = optional boolean argument (default: false)
%   > if true, notes the event number on the plot
% `savepath` is an optional string argument (default: empty)
%   > if not empty, figure saved to a file at the provided path
%   > must be a supported filetype for `exportgraphics` (not .fig)

function fighandle = plot_event_detections(detected_events, timedata, data, model_baseline, plotnums, savepath)
    
    % parse inputs
    if nargin<6
        savepath = '';
        if nargin<5
            plotnums=false;
            if nargin<4 % just plot w/o baseline subtracted
                model_baseline = [];
            end
        end
    end
    
    %% misc
    
    % full-screen figure (normalized units)
    figpos_m1 = [0,0,1,1];
    
    % get cell info
    has_cells = ~isempty(detected_events) && ~isempty(fieldnames(detected_events));
    if ~has_cells
        num_events = 0;
    else
        num_events = numel(detected_events);
    end
    
    %% plot the data
    disp('plotting results with event patches...');
    
    xdata = timedata;
    if isempty(model_baseline)
        ydata = data;
    else
        ydata = data - model_baseline;
    end

    fighandle = figure('units','normalized', 'position', figpos_m1);
    plot(xdata, ydata, 'k', 'linewidth', 1);
    hold on;

    xlim([xdata(1), xdata(end)]);

    title('all detected events');
    xlabel('time [sec]');
    ylabel('resistance [\Omega]');
    set(gca,'fontsize',18);

    if has_cells
        
        % get seq info
        seqIDs = unique([detected_events.sequence_num], 'sorted');
        
        % determine # of colors needed
        if length(seqIDs)>1
            numcolors = length(seqIDs);
        else
            numcolors = num_events;
        end

        % get a list of colors that's long enough
        S_col = load('plotcolors.mat');
        colorlist = S_col.colorlist;
        while numcolors > size(colorlist,1)
            colorlist = [colorlist; colorlist];
        end
        
        %% plot each event
        for ii=1:num_events

            event = detected_events(ii);
            starttime = event.start_time; % [sec]
            endtime = event.end_time; % [sec]
            amp = event.amplitude;
            seq = event.sequence_num;
            
            if length(seqIDs)>1
                colorix = find(seqIDs==seq);
            else
                colorix = ii;
            end

            patchx = [starttime, starttime, endtime, endtime];
            patchy = [0,         amp,       amp,     0      ];
            
            if isempty(model_baseline)
                patchy = patchy + event.baseline_amp_center;
            end

            patch(patchx, patchy, colorlist(colorix,:), ...
                'facealpha', 0.5, 'edgecolor', colorlist(colorix,:), 'linewidth',2);
            
            if plotnums
                text(starttime, 0, num2str(ii), 'verticalalignment','top', ...
                    'backgroundcolor',colorlist(colorix,:), 'edgecolor','k', ...
                    'fontunits','normalized', 'fontsize',0.015);
            end

        end
        
        %% use invisible patches to create legend
        
        lgdpatches = [];
        
        if length(seqIDs)>1
            % cycle through each sequence ID and create an invisible patch
            for seqix = 1:length(seqIDs)
                seqname = sprintf('sequence #%u', seqIDs(seqix));
                p = patch(NaN, NaN, colorlist(seqix,:), 'facealpha',0.5, 'edgecolor','none', ...
                    'displayname', seqname);
                lgdpatches = [lgdpatches, p];
            end
            
        else
            % just create a white patch to indicate that all events are the same seqID
            seqname = sprintf('sequence #%u', seqIDs);
            lgdpatches = patch(NaN, NaN, 'w', 'facealpha',0.5, 'edgecolor','k', ...
                'displayname', seqname);
        end
        
        legend(lgdpatches);
            
    end

    disp('plotting finished');

    %% save the figure, if desired
    if ~isempty(savepath)
        disp('saving plot...');
        exportgraphics(fighandle, savepath);
        disp('concatenated plot saved');
    end
    
    %%
    disp('finished plotting full results');
end