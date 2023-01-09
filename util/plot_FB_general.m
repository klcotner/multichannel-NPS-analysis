% generate plots of the matched filterbank & system model,
%   when the filter bank includes templates from one or more codes

% seq_row_starts = vector of starting row indices for each unique signal ID
function [fig, ax1, ax2] = plot_FB_general(bFB, uFB, fs, seq_row_starts)

    if nargin<4
        seq_row_starts = [];
    end

    [fig, ax1, ax2] = plot_FB(bFB, uFB, fs);

    if length(seq_row_starts)>1

        xlim1 = get(ax1,'xlim');
        xlim2 = get(ax2,'xlim');

        for seqstart = seq_row_starts

            % ax1: add horizontal lines that separate templates from different signal IDs
            hold(ax1, 'on');
            plot(ax1, xlim1, [seqstart,seqstart], 'k','linewidth',1);

            % ax1: add horizontal lines that separates templates from different signal IDs
            hold(ax2, 'on');
            plot(ax2, xlim2, [seqstart,seqstart], 'k','linewidth',1);

        end

    end
    
end