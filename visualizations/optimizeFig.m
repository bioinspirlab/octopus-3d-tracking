function [] = optimizeFig(h,figsize)
% This function optimizes a figure for legibility

% Input arguments (all optional):
% h is the figure handle
% figsize is an indication for final figure size: 'S','M','L'. this affects
% scaling of font sizes for legibility

%% Guess input arguments if not provided
if nargin < 2 || ~ischar(figsize)
    if nargin < 1
        h = gcf;
    end
    h.Units = 'pixels';
    % Auto-assign size assessment
    if max(h.OuterPosition(3:4))>1100
        figsize = 'L';
    elseif max(h.OuterPosition(3:4))>450
        figsize = 'M';
    else
        figsize = 'S';
    end
end

%% Set size-dependent settings
switch upper(figsize)
    case 'S'
        fontsize = 11;
    case 'M'
        fontsize = 14;
    case 'L'
        fontsize = 14;
    case 'XL'
        fontsize = 19;
    case 'XXL'
        fontsize = 25;
end

%% Initialize parameters
ax = h.CurrentAxes;

%% Set font size
ax.FontSize = fontsize;

%% Reduce borders
% outerpos = ax.OuterPosition
outerpos = [0 0 1 1];
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width*0.99 ax_height*0.99];

%% Change axis labels if it appears it is in radians

% Determine if X and Y data could be in radian, by looking at the maximum
% value of the children. If it is within a small distance from an integer
% times baseint (default: pi/2), assume radian
radianX = 'N';
radianY = 'N';
baseint = pi/2;
for i = 1:length(ax.Children)
    if sum(strcmp(ax.Children(i).Type,{'text','rectangle'}))==0
        if isfield(ax.Children(i),'XData')
        if (mod(max(ax.Children(i).XData)+baseint/100,baseint)<baseint/50) && ...
                (abs(max(ax.Children(i).XData))>baseint/3) && (max(ax.Children(i).XData)/baseint<9)
            radianX = 'Y';
        end
        if (mod(max(ax.Children(i).YData)+baseint/100,baseint)<baseint/50) && ...
                (abs(max(ax.Children(i).YData))>baseint/3) && (max(ax.Children(i).YData)/baseint<9)
            radianY = 'Y';
        end
        end
    end
end

if radianX == 'Y'
    ax.XTickMode = 'auto';
    numtick = length(ax.XTick);
    mintick = round(ax.XLim(1)/baseint)*baseint;
    maxtick = round(ax.XLim(2)/baseint)*baseint;
    int = baseint*round(((maxtick-mintick)/(numtick-1))/baseint);
    ticks = mintick:int:maxtick;
    ticklabels = {};
    for i = 1:length(ticks)
        ticklabels(i) = {[num2str(ticks(i)/pi) '\pi']};
        if ticks(i) == 0
            ticklabels(i) = {'0'};
        end
    end
    ax.XTick = ticks;
    ax.XTickLabels = ticklabels;
end

if radianY == 'Y'
    ax.YTickMode = 'auto';
    numtick = length(ax.YTick);
    mintick = round(ax.YLim(1)/baseint)*baseint;
    maxtick = round(ax.YLim(2)/baseint)*baseint;
    int = baseint*round(((maxtick-mintick)/(numtick-1))/baseint);
    ticks = mintick:int:maxtick;
    ticklabels = {};
    for i = 1:length(ticks)
        ticklabels(i) = {[num2str(ticks(i)/pi) '\pi']};
        if ticks(i) == 0
            ticklabels(i) = {'0'};
        end
    end
    ax.YTick = ticks;
    ax.YTickLabels = ticklabels;
end

h.SizeChangedFcn = @optimizeFig;

end


