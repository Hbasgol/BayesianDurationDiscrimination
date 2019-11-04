%% For drawing behavioral data for subjects

clear;

%% Codes of subjects that will be drawn from data
subjectList = [7 9 12 5 6 3];
% subjectList = 1:6;
nRows = 3; nCols = 2;

%% Data
load vanDrielData2015 d
standard = d.auditoryStandard;

%% Display
% Graphics constants
%load HandbookColors
range = 400; % range of stimulus values around standard
offset = 0.1; % offset on y axis for visibility
nDraw = 500; % number of posterior samples to draw
x = standard - range : 1 : standard + range;
scale = 1.25; % magnifier for areas in behavioral data points
alphaRange = [-80 80]; % axis limits for parameter space
betaRange = [-80 80];
nBins = 40; % number of bins for marginals
rng(1); % set random number seed so get same posterior samples

% Setup figure
figure(1); clf; hold on;
set(gcf,'color', 'w', 'units', 'normalized', ...
    'position', [0.2 0.1 0.45/1.0 0.8/1.0], 'paperpositionmode','auto');

for subjectIndex = 1:length(subjectList)
    subplot(nRows, nCols, subjectIndex); cla; hold on;
    axis([standard - range standard + range -offset 1]);
    set(gca, 'xtick', standard - range : 200 : standard + range, ...
        'ytick', 0 : 0.2 : 1, 'ticklength', [0.02 0], ...
        'fontsize', 12, 'box', 'off', 'tickdir', 'out',...
        'position', get(gca, 'position') + [0 0.025 0 0]);
    % legend on first axis
    if subjectIndex == 1
        plot(-100, -100, 'k+');
        plot(-100, -100, 'ko');
        [L, O] = legend('Auditory', 'Visual');
        set(L, 'fontsize', 12, 'box', 'off', 'location', 'northwest');
        set(O(1), 'position', get(O(1), 'position') + [-0.4 0.3 0]);
        set(O(2), 'position', get(O(2), 'position') + [-0.4 0.3 0]);
        set(O(4),  'markerfacecolor', 'b', ...
            'markeredgecolor', 'b', ...
            'XData', get(O(4), 'XData') -0.25, ...
            'YData', get(O(4), 'YData') + 0.3);
        set(O(6), 'markerfacecolor', 'r', ...
            'markeredgecolor', 'r', ...
            'XData', get(O(6), 'XData') -0.25, ...
            'YData', get(O(6), 'YData') + 0.3);
    end
    % white space to separate axes
    H = plot(min(get(gca, 'xlim'))*ones(1,2), [0 -offset], 'w-', 'linewidth', 2);
    % draw standard
    H = plot(standard*ones(1,2), get(gca, 'ylim'), 'k--');
    set(H, 'color', 'g', 'linewidth', 2);
    
    % subject and label
    subject = subjectList(subjectIndex);
    T = text(900, 0, char(64 + subjectIndex));
    set(T, 'horizontalalignment', 'right', 'fontsize', 13, 'fontweight', 'bold');
    
    % Draw visual behavioral data
    uniqueStimuli = unique(d.visualStimulus(subject, :));
    for i = 1:length(uniqueStimuli)
        match = find(d.visualStimulus(subject, :) == uniqueStimuli(i));
        seen = length(match);
        above = length(find(d.visualDecision(subject, match) == 1));
        H = plot(uniqueStimuli(i), above/seen, 'ko');
        set(H, 'markersize', scale*sqrt(seen), ...
            'markerfacecolor', 'r', ...
            'markeredgecolor',  'r', 'linewidth', 0.25);
    end
    
    % Draw auditory behavioral data
    uniqueStimuli = unique(d.auditoryStimulus(subject, :));
    for i = 1:length(uniqueStimuli)
        match = find(d.auditoryStimulus(subject, :) == uniqueStimuli(i));
        seen = length(match);
        above = length(find(d.auditoryDecision(subject, match) == 1));
        H = plot(uniqueStimuli(i), above/seen, 'k+');
        set(H, 'markersize', scale*sqrt(seen), ...
            'markerfacecolor', 'b', ...
            'markeredgecolor', 'b', 'linewidth', 1.5);
    end
    
    
end

%% Print
print('Images/BehavioralData_1.eps', '-depsc');
print('Images/BehavioralData_1.png', '-dpng');

