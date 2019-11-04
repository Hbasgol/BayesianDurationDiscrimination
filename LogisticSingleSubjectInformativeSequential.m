%% To check whether there is a sequential effect between trials for a subject

clear;

%% User input

% Data
load vanDrielData2015 d
% Subject and condition subjectList = [7 9 12 5 6 3]
% Subjects encoded as [A B C D E F]
subject = 12;
condition = 'auditory';

% JAGS parameters
% Sampling settings
doParallel = 1; % parallization
nThin = 1;
nChains = 3;
nSamples = 5e3;
nBurnin = 1e3;

%% Apply graphical model

% Graphical model script
modelName = 'LogisticSingleSubjectInformativeSequential';

% Data for single subject and condition to supply to graphical model
switch condition
    case 'auditory', dataStruct = struct('y', d.auditoryDecision(subject, :), ...
            'stimulus',d.auditoryStimulus(subject, :), ...
            'standard', d.auditoryStandard, ...
            'nTrials', d.nTrials);
    case 'visual', dataStruct = struct('y', d.visualDecision(subject, :), ...
            'stimulus',d.visualStimulus(subject, :), ...
            'standard', d.visualStandard, ...
            'nTrials', d.nTrials);
end

% Parameters to monitor and initial values
monitorParameters = {'alpha', 'beta', 'epsilon', 'epsilonPrior'};
for i = 1:nChains
    S.alpha = 0;
    S.beta = 10;
    S.epsilon = 0;
    init0(i) = S;
end

% Use JAGS to sample
if exist(['MCMCResults/' modelName '_' int2str(subject) '_' condition '.mat'])
    load(['MCMCResults/' modelName '_' int2str(subject) '_' condition '.mat'], 'samples', 'stats');
else
    tic % timing
    fprintf( 'Running JAGS ...\n' ); % display
    [samples, stats] = matjags( ...
        dataStruct, ...
        fullfile(pwd, [modelName '.txt']), ...
        init0, ...
        'doparallel', doParallel, ...
        'nchains', nChains,...
        'nburnin', nBurnin,...
        'nsamples', nSamples, ...
        'thin', nThin, ...
        'monitorparams', monitorParameters,...
        'savejagsoutput', 1 , ...
        'verbosity', 1 , ...
        'cleanup' , 0, ...
        'workingdir', 'C:/tmpjags' );
    toc % timing
    save(['MCMCResults/' modelName '_' int2str(subject) '_' condition '.mat'], 'samples', 'stats');
end
disp(stats.Rhat); % convergence check

%% Display

% Graphics constants
%load HandbookColors
range = 400; % range of stimulus values around standard
offset = 0.05; % offset on y axis for visibility
nDraw = 50; % number of posterior samples to draw
x = dataStruct.standard - range : 1 : dataStruct.standard + range;
scale = 2; % magnifier for areas in behavioral data points
alphaRange = [-150 150]; % axis limits for parameter space
betaRange = [0 150];
epsilonRange = [-0.2 0.2];
nBins = 40; % number of bins for marginals
nEpsilonBins = 20;
criticalEpsilon = 0;
rng(1); % set random number seed so get same posterior samples

% Setup figure
figure(subject); clf; hold on;
set(gcf,'color', 'w', 'units', 'normalized', ...
    'position', [0.2 0.2 0.6 0.6], 'paperpositionmode','auto');
axis([dataStruct.standard - range dataStruct.standard + range -offset 1]);
set(gca, 'units', 'normalized', 'position', [0.15 0.15 0.7 0.7], ...
    'xtick', dataStruct.standard - range : 100 : dataStruct.standard + range, ...
    'ytick', 0 : 0.2 : 1, ...
    'fontsize', 12, 'box', 'off', 'tickdir', 'out');
xlabel('Stimulus', 'fontsize', 15);
ylabel('Longer Response Probability', 'fontsize', 15);
H = plot(min(get(gca, 'xlim'))*ones(1,2), [0 -offset], 'w-', 'linewidth', 2);

% Draw posterior psychometric function and standard
permSamples = randperm(nSamples);
for i = 1:nDraw
    alpha = samples.alpha(1, permSamples(i));
    beta = samples.beta(1, permSamples(i));
    H = plot(x, 1./(1 + exp(-(x - dataStruct.standard -alpha)/beta)), 'k-');
    set(H, 'color', 'b', 'linewidth', 0.5);
end
H = plot(dataStruct.standard*ones(1,2), get(gca, 'ylim'), 'k--');
set(H, 'color', 'g', 'linewidth', 2);

% Draw behavioral data
uniqueStimuli = unique(dataStruct.stimulus);
for i = 1:length(uniqueStimuli)
    match = find(dataStruct.stimulus == uniqueStimuli(i));
    seen = length(match);
    above = length(find(dataStruct.y(match) == 1));
    H = plot(uniqueStimuli(i), above/seen, 'ko');
    set(H, 'markersize', scale*sqrt(seen), ...
        'markerfacecolor', 'r', ...
        'markeredgecolor', 'r', 'linewidth', 1.5);
    if strmatch(condition, 'visual', 'exact');
        set(H, 'marker', 'o');
        set(H, 'markersize', scale*sqrt(seen), ...
            'markerfacecolor', 'r', ...
            'markeredgecolor',  'w', 'linewidth', 0.25);
    end
end

% Draw inset parameter space
% joint
AX = axes; hold on;
set(AX, 'units', 'normalized', 'position', [0.225 0.7 0.15 0.25], ...
    'xtick', alphaRange, ...
    'ytick', betaRange, ...
    'fontsize', 11, 'box', 'on', 'ticklength', [0 0]);
axis([alphaRange betaRange]);
xlabel('\alpha', 'fontsize', 16,'verticalalignment', 'bottom');
ylabel('\beta', 'fontsize', 16, 'verticalalignment', 'middle');
H = plot(samples.alpha(1, permSamples(1:nDraw)), samples.beta(1, permSamples(1:nDraw)), 'k+');
set(H, 'color', 'g', 'markersize', 4);
% marginal alpha
AXT = axes; hold on;
set(AXT, 'units', 'normalized', 'position', [0.225 0.95 0.15 0.025], ...
    'xlim', alphaRange, 'ylim', [0 1]);
count = hist(samples.alpha(:), alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2));
H = bar(alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2), count/max(count), 0.7);
set(H, 'facecolor', 'g', 'edgecolor', 'g');
axis off;
% marginal beta
AXR = axes; hold on;
set(AXR, 'units', 'normalized', 'position', [0.375 0.7 0.015 0.25], ...
    'ylim', betaRange, 'xlim', [0 1]);
count = hist(samples.beta(:), betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2));
H = barh(betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2), count/max(count), 0.7);
set(H, 'facecolor', 'g', 'edgecolor','g');
axis off;

% Draw posterior and prior for epilson sequential dependency
AXE = axes; hold on;
set(AXE, 'units', 'normalized', 'position', [0.225 0.45 0.15 0.175], ...
    'xlim', epsilonRange, 'ycolor', 'w');
epsilonBins = epsilonRange(1):(epsilonRange(2) - epsilonRange(1))/nEpsilonBins:epsilonRange(2);
count = hist(samples.epsilon(:), epsilonBins);
H = bar(epsilonBins, count, 0.7);
set(H, 'facecolor', 'g', 'edgecolor', 'g', ...
    'basevalue', 0.01*max(count));
countPrior = hist(samples.epsilonPrior(:), epsilonRange(1):(epsilonRange(2) - epsilonRange(1))/nEpsilonBins:epsilonRange(2));
H = stairs(epsilonRange(1):(epsilonRange(2) - epsilonRange(1))/nEpsilonBins:epsilonRange(2), countPrior);
set(H, 'color', 'k');
set(AXE, 'ylim', [0 max(max(count), max(countPrior))]);
xlabel('\epsilon', 'fontsize', 16,'verticalalignment', 'middle');
% Savage-Dickey Bayes factor
[~, match] = min(abs(epsilonBins - criticalEpsilon));
BF = count(match)/countPrior(match);
fprintf('Bayes factor in favor of null is %1.2f\n', BF);

%% Print
print(['Images/' modelName '_' int2str(subject) '_' condition '.eps'],'-depsc');
print(['Images/' modelName '_' int2str(subject) '_' condition '.png'],'-dpng');

