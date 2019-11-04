%% Logistic with informative prior, joint model of visual and auditory
%% It assumes that two modalities use same psychometric function, alpha and beta parameters

clear;

%% Specifying parameters regarding data

% Data to be used
load vanDrielData2015 d

% Subject and condition subjectList = [7 9 12 5 6 3]
% Subjects encoded as [A B C D E F]
subject = 12;
condition = 'commoncause';

% JAGS parameters
% Sampling settings
doParallel = 1; % parallization
nThin = 1;
nChains = 3;
nSamples = 5e3;
nBurnin = 1e3;

%% Apply graphical model

% Graphical model script
modelName = 'LogisticSingleSubjectInformativeCommonCause';

% Data for single subject, both conditions, to supply to graphical model
dataStruct = struct('yA', d.auditoryDecision(subject, :), ...
            'stimulusA',d.auditoryStimulus(subject, :), ...
            'standardA', d.auditoryStandard, ...
            'nTrialsA', d.nTrials, ...
            'yV', d.visualDecision(subject, :), ...
            'stimulusV',d.visualStimulus(subject, :), ...
            'standardV', d.visualStandard, ...
            'nTrialsV', d.nTrials);

% Parameters to monitor and initial values
monitorParameters = {'alpha','beta'};
for i = 1:nChains
    S.alpha = 0;
    S.beta = 10;
    init0(i) = S;
end

% Use JAGS to sample
if exist(['MCMCResults/' modelName '_' int2str(subject) '_joint.mat'])
    load(['MCMCResults/' modelName '_' int2str(subject) '_joint.mat'], 'samples', 'stats');
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
    save(['MCMCResults/' modelName '_' int2str(subject) '_joint.mat'], 'samples', 'stats');
end
disp(stats.Rhat); % convergence check

%% Display

% Graphics constants
%load HandbookColors
range = 400; % range of stimulus values around standard
offset = 0.05; % offset on y axis for visibility
nDraw = 25; % number of posterior samples to draw
x = dataStruct.standardA - range : 1 : dataStruct.standardA + range;
scale = 2; % magnifier for areas in behavioral data points
alphaRange = [-50 50]; % axis limits for parameter space
betaRange = [0 70];
nBins = 40; % number of bins for marginals
rng(3); % set random number seed so get same posterior samples

% Setup figure
figure(subject); clf; hold on;
set(gcf,'color', 'w', 'units', 'normalized', ...
    'position', [0.2 0.2 0.6 0.6], 'paperpositionmode','auto');
[~, hostName] = system('hostname');
if strcmp(deblank(hostName), 'C16050500')
set(gcf,'color', 'w', 'units', 'normalized', ...
    'position', [0.2 0.2 0.6/1.2 0.6/1.2], 'paperpositionmode','auto');
end
axis([dataStruct.standardA - range dataStruct.standardA + range -offset 1]);
set(gca, 'units', 'normalized', 'position', [0.15 0.15 0.7 0.7], ...
    'xtick', dataStruct.standardA - range : 100 : dataStruct.standardA + range, ...
    'ytick', 0 : 0.2 : 1, ...
    'fontsize', 12, 'box', 'off', 'tickdir', 'out');
xlabel('Stimulus Duration', 'fontsize', 15);
ylabel('Longer Response Probability', 'fontsize', 15);
H = plot(min(get(gca, 'xlim'))*ones(1,2), [0 -offset], 'w-', 'linewidth', 2);

% Draw posterior psychometric function and standard
permSamples = randperm(nSamples);
for i = 1:2*nDraw % twice as many to keep same count as independent scatter
    alpha = samples.alpha(1, permSamples(i));
    beta = samples.beta(1, permSamples(i));
    H = plot(x, 1./(1 + exp(-(x - dataStruct.standardA - alpha)/beta)), 'k-');
    set(H, 'color', [0,0.7,0.9] , 'linewidth', 0.5);
end
H = plot(dataStruct.standardA*ones(1,2), get(gca, 'ylim'), 'k--');
set(H, 'color', 'g', 'linewidth', 2);

% Draw auditory behavioral data
uniqueStimuli = unique(dataStruct.stimulusA);
for i = 1:length(uniqueStimuli)
    match = find(dataStruct.stimulusA == uniqueStimuli(i));
    seen = length(match);
    above = length(find(dataStruct.yA(match) == 1));
    H = plot(uniqueStimuli(i), above/seen, 'k+');
    set(H, 'markersize', scale*sqrt(seen), ...
        'markerfacecolor', 'r', ...
        'markeredgecolor', 'r', 'linewidth', 1.5);
    if strmatch(condition, 'visual', 'exact');
        set(H, 'marker', 'o');
        set(H, 'markersize', scale*sqrt(seen), ...
            'markerfacecolor', 'b', ...
            'markeredgecolor',  'w', 'linewidth', 0.25);
    end
end

% Draw visual behavioral data
uniqueStimuli = unique(dataStruct.stimulusV);
for i = 1:length(uniqueStimuli)
    match = find(dataStruct.stimulusV == uniqueStimuli(i));
    seen = length(match);
    above = length(find(dataStruct.yV(match) == 1));
    H = plot(uniqueStimuli(i), above/seen, 'k+');
        set(H, 'marker', 'o');
        set(H, 'markersize', scale*sqrt(seen), ...
            'markerfacecolor', 'b', ...
            'markeredgecolor',  'w', 'linewidth', 0.25);
end

% Draw inset parameter space
% joint
AX = axes; hold on;
set(gca, 'units', 'normalized', 'position', [0.225 0.7 0.15 0.25], ...
    'xtick', alphaRange, ...
    'ytick', betaRange, ...
    'fontsize', 11, 'box', 'on', 'ticklength', [0 0]);
axis([alphaRange betaRange]);
xlabel('\alpha', 'fontsize', 16,'verticalalignment', 'bottom');
ylabel('\beta', 'fontsize', 16, 'verticalalignment', 'top');
H = plot(samples.alpha(1, permSamples(1:nDraw)), samples.beta(1, permSamples(1:nDraw)), 'k+');
set(H, 'color', 'g', 'markersize', 4);
% marginal alpha
AXT = axes; hold on;
set(gca, 'units', 'normalized', 'position', [0.225 0.95 0.15 0.025], ...
    'xlim', alphaRange);
count = hist(samples.alpha(:), alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2));
H = bar(alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2), count/ max(count), 0.7);
set(H, 'facecolor', 'g', 'edgecolor', 'g');
axis off;
% marginal beta
AXR = axes; hold on;
set(gca, 'units', 'normalized', 'position', [0.375 0.7 0.015 0.25], ...
    'ylim', betaRange);
count = hist(samples.beta(:), betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2));
H = barh(betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2), count/ max(count), 0.7);
set(H, 'facecolor', 'g', 'edgecolor', 'g');
axis off;

% Draw joint space showing visual and auditory independently inferred
% load relevant samples
modelName = 'LogisticSingleSubjectInformative';
load(['MCMCResults/' modelName '_' int2str(subject) '_visual.mat'], 'samples', 'stats');
samplesV = samples; statsV = stats;
load(['MCMCResults/' modelName '_' int2str(subject) '_auditory.mat'], 'samples', 'stats');
samplesA = samples; statsA = stats;
% axes
AX2 = axes; hold on;
set(AX2, 'units', 'normalized', 'position', [0.225 0.7-0.35 0.15 0.25], ...
    'xtick', alphaRange, ...
    'ytick', betaRange, ...
    'fontsize', 11, 'box', 'on', 'ticklength', [0 0]);
axis([alphaRange betaRange]);
xlabel('\alpha', 'fontsize', 16,'verticalalignment', 'bottom');
ylabel('\beta', 'fontsize', 16, 'verticalalignment', 'top');
H = plot(samplesV.alpha(1, permSamples(1:nDraw)), samplesV.beta(1, permSamples(1:nDraw)), 'k+');
set(H, 'color', 'g', 'markersize', 4);
H = plot(samplesA.alpha(1, permSamples(1:nDraw)), samplesA.beta(1, permSamples(1:nDraw)), 'k+');
set(H, 'color', 'g', 'markersize', 4);
% marginal alpha
AXT2 = axes; hold on;
set(gca, 'units', 'normalized', 'position', [0.225 0.7-0.35+.25 0.15 0.025], ...
    'xlim', alphaRange, 'ylim', [0 1]);
count = hist([samplesV.alpha(:);samplesA.alpha(:)], alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2));
H = bar(alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2), count/ max(count), 0.7);
set(H, 'facecolor', 'g', 'edgecolor', 'g');
axis off;
% marginal beta
AXR2 = axes; hold on;
set(gca, 'units', 'normalized', 'position', [0.375 0.7-0.35 0.015 0.25], ...
    'ylim', betaRange, 'xlim', [0 1]);
count = hist([samplesV.beta(:);samplesA.beta(:)], betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2));
H = barh(betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2), count/ max(count), 0.7);
set(H, 'facecolor', 'g', 'edgecolor', 'g');
axis off;
% change model name back
modelName = 'LogisticSingleSubjectInformativeCommobnCause';

%% Print
print(['Images/' modelName '_' int2str(subject) '.eps'],'-depsc');
print(['Images/' modelName '_' int2str(subject) '.png'],'-dpng');



