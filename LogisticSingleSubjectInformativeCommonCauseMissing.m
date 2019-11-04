%% Logistic with informative prior, 
%% joint model of visual and auditory
%% giving a number of trials, this model makes predictions regarding other behaviors

clear;

%% Specifying parameters regarding data

% Data to be used
load vanDrielData2015 d

% Subject and condition subjectList = [7 9 12 5 6 3]
% Subjects encoded as [A B C D E F]
subject = 12; % missing data at trial 44 of auditory
observedTrials = 60;
condition = 'joint';

% JAGS parameters
% Sampling settings

doParallel = 1; % parallization
nThin = 1;
nChains = 3;
nSamples = 5e3;
nBurnin = 1e3;

%% Apply graphical model

% Graphical model script
modelName = 'LogisticSingleSubjectInformativeCommonCauseMissing';

% Data for single subject, both conditions, to supply to graphical model
dataStruct = struct('yA', [d.auditoryDecision(subject, 1:observedTrials) nan*ones(1, d.nTrials - observedTrials + 1)], ...
            'stimulusA',d.auditoryStimulus(subject, :), ...
            'standardA', d.auditoryStandard, ...
            'nTrialsA', d.nTrials, ...
            'yV', nan*ones(1, d.nTrials), ...
            'stimulusV', d.visualStimulus(subject, :), ...
            'standardV', d.visualStandard, ...
            'nTrialsV', d.nTrials);

% Parameters to monitor and initial values
monitorParameters = {'alpha', 'beta', 'yAPredicted', 'yVPredicted'};
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
nDraw = 50; % number of posterior samples to draw
x = dataStruct.standardA - range : 1 : dataStruct.standardA + range;
scale = 2; % magnifier for areas in behavioral data points
alphaRange = [-150 150]; % axis limits for parameter space
betaRange = [0 150];
nBins = 40; % number of bins for marginals
rng(1); % set random number seed so get same posterior samples
permSamples = randperm(nSamples);

% Setup figure
figure(subject); clf; hold on;
set(gcf,'color', 'w', 'units', 'normalized', ...
    'position', [0.2 0.2 0.6 0.6], 'paperpositionmode','auto');
[~, hostName] = system('hostname');
if strcmp(deblank(hostName), 'C16050500')
set(gcf,'color', 'w', 'units', 'normalized', ...
    'position', [0.2 0.2 0.6/1.2 0.6/1.2], 'paperpositionmode','auto');
end

% Draw parameter space
AX = gca; hold on;
set(gca, 'units', 'normalized', 'position', [0.075 0.225+.1 0.35-.2/1.5 0.575-.2], ...
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
set(gca, 'units', 'normalized', 'position', [0.075 0.325+0.375 0.2167 0.1], ...
    'xlim', alphaRange);
count = hist(samples.alpha(:), alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2));
H = bar(alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2), count/ max(count), 0.7);
set(H, 'facecolor', 'g', 'edgecolor', 'g');
axis off;
% marginal beta
AXR = axes; hold on;
set(gca, 'units', 'normalized', 'position', [0.075+0.2167 0.325 0.1/1.5 0.575-.2], ...
    'ylim', betaRange);
count = hist(samples.beta(:), betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2));
H = barh(betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2), count/ max(count), 0.7);
set(H, 'facecolor', 'g', 'edgecolor', 'g');
axis off;

% Auditory partially observed (including missing) data and predictions
AX = axes; hold on;
axis([0 d.nTrials+1 0 1]);
set(gca, 'units', 'normalized', 'position', [0.425 0.6 0.5 0.3], ...
    'xtick', [1 observedTrials 240], 'ytick', 0:0.5:1, ...
    'fontsize', 11, 'box', 'off', 'ticklength', [0.02 0], 'tickdir', 'out');
correct = stats.mean.yAPredicted;
correct(find(d.auditoryDecision(subject, :) == 0 )) = 1 - correct(find(d.auditoryDecision(subject, :) == 0 ));
H = plot(1:d.nTrials, correct,' kx');
set(H, 'markersize', 5, 'color', 'b', 'linewidth', 0.25);
H = plot([1 observedTrials], ones(1,2)*mean(correct(1:observedTrials)), 'k--');
fprintf('Observed accuracy is %1.2f\n', mean(correct(1:observedTrials)));
H = plot([observedTrials+1 d.nTrials], ones(1,2)*mean(correct(observedTrials+1:d.nTrials)), 'k--');
fprintf('Prediction accuracy is %1.2f\n', mean(correct(observedTrials+1:d.nTrials)));

% Visual generalization data and predictions
AX = axes; hold on;
axis([0 d.nTrials+1 0 1]);
set(gca, 'units', 'normalized', 'position', [0.425 0.1 0.5 0.3], ...
    'xtick', [1 240], 'ytick', 0:0.5:1, ...
    'fontsize', 11, 'box', 'off', 'ticklength', [0.02 0], 'tickdir', 'out');
correct = stats.mean.yVPredicted;
correct(find(d.visualDecision(subject, :) == 0 )) = 1 - correct(find(d.visualDecision(subject, :) == 0 ));
H = plot(1:d.nTrials, correct,' kx');
set(H, 'markersize', 5, 'color', 'b', 'linewidth', 0.25);
H = plot([1 d.nTrials], ones(1,2)*mean(correct(1:d.nTrials)), 'k--');
fprintf('Generalization accuracy is %1.2f\n', mean(correct(1:d.nTrials)));


%% Print
print(['Images/' modelName '_' int2str(subject) '.eps'],'-depsc');
print(['Images/' modelName '_' int2str(subject) '.png'],'-dpng');



