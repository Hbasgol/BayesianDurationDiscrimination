%% Logistic with various reasonable informative priors
%% To test the behavior of the model in varios informative priors
%% Specifying parameters regarding data

% Data to be used
load vanDrielData2015 d

% JAGS parameters
% Sampling settings
doParallel = 1; % parallization
nThin = 1;
nChains = 3;
nSamples = 5e3;
nBurnin = 1e3;

%% Apply graphical model

% Graphical model script
modelName = 'LogisticInformativePrior';

% Data for prior predictive part of graphical model
subject = 7; % doesn't matter which one; all the same stimuli
standard = d.auditoryStandard;
stimulusDuration = unique(d.auditoryStimulus(subject, :));
nUniqueStimuli = length(stimulusDuration);
nStimulusTrials = zeros(nUniqueStimuli, 1);
for i = 1:nUniqueStimuli
    nStimulusTrials(i) = length(find(d.auditoryStimulus(subject, :) == stimulusDuration(i)));
end
dataStruct = struct('standard', standard, ...
    'nUniqueStimuli', nUniqueStimuli, ...
    'stimulusDuration', stimulusDuration, ...
    'nStimulusTrials', nStimulusTrials);

% Parameters to monitor and initial values
monitorParameters = {'alpha', 'beta', 'yPP'};
for i = 1:nChains
    S.alpha = 0;
    S.beta = 1;
    init0(i) = S;
end

% Use JAGS to sample
if exist(['MCMCResults/' modelName '_Prior.mat'])
    load(['MCMCResults/' modelName '_Prior.mat'], 'samples', 'stats');
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
        'dic', 0, ... % Normally, dic parameter is not given. However, should be zero
        'monitorparams', monitorParameters,...
        'savejagsoutput', 1 , ...
        'verbosity', 1 , ...
        'cleanup' , 0, ...
        'workingdir', 'tmpjags' );
    toc % timing
    save(['MCMCResults/' modelName '_Prior.mat'], 'samples', 'stats');
end
disp(stats.Rhat); % convergence check

%% Display
% Graphics constants

functionColor = 'b'; 
standardColor = 'g'; 
parameterColor = 'g'; 
% Constants
nDraw = 50; % number of posterior samples to draw
rng(24); % set random number seed so get same posterior samples
permSamples = randperm(nSamples);
scale = 20; % scale on predictive markers
range = 400; % range of stimulus values around standard
offset = 0.05; % offset on y axis for visibility
x = standard - range : 1 : standard + range;
alphaRange = [-150 150]; % axis limits for parameter space
betaRange = [0 300];
nBins = 40; % number of bins for marginals

%% Prior on parameters and psychophysical function figure

% Setup figure
figure(1); clf; hold on;
set(gcf,'color', 'w', 'units', 'normalized', ...
   'position', [0.2 0.2 0.6 0.6], 'paperpositionmode','auto');
axis([standard - range standard + range -offset 1]);
set(gca, 'units', 'normalized', 'position', [0.15 0.15 0.7 0.7], ...
   'xtick', standard - range : 100 : standard + range, ...
   'ytick', 0 : 0.2 : 1, ...
   'fontsize', 12, 'box', 'off', 'tickdir', 'out');
xlabel('Stimulus', 'fontsize', 15);
ylabel('Longer Response Probability', 'fontsize', 15);
H = plot(min(get(gca, 'xlim'))*ones(1,2), [0 -offset], 'w-', 'linewidth', 2);

% Draw posterior psychometric function and standard
for i = 1:nDraw
   alpha = samples.alpha(1, permSamples(i));
   beta = samples.beta(1, permSamples(i));
   H = plot(x, 1./(1 + exp(-(x-standard-alpha)/beta)), 'k-');
   set(H, 'color', functionColor, 'linewidth', 0.5);
end
H = plot(dataStruct.standard*ones(1,2), get(gca, 'ylim'), 'k--');
set(H, 'color', standardColor, 'linewidth', 2);

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
% plot joint samples
H = plot(samples.alpha(1, permSamples(1:nDraw)), samples.beta(1, permSamples(1:nDraw)), 'k+');
set(H, 'color', parameterColor, 'linewidth', 0.25, 'markersize', 4);
% marginal alpha
AXT = axes; hold on;
set(gca, 'units', 'normalized', 'position', [0.225 0.95 0.15 0.025], ...
   'xlim', alphaRange);
count = hist(samples.alpha(:), alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2));
H = bar(alphaRange(1):(alphaRange(2) - alphaRange(1))/nBins:alphaRange(2), count/ max(count), 0.7);
set(H, 'facecolor', parameterColor, 'edgecolor', parameterColor);
axis off;
% marginal beta
AXR = axes; hold on;
set(gca, 'units', 'normalized', 'position', [0.225+0.15 0.7 0.015 0.25], ...
   'ylim', betaRange);
count = hist(samples.beta(:), betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2));
H = barh(betaRange(1):(betaRange(2) - betaRange(1))/nBins:betaRange(2), count/ max(count), 0.7);
set(H, 'facecolor', parameterColor, 'edgecolor', parameterColor);
axis off;

% Print
print(['Images/' modelName '_Prior.eps'],'-depsc');
print(['Images/' modelName '_Prior.png'],'-dpng');

%% Prior Predictive figure

% Redo constants
offset = 0.5; % offset on y axis for visibility

% Use JAGS to sample
if exist(['MCMCResults/' modelName '_Prior.mat'])
    load(['MCMCResults/' modelName '_Prior.mat'], 'samples', 'stats');
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
        'dic', 0, ... 
        'monitorparams', monitorParameters,...
        'savejagsoutput', 1 , ...
        'verbosity', 1 , ...
        'cleanup' , 0, ...
        'workingdir', 'tmpjags' );
    toc % timing
    save(['MCMCResults/' modelName '_Prior.mat'], 'samples', 'stats');
end

% Setup figure
figure(2); clf; hold on;
set(gcf,'color', 'w', 'units', 'normalized', ...
   'position', [0.2 0.2 0.6 0.6], 'paperpositionmode','auto');
axis([standard - range standard + range -offset max(nStimulusTrials)+offset]);
set(gca, 'units', 'normalized', 'position', [0.15 0.15 0.7 0.7], ...
   'xtick', standard - range : 100 : standard + range, ...
   'ytick', 0 : max(nStimulusTrials), ...
   'fontsize', 12, 'box', 'off', 'tickdir', 'out');
xlabel('Stimulus', 'fontsize', 15);
ylabel('Longer Responses', 'fontsize', 15);
H = plot(min(get(gca, 'xlim'))*ones(1,2), [0 -offset], 'w-', 'linewidth', 2);

% Posterior predictive distribution of counts at each unique stimulus duration
for i = 1:nUniqueStimuli
   count = hist(samples.yPP(:, :, i), 0:nStimulusTrials(i));
   count = count/sum(count);
   for j = 0:nStimulusTrials(i)
      if count(j+1) > 0
         H = plot(stimulusDuration(i), j, 'ks');
         set(H, 'markerfacecolor', functionColor, ...
             'markeredgecolor', 'w', ...
             'linewidth', 0.5, ...
             'markersize', scale*sqrt(count(j+1)));
      end
   end
end

% Print
print(['Images/' modelName '_PriorPredictive.eps'],'-depsc');
print(['Images/' modelName '_PriorPredictive.png'],'-dpng');

