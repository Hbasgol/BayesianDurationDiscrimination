%% Logistic with various reasonable informative priors
%% To test sensitivity of the model to priors having been selected for a subject

%% Specifying parameters regarding data
clear;

% Data to be used
load vanDrielData2015 d

% Subject and condition subjectList = [7 9 12 5 6 3]
% Subjects encoded as [A B C D E F]
subject = 12;

% Condition in the experiment
condition = 'visual';

% JAGS parameters
% Sampling settings

doParallel = 1; % parallization
nThin = 1;
nChains = 3;
nSamples = 5e3;
nBurnin = 1e3;

%% Apply graphical model

% Graphical model script, which comes from a txt file.
modelName = 'LogisticSingleSubjectPriorSensitivity';

% Data for single subject and condition to supply to graphical model
switch condition
    case 'auditory', dataStruct = struct('y1', d.auditoryDecision(subject, :), ...
          'y2', d.auditoryDecision(subject, :), ... % Decision for two subsequent trials
          'y3', d.auditoryDecision(subject, :), ... 
            'stimulus',d.auditoryStimulus(subject, :), ...
            'standard', d.auditoryStandard, ...
            'nTrials', d.nTrials);
    case 'visual', dataStruct = struct('y1', d.visualDecision(subject, :), ...
          'y2', d.visualDecision(subject, :), ... % Decision for two subsequent trials
          'y3', d.visualDecision(subject, :), ...
            'stimulus',d.visualStimulus(subject, :), ...
            'standard', d.visualStandard, ...
            'nTrials', d.nTrials);
end

% Parameters to monitor and initial values
monitorParameters = {'alpha', 'beta'};
for i = 1:nChains
    S.alpha = 0*ones(1,3);
    S.beta = 10*ones(1,3);
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

%% Graphics Settings and Display

% Graphics constants
range = 400; % range of stimulus values around standard
offset = 0.05; % offset on y axis for visibility
nDraw = 200; % number of posterior samples to draw
x = dataStruct.standard - range : 1 : dataStruct.standard + range;
scale = 2; % magnifier for areas in behavioral data points
alphaRange = [-150 150]; % axis limits for parameter space
betaRange = [0 150];
nBins = 40; % number of bins for marginals
rng(5); % set random number seed so get same posterior samples
permSamples = randperm(nSamples);

%% Posterior figure
% Setup figure
figure(subject); clf; hold on;
set(gcf,'color', 'w', 'units', 'normalized', ...
   'position', [0.2 0.2 0.7/1.3 0.35/1.3], 'paperpositionmode','auto');
for i = 1:3
   subplot(1, 3, i); cla ;hold on;
set(gca, 'xtick', alphaRange, ...
    'ytick', betaRange, ...
    'fontsize', 11, 'box', 'on', 'ticklength', [0 0]);
axis([alphaRange betaRange]);
xlabel('\alpha', 'fontsize', 16,'verticalalignment', 'bottom');
ylabel('\beta', 'fontsize', 16, 'verticalalignment', 'top');
H = plot(samples.alpha(1, permSamples(1:nDraw), i), samples.beta(1, permSamples(1:nDraw), i), 'k+');
set(H, 'color', 'g', 'markersize', 4);
end

% Print
print(['Images/' modelName '_' int2str(subject) '_' condition '.eps'],'-depsc');
print(['Images/' modelName '_' int2str(subject) '_' condition '.png'],'-dpng');




