%% Hierarchical modeling for individual differences
%% Finding parameters for each subject and predicting behaviors of a new subject
%% Specifying parameters regarding data

clear;

% Data to be used
load vanDrielData2015 d

% Subject and condition subjectList = [7 9 12 5 6 3]
% Subjects encoded as [A B C D E F]
subjectList = [7 9 12 5 6 3];

% One of the condition in the experiment
condition = 'visual';

drawBothConditions = 0;

% JAGS parameters
% Sampling settings
doParallel = 1; % parallization
nThin = 1;
nChains = 3;
nSamples = 5e3;
nBurnin = 1e3;

%% Apply graphical model

% Graphical model script
modelName = 'LogisticHierarchicalInformative';

% Data for single subject and condition to supply to graphical model
switch condition
    case 'auditory', dataStruct = struct('y', d.auditoryDecision(subjectList, :), ...
            'stimulus', d.auditoryStimulus(subjectList, :), ...
            'standard', d.auditoryStandard, ...
            'nTrials', d.nTrials,...
            'nSubjects', length(subjectList));
    case 'visual', dataStruct = struct('y', d.visualDecision(subjectList, :), ...
            'stimulus', d.visualStimulus(subjectList, :), ...
            'standard', d.visualStandard, ...
            'nTrials', d.nTrials, ...
            'nSubjects', length(subjectList));
end

% Parameters to monitor and initial values
monitorParameters = {'alpha', 'beta', 'alphaPredicted', 'betaPredicted', ...
    'mualpha', 'mubeta', 'sigmaalpha', 'sigmabeta'};
for i = 1:nChains
    S.mualpha = 0;
    S.mubeta = 10;
    S.sigmaalpha = 1;
    S.sigmabeta = 1;
    init0(i) = S;
end

% Use JAGS to sample
if exist(['MCMCResults/' modelName '_' condition '.mat'])
    load(['MCMCResults/' modelName '_' condition '.mat'], 'samples', 'stats');
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
    save(['MCMCResults/' modelName '_' condition '.mat'], 'samples', 'stats');
end
disp(stats.Rhat); % convergence check

%% Display
% Graphics constants
%load HandbookColors
standard = d.auditoryStandard;
range = 400; % range of stimulus values around standard
offset = 0.1; % offset on y axis for visibility
nDraw = 30; % number of posterior samples to draw
x = standard - range : 1 : standard + range;
scale = 1.5; % magnifier for areas in behavioral data points
markerStyle = {'o','^','+','x','d','v'};
alphaRange = [-150 150]; % axis limits for parameter space
betaRange = [0 150];
nBins = 40; % number of bins for marginals
rng(3); % set random number seed so get same posterior samples
permSamples = randperm(nSamples);


%% Posterior, predicted subject, JND (JUST NOTICABLE DIFFERENCE)
figure(2); clf; hold on;
set(gcf,'color', 'w', 'units', 'normalized', ...
    'position', [0.2 0.1 0.5/1.05 0.8/1.05], 'paperpositionmode','auto');

% Individual-level parameter space
subplot(2, 2, 1); hold on;
set(gca, 'xtick', alphaRange, ...
   'ytick', betaRange, ...
   'fontsize', 12, 'box', 'on', 'ticklength', [0 0]);
axis([alphaRange betaRange]);
xlabel('\mu_\alpha', 'fontsize', 16,'verticalalignment', 'bottom');
ylabel('\mu_\beta', 'fontsize', 16, 'verticalalignment', 'top');
% plot joint samples
H = plot(samples.alpha(1, permSamples(1:nDraw)), samples.beta(1, permSamples(1:nDraw)), 'k+');
set(H, 'color', 'r', 'linewidth', 0.25, 'markersize', 4);

% draw samples for auditory condition if flag set
if drawBothConditions
    load(['MCMCResults/' modelName '_auditory.mat'], 'samples', 'stats');
H = plot(samples.alpha(1, permSamples(1:nDraw)), samples.beta(1, permSamples(1:nDraw)), 'k+');
set(H, 'color', 'r', 'linewidth', 0.25, 'markersize', 4);
    load(['MCMCResults/' modelName '_visual.mat'], 'samples', 'stats');
end

% Individual-level parameter space
subplot(2, 2, 2); hold on;
axis([standard - range standard + range -offset 1]);
set(gca, 'xtick', alphaRange, 'ytick', betaRange, ...
    'fontsize', 12, 'box', 'on', 'ticklength', [0 0]);
axis([alphaRange betaRange]);
xlabel('\alpha', 'fontsize', 18, 'verticalalignment', 'bottom');
ylabel('\beta', 'fontsize', 18, 'verticalalignment', 'top');

% observed subjects
for subjectIndex = 1:length(subjectList)
H(1) = plot([stats.ci_low.alpha(subjectIndex) stats.ci_high.alpha(subjectIndex)], ones(1,2)*stats.mean.beta(subjectIndex), 'k-');
H(2) = plot(ones(1,2)*stats.mean.alpha(subjectIndex), [stats.ci_low.beta(subjectIndex) stats.ci_high.beta(subjectIndex)], 'k-');
HC = plot(stats.mean.alpha(subjectIndex), stats.mean.beta(subjectIndex), 'ko');
set(HC, 'markersize', 10, 'markerfacecolor', 'w', 'markeredgecolor', 'w', 'linewidth', 0.5);
T = text(stats.mean.alpha(subjectIndex), stats.mean.beta(subjectIndex), char(64+subjectIndex));
set(T, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontsize', 10);
set(H(:), 'color', 'r', 'linewidth', 1);
end
% predicted subject
H(1) = plot([stats.ci_low.alphaPredicted stats.ci_high.alphaPredicted], ones(1,2)*stats.mean.betaPredicted, 'k-');
H(2) = plot(ones(1,2)*stats.mean.alphaPredicted, [stats.ci_low.betaPredicted stats.ci_high.betaPredicted], 'k-');
HC = plot(stats.mean.alphaPredicted, stats.mean.betaPredicted, 'ko');
set(HC, 'markersize', 10, 'markerfacecolor', 'w', 'markeredgecolor', 'w', 'linewidth', 0.5);
T = text(stats.mean.alphaPredicted, stats.mean.betaPredicted, 'N');
set(T, 'fontweight', 'bold', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontsize', 11);
set(H(:), 'color','r', 'linewidth', 1);

% Posterior psychophysical function
subplot(2, 2, 3); hold on;
axis([standard - range standard + range -offset 1]);
set(gca, 'xtick', standard - range : 400 : standard + range, ...
    'ytick', 0 :  1, ...
    'fontsize', 12, 'box', 'off', 'tickdir', 'out');
T = xlabel('Stimulus','fontsize', 14);
T = ylabel({'Longer Response','Probability'},'fontsize', 14);
% white space to separate axes
H = plot(min(get(gca, 'xlim'))*ones(1,2), [0 -offset], 'w-', 'linewidth', 2);
% draw standard
H = plot(standard*ones(1, 2), get(gca, 'ylim'), 'k--');
set(H, 'color', 'g', 'linewidth', 2);

% Draw posterior psychometric function and standard
for i = 1:nDraw
    alpha = samples.alphaPredicted(1, permSamples(i));
    beta = samples.betaPredicted(1, permSamples(i));
    H = plot(x, 1./(1 + exp(-(x - standard - alpha)/beta)), 'k-');
    set(H, 'color', 'b', 'linewidth', 0.5);
end
H = plot(standard*ones(1,2), get(gca, 'ylim'), 'k--');
set(H, 'color', 'g', 'linewidth', 2);

% JND
subplot(2, 2, 4); hold on;
lo = 0; hi = 220; step = 3; binsc = lo+step/2:step:hi-step/2; binse = lo:step:hi;
set(gca, 'xtick', lo:100:hi, 'xlim', [lo hi], ...
    'ycolor', 'none', ...
    'fontsize', 12, 'box', 'off', 'tickdir', 'out');
T = xlabel('JND','fontsize', 14);

% Construct JNDs
for subjectIndex = 1:length(subjectList)
        allAlpha = reshape(samples.alpha(:,:,subjectIndex), 1, []);
        allBeta = reshape(samples.beta(:,:,subjectIndex), 1, []);
        JND = LogisticPsychophysicalFunctionInverse(0.84, standard, allAlpha(:), allBeta(:)) - LogisticPsychophysicalFunctionInverse(0.5, standard, allAlpha(:), allBeta(:));
        count = histc(JND, binse);
        count = count(1:end-1);
        count = count/sum(count);
        H = plot(binsc, count,'k-');
        set(H, 'color', 'b');
        [val, ind] = max(count);
        T = text(binsc(ind), val, char(64+subjectIndex));
        set(T, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
    %end;
end
% predicted subject
JND = LogisticPsychophysicalFunctionInverse(0.84, standard, samples.alphaPredicted(:),  samples.betaPredicted(:)) - LogisticPsychophysicalFunctionInverse(0.5, standard,  samples.alphaPredicted(:),  samples.betaPredicted(:));
count = histc(JND,binse);
count = count(1:end-1);
count = count/sum(count);
H = plot(binsc, count, 'k--');
set(H, 'color', 'b', 'linewidth', 3);

% % Print
print(['Images/' modelName '_Posterior.eps'],'-depsc');
print(['images/' modelName '_Posterior.png'],'-dpng');



