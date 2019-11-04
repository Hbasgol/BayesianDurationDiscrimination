%% Logistic with latent mixture
%% For finding the function that explains better observed behaviors

clear;

% Data to be used
load vanDrielData2015 d

% Subject and condition subjectList = [7 9 12 5 6 3]
% Subjects encoded as [A B C D E F]
subjectList = [7 9 12 5 6 3];
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
modelName = 'LogisticCauchyLatentMixture';

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
monitorParameters = {'alpha', 'beta', 'z', 'phi'};
for i = 1:nChains
    S.alpha = zeros(length(subjectList), 2);
    S.beta = 10*ones(length(subjectList), 2);
    S.z = round(rand(length(subjectList),1));
    S.phi = 0.5;
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
eps = 0.01; binsc = eps/2:eps:1-eps/2; binse = 0:eps:1;
standard = d.auditoryStandard;
subjectStr = int2str(subjectList(1));
for i = 2:length(subjectList)
    subjectStr = char(subjectStr, int2str(subjectList(i)));
end

% Setup figure
figure(1); clf; hold on;
set(gcf,'color', 'w', 'units', 'normalized', ...
    'position', [0.2 0.2 0.6 0.35], 'paperpositionmode','auto');

% Posterior for model indicators
subplot(1, 2, 1); cla; hold on;
H = bar(1:6, stats.mean.z, 0.6);
set(H, 'facecolor', 'g', 'edgecolor', 'g');
xlabel('Subject', 'fontsize',16);
ylabel('Probability Logistic', 'fontsize',16);
set(gca, 'xlim', [0.5 6.5], 'xtick', 1:6, 'xticklabel', {'A','B','C','D','E','F'}, ...
    'ytick', 0:.2:1, 'fontsize', 12, 'tickdir', 'out', 'ticklength', [0.02 0]);

% Posterior for base-rate
subplot(1, 2, 2); cla; hold on;
count = histc(samples.phi(:), binse);
count = count(1:end-1);
H = bar(binsc, count/sum(count(:)), 1);
set(H, 'facecolor', 'g', 'edgecolor', 'g');
xlabel('Base-rate Logistic');
set(gca, 'xlim', [0 1], 'xtick', 0:.2:1, 'ycolor', 'w', 'fontsize', 12, ...
    'tickdir', 'out', 'ticklength', [0.02 0]);

%% Print
print(['Images/' modelName '_' condition '.eps'],'-depsc');
print(['Images/' modelName '_' condition '.png'],'-dpng');
