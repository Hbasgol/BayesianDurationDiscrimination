%% For finding the posterior predictive distribution for a subject with informative priors

clear;

%% Specifying parameters regarding data

% Data to be used
load vanDrielData2015 d

% Subject and condition subjectList = [7 9 12 5 6 3]
% Subjects encoded as [A B C D E F]
subject = 12; % Subject C

% Condition in the experiment
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
modelName = 'LogisticSingleSubjectInformativePosteriorPredictive';

% Data for single subject and condition to supply to graphical model
switch condition
   case 'auditory',
      stimulusDuration = unique(d.auditoryStimulus(subject, :));
      nUniqueStimuli = length(stimulusDuration);
      nStimulusTrials = zeros(nUniqueStimuli, 1);
      for i = 1:nUniqueStimuli
         nStimulusTrials(i) = length(find(d.auditoryStimulus(subject, :) == stimulusDuration(i)));
      end
      dataStruct = struct('y', d.auditoryDecision(subject, :), ...
         'stimulus',d.auditoryStimulus(subject, :), ...
         'standard', d.auditoryStandard, ...
         'nTrials', d.nTrials, ...
         'nUniqueStimuli', nUniqueStimuli, ...
         'stimulusDuration', stimulusDuration, ...
         'nStimulusTrials', nStimulusTrials);
   case 'visual',
      stimulusDuration = unique(d.visualStimulus(subject, :));
      nUniqueStimuli = length(stimulusDuration);
      nStimulusTrials = zeros(nUniqueStimuli, 1);
      for i = 1:nUniqueStimuli
         nStimulusTrials(i) = length(find(d.visualStimulus(subject, :) == stimulusDuration(i)));
      end;dataStruct = struct('y', d.visualDecision(subject, :), ...
         'stimulus',d.visualStimulus(subject, :), ...
         'standard', d.visualStandard, ...
         'nTrials', d.nTrials, ...
         'nUniqueStimuli', nUniqueStimuli, ...
         'stimulusDuration', stimulusDuration, ...
         'nStimulusTrials', nStimulusTrials);
end

% Parameters to monitor and initial values
monitorParameters = {'alpha', 'beta', 'yPP'};
for i = 1:nChains
   S.alpha = 0;
   S.beta = 10;
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
standard = 500;
range = 400; % range of stimulus values around standard
offset = 0.5; % offset on y axis for visibility
scale = 20; % magnifier for posterior predictive

%% Posterior predictive figure
% Setup figure
figure(1); clf; hold on;
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

% Subject behavior
yObserved = zeros(nUniqueStimuli, 1);
for i = 1:nUniqueStimuli
   switch condition
      case 'visual'
         match = find(d.visualStimulus(subject, :) == stimulusDuration(i));
         above = length(find(d.visualDecision(subject, match) == 1));
      case 'auditory'
         match = find(d.auditoryStimulus(subject, :) == stimulusDuration(i));
         above = length(find(d.auditoryDecision(subject, match) == 1));
   end
   yObserved(i) = above;
end

% Posterior predictive distribution of counts at each unique stimulus duration
% whole distribution
for i = 1:nUniqueStimuli
   count = hist(samples.yPP(:, :, i), 0:nStimulusTrials(i));
   count = count/sum(count);
   for j = 0:nStimulusTrials(i)
      if count(j+1) > 0
         H = plot(stimulusDuration(i), j, 'ks');
         set(H, 'markerfacecolor', 'b', ...
             'markeredgecolor', 'w', ...
             'linewidth', 0.4, ...
             'markersize', scale*sqrt(count(j+1)));
         % check whether matches observed
         if yObserved(i) == j
            set(H, 'markeredgecolor', 'k');
         end
      end
   end
end

% just where data are observed
for i = 1:nUniqueStimuli
   count = hist(samples.yPP(:, :, i), 0:nStimulusTrials(i));
   count = count/sum(count);
   for j = 0:nStimulusTrials(i)
      if count(j+1) > 0
                  if yObserved(i) == j
H = plot(stimulusDuration(i), j, 'ks');
         set(H, 'markerfacecolor', 'b', ...
             'markeredgecolor', 'k', ...
             'linewidth', 1.5, ...
             'markersize', scale*sqrt(count(j+1)));
         end
      end
   end
end

% Overlay line
H = plot(stimulusDuration, yObserved, 'k-');
set(H, 'color', 'b');

% Print
print(['Images/' modelName '_' int2str(subject) '_' condition '_PosteriorPredictive.eps'],'-depsc');
print(['Images/' modelName '_' int2str(subject) '_' condition '_PosteriorPredictive.png'],'-dpng');



