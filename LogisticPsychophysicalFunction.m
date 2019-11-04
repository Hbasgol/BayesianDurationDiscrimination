% Logistic psychophysical function
% theta = LogisticPsychophysicalFunction(x, standard, alpha, beta)

function theta = LogisticPsychophysicalFunction(x, standard, alpha, beta)

theta = 1./(1 + exp(-(x-standard-alpha)/beta));
