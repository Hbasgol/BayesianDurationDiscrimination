% Inverse logistic psychophysical function
% x = LogisticPsychophysicalFunctionInverse(theta, standard, alpha, beta)
function x = LogisticPsychophysicalFunctionInverse(theta, standard, alpha, beta)
x = standard+alpha-beta*log((1-theta)./theta);