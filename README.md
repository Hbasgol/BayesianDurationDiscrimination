# Hierarchical Bayesian Modeling for Duration Discrimination

This project is conducted for Bayesian Statistics and Machine Learning course given by Ali Taylan Cemgil in Boğaziçi University and re-implemented as a group project.

The reference paper of this project is "Bayesian Methods in Cognitive Modeling" written by Michael Lee, who is from University of California Irvine. The code and data of Michael Lee can be found at https://osf.io/zur8m/

In this project, MATJAGS (MATLAB JAGS interface) is used. As an inference method, Gibbs sampling algorithm in MCMC family is used in JAGS. To select the better model, Savage-Dickey Method is implemented, which is often used in Cognitive Psychology literature.

Normally, the project is implemented by Michael Lee has more models than that we have implemented. The models that we have implemented are for:

- Assessing informative prior and vague prior for likelihood function
- Finding sequential effects between trials
- Finding contaminant trials
- Selecting likelihood function (with latent-mixture method)
- Assessing individual differences (with hierarchical bayesian modeling) 
- Predicting behaviors of subjects (with a common cause model)

Although we used the same dataset with Michael Lee, we considered different subjects. Additionally, Michael Lee did not compute models that we have mentioned above for each modality and for six different subjects. Moreover, the code presented by Michael Lee made convenient to run again by deciding missing parameters in the function and modules for graphical representations, we have solved these problems.

Reference paper: Lee, M. D. (2018). Bayesian Methods in Cognitive Modeling. Stevens Handbook of Experimental Psychology and Cognitive Neuroscience, 1-48. doi:10.1002/9781119170174.epcn502
