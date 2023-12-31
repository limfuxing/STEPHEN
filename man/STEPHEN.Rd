% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STEPHEN.R
\name{STEPHEN}
\alias{STEPHEN}
\title{This function takes Fitbit's Step Count and Heart Rate (HR) data to estimate physical activity levels (Sedentary, LowPA, ModeratePA and VigorousPA) using an ensemble of Hidden semi-Markov models. Users need to supply the Step Count and Heart Rate Data and has the option to also supply the pre-trained set of models that will be used to predict physical activity. We recommend that this set of models have been trained on individuals with similar lifestyles to that we will predict. If not supplied, the default set based on 11 individuals from OPTIMISE trial (Brakenridge et al., 2022, BMC Public Health) will be used.}
\usage{
STEPHEN(steps.data, HR.data, preTrainedSet = NULL)
}
\arguments{
\item{steps.data}{A Data Frame containing step count data at 1-min resolution. Run the examples to understand the required structure.}

\item{HR.data}{A Data Frame containing heart rate data at 1-min resolution. Run the examples to understand the required structure.}

\item{preTrainedSet}{A list with each element containing a pre-trained Hidden semi-Markov model to predict physical activity using step count and heart rate data.}
}
\value{
A Data Frame containing the date, step count, heart rate and predicted physical activity levels (1=Sedentary, 2= LowPA, 3=ModeratePA, 4=VigorousPA).
}
\description{
This function takes Fitbit's Step Count and Heart Rate (HR) data to estimate physical activity levels (Sedentary, LowPA, ModeratePA and VigorousPA) using an ensemble of Hidden semi-Markov models. Users need to supply the Step Count and Heart Rate Data and has the option to also supply the pre-trained set of models that will be used to predict physical activity. We recommend that this set of models have been trained on individuals with similar lifestyles to that we will predict. If not supplied, the default set based on 11 individuals from OPTIMISE trial (Brakenridge et al., 2022, BMC Public Health) will be used.
}

\examples{
# load step count data
data(steps.ex)
# load heart rate data
data(hr.ex)
res <- STEPHEN(steps.data=steps.ex,HR.data=hr.ex)
head(res)
# plot data
plot_data(res)
# plot proportion of PA
plot_propPA(res)
}
