#' This function takes Fitbit's Step Count and Heart Rate (HR) data to estimate physical activity levels (Sedentary, LowPA, ModeratePA and VigorousPA) using an ensemble of Hidden semi-Markov models. Users need to supply the Step Count and Heart Rate Data and has the option to also supply the pre-trained set of models that will be used to predict physical activity. We recommend that this set of models have been trained on individuals with similar lifestyles to that we will predict. If not supplied, the default set based on 11 individuals from OPTIMISE trial (Brakenridge et al., 2022, BMC Public Health) will be used.
#'
#' @param steps.data A Data Frame containing step count data at 1-min resolution. Run the examples to understand the required structure.
#' @param HR.data A Data Frame containing heart rate data at 1-min resolution. Run the examples to understand the required structure.
#' @param preTrainedSet  A list with each element containing a pre-trained Hidden semi-Markov model to predict physical activity using step count and heart rate data.
  
#' @return  A Data Frame containing the date, step count, heart rate and predicted physical activity levels (1=Sedentary, 2= LowPA, 3=ModeratePA, 4=VigorousPA).
#' @export

STEPHEN <- function (steps.data, HR.data, preTrainedSet = NULL) {
    if (nrow(steps.data) != nrow(HR.data)) 
        stop("Step and HR data do not have the same number of rows! Please fix.")
    if (any(steps.data[, 1] != HR.data[, 1])) 
        stop("Timestamp in steps. data and HR.data do not match. We need to have matched data")
    if (is.null(preTrainedSet)) {
        data(optimise1019)
        data(optimise1022)
        data(optimise1050)
        data(optimise1061)
        data(optimise1083)
        data(optimise1088)
        data(optimise1104)
        data(optimise1112)
        data(optimise1161)
        data(optimise1162)
        data(optimise1167)
        preTrainedSet <- list(optimise1019,optimise1022,optimise1050,optimise1061,
				optimise1083,optimise1088,optimise1104,optimise1112,
				optimise1161,optimise1162,optimise1167)
        names(preTrainedSet) <- paste0('subject',1:length(preTrainedSet))
    }

    Trained.subjects <- names(preTrainedSet)
    if (!is.null(names(preTrainedSet))) 
        Trained.subjects <- 1:length(preTrainedSet)
    date <- as.Date(steps.data[, 1], format = "%m/%d/%Y")
    time = lubridate::mdy_hms(HR.data[, 1])
    # check if assumed time format is correct
    if(! as.numeric(diff(time[1:2],unit='secs'))==60) 
      time = lubridate::mdy_hm(HR.data[, 1])
    hour = lubridate::hour(time)
    mins = lubridate::minute(time)
    secs = lubridate::second(time)

    if(! min(hour) == 0 & max(hour)==23 ) 
     print('Warnings: please check your timestamp date variable. It must in MDY HMS or MDY HM format')
    if(! min(mins) == 0 & max(mins)==59 ) 
     print('Warnings: please check your timestamp date variable. It must in MDY HMS or MDY HM format')

    time.24hr <- hour + mins/60 + secs/3600
    date.cont <- (as.numeric(date) - min(as.numeric(date))) *24 + time.24hr
    rownames(steps.data) <- paste0(date, "-", time.24hr)
    rownames(HR.data) <- paste0(date, "-", time.24hr)
    HR.data[hour >= 0 & hour < 6, 2] = 0
    non.usage = HR.data[, 2] == 0
    chpts = which(as.logical(abs(diff(non.usage)))) + 1
    idx = data.frame(st = chpts[-length(chpts)], en = chpts[-1] - 1)
    # if first segment does not start at 1
    if(idx$st[1]>1) 
      idx <- rbind(data.frame(st=1,en=idx$st[1]-1),idx)
    # if last segment does not end on the last row of data
    if(idx$en[nrow(idx)]<nrow(HR.data)) 
      idx <- rbind(idx,data.frame(st=idx$en[nrow(idx)]+1,en=nrow(HR.data)))

    idx$len = idx$en - idx$st + 1
    getSubset <- function(idx, x) {
        x[idx[1]:idx[2]]
    }
    # when the first segment is usage
    if(HR.data[idx$st[1],2]!=0) 
      idx.use <- idx[seq(1, nrow(idx), 2), ]
    # when the first segment is non-usage
    if(HR.data[idx$st[1],2]==0) 
      idx.use <- idx[seq(2, nrow(idx), 2), ]

    idx.use <- idx.use[idx.use$len >= 5, ]
    steps.use <- unlist(apply(idx.use, 1, getSubset, x = steps.data[,2]))
    hr.use <- unlist(apply(idx.use, 1, getSubset, x = HR.data[,2]))
    cum.idx <- NULL
    for (i in 1:nrow(idx.use)) 
	cum.idx <- c(cum.idx, idx.use[i,1]:idx.use[i,2])

    date.use <- date[cum.idx]
    time.use <- time[cum.idx]

    data <- list(x = data.frame(steps.use, hr.use), N = idx.use$len)
    class(data) <- "hsmm.data"
    class.HSMM1 <- NULL
    for (set in Trained.subjects) {
        shmm.int <- rep(NA, length(hr.use))
        out <- preTrainedSet[[set]]
        base.HR <- mean(hr.use[steps.use == 0])
        bias.HR <- base.HR - min(out$model$parms.emission$mu2[out$model$parms.emission$mu1 < 
            1])
        out$model$parms.emission$mu2 <- out$model$parms.emission$mu2 + 
            bias.HR
        out$model$parms.emission$size1 <- ifelse(out$model$parms.emission$size1 < 
            0.2, 0.2, out$model$parms.emission$size1)

        # ensure class 1-4 have increasing mean HR
        HR.order <- order(out$model$parms.emission$mu2)
        out$model$parms.emission$mu2 <- out$model$parms.emission$mu2[HR.order]
        out$model$parms.emission$mu1 <- out$model$parms.emission$mu1[HR.order]

        smoNB5.class <- predict(out, newdata = data, trace = FALSE)
        print(paste0("Finished predicting using ", match(set, 
            Trained.subjects), " out of ", length(Trained.subjects), 
            " preTrained models"))
        shmm.int <- smoNB5.class$s
        class.HSMM1 <- cbind(class.HSMM1, shmm.int)
    }
    getmode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    class.HSMM1 <- apply(class.HSMM1, 1, getmode)
    PA_class <- c("SED", "LPA", "MPA", "VPA")
    output <- data.frame(date=date.use,time=time.use,steps = steps.use, HR = hr.use, predicted_PA = factor(PA_class[class.HSMM1], levels=PA_class))
    output
}

