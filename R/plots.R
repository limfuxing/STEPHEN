#' Plot Step Count and Heart Rate Data coloured by STEPHEN's Physical Activity Class.
#'
#' @param object Output of STEPHEN with steps, HR and predicted_PA as columns. 
  
#' @return  A ggplot object
#' @export
plot_data <- function(object) {
  require(ggplot2)
  ggplot(object, aes(x=steps,y=HR,col=predicted_PA)) + geom_point(size=0.3) + xlab('Step Counts') + ylab('Heart Rate')
}


#' Plot proportion of different class of physical activity (PA) across time  
#'
#' @param object Output of STEPHEN with steps, HR and predicted_PA as columns. 
  
#' @return  A ggplot object
#' @export
plot_propPA <- function(object) {
  out <- aggregate(object$predicted_PA,by=list(object$date),length)
  n   <- out$x
  names(n) <- out$Group.1
  out2 <- aggregate(object$predicted_PA,by=list(Date=object$date,PA=object$predicted_PA),length)
  out2$n <- n[as.character(out2$Date)]
  out2$prop <- out2$x/out2$n
  require(ggplot2)
  ggplot(out2, aes(x=Date,y=prop,col=PA)) + geom_line()
}



