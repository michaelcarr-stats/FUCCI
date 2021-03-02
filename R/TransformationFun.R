TransformationFun <- function(theta,b,inverse) {

#This function transforms parameters using the logit function and specifies
#the upper bound b and the direction of the transformation (inverse or not)
  
  if (inverse) {
    theta = b/(1+exp(-theta));
  } else {
    theta = log(theta/(b-theta));
  }
  return(theta)
}