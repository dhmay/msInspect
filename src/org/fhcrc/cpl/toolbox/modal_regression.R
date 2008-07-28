##################################################################
#Modal regression for Damon's AMT project 01/03/06               #
#Function modal_regress runs quantile regression a few           #
#times till the mode of residual converge to zero                #
#function inputs: time and hydrophbicity                         #
#function output: a vector contains two components:              #
#                  regression intercept and slope                #
#                                                                #
# NOTE: this code depends on the 'quantreg' library, which must  #
# be loaded prior to running                                     #
#                                                                #
##################################################################

###############
#the function#
###############
 
modal_regress<-function(time, H)
{
  t<-0.5
  for (i in 1:20) 
  {
    ftrq<-rq(H~time, tau=t)
    dr<-density(ftrq$residual)
    drx<-dr$x[dr$y==max(dr$y)]
    
    if (drx> 0.01 )
       t <- t+0.01
    else if (drx< -0.01 ) 
       t <- t-0.01
    else break 
  }

  if (0)
  {
    plot(time, H)
    abline(ftrq$coef, col=3) 
    plot(density(ftrq$residual), main="Density of residual")
    abline(v=0, col=3)
  }

  return(ftrq$coef)

}

