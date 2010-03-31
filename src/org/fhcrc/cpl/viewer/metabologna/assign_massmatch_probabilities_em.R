#test Marty's program using simulated data##

#proportion = initial proportion estimate
#numiterations = number of iterations
#showcharts=TRUE
#max_deltap_for_stable <- .05;
#iters_stable_for_converg <- 5;

emtest <- function(x,initial_proportion, miniterations, maxiterations, showcharts,
                   max_deltap_proportion_for_stable, iters_stable_for_converg, area)
{
  sigx<-sd(x)
  mux<-mean(x)
  proportion <- initial_proportion

  iteration_proportions=proportion;
  iteration_muxs=mux;
  iteration_msigxs=sigx;
  iteration_max_prob_diff=0;

  converged=FALSE;

  z=rep(0,length(x));

  num_iterations_stable = 0;

  ### e step
  for(i in 1:maxiterations)
  {

    proportion_times_normdist <- proportion *  dnorm(x,mean=mux,sd=sigx);
    znew <- proportion_times_normdist / (proportion_times_normdist + ((1-proportion) / area) )


    #keep track of the /biggest percent difference/ between probabilities from one iteration to next.
    #have to be careful about division by 0.  We need a dummy value for the first iteration, and it can't
    #be 0 because then if only 1 iteration is required, we'll consider ourselves converged
    if (i == 1)
    {
        max_deltap_proportion <- max_deltap_proportion_for_stable + 0.001;
    }
    else
    {
        max_deltap_proportion <- max(abs(znew-z) / max(z,.000001));
    }

    iteration_max_prob_diff[i] <- max_deltap_proportion;

    z <- znew;

    if (max_deltap_proportion <= max_deltap_proportion_for_stable)
    {
        num_iterations_stable = num_iterations_stable + 1;
    }
    else
    {
        num_iterations_stable = 0;
    }

    if (num_iterations_stable >= iters_stable_for_converg && i >= miniterations)

    {
        converged=TRUE;
        break;
    }

### mstep
    mux <- sum(z*x)/sum(z)
    muxx <- sum(z*x*x)/sum(z)
    sigx <- sqrt(muxx-mux*mux)

    proportion <- mean(z)

#print( c(mux,muy,sigx,sigy,proportion) )
    iteration_proportions <- c(iteration_proportions, proportion);
    iteration_muxs <- c(iteration_muxs, mux);
    iteration_msigxs <- c(iteration_msigxs, sigx);
  }

  if (showcharts)
  {
    jpeg(chart_out_file, 800, 800);

    par(mfrow=c(2,2));
    plot(iteration_muxs, pch=20, xlab="Iteration", ylab="Mu(x)")
    plot(iteration_msigxs, pch=20, xlab="Iteration", ylab="Sigma(x)")
    plot(iteration_proportions, pch=20, xlab="Iteration", ylab="Proportion True")
    plot(iteration_max_prob_diff, pch=20, xlab="Iteration", ylab="Max % Probability Change")

    dev.off();
  }

  num_iterations_performed = length(iteration_muxs);

  return(list(z=z, para=c(mux,sigx,proportion), num_iterations=num_iterations_performed,
         converged=converged));
}


#print(proportion);
est<-emtest(x=targetx, initial_proportion=proportion,
            miniterations=miniterations, maxiterations=maxiterations,
            showcharts=showcharts,
            max_deltap_proportion_for_stable=max_deltap_proportion_for_stable,
            iters_stable_for_converg=iters_stable_for_converg, area=area)

#test normality#
set.seed(100)
num_sample_points=10000;
estdist.x<-c(rnorm(est$para[3]*num_sample_points,est$para[1], est$para[2]),
             runif((1-est$para[3])*num_sample_points, min(targetx), max(targetx)))

ksx <- ks.test(targetx, estdist.x)$p.value

qtargetx<-quantile(targetx, probs=seq(0,1,0.001))
qestx<-quantile(estdist.x, probs=seq(0,1,0.001))

corx<-cor(qtargetx,qestx)

qb1x<-summary(lm(qtargetx ~ (-1) + qestx))$coef[1,1]

#Write out charts indicating how well our distribution matches the theoretical distribution
if (showcharts)
{
    jpeg(normtest_chart_out_file, 800, 800);

    par(mfrow=c(2,1))

    #x dimension
    plot(density(targetx), type="l", lty=1)
    points(density(estdist.x), type="l", lty=2, col=2)

#    plot(density(rnorm(est$para[3]*num_sample_points,est$para[1], est$para[2])), type="l", lty=1);
#    points(density((1-est$para[3])*num_sample_points, min(targetx), max(targetx)), type="l",lty=2, col=2);

    qqplot(targetx, estdist.x, main="QQplot of estimated mixed x dist", pch=20, cex=0.5)
    abline(0,1, col=2)
 
    dev.off();
}

list(probs=est$z, ksresults=ksx, qbeta=qb1x, corresults=corx, converged=est$converged, num_iterations=est$num_iterations, dist_params=est$para)


