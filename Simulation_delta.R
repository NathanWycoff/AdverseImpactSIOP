#######Params You Can Modify
#Iterable Params#
#These will modify the params below when necessary.
SR <- seq(0.1,0.4,0.1)
percent <- seq(0.05,0.2,0.05)
N <- c(100,200,300,400,500,1000)
difference <- seq(0.1,1,0.1)

majMean <- 1#What is the mean of the Majorities' test scores?
majVar <- 1#What is the variance of the Majorities' test scores?
minMean <- 0#What is the mean of the Minorities' test scores?
minVar <- 1#What is the mean of the Minorities' test scores?

majCount <- 5000#How many majorities apply?
minCount <- 5000#How many minorities apply?

hireCount <- 200#How many people are hired?

alpha <- 0.05#What is the confidence level of the interval estimate?

iters <- 200#How many times do we run the cycle? Higher is slower, but gives\
#better accuracy.

########Some insight into what the sampling distribution looks like.
#This code is run automatically by the loops.

analyzeSamp <- function(AIR)	{
	AIR = AIR[1:(length(AIR)-1)]
	histdata <- hist(AIR, breaks = min(500000,1000), plot = FALSE)

	#Point estimate
	result <- histdata$mids[which.max(histdata$density)]

	#Interval Estimate
	result <- c(result,quantile(AIR,c(0.05,0.95)))

	#Other interesting statistics
	result <- c(result, mean(AIR))
	result <- c(result,median(AIR))
	result <- c(result,sd(AIR))

	#Percent with same value as upper bound
	result <- c(result, sum(AIR==result[3])/length(AIR))

	#Probability AIR < 1
	result <- c(result, sum(AIR < 1)/length(AIR))

	return(result)
}

##########These are functions which do hypothesis tests not built into R
#Z test function declaration, modified from http://www.r-bloggers.com/comparison-of-two-proportions-parametric-z-test-and-non-parametric-chi-squared-methods/
z.prop = function(x1,x2,n1,n2){
  numerator <- (x1/n1) - (x2/n2)
  p.common <- (x1+x2) / (n1+n2)
  denominator <- sqrt(p.common * (1-p.common) * (1/n1 + 1/n2))
  z.prop.ris <- numerator / denominator
  p.val <- pnorm(z.prop.ris, lower.tail = FALSE)
  return(p.val)
}

#Conduct a Zir test (a log-odds proportion comparison) as described by http://adverseimpact.org/CalculatingAdverseImpact/ZIR.htm
zir <- function(x1,x2,n1,n2, alternative = 'two.sided'){
  if (alternative == 'two.sided'){
    #Store some useful values
    SRmin <- min(c(x1/n1,x2/n2))
    SRmaj <- max(c(x1/n1,x2/n2))
    SRt <- (x1+x2)/(n1+n2)
    Pmin <- c(n1/(n1+n2),n2/(n1+n2))[which.min(c(x1/n1,x2/n2))]
    
    #Calculation of test statistic is split for readability
    numerator <- log(SRmin/SRmaj)
    denom1 <- 1-SRt
    denom2 <- SRt * (n1 + n2) * Pmin * (1-Pmin)
    Z <- numerator / sqrt(denom1/denom2)
    
    #Get a p-value from the test statistic
    p.val <- pnorm(Z)
    return(p.val)
  }
  else if (alternative == 'greater'){
   
    #
    if(x2==0)	{
        return(NA)
    }
    #Store some useful values
    SRmaj <- x1/n1
    SRmin <- x2/n2
    SRt <- (x1+x2)/(n1+n2)
    Pmin <- n1/(n1+n2)
    
    #Calculation of test statistic is split for readability
    numerator <- log(SRmin/SRmaj)
    denom1 <- 1-SRt
    denom2 <- SRt * (n1 + n2) * Pmin * (1-Pmin)
    Z <- numerator / sqrt(denom1/denom2)
    
    #Get p-value from test statistic
    p.val <- pnorm(Z, lower.tail = TRUE)
    return(p.val/2)
  }
  else if (alternative == 'less'){
      return(zir(x2,x1,n2,n1, alternative = 'greater'))
  }
  else{
    print("\'alternative\' must be in {two.sided,less,greater}")
  }
}
		

#######Actual Sampling Iterations

output <- list(c())
i <- 0
#Variables to store run parameters
scol <- 1:960
pcol <- 1:960
ncol <- 1:960
dcol <- 1:960

majCounts <- c()
minCounts <- c()
for (s in SR)	{
	for (p in percent)	{
		for (n in N)	{
			for (d in difference)	{
            i <- i + 1
            print(i)
            #if (i > 20)	{break}
            #Store configuration
            scol[i] <- s
            pcol[i] <- p
            ncol[i] <- n
            dcol[i] <- d

            #Set the values for this iteration
            majMean <- d
            #floating point arithmetic doesn't work correctly without this adjustment.
            majCount <- ceiling(n*(1-p) - 10e-5)
            minCount <- ceiling(n*p - 10e-5)
            hireCount <- ceiling(s * n - 10e-5)

			#Create populations
			majs = sapply(1:iters, function(x) rnorm(majCount, majMean, sqrt(majVar)))
			mins = sapply(1:iters, function(x) rnorm(minCount, minMean, sqrt(minVar)))
			minHired <- sapply(1:iters, function(x) sum(head(order(-c(majs[,x],mins[,x])),hireCount)>majCount)/minCount);
			majHired <- sapply(1:iters, function(x) sum(head(order(-c(majs[,x],mins[,x])),hireCount)<=majCount)/majCount);
			AIR <- minHired/majHired
			
			#Analyze posteriors
			ret <- analyzeSamp(AIR)

            #Record the number of minorities hired at the upper CI.
            new_ret <- unique(minHired[AIR == ret[3]] * minCount)

            #How many minorities would need to be hired before we don't reject?
            for (min_hired in 0:minCount) {
                maj_hired <- hireCount - min_hired
                mat <- t(matrix(c(majCount - maj_hired, maj_hired, minCount - min_hired, min_hired), nrow = 2))
                p.val <- z.prop(mat[1,2],mat[2,2],mat[1,1] +  mat[1,2], mat[2,1] +  mat[2,2])
                z.val <- qnorm(p.val, lower.tail = TRUE)
                if (z.val > -1.64) {
                    if (min_hired == 0) {
                        crit_min_hired <- "ND"
                    }
                    else {
                        crit_min_hired <- min_hired - 1
                    }
                    break
                }
            }
            new_ret <- c(new_ret, crit_min_hired)
			
			#Conduct hypothesis testing
			#4/5
			fourFifths <- sum(AIR<0.8) / length(AIR)
			ret <- c(ret, fourFifths)

			#Fisher Exact
			mat <- lapply(1:iters, function(x) t(matrix(c(majCount - majHired[x]*majCount , majHired[x]*majCount, minCount - minHired[x]*minCount,minHired[x]*minCount), nrow = 2)))
			fisherResults <- sapply(1:iters, function(x) fisher.test(mat[[x]],alternative='less')$p.value)
  			good <- !is.na(fisherResults)
 		 	
			ret <- c(ret,sum(fisherResults[good] < alpha)/sum(good))
			
			#Lancaster Adjustment to Fisher Test.
			lancaster <- sapply(1:iters, function(x) fisherResults[x] - (phyper(mat[[x]][2,2],sum(mat[[x]][2,1:2]),sum(mat[[x]][1,1:2]),hireCount))/2)
			ret <- c(ret,sum(lancaster[good] < alpha) / sum(good))
			
			#Chi-2
			chi2 <- sapply(1:iters, function(x) chisq.test(mat[[x]],simulate.p.value = TRUE)$p.value)
			ret <- c(ret,sum(chi2 < alpha) / length(chi2))

			#2 sample z-test
			zProp <- sapply(1:iters, function(x) z.prop(mat[[x]][1,2],mat[[x]][2,2],mat[[x]][1,1] +  mat[[x]][1,2], mat[[x]][2,1] +  mat[[x]][2,2]))
			ret <- c(ret,sum(zProp < alpha) / length(zProp ))

			#ZiR
			zirs <- sapply(1:iters, function(x) zir(mat[[x]][1,2],mat[[x]][2,2],mat[[x]][1,1] +  mat[[x]][1,2], mat[[x]][2,1] +  mat[[x]][2,2], alternative = 'greater'))
		  	good <- !is.na(zirs)
			ret <- c(ret,sum(zirs[good] < alpha) / sum(good ))

			#How many zirs were good?
			

			ret <- c(ret,sum(good))
            output[[i]] <- c(ret, new_ret)
			#ret})
			}
		}
	}
}

#Store results in DF
modes <- sapply(1:length(output), function(x) output[[x]][1])
lower <- sapply(1:length(output), function(x) output[[x]][2])
upper <- sapply(1:length(output), function(x) output[[x]][3])
mean <- sapply(1:length(output), function(x) output[[x]][4])
median <- sapply(1:length(output), function(x) output[[x]][5])
stddev <- sapply(1:length(output), function(x) output[[x]][6])
sameAsUpper <- sapply(1:length(output), function(x) output[[x]][7])
lessThanOne <- sapply(1:length(output), function(x) output[[x]][8])
four5 <- sapply(1:length(output), function(x) output[[x]][9])
fisher <- sapply(1:length(output), function(x) output[[x]][10])
lancaster <- sapply(1:length(output), function(x) output[[x]][11])
chi2 <- sapply(1:length(output), function(x) output[[x]][12])
zd <- sapply(1:length(output), function(x) output[[x]][13])
zirRes <- sapply(1:length(output), function(x) output[[x]][14])
zirCount <- sapply(1:length(output), function(x) output[[x]][15])
minHired <- sapply(1:length(output), function(x) output[[x]][16])
minZd <- sapply(1:length(output), function(x) output[[x]][17])


#Store the output in a dataframe
df <- data.frame(SR = scol, percent = pcol, n = ncol, difference = dcol, mode = modes, 
lowerCI = lower, upperCI = upper, mean = mean, median = median, stddev = stddev,
fourFifths = four5, fisher = fisher, lancaster = lancaster, chi2 = chi2, zd = zd,
zir = zirRes, zirCount = zirCount, sameAsUpper = sameAsUpper, lessThanOne = lessThanOne, minHired = minHired, minZd = minZd)

#If you want a file with your results, uncomment the next line (delete the "#")
#write.csv(file='sim_delta_output.csv', x = df)
