#######Params You Can Modify
#Iterable Params#
#These will modify the params below when necessary.
SR <- seq(0.1,0.4,0.1)
percent <- seq(0.05,0.2,0.05)
N <- c(100,200,300,400,500,1000)

majMean <- 1#What is the mean of the Majorities' test scores?
majVar <- 1#What is the variance of the Majorities' test scores?
minMean <- 0#What is the mean of the Minorities' test scores?
minVar <- 1#What is the mean of the Minorities' test scores?

majCount <- 5000#How many majorities apply?
minCount <- 5000#How many minorities apply?

hireCount <- 200#How many people are hired?

alpha <- 0.05#What is the confidence level of the interval estimate?

iters <- 20000#How many times do we run the cycle? Higher is slower, but gives\
#better accuracy.


#Note: When s = 0.4, n = 100, p = 0.2, for example, we get a case where our discrete CI has to cover almost 60% of the mass, even though we only want a 40% interval.

########Some insight into what the sampling dist looks like.
#This code is run automatically by the loops.

analyzeSamp <- function(AIR)	{
	AIR = AIR[1:(length(AIR)-1)]
	histdata <- hist(AIR, breaks = min(500000,1000), plot = FALSE)

    result <- c()

	#Interval Estimate with 0.6 CI
	result <- c(result,quantile(AIR,c(0.2, 0.8)))

	#Percent with same value as LOWER (as opposed to upper) bound
	result <- c(result, sum(AIR==result[1])/length(AIR))

    #Percent below the lower bound.
    result <- c(result, sum(AIR < result[1]) / length(AIR))

	return(result)
}


#######Multicore iterations. This code block iterates over the "Iterable Params"
#splitting the work accross multiple cores.
output <- list(c())
i <- 0
#Variables to store run parameters
run_length <- length(SR) * length(percent) * length(N)
scol <- rep(0, run_length)
pcol <- rep(0, run_length)
ncol <- rep(0, run_length)
dcol <- rep(0, run_length)

d <- 0#Difference will always be zero in this sim

majCounts <- c()
minCounts <- c()
for (s in SR)	{
	for (p in percent)	{
		for (n in N)	{
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

            output[[i]] <- ret
			#ret})
		}
	}
}

#output <- list(c())
#for (i in 1:length(fs))	{
#	output[[i]] <- value(fs[[i]])
#}

#Store results in DF
lower <- sapply(1:length(output), function(x) output[[x]][1])
upper <- sapply(1:length(output), function(x) output[[x]][2])
sameAsLower <- sapply(1:length(output), function(x) output[[x]][3])
belowLower <- sapply(1:length(output), function(x) output[[x]][4])


#Store the output in a dataframe
df <- data.frame(SR = scol, percent = pcol, n = ncol, difference = dcol, lowerCI = lower, upperCI = upper, sameAsLower = sameAsLower, belowLower = belowLower)

#If you want a file with your results, uncomment the next line (delete the "#")
#write.csv(file='sim_no_delta_output.csv', x = df)
