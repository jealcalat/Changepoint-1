##################################################################
##	Bayesian change point detection	for poisson data			##
##	baed on Chib (1997) by Lei Gong				  		  		##
##	fit model with pre-assigned number of change points	  		##
##	use Bayes factor/marginal liklelihood to compare models		##	
##	relative fixed-width rule as stopping criteria				## 
##	input arguments:											##
##		m:	number of change points								##
##		y:	original dataset									##
##	output files:												##
##		P.txt, S.txt, theta.txt, cond.txt						##
##################################################################


library(mcmcse)


## input m ##
args <- commandArgs(TRUE)
m<- as.numeric(args[1])


## input y ##
raw<- read.table("change-point-2.txt")
y<- raw[,1]
n<- length(y)	# sample size


## functions ##
# equation 6
eq6<- function(t, k, P, mass){
	if(k == 0){
		output<- 0
	}
	else if(k == 1){
		output<- P[1]*mass[t-1, 1]
	}
	else{
		output<- (1-P[k-1])*mass[t-1, k-1] + P[k]*mass[t-1, k]
	}
	return(output)
}


## initial setups ##
# hyper parameters
a<- round(0.1*n/(m+1)-0.1)
b<- 0.1

# transition matrix P
P<- array(0.5, c(1, m+1))
P[1, m+1]<- 1

# latent variable S
S<- array(1, c(1, n))
S[1, n]<- m+1

# model parameters theta
theta<- array(n/(m+1), c(1, m+1))

comb<- c(P[-c(m+1)], S[2:(n-1)], theta)	# all combined parameters exclude P[m+1] = 1 & S[1] = 1 & S[n] = m+1
iter<- 0	# iteration counts
thresh<- 5000	# threshold for stopping rule


## mcmc simulation ##
while(1){
	iter<- iter+1
	P_cur<- P[iter, ]
	theta_cur<- theta[iter, ]
	S_cur<- S[iter, ]
	
	# update mass function
	mass<- array(0, c(n, m+1))
	mass[1, 1]<- 1
	for(t in 2:n){
		mass[t, 1]<- 1	# if k == 1 then k-1 == 0
		for(k in 2:(m+1)){
			mass[t, k]<- 1/(1+eq6(t, k-1, P_cur, mass)*ppois(y[t], theta_cur[k-1])/(eq6(t, k, P_cur, mass)*ppois(y[t], theta_cur[k])))
		}
		mass[t,]<- mass[t,]/sum(mass[t,])
	}
	
	# update S
	S_cur[n]<- m+1
	S_cur[1]<- 1
	for(t in (n-1):2){
		prob<- c(prod(mass[t, S_cur[t+1]-1])*(1-prod(P_cur[S_cur[t+1]-1])), mass[t, S_cur[t+1]]*P_cur[S_cur[t+1]])
		prob<- prob/sum(prob)	# probability for pick from S_cur[t+1]-1:S_cur[t+1]
		print(prob)
		S_cur[t]<- sample((S_cur[t+1]-1):S_cur[t+1], 1, prob = prob)
	}
	
	# update P
	count<- ftable(S_cur)-1	#n_ii
	for(i in 1:m){
		x1<- rgamma(1, a+count[i], 1)
		x2<- rgamma(1, b+1, 1)
		P_cur[i]<- x1/(x1+x2)
	}
	
	# update theta
	for(i in 1:(m+1)){
		theta_cur[i]<- rgamma(1, m+1+sum(y[S_cur == i]), count[i]+2)
	}
	
	P<- rbind(P, P_cur)
	S<- rbind(S, S_cur)
	theta<- rbind(theta, theta_cur)
	write.table(P, file = paste(m, "P.txt", sep=''))
	write.table(S, file = paste(m, "S.txt", sep=''))
	write.table(theta, file = paste(m, "theta.txt", sep=''))
	
	# store all parameters together
	comb_cur<- c(P_cur[-c(m+1)], S_cur[2:(n-1)], theta_cur)
	comb<- rbind(comb, comb_cur)
	
	# check stopping rule
	if(iter > thresh){
		comb_mcse<- mcse.mat(comb)[,2]
		comb_sd<- apply(comb, 2, sd)
		thresh<- thresh+1000
		cond<- comb_mcse*1.645+1/iter < 0.02*comb_sd
		write.table(cond, file = paste(m, "cond.txt", sep=''), append = T)
		if(prod(cond)){
			break
		}
	}
}

	
