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
		S_cur[t]<- sample((S_cur[t+1]-1):S_cur[t+1], 1, prob = prob)
	}
	
	# update P
	count<- table(S_cur)-1	#n_ii
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


## posterior analysis ##
# total simulation effort
G<- dim(P)[1]

# posterior mean
P_star<- apply(P, 2, mean)
theta_star<- apply(theta, 2, mean)

# likelihood at star
y_like<- array(0, n)
mass<- array(0, c(n, m+1))
mass[1, 1]<- 1
y_like[1]<- ppois(y[1], theta_star[1])
for(t in 2:n){
	mass[t, 1]<- 1	# if k == 1 then k-1 == 0
	for(k in 2:(m+1)){
		mass[t, k]<- 1/(1+eq6(t, k-1, P_star, mass)*ppois(y[t], theta_star[k-1])/(eq6(t, k, P_star, mass)*ppois(y[t], theta_star[k])))
	}
	mass[t,]<- mass[t,]/sum(mass[t,])
	for(k in 1:(m+1)){
		y_like[t]<- y_like[t]+ppois(y[t], theta_star[k])*eq6(t, k, P_star, mass)
	}
}
ln_y_like<- sum(log(y_like))	# likelihood function

#marginal likelihood
time<- array(NA, c(G, m))	# change-points
theta_post<- array(NA, c(G, m+1))
S_post<- array(NA, c(G, n))
P_post<- array(NA, c(G, m))
P_update<- P_star
S_plot<- array(0, c(G, n, m+1))
for(i in 1:G){
	# update time
	time[i,]<- cumsum(table(S[i,]))[-(m+1)]
	
	# update theta_post
	for(k in 1:(m+1)){
		theta_post[i, k]<- pgamma(theta_star[k], m+1+sum(y[S[i,] == k]), table(S[i,])[k]+1)
	}
	
	# update S_post
	mass<- array(0, c(n, m+1))
	mass[1, 1]<- 1
	for(t in 2:n){
		mass[t, 1]<- 1	# if k == 1 then k-1 == 0
		for(k in 2:(m+1)){
			mass[t, k]<- 1/(1+eq6(t, k-1, P_update, mass)*ppois(y[t], theta_star[k-1])/(eq6(t, k, P_update, mass)*ppois(y[t], theta_star[k])))
		}
		mass[t,]<- mass[t,]/sum(mass[t,])
	}
	S_post[i, n]<- m+1
	S_post[i, 1]<- 1
	for(t in (n-1):2){
		prob<- c(prod(mass[t, S_post[i, t+1]-1])*(1-prod(P_update[S_post[i, t+1]-1])), mass[t, S_post[i, t+1]]*P_update[S_post[i, t+1]])
		prob<- prob/sum(prob)	# probability for pick from S_cur[t+1]-1:S_cur[t+1]
		S_post[i, t]<- sample((S_post[i, t+1]-1):S_post[i, t+1], 1, prob = prob)
	}
	count_post<- table(S_post[i,])-1	#n_ii
	for(k in 1:m){
		x1<- rgamma(1, a+count_post[k], 1)
		x2<- rgamma(1, b+1, 1)
		P_update[k]<- x1/(x1+x2)
		P_post[i, k]<- gamma(a+b+count_post[k]+1)/(gamma(a+count_post[k])*gamma(b+1))*P_star[k]^(a+count_post[k]-1)*(1-P_star[k])^(b+1-1)
	}
	
	# update S_plot
	mass<- array(0, c(n, m+1))
	mass[1, 1]<- 1
	for(t in 2:n){
		mass[t, 1]<- 1	# if k == 1 then k-1 == 0
		for(k in 2:(m+1)){
			mass[t, k]<- 1/(1+eq6(t, k-1, P[i,], mass)*ppois(y[t], theta[i, k-1])/(eq6(t, k, P[i,], mass)*ppois(y[t], theta[i, k])))
		}
		mass[t,]<- mass[t,]/sum(mass[t,])
	}
	for(t in 1:n){
		for(k in 1:(m+1)){
			S_plot[i, t, k]<- eq6(t, k, P[i,], mass)
		}
	}
}

ln_theta_post_den<- log(mean(apply(theta_post, 1, prod)))
ln_P_post_den<- log(mean(apply(P_post, 1, prod)))
ln_theta_den<- sum(log(apply(as.array(theta_star), 1, pgamma, m+1, 1)))
ln_P_den<- sum(log(apply(as.array(P_star), 1, pbeta, a, b)))

# bayes factor
ln_bayes<- ln_y_den+ln_theta_den+ln_P_den - ln_theta_post_den-ln_P_post_den
write.table(ln_bayes, file = paste(m, "bayes.txt", sep = ''))

# plots
S_plots<- apply(S_plot, c(2,3), mean)
write.table(S_plots, file = paste(m, "Splots.txt", sep = ''))




	
