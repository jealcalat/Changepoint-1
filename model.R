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
##		P.txt, S.txt, theta.txt, cond.txt, bayes.txt, Rplot.pdf	##
##################################################################


library(mcmcse)


## input m ##
args <- commandArgs(TRUE)
m<- as.numeric(args[1])


## input y ##
raw<- read.table("change-point-num.txt", sep='')
y<- raw[,1]
n<- length(y)	# sample size


## initial setups ##
# hyper parameters
a<- round(0.1*n/(m+1)-0.1)
b<- 0.1

# transition matrix P
P<- array(0.5, c(1, m+1))
P[1, m+1]<- 1

# latent variable S
S<- array(m+1, c(1, n))
for(t in 1:(m+1)){
	S[t]<- t
}

# model parameters theta
theta<- array(mean(y), c(1, m+1))

comb<- c(P[-c(m+1)], S[2:(n-1)], theta)	# all combined parameters exclude P[m+1] = 1 & S[1] = 1 & S[n] = m+1
iter<- 0	# iteration counts
thresh<- 5000	# threshold for stopping rule


## mcmc simulation ##
while(1){
	iter<- iter+1
	P_cur<- P[iter, ]
	theta_cur<- theta[iter, ]
	S_cur<- S[iter, ]
	
	# update S
	eq6<- array(0, c(n, m+1))
	mass<- array(0, c(n, m+1))
	eq6[1, 1]<- 1
	mass[1, 1]<- 1
	for(tt in 2:n){
		eq6[tt, 1]<- P_cur[1]*mass[tt-1, 1]
		for(kk in 2:(m+1)){
			eq6[tt, kk]<- (1-P_cur[kk-1])*mass[tt-1, kk-1]+P_cur[kk]*mass[tt-1, kk]
		}
		for(kk in 1:(m+1)){
			mass[tt, kk]<- eq6[tt, kk]*dpois(y[tt], theta_cur[kk])
		}
		mass[tt,]<- mass[tt,]/sum(mass[tt,])
	}
	S_cur[n]<- m+1
	for(t in (n-1):1){
		k<- S_cur[t+1]
		if(k > 1){
			prob<- c(mass[t, k-1]*(1-P_cur[k-1]), mass[t, k]*P_cur[k])
			prob<- prob/sum(prob)	# probability for pick from S_cur[t+1]-1:S_cur[t+1]
			S_cur[t]<- sample((k-1):k, 1, prob = prob)
		}
		if(k == 1){
			S_cur[t]<- 1	
		}
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
		cond<- comb_mcse*1.96+1/iter < 0.05*comb_sd
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
eq6<- array(0, c(n, m+1))
mass<- array(0, c(n, m+1))
eq6[1, 1]<- 1
mass[1, 1]<- 1
y_like[1]<- dpois(y[1], theta_star[1])
for(tt in 2:n){
	eq6[tt, 1]<- P_star[1]*mass[tt-1, 1]
	for(kk in 2:(m+1)){
		eq6[tt, kk]<- (1-P_star[kk-1])*mass[tt-1, kk-1]+P_star[kk]*mass[tt-1, kk]
	}
	for(kk in 1:(m+1)){
		mass[tt, kk]<- eq6[tt, kk]*dpois(y[tt], theta_star[kk])
	}
	mass[tt,]<- mass[tt,]/sum(mass[tt,])
	for(kk in 1:(m+1)){
		y_like[tt]<- y_like[tt]+dpois(y[tt], theta_star[kk])*eq6[tt, kk]
	}
}
ln_y_like<- sum(log(y_like))	# likelihood function
write.table(ln_y_like, file = paste(m, "ln_y_like.txt", sep = ''))

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
		theta_post[i, k]<- dgamma(theta_star[k], m+1+sum(y[S[i,] == k]), table(S[i,])[k]+1)
	}
	
	# update S_post
	eq6<- array(0, c(n, m+1))
	mass<- array(0, c(n, m+1))
	eq6[1, 1]<- 1
	mass[1, 1]<- 1
	for(tt in 2:n){
		eq6[tt, 1]<- P_update[1]*mass[tt-1, 1]
		for(kk in 2:(m+1)){
			eq6[tt, kk]<- (1-P_update[kk-1])*mass[tt-1, kk-1]+P_update[kk]*mass[tt-1, kk]
		}
		for(kk in 1:(m+1)){
			mass[tt, kk]<- eq6[tt, kk]*dpois(y[tt], theta_star[kk])
		}
		mass[tt,]<- mass[tt,]/sum(mass[tt,])
	}
	S_post[i, n]<- m+1
	for(t in (n-1):1){
			k<- S_post[i, t+1]
		if(k > 1){
			prob<- c(mass[t, k-1]*(1-P_update[k-1]), mass[t, k]*P_update[k])
			prob<- prob/sum(prob)	# probability for pick from S_cur[t+1]-1:S_cur[t+1]
			S_post[i, t]<- sample((k-1):k, 1, prob = prob)
		}
		if(k == 1){
			S_post[i, t]<- 1	
		}
	}

	count_post<- table(S_post[i,])-1	#n_ii
	for(k in 1:m){
		x1<- rgamma(1, a+count_post[k], 1)
		x2<- rgamma(1, b+1, 1)
		P_update[k]<- x1/(x1+x2)
		P_post[i, k]<- gamma(a+b+count_post[k]+1)/(gamma(a+count_post[k])*gamma(b+1))*P_star[k]^(a+count_post[k]-1)*(1-P_star[k])^(b+1-1)
	}
	
	# update S_plot
	eq6<- array(0, c(n, m+1))
	mass<- array(0, c(n, m+1))
	eq6[1, 1]<- 1
	mass[1, 1]<- 1
	for(tt in 2:n){
		eq6[tt, 1]<- P[i, 1]*mass[tt-1, 1]
		for(kk in 2:(m+1)){
			eq6[tt, kk]<- (1-P[i, kk-1])*mass[tt-1, kk-1]+P[i, kk]*mass[tt-1, kk]
		}
		for(kk in 1:(m+1)){
			mass[tt, kk]<- eq6[tt, kk]*dpois(y[tt], theta[i, kk])
		}
		mass[tt,]<- mass[tt,]/sum(mass[tt,])
	}
	S_plot[i, ,]<- eq6
}

ln_theta_post_den<- log(mean(apply(theta_post, 1, prod)))
write.table(ln_theta_post_den, file = paste(m, "ln_theta_post_den.txt", sep = ''))
ln_P_post_den<- log(mean(apply(P_post, 1, prod)))
write.table(ln_P_post_den, file = paste(m, "ln_P_post_den.txt", sep = ''))
ln_theta_den<- sum(log(apply(as.array(theta_star), 1, dgamma, m+1, 1)))
write.table(ln_theta_den, file = paste(m, "ln_theta_den.txt", sep = ''))
ln_P_den<- sum(log(apply(as.array(P_star), 1, dbeta, a, b)))
write.table(ln_P_den, file = paste(m, "ln_P_den.txt", sep = ''))

# bayes factor
ln_bayes<- ln_y_like+ln_theta_den+ln_P_den - ln_theta_post_den-ln_P_post_den
write.table(ln_bayes, file = paste(m, "bayes.txt", sep = ''))

# plots
S_plots<- apply(S_plot, c(2,3), mean)
x<- seq(1, n, 1)

# device
pdf(file = paste(m, "Rplot.pdf", sep=''))
par(mfrow = c(ceiling(m/2)+1, 2))	# split plot
# for(k in 1:m){
	# hist(time[,k], main = paste(k, "th change-point", sep=''), xlab = "Time")
# }
for(k in 1:(m+1)){
	plot(density(theta[, k]), main = paste("lambda", k))
}
plot(x, S_plots[,1], "l", ylim = c(0, 1), xlab = "Time", ylab = "Pr(S|Y)", main = "Prob for change points")
for(k in 2:(m+1)){
	lines(x, S_plots[,k], lty = 2)
}
dev.off()




	
