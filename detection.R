#### change point detection ####
library(mcmcse)

args <- commandArgs(TRUE)
m<- as.numeric(args[1])	# number of change point

raw<- read.table("change-point-2.txt")
y<- raw[,1]

n<- length(y)	# sample size

#prior information for updating P
a<- 0.1*n/(m+1)-0.1
b<- 0.1

P<- array(0.5, c(1, m+1))	# initial transition matrix
P[1, m+1]<- 1
theta<- array(n/(m+1), c(1, m+1))	# initial theta
S<- array(0, c(1, n))	# initial S
comb<- c(P[-c(m+1)], S[2:(n-1)], theta)	# exclude P[m+1] = 1 & S[1] = 1 & S[n] = m+1
iter<- 0

p6<- function(t, k, theta, P, mass){	# equition (6)
	if(k == 1){
		output<- P[k]*mass[t-1, k]
	}
	if(k > 1){
		output<- (1-P[k-1])*mass[t-1, k-1] + P[k]*mass[t-1, k]  
	}
	return(output)
}

thresh<- 1000	# threshold for checking stopping criteria

while(1){
	iter<- iter+1
	P_cur<- P[iter, ]
	theta_cur<- theta[iter, ]
	S_cur<- S[iter, ]
	
	# update mass distribution
	mass<- array(0, c(n, m+1))	# mass distribution for S
	for(t in 1:n){
		if(t == 1){
			mass[1, 1]<- 1
		}	
		if(t > 1){
			for(k in 1:(m+1)){
				#temp<- p6(t, k, theta_cur, P_cur, mass)*ppois(y[t], theta_cur[k])+p6(t, k-1, theta_cur, P_cur, mass)*ppois(y[t], theta_cur[k-1])
				mass[t, k]<- p6(t, k, theta_cur, P_cur, mass)*ppois(y[t], theta_cur[k])
			}
			#mass[t, 1]<- 1-sum(mass[t, -1])
			mass[t,]<- mass[t,]/sum(mass[t,]) 
		}
	}
	
	# update S
	S_cur[n]<- m+1
	S_cur[1]<- 1
	for(t in (n-1):2){
		if(S_cur[t+1] == 1){
			S_cur[t] = 1
		}
		else{
			prob<- c(mass[t, (S_cur[t+1]-1)]*(1-P_cur[S_cur[t+1]-1]), mass[t, S_cur[t+1]]*P_cur[S_cur[t+1]])
			S_cur[t]<- sample((S_cur[t+1]-1):S_cur[t+1], 1, prob = prob)
		}
	}
	
	# update P
	count<- ftable(S_cur)-1	# counts for n_ii
	for(i in 1:m){
		P_cur[i]<- rbeta(1, a+count[i], b+1)
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
	
	#store all parameters together
	comb_cur<- c(P_cur[-c(m+1)], S_cur[2:(n-1)], theta_cur)
	comb<- rbind(comb, comb_cur)
	
	if(iter > thresh){
		comb_mcse<- mcse.mat(comb)[,2]
		comb_sd<- apply(comb, 2, sd)
		thresh<- thresh+500
		cond<- comb_mcse*1.645+1/iter < 0.1*comb_sd
		write.table(cond, file = paste(m, "cond.txt", sep=''), append = T)
		if(prod(cond)){
			break
		}
	}
	
}
