#### change point detection ####
library(mcmcse)


m<- 1	# number of change point

raw<- read.table("2change-point.txt")
y<- raw[,1]

#prior information for updating P
a<- 1
b<- 0.1

n<- length(y)	# sample size

P<- array(0.5, c(1, m+1))	# initial transition matrix
P[1, m+1]<- 1
theta<- array(n/(m+1), c(1, m+1))	# initial theta
S<- array(NA, c(1, n))	# initial S
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

thresh<- 1000

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
		theta_cur[i]<- rgamma(1, m+1+sum(y[S_cur == i]), count[i]+1+1)
	}
	
	P<- rbind(P, P_cur)
	S<- rbind(S, S_cur)
	theta<- rbind(theta, theta_cur)
	write.table(P, "P.txt")
	write.table(S, "S.txt")
	write.table(theta, "theta.txt")
	
	#store all parameters together
	if(iter == 1){
		comb<- c(P_cur, S_cur, theta_cur)
	}
	else{
		comb_cur<- c(P_cur, S_cur, theta_cur)
		comb<- rbind(comb, comb_cur)
	}
	print(iter)
	if(iter > thresh){
		comb_mcse<- mcse.mat(comb)[,2]
		comb_sd<- apply(comb, 2, sd)
		thresh<- thresh+100
		if(prod((comb_mcse*1.645+1/iter) < 0.5*comb_sd)){
			break
		}
	}
	
}
