#### model comparision ####

args <- commandArgs(TRUE)
m<- as.numeric(args[1])	# number of change point

# read mcmc results
P<- read.table(paste(m, "P.txt", sep=''), row.names = NULL)
P<- as.matrix(P[,-1])
P<- P[-1,]
S<- read.table(paste(m, "S.txt", sep=''), row.names = NULL)
S<- as.matrix(S[,-1])
S<- S[-1,]
theta<- read.table(paste(m, "theta.txt", sep=''), row.names = NULL)
theta<- as.matrix(theta[,-1])
theta<- theta[-1,]

G<- dim(P)[1]	# number of total simulation effort G

raw<- read.table("change-point-2.txt")
y<- raw[,1]
n<- length(y)	# sample size

a<- 0.1*n/(m+1)-0.1
b<- 0.1

# posterior mean estimates
P_star<- apply(P, 2, mean)
theta_star<- apply(theta, 2, mean)

# function
p6<- function(t, k, theta, P, mass){	# equition (6)
	if(k == 0){
		output<- 0
	}
	if(k == 1){
		output<- P[k]*mass[t-1, k]
	}
	if(k > 1){
		output<- (1-P[k-1])*mass[t-1, k-1] + P[k]*mass[t-1, k]  
	}
	return(output)
}

# posterior
theta_post<- array(NA, c(G, m+1))
S_post<- array(NA, c(G, n))
P_post<- array(NA, c(G, m))
P_update<- P_star
S_plot<- array(0, c(G, n, m+1))

for(i in 1:G){
	# update theta_post
	S_cur<- S[i,]
	count<- ftable(S_cur)-1
	for(j in 1:(m+1)){
		theta_post[i, j]<- pgamma(theta_star[j], m+1+sum(y[S_cur == j]), count[j]+2)
	}
	
	# update mass distribution
	mass<- array(0, c(n, m+1))	# mass distribution for S
	mass_plot<- array(0, c(n, m+1))
	for(t in 1:n){
		if(t == 1){
			mass[1, 1]<- 1
			mass_plot[1, 1]<- 1
			S_plot[i, t, 1]<- 1
		}
		else{
			for(k in 1:(m+1)){
				temp<- p6(t, k, theta_star, P_update, mass)*ppois(y[t], theta_star[k])+p6(t, k-1, theta_star, P_update, mass)*ppois(y[t], prod(theta_star[k-1]))
				mass[t, k]<- p6(t, k, theta_star, P_update, mass)*ppois(y[t], theta_star[k])/temp
				
				temp_plot<- p6(t, k, theta[i,], P[i,], mass_plot)*ppois(y[t], theta[i,k])+p6(t, k-1, theta[i,], P[i,], mass_plot)*ppois(y[t], prod(theta[i,k-1]))
				mass_plot[t, k]<- p6(t, k, theta[i,], P[i,], mass_plot)*ppois(y[t], theta[i,k])/temp_plot
			}
			}
		mass_plot[t,]<- mass_plot[t,]/sum(mass_plot[t,])
	}
	for(t in 2:n){
		for(k in 1:(m+1)){
			S_plot[i, t, k]<- p6(t, k, theta[i,], P[i,], mass_plot)
		}
	}
	
	# update S_post
	S_post[i, n]<- m+1
	S_post[i, 1]<- 1
	for(t in (n-1):2){
		if(S_post[i, t+1] == 1){
			S_post[i, t] = 1
		}
		else{
			prob<- c(mass[t, (S_post[i, t+1]-1)]*(1-P[i, S_post[i, t+1]-1]), mass[t, S_post[i, t+1]]*P[i, S_post[i, t+1]])
			S_post[i, t]<- sample((S_post[i, t+1]-1):S_post[i, t+1], 1, prob = prob)
		}
	}
	
	# update P_post
	count_post<- ftable(S_post[i,])-1
	for(j in 1:m){
		P_update[j]<- rbeta(1, a+count_post[j], b+1)
		P_post[i, j]<- pbeta(P_star[j], a+count_post[j], b+1)
	}
}

y_post<- array(0, n)
mass<- array(0, c(n, m+1))	# mass distribution for S
for(t in 1:n){
	if(t == 1){
		mass[1, 1]<- 1
		y_post[t]<- y_post[t]+ppois(y[t], theta_star[k])
	}
	else{
		for(k in 1:(m+1)){
			temp<- p6(t, k, theta_star, P_star, mass)*ppois(y[t], theta_star[k])+p6(t, k-1, theta_star, P_star, mass)*ppois(y[t], prod(theta_star[k-1]))
			mass[t, k]<- p6(t, k, theta_star, P_star, mass)*ppois(y[t], theta_star[k])/temp
		}
		mass[t,]<- mass[t,]/sum(mass[t,])
		for(k in 1:(m+1)){
			y_post[t] = y_post[t]+ppois(y[t], theta_star[k])*p6(t, k, theta_star, P_star, mass)
		}
	}

}

theta_post_den<- mean(apply(theta_post, 1, prod))
P_post_den<- mean(apply(P_post, 1, prod))
ln_theta_den<- sum(log(apply(as.array(theta_star), 1, pgamma, m+1, 1)))
ln_P_den<- sum(log(apply(as.array(P_star), 1, pbeta, a, b)))
ln_y_den<- sum(log(y_post))

ln_bayes<- ln_y_den+ln_theta_den+ln_P_den - log(theta_post_den)-log(P_post_den)
write.table(ln_bayes, file = paste(m, "bayes.txt", sep = ''))
S_plots<- apply(S_plot, c(2,3), mean)
write.table(S_plots, file = paste(m, "Splots.txt", sep = ''))
