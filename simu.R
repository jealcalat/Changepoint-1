#### Change point simulator for poisson distribution ####

## 2 change-point case ##
n1<- 8
n2<- 12
lambda1<- 2
lambda2<- 5

x1<- rpois(n1, lambda1)
x2<- rpois(n2, lambda2)
x<- c(x1, x2)	# combined sample
write.table(x, "2change-point.txt")