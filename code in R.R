#----Question 1-------
IterationsAlg <- function(B, z, x0, delta=10^-5){
  gap <- Inf
  temp <- x0
  iter <- 1
  iterations<-matrix(data = NA, ncol = 2, nrow = 100000)
  iterations[iter,] <- temp
  while(gap>delta & iter<100000){
    iter <- iter+1
    iterations[iter,] <- B%*%temp+z
    gap <- max(abs(iterations[iter,]-temp))
    temp <- iterations[iter,]
  }
  iterations <- na.omit(iterations)
  return(iterations)}
#-----Gauss-Seidel-----
initGauSei <- function(A, b){
  L <- lower.tri(A, diag = TRUE)*A
  U <- A-L
  Linv <- inv(L)
  B <- -Linv%*%U
  z <- Linv%*%b
  res <- list(B,z)
  return(res)
}

A <- matrix(data = c(2,3,3,5), nrow = 2)
b <- c(13,21)
x0 <- c(1,1)
temp <- initGauSei(A, b)
B <- temp[[1]]
z <- temp[[2]]
res<-IterationsAlg(B,z,x0)
par(mfrow=c(1,2))
plot(res, main = '2D',
     xlab = 'x1', ylab = 'x2', type = "b")
scatterplot3d(z= 1:101, x=res[,1], y = res[,2], type = "b", main='3D',
              zlab = "iteration", xlab = 'x1', ylab = 'x2')
#-----Jacobi-----
initJacobi <- function(A, b){
  D <- diag(diag(A))
  E <- A-D
  Dinv <- inv(D)
  B <- -Dinv%*%E
  z <- Dinv%*%b
  res <- list(B,z)
  return(res)}

A <- matrix(data = c(2,3,3,5), nrow = 2)
b <- c(13,21)
x0 <- c(1,1)
temp <- initJacobi(A, b)
B <- temp[[1]]
z <- temp[[2]]
res<-IterationsAlg(B,z,x0)
par(mfrow=c(1,2))
plot(res, main = '2D', xlab = 'x1', ylab = 'x2', type = "b",)
scatterplot3d(z= 1:248, x=res[,1], y = res[,2], type = "b",
              main='3D',zlab = "iteration", xlab = 'x1', ylab = 'x2')

#----Gradient Descent----
alpha_i <- function(A, x, b){
  temp <- A%*%x-b
  return((t(temp)%*%temp)/(t(temp)%*%A%*%temp))
}
A <- matrix(data = c(2,3,3,5), nrow = 2)
b <- c(13,21)
x0 <- c(1,1)
eigen(A)
gap <- Inf
delta <- 10^-5
iter <- 1
xi <- x0
iterations<-matrix(data = NA, ncol = 2, nrow = 100000)
iterations[iter,] <- xi
while(gap>delta & iter<100000){
  iter <- iter+1
  alpha <- alpha_i(A, xi, b) 
  iterations[iter,] <- xi-alpha%*%t((A%*%xi-b)) 
  gap <- max(abs(iterations[iter,]-xi))
  xi <- iterations[iter,]
}
iterations <- na.omit(iterations)
par(mfrow=c(1,2))
plot(iterations, xlab = 'x1', ylab = 'x2', type = "b", main = '2D')
scatterplot3d(z= 1:6, x=iterations[,1], y = iterations[,2], type = "b",
              zlab = "iteration", xlab = 'x1', ylab = 'x2', main = '3D')

#---Coordinate Descent----
A <- matrix(data = c(2,3,3,5), nrow = 2)
b <- c(13,21)
x0 <- c(1,1)
gap <- Inf
delta <- 10^-5
iter <- 1
x1i <- x0[1]
x2i <- x0[2]
iterations<-matrix(data = NA, ncol = 2, nrow = 100000)
iterations[iter,] <- c(x1i,x2i)
while(gap>delta & iter<100000){
  iter <- iter+1
  x1i <- (13-3*x2i)/2
  iterations[iter,] <- c(x1i,x2i)
  iter <- iter+1
  x2i <- (21-3*x1i)/5
  iterations[iter,] <- c(x1i,x2i) 
  gap <- max(abs(iterations[iter,]-xi))
  xi <- iterations[iter,]
}
iterations <- na.omit(iterations)
par(mfrow=c(1,2))
plot(iterations, xlab = 'x1', ylab = 'x2', type = "b", main = '2D')
scatterplot3d(z= 1:201, x=iterations[,1], y = iterations[,2], type = "b",
              zlab = "iteration", xlab = 'x1', ylab = 'x2', main = '3D')



#----Question 2-------
#--- inverse probability transform-----
samples <- runif(1000)
results<- qbinom(samples, 10, 0.5)
hist(results, col='green', main = 'Histogram of the samples', breaks = -1:10)

probs <- rep(NA,11)
quantiales <- rep(NA,12)
binomProbs <- function(k){return(0.5^10*nchoosek(10,k))}
temp <- 0
for (i in 1:11) {
  probs[i] <- binomProbs(i-1)
  quantiales[i] <- probs[i]+temp
  temp <- quantiales[i]
}
samples <- runif(1000)
results <- rep(NA,1000)
for (i in 1:1000) {
  quantiales[12] <- samples[i]
  results[i] <- rank(quantiales)[12]-1}
hist(results, col='green', main = 'Histogram of the samples', breaks = -1:10)

#---Accept reject----
probs <- rep(NA,11)
binomProbs <- function(k){return(0.5^10*nchoosek(10,k))}
for (i in 1:11) {
  probs[i] <- binomProbs(i-1)
}
m <- max(probs)*11
results <- rep(NA,1000)
for (i in 1:1000) {
  flag <- TRUE
  while (flag) {
    x <- floor(runif(1)*11)
    u <- runif(1)
    if(probs[x+1]*11>=m*u){
      flag <- FALSE}}
  results[i] <- x
}
hist(results, col='green', main = 'Histogram of the samples', breaks = -1:10)

#----Metropolis-Heistings----
binomProbs <- function(k){
  if(k<0 | k>10){ans=0}
  else{ans=0.5^10*nchoosek(10,k)}
  return(ans)}
proposal <- function(x){x+floor(5*runif(1))-2}
results <- rep(NA,1000)
xt <- 5
for (i in 1:1000) {
  proposedx <- proposal(xt)
  a <- binomProbs(proposedx)/binomProbs(xt)
  if(runif(1)<=a){xt <- proposedx}
  results[i] <- xt}
hist(results, col='green', main = 'Histogram of the samples', breaks = -1:10)
#----Question 3-------
LC_Generator<- function(a, m, c, x0, num){
  res <- rep(NA, num)
  temp <- x0
  for (i in 1:num) {
    temp <- (a*temp+c) %% m
    res[i] <- temp}
  print(res)}

LC_Generator(1664525,2^32,0,3,10)


#----Question 4-----
nonconvex <- function(lambda, X){
    return(lambda*exp(-lambda*X))
}
samples <- seq(0.1,10, 0.1)
likelihood<-nonconvex(samples,2)
plot(x= samples, xlab ='lambda', likelihood, type = 'l', main = 'Likelihood of Lambda ')
