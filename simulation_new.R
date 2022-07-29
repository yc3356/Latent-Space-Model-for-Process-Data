library(latentnet)
# for each student generate their transition matrix
# assume the latent space is in 2 dimension
n.correct = 5
n.incorrect = 10
N <- n.correct + n.incorrect
I = 1000
latent.position.correct = mvtnorm::rmvnorm(n = n.correct,mean =c(0,0))

theta <- sample(x = c(2, 10),size = I,replace = T)
TransitionMatrixList <- list()
#PersonalIntensity <- runif(n = I,0,10)
PersonalIntensity = rep(0,I)
for (i in 1:I){
  latent.position.incorrect = mvtnorm::rmvnorm(n=n.incorrect,mean = c(theta[i],theta[i]))
  latent.position = rbind(latent.position.correct,latent.position.incorrect)
  TempMatrix <- matrix(NA,ncol = N, nrow = N)
  
  for (a in 1:(N-1)){
    for (b in (a+1):N){
      #TempMatrix[a,b] <- rpois(n=1, lambda = exp(PersonalIntensity[i] - dist(rbind(latent.position[a,],latent.position[b,]))))
      p <- boot::inv.logit(PersonalIntensity[i] - dist(rbind(latent.position[a,],latent.position[b,])))
      TempMatrix[a,b] <- sample(x = c(0,1),size = 1,prob = c(1-p,p))
      TempMatrix[b,a] <- TempMatrix[a,b]
    }
  }
  diag(TempMatrix) <- 0 
  TransitionMatrixList[[i]] <- TempMatrix
}

AggMatrix <- matrix(0,ncol = N,nrow = N)
for (i in 1:I) {
  AggMatrix <- AggMatrix + TransitionMatrixList[[i]]
}


G <- as.network(AggMatrix,directed = FALSE,loops = FALSE,multiple = FALSE,hyper = FALSE,matrix.type = 'adjacency')
G %e% 'weight' <- AggMatrix
m1 <- ergmm(G ~ euclidean(d = 2, G =2),family = "Poisson",response = 'weight',tofit = c("mle"))

# calculate the average Linkage
positions <- m1$mle$Z

DistanceMatrix <- matrix(NA,ncol=N,nrow=N)
for (i in 1:N){
  for (j in 1:N){
    DistanceMatrix[i,j] = dist(rbind(positions[i,],positions[j,]))
  }
}
DistanceMatrix <- round(DistanceMatrix,2)


optimal.positions <- positions[1:n.correct,]
PartialScore <- c()
for (i in 1:I){
  temp <- apply(TransitionMatrixList[[i]],1,sum)
  TempDist <- rep(DistanceMatrix[which(temp!=0),n.correct],temp[temp!=0])
  PartialScore <- c(PartialScore, mean(TempDist))
}

t.test(PartialScore[theta==2],PartialScore[theta==10])



which(is.na(PartialScore))




# based on idea transition matrix
n = 26
nkey = 5
optimal.strategy <- c("a",sample(x = letters[2:25],size = nkey,replace = F),"z")

optimal.matrix <- matrix(0,ncol = n,nrow = n)
colnames(optimal.matrix) <- letters
rownames(optimal.matrix) <- letters
for (i in 1:(length(optimal.strategy)-1)){
  optimal.matrix[optimal.strategy[i],optimal.strategy[i+1]] = 1
}
I = 1000
theta <- runif(I,0,10)
response_sequence <- list()
for (i in 1:I){
  noise <- matrix(runif(n*n,0,theta[i]),nrow = n,ncol = n)
  temp <- noise + optimal.matrix
  for (a in 1:(n-1)){
    temp[a,] <- exp(temp[a,]) / sum(exp(temp[a,]))
  }
  temp[26,] <- 0
  diag(temp) <- 0
  tempseq <- c("a")
  while (tail(tempseq,1) != "z"){
    tempseq <- c(tempseq, sample(letters,size=1,prob = temp[tail(tempseq,1),]))
  }
  response_sequence[[i]] <- tempseq
}


adjacent <- matrix(0,nrow = n,ncol = n)
colnames(adjacent) <- letters
rownames(adjacent) <- letters
for (curseq in response_sequence){
  for (j in 1:(length(curseq)-1)){
    adjacent[curseq[j],curseq[j+1]] = adjacent[curseq[j],curseq[j+1]] + 1
  }
}


G <- as.network(adjacent,directed = TRUE,loops = FALSE,multiple = FALSE,hyper = FALSE,matrix.type = 'adjacency')
G %e% 'weight' <- adjacent
m1 <- ergmm(G ~ euclidean(d = 1),family = "Poisson",response = 'weight',tofit = c("mle"))


positions <- m1$mle$Z
rownames(positions) <- letters

PartialScore <- c()
for (i in 1:I){
  optimal.positions <- positions[optimal.strategy,]
  curseq <- response_sequence[[i]]
  tempdist <- 0
  for (j in curseq){
    tempdist = tempdist + sum(abs(positions[j,] - optimal.positions))
  }
  PartialScore <- c(PartialScore, tempdist/(length(curseq)*(nkey+2)))
}

cor(PartialScore,theta)



  








