library(ggplot2)
library(latentnet)
library(RSiena)  
library(network)
set.seed(3356)



# assume that there are 2 response pattern
n = 10
GroupOneMatrix <- matrix(0,nrow = n,ncol = n)
GroupTwoMatrix <- matrix(0,nrow = n,ncol = n)
GroupOneSequence <- c(1,2,3,4,5,10)
GroupTwoSequence <- c(1,6,7,8,9,10)
m = length(GroupOneSequence)
for (i in 1:(m-1)){
    GroupOneMatrix[GroupOneSequence[i],GroupOneSequence[i+1]] <- GroupOneMatrix[GroupOneSequence[i],GroupOneSequence[i+1]] + 1
    GroupTwoMatrix[GroupTwoSequence[i],GroupTwoSequence[i+1]] <- GroupTwoMatrix[GroupTwoSequence[i],GroupTwoSequence[i+1]] + 1
}

Nstudent <- c(200,500,1000)
prop <- c(0.3,0.5,0.7)
R <- c(0.1,0.2)


grouponedist <- c()
grouptwodist <- c()
tstatistics <- c()
pvalue <- c()

for (m in Nstudent){
  for (r in R){
    for (l in prop){
      response_sequence <- list()
      groupindex <- sample(x = c(0,1),size = m,prob = c(1-l,l),replace = T)
      for (i in 1:m){
        if (groupindex[i]==0){
          temp_matrix <- GroupOneMatrix
        }else{
          temp_matrix <- GroupTwoMatrix
        }
        
        for (a in 1:n){
          for (b in 1:n){
            temp_matrix[a,b] <- temp_matrix[a,b] + runif(1,0,r)
          }
        }
        for (a in 1:n){
          temp_matrix[a,] <- temp_matrix[a,] / sum(temp_matrix[a,])
        }
        temp_matrix[,1] <- 0
        temp_matrix[10,] <- 0
        
        sequences <- c(1)
        while (tail(sequences,1)!=10){
          sequences <- append(sequences,sample(1:10,size=1,prob = temp_matrix[tail(sequences,1),]))
        }  
        response_sequence[[i]] <-sequences
      }
      
      
      transionmatrix <- matrix(0,nrow = n,ncol = n)
      for (i in 1:m){
        cursequence <- response_sequence[[i]]
        for (j in 1:(length(cursequence)-1)){
          transionmatrix[cursequence[j],cursequence[j+1]] = transionmatrix[cursequence[j],cursequence[j+1]] + 1
        }
      }
      rownames(transionmatrix) <- 1:10
      colnames(transionmatrix) <- 1:10
      
      transionmatrix <- transionmatrix[2:(n-1),2:(n-1)]
      
      G <- as.network(transionmatrix,directed = TRUE,loops = TRUE)
      ##### add weight of edge
      G %e% 'weight' <- transionmatrix 
      fit <- ergmm(G ~ euclidean(d = 2,G=2),family = "Poisson",response = 'weight',tofit =c("mle"),seed=3356)
      fname <- paste("simulation","NS",m,"Prop",l,"R",r,".rds",collapse = "",sep = "")
      saveRDS(fit,fname)  
      
      #print(fit$mle$Z.K)
      latent_position <- fit$mle$Z
      key_action <- c(2,3,4,5) - 1
      d1 <- c()
      d2 <- c()
      ind <- 1
      for (sequence in response_sequence){
        if (length(sequence)>2){
          cur_dist <- c()
          for (a in (sequence[2:(length(sequence)-1)])){
            for (b in key_action){
              cur_dist <- c(cur_dist,sqrt(sum((latent_position[a-1,] - latent_position[b,])^2)))
            }
          }
          if (groupindex[ind]==0){
            d1 <- c(d1,mean(cur_dist))
          }else{
            d2 <- c(d2,mean(cur_dist))
          }
        }
        ind = ind + 1
      }
      
      check <- t.test(d1,d2)
      
      
      grouponedist <- c(grouponedist,round(check$estimate[1],3))
      grouptwodist <- c(grouptwodist,round(check$estimate[2],3))
      tstatistics <- c(tstatistics,round(check$statistic,3))
      pvalue <- c(pvalue,round(check$p.value,3))
      
      
    }
  }
}

result <- data.frame(a1=grouponedist,
                     a2=grouptwodist,
                     t=tstatistics,
                     p=pvalue)



write.csv(result,"resultnew.csv")














plots <- list()
cor_r <- c()
sim <- c(100,200,500,1000,5000)
inclusive_cor <-c()
regression <- list()
for (z in 1:5){
  I = sim[z]
  n = 26
  abilities <- rnorm(I,0,3)
  threshold <- (1 - (1/(1+exp(-abilities))))/5
  
  # set ideal matrix
  ideal_matrix <- matrix(0,nrow = n,ncol = n)
  key_minum_sequence <- c(1,sample(2:25,5,replace = F),26)
  for (i in key_minum_sequence[1:(length(key_minum_sequence)-1)]){
    for (j in key_minum_sequence[2:(length(key_minum_sequence))]){
      ideal_matrix[i,j] <- ideal_matrix[i,j] + 1
    }
  }
  # add random noise and normalize
  tranistion_matrix <- list()
  for (i in 1:I){
    tranistion_matrix[[i]] = ideal_matrix
    for (a in 1:n){
      for (b in 1:n){
        noise <- runif(n = 1,min = 0,max = threshold[i])
        tranistion_matrix[[i]][a,b] = tranistion_matrix[[i]][a,b] + noise
      }
    }
    tranistion_matrix[[i]][,1] <- 0
    tranistion_matrix[[i]][n,] <- 0
    for (a in 1:n){
      tranistion_matrix[[i]][a,] <- tranistion_matrix[[i]][a,] / sum(tranistion_matrix[[i]][a,])
    }
    tranistion_matrix[[i]][is.na(tranistion_matrix[[i]])]  <- 0
    #diag(tranistion_matrix[[i]]) <- 0
  }
  
  
  # generate sequence
  response_sequence <- list()
  len <- c()
  for (i in 1:I){
    sequences <- c(1)
    while (length(sequences) < 5) {
      sequences <- c(1)
      while (tail(sequences,1)!=26){
        sequences <- append(sequences,sample(1:26,size=1,prob = tranistion_matrix[[i]][tail(sequences,1),]))
      }  
    }
    response_sequence[[i]] <-sequences
    len <- c(len,length(sequences))
  }

  
  
  
  ## scoring rubric
  ### inclusive approach
  inclusive_scores <- c()
  for (i in 1:I){
    inclusive_scores <- c(inclusive_scores,all(key_minum_sequence %in% response_sequence[[i]]))
  }
  
  ## latent space model for partial scoring
  
  label <- LETTERS
  adjacent_matrix <- matrix(0,n,n)
  rownames(adjacent_matrix) <- label
  colnames(adjacent_matrix) <- label
  
  for (i in 1:I){
    sequences <- response_sequence[[i]]
    for (a in 1:(length(sequences)-1)){
      adjacent_matrix[sequences[a],sequences[a+1]] <- adjacent_matrix[sequences[a],sequences[a+1]] + 1
    }
  }
  cor(as.vector(adjacent_matrix),as.vector(ideal_matrix))
  
  graph_a <- list(adj = adjacent_matrix)
  G <- as.network(graph_a$adj,directed = TRUE,loops = TRUE)
  ##### add weight of edge
  G %e% 'weight' <- graph_a$adj 
  
  
  
  D = 10
  m1 <- ergmm(G ~ euclidean(d = D),family = "Poisson",response = 'weight',tofit = c("mle"))
  dist_matrx <- matrix(0,26,26)
  for (i in 1:26){
    for (j in 1:26){
      dist_matrx[i,j] <- sqrt(sum((m1$mle$Z[i,] - m1$mle$Z[j,])^2))
    }
  }
  
  cor(as.vector(dist_matrx),as.vector(ideal_matrix))
  
  
  parital_distance <- c()
  for (i in 1:I){
    dist <- c()
    for (k in key_minum_sequence){
      for (l in response_sequence[[i]]){
        dist <- c(dist, sqrt(sum((m1$mle$Z[k,] - m1$mle$Z[l,])^2)))
      }
    }
    parital_distance <- c(parital_distance,mean(dist))
  }
  cor_r <- c(cor_r,cor(abilities,parital_distance))
  
  plot_data <- data.frame(
    abilities = abilities,
    parital_distance = parital_distance,
    len = len,
    inclusive_scores=inclusive_scores
  )
  
  regression[[z]] <- plot_data
  a <- ggplot(plot_data,aes(x = parital_distance, y = abilities,col=inclusive_scores)) +
    geom_point(size =1)+
               #aes(col=factor(inclusive_scores))) + 
    geom_smooth(method='lm', formula= y~x)+
    scale_colour_grey()+
    xlab("Partial Scoing") + 
    ylab("Latent Ability") + 
    labs(col="Binary Scoring") +
    theme_minimal()+
    theme(axis.text=element_text(size=10,face="bold"),
          axis.title=element_text(size=10,face="bold"),
          legend.text =element_text(size=10),
          legend.title =element_text(size=10),
          legend.key.size = unit(1, "cm"))
  plots[[z]]<- a
}


library(ggpubr)
ggarrange(plots[[1]], plots[[2]], plots[[3]],plots[[4]],plots[[5]], ncol=2, nrow=3, common.legend = TRUE, legend="bottom",label.x = seq(0,3,0.5),align = "hv")

round(cor_r,3)
#colnames(plot_data)
fit <- lm(abilities~parital_distance+len+inclusive_scores,data=regression[[1]])
summary(fit)

fit <- lm(abilities~parital_distance+len+inclusive_scores,data=regression[[2]])
summary(fit)

fit <- lm(abilities~parital_distance+len+inclusive_scores,data=regression[[3]])
summary(fit)

fit <- lm(abilities~parital_distance+len+inclusive_scores,data=regression[[4]])
summary(fit)

fit <- lm(abilities~parital_distance+len+inclusive_scores,data=regression[[5]])
summary(fit)


table(regression[[5]]$inclusive_scores)




a <- c(1,1,0,0,0,0)
b<- c(0,0,0,0,1,0)
