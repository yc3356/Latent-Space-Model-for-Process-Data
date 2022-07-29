fit <- readRDS("~/Desktop/latent space model/code/D3new.rds")
G <- fit$model$Yg


require(sna)
library(latentnet)
library(RSiena)  
library(network)

adjacent_RT <- as.matrix(G,matrix.type="adjacency",attrname="RT")
adjacent_RA <- as.matrix(G,matrix.type="adjacency",attrname="RA")
averageRT <- adjacent_RT / adjacent_RA
averageRT[is.na(averageRT)] <- 0
G %v% 'sender_RT' <- apply(averageRT, 1, mean)
G %v% 'receiver_RT' <- apply(averageRT, 2, mean)

number_of_dimension = 2
number_of_group = 3
fit <- ergmm(G ~ euclidean(d = number_of_dimension, G = number_of_group) + edgecov(averageRT),response = "RA",family = "Poisson")
#fit <- ergmm(G ~ euclidean(d=number_of_dimension,G=number_of_group)+nodecov("sender_RT") + nodecov("receiver_RT"),response = "RA",family = "Poisson")


summary(fit)
label <- colnames(adjacent_RA)
distance_matrix <- matrix(0,length(label),length(label))
rownames(distance_matrix) <- label
colnames(distance_matrix) <- label

for (i in 1:length(label)){
  for (j in 1:length(label)){
    a <- which(colnames(distance_matrix) == label[i])
    b <- which(colnames(distance_matrix) == label[j])
    distance_matrix[i,j] <- round(sqrt(sum((fit$mkl$Z[a,] - fit$mkl$Z[b,])^2)),2)
  }
}


rownames(distance_matrix) <- str_c('p',as.character(1:23))
colnames(distance_matrix)<- str_c('p',as.character(1:23))

library(reshape2)
library(ggplot2)
plot_data <- melt(distance_matrix)
colnames(plot_data) <- c("action 1", "action 2", "value")
ggplot(data = plot_data, aes(x=`action 1`, y=`action 2`, fill=value)) + 
  geom_tile(aes(fill = value))+
  scale_colour_grey()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))+    
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "blue") 




library(dplyr)
x <- fit$mkl$Z[,1]
y <- fit$mkl$Z[,2]
z <- as.factor(fit$mkl$Z.K)
plotdata <- data.frame(dimension1=x,dimenion2=y,cluster=z)

plottext <- substring(label,5)
rownames(plotdata) <- plottext

ggplot(plotdata,aes(x=dimension1, y=dimenion2,col=cluster)) +
  #geom_point(aes(x=fit$mkl$Z.mean[1,1], y=fit$mkl$Z.mean[1,2]),col='black') +
  #geom_point(aes(x=fit$mkl$Z.mean[2,1], y=fit$mkl$Z.mean[2,2]),col='black') +
  geom_text(aes(label = rownames(plotdata)),size=4,fontface = "bold") + 
  xlab("Dimension 1") + 
  ylab("Dimension 2") + 
  labs(col="cluster") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))















latent_position <- fit$mcmc.pmode$Z



all_sequence <- list()
for (stdid in unique(US$StIDStd)){
  temp <- pull(
    US %>% filter(StIDStd==stdid) %>%
      select(Event_value))
  ans <- c()
  for (t in temp){
    ans <- c(ans,which(t==colnames(adjacent_matrix)))
  }
  all_sequence[[stdid]] <- ans
}

key_action <- c(6,4,5,2,7,10)


all_centeriod <- c()
for (s in all_sequence){
  dist <- c()
  for (k in key_action){
    for (i in s){
      dist <- c(dist,sqrt(sum((fit$mkl$Z[k,] - fit$mkl$Z[i,])^2))) 
    }
  }
  all_centeriod <- c(all_centeriod,mean(dist))
}


student <- unique(US$StIDStd) 

response_accuracy <-table %>% filter(Event == "ACER_EVENT" & cnt=="USA") %>%
  filter(StIDStd %in% student) %>%
  mutate(Event_value = str_replace_all(Event_value, pattern = "'",replacement = "")) %>%
  group_by(StIDStd) %>%
  mutate(ranking = min_rank(desc(Event_number))) %>%
  filter(ranking == 1) %>%
  mutate(RA = (Event_value == "10001011000010001000000")) %>%
  select(StIDStd,RA)



all_dist <- data.frame(
  avg_link = all_centeriod,
  RA = response_accuracy$RA
)
rownames(all_dist) <- names(all_sequence)
colnames(all_dist) <- c('avg_link','RA')


correct_distance <- all_dist %>% filter(RA==TRUE) %>% select(avg_link) %>% pull()
sd(correct_distance)
incorrect_distance <- all_dist %>% filter(RA==FALSE) %>% select(avg_link) %>% pull()
sd(incorrect_distance)
t.test(correct_distance,incorrect_distance)





student_list <- sample(unique(US$StIDStd),8,replace = F)
student_list <- c(student_list,'04648')
student_list <- c(student_list,'03639')
student_list <- c(student_list,'00427')
student_list <- c(student_list,'00849')
student_list <- c(student_list,'04562')
student_list <- c(student_list,'00128')
student_list <- c(student_list,'00849')



plot_data2 <- data_frame(stdid = student_list,
                         avg_link= round(all_dist[student_list,]$avg_link,3),
                         RA = all_dist[student_list,]$RA)


plot_data2 <-plot_data2 %>% mutate(ranking = min_rank(avg_link))
plot_data2$RA <- as.factor(plot_data2$RA)
ggplot(plot_data2,aes(x = ranking, y = avg_link)) +
  geom_point(aes(col=factor(RA)),size =5) + 
  geom_label(aes(x = ranking, y = avg_link+0.1),label = plot_data2$stdid,size=5)+
  geom_hline(yintercept=1.12)+
  scale_colour_grey()+
  #xlim(0,16) + 
  xlab("Rank of Task-takers") + 
  ylab("Average Linkage") + 
  labs(col="Correct Or Not") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=30),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))




all_dist <-all_dist %>% mutate(ranking = min_rank(avg_link))
all_dist$RA <- as.factor(all_dist$RA)
ggplot(all_dist,aes(x = ranking, y = avg_link)) +
  geom_point(aes(col=factor(RA)),size =1) + 
  #geom_label(aes(x = ranking, y = avg_link+0.1),label = plot_data2$stdid,size=5)+
  scale_colour_grey()+
  #xlim(0,16) + 
  xlab("Rank of Task-takers") + 
  ylab("Average Linkage") + 
  labs(col="Correct Or Not") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=30),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))

