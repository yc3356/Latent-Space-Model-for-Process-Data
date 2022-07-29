library(igraph)
library(foreign)
library(ggraph)
library(tidyverse)
### read and clearn the log data

data <- read.spss("CBA_cp007q02_logs12_SPSS.sav")
table <- tibble("StIDStd" = data$StIDStd,
                "cnt" = data$cnt,
                "school_id" = data$schoolid,
                "Time" = data$time,
                "Event_number" = data$event_number,
                "Event" = str_replace_all(data$event,pattern = " ",replacement = ""),
                "Event_value" = str_replace_all(data$event_value,pattern = " ",replacement = "")) 



US <- table %>% filter(cnt=="USA") %>%  ## only US data are used in this study
  filter(Event != "ACER_EVENT") %>% 
  filter(Event =="click") %>%
  filter(Event_value != "map") %>%
  filter(Event_value != "paragraph01") %>%
  filter(Event_value != "paragraph02") %>%
  filter(Event_value != "Silver") %>%
  filter(Event_value != "Diamond") %>%
  filter(Event_value != "ItemContent") %>%
  filter(Event_value != "StimulusContent") %>%
  filter(Event_value != "HelpButtonPForm") %>%
  filter(Event_value != "toolbarCopy") %>%
  filter(Event_value != "section05_menu") %>%
  filter(Event_value != "ToolbarInnerFrame") %>%
  filter(Event_value != "Unity") %>%
  filter(Event_value != "mapPos") %>%
  filter(Event_value != "timeDisplay") %>%
  filter(Event_value != "toolbarButtonList") %>%
  filter(Event_value != "Sakharov") %>%
  filter(Event_value != "Market") %>%
  filter(Event_value != "help_28") %>%
  filter(Event_value != "Emerald") %>%
  filter(Event_value != "Park") %>%
  filter(Event_value != "timeMinutes") %>%
  filter(Event_value != "Einstein") %>%
  filter(Event_value != "ItemUnitTitle") %>%
  filter(Event_value != "toolbarPaste") %>%
  filter(Event_value != "Diamondnowhere") %>%
  filter(Event_value != "ItemBottom") %>%
  filter(Event_value != "item1") %>%
  filter(Event_value != "StimulusBottom") %>%
  filter(Event_value != "Mandela") %>%
  filter(Event_value != "minutesWord") %>%
  filter(Event_value != "Nobel") %>%
  filter(Event_value != "sec18_sub02_btn_Copy_on_png") %>%
  filter(Event_value != "ItemQuestionHead") %>%
  filter(Event_value != "tTimeText") %>%
  filter(Event_value != "StimulusHeader") %>%
  filter(Event_value != "ToolbarOutterFrame") %>%
  filter(Event_value != "resetButton") %>%
  select(StIDStd,Event_number, Event_value,Time)


#ignore the students with too few actions
not_quit_student  <- US %>% filter(str_detect(Event_value, 'hit_')) %>% group_by(StIDStd) %>%
  summarise(n_action = n()) %>%
  filter(n_action >=2) %>%
  distinct(StIDStd) %>%
  pull()


US <- US %>% filter(StIDStd %in% not_quit_student)



person_adjmatrix <-list()
person_adjmatrix_RT <-list()
label <- unique(US$Event_value)
for (i in unique(US$StIDStd)){
  cur <- US[which(US$StIDStd==i),]
  cur <- cur %>% arrange(Time)
  tempmatrix <- matrix(0,length(label),length(label))
  rownames(tempmatrix) <- label
  colnames(tempmatrix) <- label
  tempmatrixRT <- matrix(0,length(label),length(label))
  rownames(tempmatrixRT) <- label
  colnames(tempmatrixRT) <- label
  for (j in 1:(nrow(cur)-1)){
    tempmatrix[cur$Event_value[j],cur$Event_value[j+1]] = tempmatrix[cur$Event_value[j],cur$Event_value[j+1]] + 1
    tempmatrixRT[cur$Event_value[j],cur$Event_value[j+1]] = tempmatrixRT[cur$Event_value[j],cur$Event_value[j+1]] + (cur$Time[j+1]-cur$Time[j])
  }
  person_adjmatrix[[i]] <- tempmatrix
  person_adjmatrix_RT[[i]] <- tempmatrixRT
}


adjacent_matrix <-  matrix(0,length(label),length(label))
adjacent_matrix_RT <-  matrix(0,length(label),length(label))
for (i in unique(US$StIDStd)){
  adjacent_matrix = adjacent_matrix + person_adjmatrix[[i]]
  adjacent_matrix_RT = adjacent_matrix_RT + person_adjmatrix_RT[[i]]
}

diag(adjacent_matrix) <- 0
diag(adjacent_matrix_RT) <- 0
 
 
### transfer the data information into the social network data
##### create the adjacent matrix
library(network)
#G <- as.network(adjacent_matrix,matrix.type = 'a',directed = TRUE,loops = TRUE)
G <- network(adjacent_matrix,matrix.type = 'a',directed = TRUE,loops = FALSE,hyper=FALSE)
##### add weight of edge
G %e% 'RA' <- adjacent_matrix
G %v% 'sender_RT' <- apply(adjacent_matrix_RT, 1, mean)
G %v% 'receiver_RT' <- apply(adjacent_matrix_RT, 2, mean)
averaRT <- adjacent_matrix_RT / adjacent_matrix
averaRT[is.na(averaRT)] <- 0

require(sna)
library(latentnet)
library(RSiena)  
number_of_dimension = 3
number_of_group = 3
fit <- ergmm(G ~ euclidean(d=number_of_dimension,G=number_of_group)+edgecov(averaRT) ,response = "RA",family = "Poisson")
#rownames(adjacent_matrix)[which(fit$mcmc.mle$Z.K==2)]

#fit <- ergmm(G ~ euclidean(d=number_of_dimension,G=number_of_group)+sendercov("sender_RT") + receivercov("receiver_RT"),response = "RA",family = "Poisson",ergmm.control(sample.size=5000,burnin=200000))
                 #sendercov("sender_RT") + receivercov("receiver_RT"), family="Poisson",control=ergmm.control(sample.size=8000,burnin=100000,interval=100))
                 #edgecov(averaRT),family = "Poisson",control=ergmm.control(sample.size=8000,burnin=100000,interval=100))
                 #nodecov("sender_RT") + nodecov("receiver_RT"), family = "Poisson")
                 #nodecov("receiver_RT"), family = "Poisson")
               #tofit = c("mcmc"),control=ergmm.control(sample.size=8000,burnin=100000,interval=100))

summary(fit)


#fit <- D3new

#D3new

distance_matrix <- matrix(0,length(label),length(label))
rownames(distance_matrix) <- rownames(adjacent_matrix)
colnames(distance_matrix) <- rownames(adjacent_matrix)


for (i in 1:length(label)){
  for (j in 1:length(label)){
    a <- which(colnames(adjacent_matrix) == label[i])
    b <- which(colnames(adjacent_matrix) == label[j])
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

