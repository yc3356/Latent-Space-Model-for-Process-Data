
#fit <- ergmm(G ~ euclidean(d=number_of_dimension,G=number_of_group)+receivercov("receiver_RT"),response = "RA",family = "Poisson",ergmm.control(sample.size=5000,burnin=200000))
#fit <- ergmm(G ~ euclidean(d=number_of_dimension,G=number_of_group)+sendercov("sender_RT") ,response = "RA",family = "Poisson",ergmm.control(sample.size=5000,burnin=200000))
#fit <- ergmm(G ~ euclidean(d=number_of_dimension,G=number_of_group)+sendercov("sender_RT") + receivercov("receiver_RT"),response = "RA",family = "Poisson",ergmm.control(sample.size=5000,burnin=200000))



best_fit = fit
summary(best_fit)
#saveRDS(best_fit,'./best_fit.rds')

## visualization
library(dplyr)
x <- best_fit$mcmc.mle$Z[,1]
y <- best_fit$mcmc.mle$Z[,2]
z <- as.factor(best_fit$mcmc.mle$Z.K)
plotdata <- data.frame(dimension1=x,dimenion2=y,cluster=z)

plottext <- substring(rownames(adjacent_matrix),5)
rownames(plotdata) <- plottext

ggplot(plotdata,aes(x=dimension1, y=dimenion2,col=cluster)) +
  geom_point(aes(x=best_fit$mcmc.mle$Z.mean[1,1], y=best_fit$mcmc.mle$Z.mean[1,2]),col='black') +
  geom_point(aes(x=best_fit$mcmc.mle$Z.mean[2,1], y=best_fit$mcmc.mle$Z.mean[2,2]),col='black') +
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











data <- read.spss("CBA_cp007q02_logs12_SPSS.sav")
table <- tibble("StIDStd" = data$StIDStd,
                "cnt" = data$cnt,
                "school_id" = data$schoolid,
                "Time" = data$time,
                "Event_number" = data$event_number,
                "Event" = str_replace_all(data$event,pattern = " ",replacement = ""),
                "Event_value" = str_replace_all(data$event_value,pattern = " ",replacement = "")) 


US <- table %>% filter(cnt=="USA") %>% filter(StIDStd!="     ")
#ignore the students with too few actions
not_quit_student  <- US %>% filter(str_detect(Event_value, 'hit_')) %>% group_by(StIDStd) %>%
  summarise(n_action = n()) %>%
  filter(n_action >=2) %>%
  distinct(StIDStd) %>%
  pull()

US <- US %>% filter(StIDStd %in% not_quit_student)
US <- US %>% filter(Event != "ACER_EVENT")


check <- c()
for (u in unique(US$StIDStd)){
  if ("END_ITEM" %in% US$Event[US$StIDStd ==u]){
    if ("START_ITEM" %in% US$Event[US$StIDStd ==u]){
      check <- c(check, u)
  }
}}

US <- US %>% filter(StIDStd %in% check)

response_length <- c()
for (u in unique(US$StIDStd)){
  response_length <- c(response_length,length(US$Event[US$StIDStd==u & US$Event=="click"]))
}

mean(response_length)
sd(response_length)



n_hit <- c()
for (u in unique(US$StIDStd)){
  temp <- US$Event_value[US$StIDStd==u]
  n <- 0
  for (t in temp){
    if (grepl('hit_', t, fixed=TRUE)){
      n <- n + 1
    }
  }
  n_hit <- c(n_hit,n)
}

mean(n_hit)
sd(n_hit)



response_squence_RT <- c()
for (u in unique(US$StIDStd)){
  begin_time <- US$Time[US$StIDStd==u & US$Event=="START_ITEM"]
  end_time <- US$Time[US$StIDStd==u & US$Event=="END_ITEM"]
  response_squence_RT <- c(response_squence_RT, end_time - begin_time)
}

mean(response_squence_RT)
sd(response_squence_RT)




RT_action <- c()
for (u in unique(US$StIDStd)){
  temp <- sort(US$Time[US$StIDStd ==u])
  RT_action <- c(RT_action,diff(temp))
}
mean(RT_action)
sd(RT_action)

mean(response_length)
sd(response_length)

mean(response_squence_RT)
sd(response_squence_RT)







