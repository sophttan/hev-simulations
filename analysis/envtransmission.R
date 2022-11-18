SEIR_environment <- function(d, res, b, inc, inf) {
  
  for (i in 1:time) {
    d$R[d$Icounter==inf] <- 1
    d$I[d$Icounter==inf] <- 0
    d$Icounter[d$Icounter==inf] <- 0
    
    d$I[d$Ecounter==inc] <- 1
    d$E[d$Ecounter==inc] <- 0
    d$Ecounter[d$Ecounter==inc] <- 0
    
    summary_data <- d %>% group_by(HH) %>% summarise(Sh=sum(S), Ih=sum(I)) %>% 
      mutate(Sc=sum(Sh)-Sh, Ic=sum(Ih)-Ih) %>% ungroup()
    res <- res %>% cbind(I=summary_data$Ih)
    
    new_inf <- b*1000
    new_inf_sus <- b*1000-b*sum(summary_data$Ih)
    print(new_inf)
    print(new_inf_sus)
    
    if(floor(new_inf_sus)>0 & sum(d$S)>0) {
      test_inf <- d %>% filter(S==1) %>% 
        sample_n(ifelse(floor(new_inf_sus)<nrow(.), new_inf_sus, nrow(.)), replace=F)
      d$E[d$No %in% test_inf$No] <- 1
    }
    
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
    
    print(head(d))
  }
  
  return(res)
}

initial_pop <- data.frame(No=1:1000, 
                          HH=rep(1:200, each=5), 
                          S=1, 
                          E=0, Ecounter=0, 
                          I=0, Icounter=0,
                          R=0)
final2 <- SEIR_environment(initial_pop, results, 0.001, num_weeks_inc, num_weeks_inf)
final2[,2:ncol(final2)] %>% sum()

compare <- expand.grid(HH=1:200, time=1:364) 
names(final) <- c("No", "HH", 1:364)
names(final2) <- c("HH", 1:364)
f2 <- final2 %>% pivot_longer(cols = 2:(time+1), names_to = "time") %>% mutate(time=as.numeric(time))

compare <- compare %>% left_join(f, by=c("HH", "time")) %>% left_join(f2, by=c("HH","time"))
compare

compare %>% group_by(time) %>% summarise_all(sum) %>% select(!HH) %>% ggplot(aes(time)) + 
  geom_line(aes(y=value.x), color="blue") + 
  geom_line(aes(y=value.y), color="black") + 
  scale_y_continuous("Number of infections") + 
  scale_x_continuous("Weeks")

f %>% group_by(time) %>% summarise_all(sum) %>% select(!HH) %>% ggplot(aes(time)) + 
  geom_line(aes(y=value), color="black") + 
  scale_y_continuous("Number of current infections") + 
  scale_x_continuous("Days")
