ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

model.valid.me <- function(data.re, fit.re,
                           data.fe, fit.fe,
                           B=5, seed=1){
  test_list <- list()
  one_test_list <- list()
  roc_list <- list()
  # Random effect
  auc.train.re <- vector(mode = "numeric", length = B)
  auc.test.re <- vector(mode = "numeric", length = B)
  data.re$pred.prob.full <- fitted(fit.re)
  auc.full.re <- roc(data.re$recruit_yoy,
                     data.re$pred.prob.full, data=data.re)$auc
  
  # Fixed effect
  auc.train.fe <- vector(mode = "numeric", length = B)
  auc.test.fe <- vector(mode = "numeric", length = B)
  data.fe$pred.prob.full <- fitted(fit.fe)
  auc.full.fe <- roc(data.fe$recruit_yoy,
                     data.fe$pred.prob.full, data=data.fe)$auc
  
  # Core loop
  set.seed(seed)
  for(i in 1:B){
    sample <- sample.split(data.fe$recruit_yoy, 
                           SplitRatio = .7)
    
    # Random effect
    train.re <- subset(data.re, sample == TRUE)
    test.re <- subset(data.re, sample == FALSE)
    
    fit.train.re <- glmer(formula(fit.re),
                          data = train.re,
                          family = "binomial",
                          nAGQ=0)
    
    train.re$pred.prob <- fitted(fit.train.re)
    roc.train.re <- roc(train.re$recruit_yoy, 
                        train.re$pred.prob,
                        data=train.re)
    auc.train.re[i] <- roc.train.re$auc
    
    test.re$pred.prob.back <- predict(fit.train.re,
                                      newdata=test.re,
                                      type="response",
                                      allow.new.levels=T)
    roc.test.re <- roc(test.re$recruit_yoy,
                       test.re$pred.prob.back,
                       data=test.re)
    auc.test.re[i] <- roc.test.re$auc
    
    # Fixed effect
    train.fe <- subset(data.fe, sample == TRUE)
    test.fe <- subset(data.fe, sample == FALSE)
    
    fit.train.fe <- glm(formula(fit.fe),
                        data = train.fe,
                        family = "binomial")
    
    train.fe$pred.prob <- fitted(fit.train.fe)
    roc.train.fe <- roc(train.fe$recruit_yoy, 
                        train.fe$pred.prob, data=train.fe)
    auc.train.fe[i] <- roc.train.fe$auc
    
    
    test.fe$pred.prob.back <- predict.glm(fit.train.fe,
                                          newdata=test.fe,
                                          type="response",
                                          allow.new.levels=T)
    roc.test.fe <- roc(test.fe$recruit_yoy,
                       test.fe$pred.prob.back, data=test.fe)
    auc.test.fe[i] <- roc.test.fe$auc
    
    ## Return coords of both ROC ##
    tr.re <- roc.train.re %>%
      coords(ret = c("sensitivity", "specificity"), 
             transpose=F) %>%
      add_column("train_or_test"="train",
                 "FE_or_ME"="ME") %>%
      mutate("point"=1:n())
    te.re <- roc.test.re %>%
      coords(ret = c("sensitivity", "specificity"),
             transpose=F) %>%
      add_column("train_or_test"="test",
                 "FE_or_ME"="ME") %>%
      mutate("point"=1:n())
    tr.fe <- roc.train.fe %>%
      coords(ret = c("sensitivity", "specificity"), 
             transpose=F) %>%
      add_column("train_or_test"="train",
                 "FE_or_ME"="FE") %>%
      mutate("point"=1:n())
    te.fe <- roc.test.fe %>%
      coords(ret = c("sensitivity", "specificity"), 
             transpose=F) %>%
      add_column("train_or_test"="test",
                 "FE_or_ME"="FE") %>%
      mutate("point"=1:n())
    
    roc_list[[i]] <- rbind(tr.re, te.re, tr.fe, te.fe) %>%
      add_column("Loop"=i)
    
    ## One-ROC test ##
    # Train
    auc.test.train.re <- roc.area(obs=as.numeric(
      as.character(train.re$recruit_yoy)),
      pred=train.re$pred.prob) %>%
      as_tibble() %>%
      add_column("train_or_test"="train",
                 "ME_or_FE"="ME")
    
    auc.test.train.fe <- roc.area(obs=as.numeric(
      as.character(train.fe$recruit_yoy)),
      pred=train.fe$pred.prob) %>%
      as_tibble() %>%
      add_column("train_or_test"="train",
                 "ME_or_FE"="FE")
    
    # Test
    auc.test.test.re <- roc.area(obs=as.numeric(
      as.character(test.re$recruit_yoy)),
      pred=test.re$pred.prob.back) %>%
      as_tibble() %>%
      add_column("train_or_test"="test",
                 "ME_or_FE"="ME")
    
    auc.test.test.fe <- roc.area(obs=as.numeric(
      as.character(test.fe$recruit_yoy)),
      pred=test.fe$pred.prob.back) %>%
      as_tibble() %>%
      add_column("train_or_test"="test",
                 "ME_or_FE"="FE")
    
    # Combine and put in list
    one_test_list[[i]] <- rbind(auc.test.train.re,
                                auc.test.train.fe,
                                auc.test.test.re,
                                auc.test.test.fe)
    
    ## Two-ROC test ##
    # Train
    htest.train <- roc.test(roc.train.re, roc.train.fe,
                            method="delong", 
                            alternative="g",
                            paired=T) %>% 
      tidy() %>%
      add_column("train_or_test"="train")
    # Test
    htest.test <- roc.test(roc.test.re, roc.test.fe,
                           method="delong", 
                           alternative="g",
                           paired=T) %>% 
      tidy() %>%
      add_column("train_or_test"="test")
    # Combine and put in list
    test_list[[i]] <- rbind(htest.train, htest.test)
    
  }
  
  ## Putting results together for return() ##
  # Tibble of all roc.test() results
  test_tibble <- bind_rows(test_list)
  one_test_tibble <- bind_rows(one_test_list)
  
  # Tibble of summary statistics for each point along ROC curve
  # Right now, only sensitivity and specificity
  # Add precision and recall?
  roc_tibble <- bind_rows(roc_list) %>% 
    group_by(train_or_test, FE_or_ME, point) %>%
    summarize(mean_sens = mean(sensitivity),
              mean_spec = mean(specificity),
              L_sens = quantile(sensitivity, 0.025),
              U_sens = quantile(sensitivity, 0.975))
  
  ## Summary statistics for AUC ##
  ## Clean this up using dplyr::summarize()
  # Fixed effects
  auc.train.min.fe <- round(min(auc.train.fe), digits=4)
  auc.train.max.fe <- round(max(auc.train.fe), digits=4)
  auc.train.q1.fe <- round(quantile(auc.train.fe, 0.25), digits=4)
  auc.train.q3.fe <- round(quantile(auc.train.fe, 0.75), digits=4)
  auc.train.median.fe <- round(median(auc.train.fe), digits=4)
  auc.train.mean.fe <- round(mean(auc.train.fe), digits=4)
  auc.train.var.fe <- round(var(auc.train.fe), digits=4)
  
  auc.test.min.fe <- round(min(auc.test.fe), digits=4)
  auc.test.max.fe <- round(max(auc.test.fe), digits=4)
  auc.test.q1.fe <- round(quantile(auc.test.fe, 0.25), digits=4)
  auc.test.q3.fe <- round(quantile(auc.test.fe, 0.75), digits=4)
  auc.test.median.fe <- round(median(auc.test.fe), digits=4)
  auc.test.mean.fe <- round(mean(auc.test.fe), digits=4)
  auc.test.var.fe <- round(var(auc.test.fe), digits=4)
  auc.n <- length(auc.train.re)
  
  # Random effects
  auc.train.min.re <- round(min(auc.train.re), digits=4)
  auc.train.max.re <- round(max(auc.train.re), digits=4)
  auc.train.q1.re <- round(quantile(auc.train.re, 0.25), digits=4)
  auc.train.q3.re <- round(quantile(auc.train.re, 0.75), digits=4)
  auc.train.median.re <- round(median(auc.train.re), digits=4)
  auc.train.mean.re <- round(mean(auc.train.re), digits=4)
  auc.train.var.re <- round(var(auc.train.re), digits=4)
  
  auc.test.min.re <- round(min(auc.test.re), digits=4)
  auc.test.max.re <- round(max(auc.test.re), digits=4)
  auc.test.q1.re <- round(quantile(auc.test.re, 0.25), digits=4)
  auc.test.q3.re <- round(quantile(auc.test.re, 0.75), digits=4)
  auc.test.median.re <- round(median(auc.test.re), digits=4)
  auc.test.mean.re <- round(mean(auc.test.re), digits=4)
  auc.test.var.re <- round(var(auc.test.re), digits=4)
  
  # Inefficiently combining AUC summary stats into tibble
  auc_tibble <- data.frame(auc.train.min = c(auc.train.min.fe,
                                             auc.train.min.re),
                           auc.train.max = c(auc.train.max.fe,
                                             auc.train.max.re),
                           auc.train.q1 = c(auc.train.q1.fe,
                                            auc.train.q1.re),
                           auc.train.q3 = c(auc.train.q3.fe,
                                            auc.train.q3.re),
                           auc.train.median = c(auc.train.median.fe,
                                                auc.train.median.re),
                           auc.train.mean = c(auc.train.mean.fe,
                                              auc.train.mean.re),
                           auc.train.var = c(auc.train.var.fe,
                                             auc.train.var.re),
                           auc.test.min = c(auc.test.min.fe,
                                            auc.test.min.re),
                           auc.test.max = c(auc.test.max.fe,
                                            auc.test.max.re),
                           auc.test.q1 = c(auc.test.q1.fe,
                                           auc.test.q1.re),
                           auc.test.q3 = c(auc.test.q3.fe,
                                           auc.test.q3.re),
                           auc.test.median = c(auc.test.median.fe,
                                               auc.test.median.re),
                           auc.test.mean = c(auc.test.mean.fe,
                                             auc.test.mean.re),
                           auc.test.var = c(auc.test.var.fe,
                                            auc.test.var.re),
                           auc.n = c(auc.n, auc.n),
                           "ME_or_FE" = c("FE", "ME")) %>% as.tibble()
  
  # Combine all tibbles into list for return()
  ret_list <- lst(test_tibble, roc_tibble, auc_tibble,
                  one_test_tibble)
  
  return(ret_list)
}