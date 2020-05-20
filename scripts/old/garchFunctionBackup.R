# get all terminal nodes of tree given subsample
collectSubsamplesTree = function(returns_df, split_variables_df, split_list){
  list_subsamples = list()
  # get subsample from every row
  for (i in 1:nrow(split_list)){
    list_subsamples[[i]] = sampleToSubsampleTree(returns_df, split_variables_df, split_list[i,])
  }
  # remove subsamples that are equivalent (because splitting ended earlier)
  final_list_subsamples = unique(list_subsamples)
  return(final_list_subsamples)
}

# get subsample given tree splits
sampleToSubsampleTree = function(returns_df, split_variables_df, split_information){
  # make sure same format for merging by date
  returns_df$date = as.Date(returns_df$date)
  split_variables_df$date = as.Date(split_variables_df$date)
  
  # return old data if no splits 
  subsample = returns_df
  
  # first split
  # check if split was done
  if(split_information$split_variable1[1]=="split_removed"){
    print("No first split")
  } else {
    
    # select column return and splitting variable
    split_data = left_join(returns_df, split_variables_df, by="date")
    split_covariate = drop_na(split_data[,colnames(split_data) %in% c("date","return",split_information$split_variable1[1])])
    split_covariate = split_covariate[,c("date","return", split_information$split_variable1[1])] # order covariates
    colnames(split_covariate) = c("date","return","selected_covariate")
    split_covariate$date = as.Date(split_covariate$date)
    
    if(split_information$higher_lower1[1] == "lower") {
      subsample = filter(split_covariate, selected_covariate < split_information$threshhold1[1])
    } else {
      subsample = filter(split_covariate, selected_covariate >= split_information$threshhold1[1])
    }
    
    # Second split
    # check if split was done
    if(split_information$split_variable2[1]=="split_removed"){
      print("No second split")
    } else {
      # select column return and splitting variable 
      split_data = left_join(subsample, split_variables_df, by="date")
      split_covariate = drop_na(split_data[,colnames(split_data) %in% c("date","return",split_information$split_variable2[1])])
      split_covariate = split_covariate[,c("date","return", split_information$split_variable2[1])] # order covariates
      colnames(split_covariate) = c("date","return","selected_covariate")
      split_covariate$date = as.Date(split_covariate$date)
      
      if(split_information$higher_lower2[1] == "lower") {
        subsample = filter(split_covariate, selected_covariate < split_information$threshhold2[1])
      } else {
        subsample = filter(split_covariate, selected_covariate >= split_information$threshhold2[1])
      }
    }
  }
  subsample = select(subsample, date, return) # drop splitting variable
  return(subsample)
}

# Build a tree with 3 splits and prunes via AIC.
  # model specification: GARCH 1,1 with t-innovations
  # also splits if splits decrease likelihood (in unlikely case)
  # simplyfication: splits each leaf again (does not search for optimal split across leaves)
buildAndPruneTree = function(returns, split_variables, list_split_variables,model_specification, max_lags)  {
  # 2) Build tree with 3 splits
  
    # first split
      # run split function for both nodes and choose best split according to likelihood improvement 
      split_1 = find_split(returns, split_variables, list_split_variables,model_specification, max_lags)
      # execute split (since only one choosen split)
      
    # second + third split
      #  Simplification: split both nodes once (do not choose best split between nodes)
      # run split function for both nodes. 
      split_2_1 = find_split(split_1$subsample_lower, split_variables, list_split_variables,model_specification, max_lags)
      split_2_2 =  find_split(split_1$subsample_higher, split_variables, list_split_variables,model_specification ,max_lags)
      
    # complete split information table
      split_order =  as.data.frame(matrix(nrow = 4, ncol = 7,))
      colnames(split_order) = c("regime","split_variable1","threshhold1","higher_lower1", "split_variable2","threshhold2","higher_lower2")
      split_order$regime =  paste0("regime", seq(1:4))
      split_order$split_variable1 =split_1$selected_split$var_name
      split_order$threshhold1 = split_1$selected_split$threshhold_split
      split_order$higher_lower1 = c(rep("lower",2),rep("higher",2))
      split_order$split_variable2 =c(rep(split_2_1$selected_split$var_name,2), rep(split_2_2$selected_split$var_name,2))
      split_order$threshhold2 = c(rep(split_2_1$selected_split$threshhold_split,2), rep(split_2_2$selected_split$threshhold_split,2))
      split_order$higher_lower2 = rep(c("lower","higher"),2)
      
  # 3) Prune the tree
    # estimate models for subtrees of tree
    # all possible models to evaluatate
    submodels = list(returns, split_1$subsample_lower, split_1$subsample_higher,
                     split_2_1$subsample_lower, split_2_1$subsample_higher,
                     split_2_2$subsample_lower, split_2_2$subsample_higher)
    loglik_submodels = rep(-10^10, length(submodels))
  
    # estimate models
      ar = 1
      ma = 1
      threshhold = F
      th_value  = 0 # not optimized within fct
      data_threshhold = 0
      distribution ="t"
      start_parms = c(0,0.1,  rep(0.1/ma,ma), rep(0.9/ar,ar)) # initialize parms. 
      if(threshhold==T){
        start_parms=  c(start_parms, 0) # set asymmetry parameter to 0 
      }
      if(distribution=="t"){
        start_parms=  c(start_parms, 6) # keep df_t > 2 due to likelihood fct
      }
      number_restrictions = length(start_parms)
      for (submodels_iter in 1:length(submodels)){
        opt_parms= nlm(garchEstimation,start_parms,
                       returns = submodels[[submodels_iter]]$return[(max_lags+1):nrow(submodels[[submodels_iter]])],  ar = ar, ma = ma,
                       threshhold = threshhold, th_value = th_value, data_threshhold = data_threshhold,
                       distribution=distribution,
                       print.level=0,iterlim=1000, check.analyticals=1)
        loglik_submodels[submodels_iter] = -opt_parms$minimum
      }
      print("loglik_submodels")
      print(loglik_submodels)
    # define all possible subtrees with number of active submodels. number is position in submodels
      possible_subtrees = list(c(1), c(2,3), c(3,4,5), c(2,6,7), c(4,5,6,7))
      samples_groups_in_subtrees = list(c(1:4), c(1:2,3:4),c(1:4),c(1:4),c(1:4),c(1:4),c(1:4))
      
    # aic for all subtrees
      aic_subtrees= rep(10^10,length(possible_subtrees))
      
      for (subtree_iter in 1:length(aic_subtrees)){
        logLiks  = loglik_submodels[possible_subtrees[[subtree_iter]]]
        aic_subtrees[subtree_iter] = my_aic(sum(logLiks), length(logLiks)*number_restrictions) # parms is estimated parm/model * number models (each model has same number parms)
      }
      print("aic_subtrees")
      print(aic_subtrees)
    # select subtree with min aic
      best_subtree_indic = which.min(aic_subtrees)
      best_subtrees = possible_subtrees[[best_subtree_indic]]
      
    # create data of split for prediction
      split_order_pruned = split_order
      
    # add NAs for splitting variable if no split is performed
      if(best_subtree_indic==1) {
        split_order_pruned$split_variable1 = "split_removed"
        split_order_pruned$split_variable2 = "split_removed"
      } else 
        if(best_subtree_indic==2) {
          split_order_pruned$split_variable2 = "split_removed"
        } else
          if(best_subtree_indic==3) {
            split_order_pruned$split_variable2[3:4] = "split_removed"
          } else
            if(best_subtree_indic==4) {
              split_order_pruned$split_variable2[1:2] = "split_removed"
            }
      
  # return split lists
  split_list = list(split_order_pruned, split_order)
  names(split_list) = c("split_order_pruned","split_order") 
  return(split_list)
}      


# function to find optimal split for GARCH 1,1 models
find_split = function(returns, split_variables, list_split_variables, model_specification,max_lags){
  # step 1) find optimal GARCH 1/1 for sample. remove first max_lags obs since they are not used by TreeGarch either
  # inputs fct
  ar = 1
  ma = 1
  threshhold = F
  th_value  = 0 # not optimized within fct
  data_threshhold = 0
  distribution ="t"
  start_parms = c(0,0.1,  rep(0.1/ma,ma), rep(0.9/ar,ar)) # initialize parms. 
  if(threshhold==T){
    start_parms=  c(start_parms, 0) # set asymmetry parameter to 0 
  }
  if(distribution=="t"){
    start_parms=  c(start_parms, 6) # keep df_t > 2 due to likelihood fct
  }
  # estimate basic Garch (1,1) for sample
  opt_parms= nlm(garchEstimation,start_parms,
                 returns = returns$return[(max_lags+1):nrow(returns)],  ar = ar, ma = ma,
                 threshhold = threshhold, th_value = th_value, data_threshhold = data_threshhold,
                 distribution=distribution,
                 print.level=0,iterlim=1000, check.analyticals=1)
  
  # step 2) split sample via reduction in log likelihood
  optimal_split = as.data.frame(matrix(nrow = length(list_split_variables), ncol = 3,))
  colnames(optimal_split) = c("var_name","threshhold_split","improvementLogLik")
  optimal_split$var_name = list_split_variables
  
  # choose variable for split
  for (split_var_iter in 1:length(list_split_variables)) {
    # assign current split variable and remove first missing obs
    split_one_covariate = drop_na(select(split_variables, c(list_split_variables[split_var_iter],"date"))) 
    
    # get return and splitting variable in same dataframe
    split_data = select(left_join(returns, split_one_covariate, by="date"), -c("date"))
    colnames(split_data) = c("return","split_var")
    split_data = split_data[(max_lags+1):nrow(split_data),]
    
    # get grid of threshholds via quantile in subsample
    split_threshholds = quantile(split_data$split_var,vector_quantiles)
    
    # starting values are parms of first GARCH
    start_parms = opt_parms$estimate  
    
    # initalize result table for one covariate
    single_split_criterion_table = as.data.frame(matrix(nrow = length(split_threshholds), ncol = 5,))
    colnames(single_split_criterion_table) = c("threshhold_split","improvementLogLik","logLikSample1","logLikSample2","logLikFullSample")
    single_split_criterion_table$value_split = split_threshholds
    single_split_criterion_table$logLikFullSample = -opt_parms$minimum
    
    
    
    for (i in 1:length(split_threshholds)) {
      # split sample starting from first obs that is not NA for splitting value
      subsample1 = split_data$return[(split_data$split_var < split_threshholds[i])]
      subsample2 = split_data$return[(split_data$split_var >=split_threshholds[i])]
      
      # estimate GARCH in subsamples and get sum of likelihood
      opt_parms1= nlm(garchEstimation,start_parms,
                      returns = subsample1,  ar = ar, ma = ma,
                      threshhold = F, th_value = 0, data_threshhold = data_threshhold,
                      distribution=distribution,
                      iterlim=1000, check.analyticals=1)
      opt_parms2= nlm(garchEstimation,start_parms,
                      returns = subsample2,  ar = ar, ma = ma,
                      threshhold = F, th_value = 0, data_threshhold = data_threshhold,
                      distribution=distribution,
                      iterlim=1000, check.analyticals=1)
      single_split_criterion_table$logLikSample1[i] = -opt_parms1$minimum
      single_split_criterion_table$logLikSample2[i] = -opt_parms2$minimum
    }
    
    # calc improvement in loglik and save best split
    single_split_criterion_table$improvementLogLik = single_split_criterion_table$logLikSample1+single_split_criterion_table$logLikSample2-single_split_criterion_table$logLikFullSample
    indic_max_improvement = which.max(single_split_criterion_table$improvementLogLik)
    optimal_split[split_var_iter,2:3] = c(single_split_criterion_table[indic_max_improvement,c("value_split","improvementLogLik")])
  }
  
  # pick optimal covariate and threshhold to split
  indic_optimal_covariate_and_split = which.max(optimal_split$improvementLogLik)
  selected_split  =optimal_split[indic_optimal_covariate_and_split,]

  # If no split is found: return old data, do n
  if(selected_split$improvementLogLik <=0) {
    print("no improving split found")
    selected_split$var_name="no_split_found"
    result = list(returns, selected_split) # return back old data and store in list
    names(result) = c("input_returns", "selected_split")
    
  } else {
    
    # If split is found: add new subsamples to output
    print("improving split found")
    
    split_data = select(left_join(returns, split_one_covariate, by="date"), -c("date"))
    colnames(split_data) = c("return","split_var")
    split_data = split_data[(max_lags+1):nrow(split_data),]
    
    final_subsample_split =  left_join(returns, split_variables, by="date") %>%
      select(return, selected_split$var_name, date)
    
    colnames(final_subsample_split) =c("return","opt_split_var","date") # set colnames because comparision doesnt work with referencing
    
    final_subsample_split1 = filter(final_subsample_split, opt_split_var < selected_split$threshhold_split) %>%
      select(return, date)# lower than threshhold
    final_subsample_split2 = filter(final_subsample_split, opt_split_var >= selected_split$threshhold_split) %>%
      select(return, date)# equal or higher than threshhold
    
    # create list and get names 
    result = list(returns, selected_split, final_subsample_split1, final_subsample_split2) # return back split sample
    print("selected_split")
    print(selected_split)
    names(result) = c("input_returns", "selected_split", "subsample_lower", "subsample_higher")
  }
  return(result)
}


# flexible garch estimation function for all orders, threshhold and distributions 

garchEstimation = function(theta, returns, ma,  ar,threshhold,th_value,data_threshhold, distribution) { 
  {
    # Inputs
      # theta: AR-coefs, MA-coefs, Treshhold, df of t-distribution
      # returns: univariate return vector, double
      # ar: order AR process for squared returns
      # ma: order MA process for squared returns
      # threshhold: T/F. if T, then # cols of data_treshhold is number of threshhold parameters. if F then inactive
      # data_treshhold: only evaluated if threshhold = T. ATM not used
      # distribution: normal / t
    
    # Outputs
      # likelihood of sample given the model specification
      # !parameters square root! (when using an optimizer)
    
    # conditions
      # distribution
      if (distribution %in% c("normal", "t") ==F ){
        print("Error: non-supported distribution type in garchEstimation")
      } 
      # number of parameters to be estimated equals dim of theta
      if (length(theta) != (1+1 + ma+ ar+ as.numeric(threshhold) + as.numeric(distribution=="t"))) { # number of parms: mean return+ constant + ar +ma + threshhold_parameter (if active) + degrees of freedom t (if active)
        print("Error: Number of input parameters does not match length of parameter vector in garchEstimation")
      } 
    
    
      # assign coefficients:
      mu_coef = theta[1]
      constant_coef = theta[2]
      ma_coef = theta[(1:ma)+2]
      ar_coef = theta[((ma+1):(ma+ar))+2]
      th_coef = ifelse(threshhold==T, theta[(ar+ma+1)+2],NA)
      df_t_coef = ifelse(distribution=="t", theta[(ar+ma+1+as.numeric(threshhold))+2],NA)
    
      # number of timesteps, starting values
      max_lags = max(ar, ma)
      n=length(returns)
      x.start= mean(returns)
      sigmasq.start= var(returns)
      
      # define returns
      data=c(x.start,returns)
      
      # initialize variance as sample variance
      sigmasq= rep(0,n+1)
      sigmasq[1]=sigmasq.start # initialize with unconditional variance
      
      my.sigma=c(sqrt(sigmasq[1]),rep(0,n))
      
      mean_ret=rep(mu_coef,n+1) # constant conditional mean
      # mean_ret[1] = 0 
      ar_part = rep(0, n+1) # initialize AR
      ma_part = rep(0, n+1) # initialize MA
      th_active_iteration =  rep(0, n+1) # initialize th_part (0 if not active)
      
      # calculate ma part. loop for all ma parts
      for (k in 1:length(ma_coef)) {
        ma_part[(1+k):(n+1)] = ma_part[(1+k):(n+1)] + ma_coef[k]^2*(data[1:(n+1-k)]-mean_ret[1:(n+1-k)])^2 # calc MA part
      }

      # threshhold: th_value as input parameter is not optimized. Parm "active" when epsilon < th_value
      if (threshhold==T)  {
        th_active_iteration[2:(n+1)] = ((data[1:n]-mean_ret[1:n])<=th_value) * ((data[1:n]-mean_ret[1:n])^2)*th_coef
      }

      # loop for time: calc AR and conditional variance
      for (i in (max_lags+1):(n+1))
      {
        # calculate AR parts. loop for all AR parts
        for (j in 1:length(ar_coef)) {
          ar_part[i] = ar_part[i] + ar_coef[j]^2*sigmasq[i-j]  # calc AR part
        }
        # ar_part[i]  = theta[3]^2*sigmasq[i-1] # iteratively add ar part
        sigmasq[i]  = constant_coef^2 + ar_part[i] + ma_part[i] +th_active_iteration[i] # conditional variance. force constant to be positive
      }
    
    # calc log likelihood of model
      if (distribution=="normal"){
        log_liklihood = (n+1-max_lags)*log(sqrt(2*pi))+sum(0.5*((data[(max_lags+1):(n+1)]-mean_ret[(max_lags+1):(n+1)])^2)/sigmasq[(max_lags+1):(n+1)]) + sum(0.5*log(sigmasq[(max_lags+1):(n+1)]))
        return(log_liklihood)
      } else  if (distribution=="t"){
        log_liklihood = -(n+1-max_lags)*log(gamma((df_t_coef+1)/2)/(gamma(df_t_coef/2)*sqrt(pi*(df_t_coef-2)))) + (df_t_coef+1)/2*sum(log(1+((data[(max_lags+1):(n+1)]-mean_ret[(max_lags+1):(n+1)])^2)/((df_t_coef-2)*sigmasq[(max_lags+1):(n+1)]))) + 0.5*sum(log(sigmasq[(max_lags+1):(n+1)]))
        return(log_liklihood)
      }
      
    # audrino implementation likelihood to compare. Delivers same results
      # normal
           # 1/2*sum(log(sigmasq[2:(n+1)])) - sum(log(dnorm((data[2:(n+1)]-mean_ret[2:(n+1)])/sqrt(sigmasq[2:(n+1)])))) # audrino code
      #tdistrib 
           # 1/2*sum(log(sigmasq[2:(n+1)]*(theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n+1)]-mean_ret[2:(n+1)])/sqrt(sigmasq[2:(n+1)]*(theta[6]-2)/theta[6]),df=theta[6])))+10^(10)*(theta[6]<2)+10^(10)*(theta[6]>10) # audrino code
  }
}


my.loglike.t=function(theta) #Estimate an asymmetric GARCH(1,1) model with Student's t innovations
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
  
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start
  
  my.sigma=c(sqrt(my.sigmasq[1]),rep(0,n))
  
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }
  
  for (i in 2:(n+1))
  {
    my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)
    
    # my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    # my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]*(data[i-1]-my.mean[i-1])^2*((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }
  
  #my.sigmasq=my.sigma^2cd
  #my.sigmasq=exp(log.sigmasq)
  
  # normdistrib, GARCH 1/1
  1/2*sum(log(my.sigmasq[2:(n+1)])) - sum(log(dnorm((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)]))))
  
  #tdistrib
  #1/2*sum(log(my.sigmasq[2:(n+1)]*(theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)]*(theta[6]-2)/theta[6]),df=theta[6])))+10^(10)*(theta[6]<2)+10^(10)*(theta[6]>10)
  # return(-n*log(gamma((theta[6]+1)/2)/(gamma(theta[6]/2)*sqrt(pi*(theta[6]-2)))) + (theta[6]+1)/2*sum(log(1+((data[2:(n+1)]-my.mean[2:(n+1)])^2)/((theta[6]-2)*my.sigmasq[2:(n+1)]))) + 0.5*sum(log(my.sigmasq[2:(n+1)])))
}
