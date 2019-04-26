bootstrap_replicate_1d<-function(data, func){
  bs_sample = sample(data, size=length(data), replace=TRUE)
  do.call(func, list(x=bs_sample))
}

generate_bootstrap_replicate<-function(data, func, n=10000){
  replicate(n, do.call(bootstrap_replicate_1d, list(data=data, func=func)))
}


calculate_bootstrap_ci<-function(data, func, n=10000, as_string=FALSE){
  values = generate_bootstrap_replicate(data, func, n)
  rtnVal = quantile(values, c(.025, .5,  .975))
  if (as_string){
    rtnVal = paste(formatC(rtnVal[2], big.mark=','), ' (', 
                   formatC(rtnVal[1], big.mark=','), 
                   '-', 
                   formatC(rtnVal[3], big.mark=','), 
                   ' 95% CI)', sep='')
  }
  return(rtnVal)
}

permutation_sample<-function(data_1, data_2){
  combined = c(data_1, data_2)
  shuffled = sample(combined)
  list(d1 = shuffled[1:length(data_1)], 
       d2=shuffled[(length(data_1)+1):(length(shuffled))])
  
}

calculate_difference<-function(data_1, data_2, func){
  do.call(func, list(x=data_1)) - do.call(func, list(x=data_2))
}

permute_difference<-function(data_1, data_2, func){
  permutation = permutation_sample(data_1, data_2)
  calculate_difference(permutation[[1]], permutation[[2]], func)
}

bootstrap_nhst_two_populations<-function(data_1, data_2, func, n=10000){
  reps = replicate(n, do.call(permute_difference, 
                       list(data_1=data_1, 
                                  data_2=data_2, 
                                  func=func)))
  
  actual_difference = calculate_difference(data_1, data_2, func)
  
  p_value = sum(reps > actual_difference) / length(reps)
  
  if(actual_difference - mean(reps) > 0){
    count = sum(reps > actual_difference)
  } else{
    count = sum(reps < actual_difference)
  }
  
  
  min_value = min(c(reps, actual_difference))
  min_value = min_value - (0.1*abs(min_value))
  
  max_value = max(c(reps, actual_difference))
  max_value = max_value + (0.1*abs(max_value))
  
  plot(density(reps), 
       xlim=c(min_value, max_value),
       main='Permuted differences between datasets\nfor null hypothesis')
  abline(v=actual_difference, col='red')
  
  print(paste("A value at least as extreme as the empirical difference between the two datasets was observed ", 
              count, " times out of ", length(reps), " (", count/length(reps)*100, "%)",
              " of the time in the null model", sep=''))
}


bootstrap_difference<-function(data_1, data_2, func, n=1){
  generate_bootstrap_replicate(data_1, func, n) - generate_bootstrap_replicate(data_2, func, n)
}

N50<-function(x, as_string=FALSE){
  tmp <- rev(sort(x))
  tmp2 <- cumsum(tmp) <= sum(x)/2
  rtnValue = tmp[sum(tmp2)]
  if (as_string){
    rtnValue=formatC(rtnValue, big.mark=',')
  }
  return(rtnValue)
}
