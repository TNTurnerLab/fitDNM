cgf_0 = function(mu,s,lambdas,weights){
	
	nloci = length(lambdas)
	##check if length of lambdas is the same as that of weights
	if(nloci != length(weights)){
		stop("check the length of weights and lambdas")
	}
	result = 0
	
	
	for(i in c(1:nloci)){
		result = result + lambdas[i] *(exp(weights[i] * mu + s) -1)		
	}	
	return(result)
		
}

cgf_1 =  function(mu,s,lambdas,weights){
	nloci = length(lambdas)
	##check if length of lambdas is the same as that of weights
	if(nloci != length(weights)){
		stop("check the length of weights and lambdas")
	}
	result_s = 0
	result_mu = 0
	
	for(i in c(1:nloci)){
		result_s = result_s + lambdas[i] *(exp(weights[i] * mu + s))
		result_mu = result_mu + weights[i] * lambdas[i] *(exp(weights[i] * mu + s))		
	}	
	return(c(result_s, result_mu))
		
}


cgf_2 =  function(mu,s,lambdas,weights){
	nloci = length(lambdas)
	##check if length of lambdas is the same as that of weights
	if(nloci != length(weights)){
		stop("check the length of weights and lambdas")
	}
	result_ss = 0
	result_mumu = 0
	result_mus = 0

	
	for(i in c(1:nloci)){
		result_ss = result_ss + lambdas[i] *(exp(weights[i] * mu + s))
		result_mumu = result_mumu + weights[i]^2 * lambdas[i] *(exp(weights[i] * mu + s))
		result_mus = result_mus + weights[i] * lambdas[i] *(exp(weights[i] * mu + s))

		
				
	}	
	#return(c(result_ss, result_mus, result_mumu))
	return(	result_ss* result_mumu - (result_mus^2))
}

solve_s_u = function(x,y,lambdas, weights,delta =10,refine=5,start=0){
	##solve mu first
    tmp = y/x
	mu.old = start
	tmp1=  tmp - (sum(lambdas *weights*exp(weights * mu.old)) / sum(lambdas*exp(weights * mu.old)))
	sign.old = tmp1/abs(tmp1)
	stop=0
	while(stop == 0){
		
		mu.new = mu.old + sign.old * delta
		tmp2=  tmp - (sum(lambdas *weights*exp(weights * mu.new)) / sum(lambdas*exp(weights * mu.new)))
		sign.new = tmp2/abs(tmp2)
		
		if(is.na(sign.new)){
			stop=2
		}else if((sign.old * sign.new) <=0){
			stop=1
		}else{
			mu.old = mu.new
			sign.old = sign.new
		}
			
	}
	
	##refine
	
	if(stop==1){
		for(i in c(1:refine)){
			delta = delta/10
			stop=0
			while(stop == 0){
			
				mu.new = mu.old + sign.old * delta
				tmp2=  tmp - (sum(lambdas *weights*exp(weights * mu.new)) / sum(lambdas*exp(weights * mu.new)))
				sign.new = tmp2/abs(tmp2)
		
				if((sign.old * sign.new) <=0){
					stop=1
				}else{
					mu.old = mu.new
					sign.old = sign.new
				}
			}	
		}

	}		
	
	mu.est = (mu.old+mu.new)/2
    s.est = log(x) - log((sum(lambdas*exp(weights * mu.est))))
    if(s.est == Inf){return(NA)}
	return(c(s.est, mu.est))
}


conditional_approximation = function(x,y,lambdas, weights){
	
    ##solve s0
	s0 = log(x/(sum(lambdas)))
	## solve s and mu
	tmp = solve_s_u(x,y,lambdas, weights)
	s = tmp[1]
	mu = tmp[2]
	
	## compute w
	w.part = 2 * (s*x+mu*y-cgf_0(mu=mu,s=s,lambdas= lambdas,weights= weights) - s0*x+cgf_0(mu=0,s=s0,lambdas= lambdas,weights= weights))
	
	## make sure w.part >=0
	refine.delta = 10
	iter = 0
	while(w.part <0 & iter<100){
		iter = iter+1
		tmp = solve_s_u(x,y,lambdas, weights,refine=refine.delta*iter)
		s = tmp[1]
		mu = tmp[2]
		w.part = 2 * (s*x+mu*y-cgf_0(mu=mu,s=s,lambdas= lambdas,weights= weights) - s0*x+cgf_0(mu=0,s=s0,lambdas= lambdas,weights= weights))
	}
	if(iter >=100){
		return(NA)
	}
	
	if(abs(mu) <=10^(-4)){
		return(NA)
	}

	w = mu/abs(mu) * sqrt(w.part)
	
	##compute K_ss(s0,0)
	K_ss_s0 = exp(s0) *sum(lambdas)
	##compute |K_2_(s,mu)|
	K_2_smu = cgf_2(mu=mu,s=s,lambdas,weights)
	
	result = 1 - pnorm(w)+dnorm(w)*(sqrt(K_ss_s0/(K_2_smu))/mu - 1/w)
	
	return(result)
}

conditional_approximation_2 = function(x,y,lambdas, weights){
	
	##solve s0
	s0 = log(x/(sum(lambdas)))
	## solve s and mu
    tmp = solve_s_u(x,y,lambdas, weights)
	s = tmp[1]
	mu = tmp[2]
	
	## compute w
	w.part = 2 * (s*x+mu*y-cgf_0(mu=mu,s=s,lambdas= lambdas,weights= weights) - s0*x+cgf_0(mu=0,s=s0,lambdas= lambdas,weights= weights))
	
    if(is.na(w.part)){return(NA)}
	## make sure w.part >=0
	refine.delta = 10
	iter = 0
    while(w.part <0 & iter<10){
        iter = iter+1
		tmp = solve_s_u(x,y,lambdas, weights,refine=refine.delta*iter)
		s = tmp[1]
		mu = tmp[2]
		w.part = 2 * (s*x+mu*y-cgf_0(mu=mu,s=s,lambdas= lambdas,weights= weights) - s0*x+cgf_0(mu=0,s=s0,lambdas= lambdas,weights= weights))
	}
	if(iter >=100){
		return(NA)
	}
	
	if(abs(w.part) <=10^(-4)){
		stop2=0
		while(stop2 ==0){
			y = y + 0.01
			tmp = solve_s_u(x,y,lambdas, weights)
			s = tmp[1]
			mu = tmp[2]
			w.part2 = 2 * (s*x+mu*y-cgf_0(mu=mu,s=s,lambdas= lambdas,weights= weights) - s0*x+cgf_0(mu=0,s=s0,lambdas= lambdas,weights= weights))
			if(abs(w.part2)>10^(-4)){
				stop2=1
			}
			
		}
		w = mu/abs(mu) * sqrt(w.part2)
	
		##compute K_ss(s0,0)
		K_ss_s0 = exp(s0) *sum(lambdas)
		##compute |K_2_(s,mu)|
		K_2_smu = cgf_2(mu=mu,s=s,lambdas,weights)
	
		result1 = 1 - pnorm(w)+dnorm(w)*(sqrt(K_ss_s0/(K_2_smu))/mu - 1/w)
		
		stop3=0
		while(stop3 ==0){
			y = y - 0.01
			tmp = solve_s_u(x,y,lambdas, weights)
			s = tmp[1]
			mu = tmp[2]
			w.part2 = 2 * (s*x+mu*y-cgf_0(mu=mu,s=s,lambdas= lambdas,weights= weights) - s0*x+cgf_0(mu=0,s=s0,lambdas= lambdas,weights= weights))
			if(abs(w.part2)>10^(-4)){
				stop3=1
			}
			
		}
		w = mu/abs(mu) * sqrt(w.part2)
	
		##compute K_ss(s0,0)
		K_ss_s0 = exp(s0) *sum(lambdas)
		##compute |K_2_(s,mu)|
		K_2_smu = cgf_2(mu=mu,s=s,lambdas,weights)
	
		result2 = 1 - pnorm(w)+dnorm(w)*(sqrt(K_ss_s0/(K_2_smu))/mu - 1/w)
		
		return((result1+result2)/2)
		
		
		
	}else{

		w = mu/abs(mu) * sqrt(w.part)
	
		##compute K_ss(s0,0)
		K_ss_s0 = exp(s0) *sum(lambdas)
		##compute |K_2_(s,mu)|
		K_2_smu = cgf_2(mu=mu,s=s,lambdas,weights)
	
		result = 1 - pnorm(w)+dnorm(w)*(sqrt(K_ss_s0/(K_2_smu))/mu - 1/w)
	
		return(result)
	}	
}

double_saddle_point_approximation = function(y,lambdas, weights){
	
    tmp_result = 0;
	tmp_result2 = 0;
	sum.lam = sum(lambdas)
	
	max_weights = max(weights,na.rm=T)
	start = ceiling(y/max_weights)
	if(start <=0){
		return(1)
#	}else if(start==1){
#		tmp_result = sum(lambdas[weights>=y])/sum.lam *dpois(start,sum(lambdas))
	}else{
        tmp_result = conditional_approximation_2(x=start,y=y,lambdas = lambdas,weights=weights) *dpois(start,sum(lambdas))
	}

    if(is.na(tmp_result)){
		return(NA)
	}
    for(i in c((start+1):100)){
			tmp_result2 = conditional_approximation_2(x=i,y=y,lambdas = lambdas,weights=weights) *dpois(i,sum(lambdas))
			
			if(	is.na(tmp_result2) || tmp_result2 == -Inf){
				return(NA)
			}else if(abs(tmp_result2/tmp_result) < 0.00001){
				tmp_result = tmp_result + tmp_result2
				break
			}else{
				tmp_result = tmp_result + tmp_result2
			}
			
	}
	
	return(tmp_result)
	
}










