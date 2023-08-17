source("pankmer_reading.R")

get_fitting = function(sp, program, start_fit=2, ...){
    ret = read_avg_new(sp, program, ...)

    N = max(ret$x)
        
    pos = ret$x >= start_fit
    retopt = heapslaw(ret$x, ret$y, ret$p0)
    retlm = lmheapslaw(ret$x[pos],ret$y[pos])  
    #retlm_norm = lmheapslaw_norm(ret$x[pos], ret$y[pos], sp)
    #return(list(retopt=retopt, retlm=retlm, retlm_norm=retlm_norm))
    return(list(retopt=retopt, retlm=retlm))
}

get_fitting_gamma = function(sp, program, start_fit, ...){
    ret = read_avg_tot(sp, program, ...)

    N = max(ret$x)
        
    pos = ret$x >= start_fit
    retopt = heapslaw_gamma(ret$x, ret$y, ret$p0)
    retlm = lmheapslaw_gamma(ret$x[pos],ret$y[pos])  
    return(list(retopt=retopt, retlm=retlm))
}

get_alphas = function(program_vec=c("kmer","pantools_mat","roary_mat","bpga_mat"), max_start_fit=10, ...) {
    results = as.data.frame(matrix(nrow=0,ncol=8))
    for (program in program_vec) {
        for (start_fit in 2:max_start_fit) { 
        #for (start_fit in 2:10) { 
            for (i in 1:nrow(species_df())) {
                sp = species_df()[i, "species"]
                cat(sp, start_fit, "...")
                N = species_df()[i, "N"]
                if (N - start_fit >= 5) {
                    fit = get_fitting(sp, program, start_fit = start_fit, ...)
                    results = rbind(results, c(program, sp, start_fit, unlist(fit)))
                } else {
                    results = rbind(results, c(program, sp, start_fit, rep(NA,ncol(results)-3)))
                }
                cat("ok\n")
            }
        }
        colnames(results) = c("program", "species", "start_fit", "opt.k1", "opt.k2", "lm.k1","lm.k2","lm.R2")
    }
    results$species = as.factor(results$species)
    results$program = as.factor(results$program)
    results$lm.k2 = sapply(results$lm.k2,as.numeric)
    results$lm.R2 = sapply(results$lm.R2,as.numeric)
    results$start_fit = sapply(results$start_fit,as.numeric)
    results = results[order(results$species,results$program), ]
    return(results)
}

get_alpha <- function (start, N, sp, type, start_fit=2, ...)  {
    x = 1:N
    res = get_fitting(sp, type, start_fit, ...)
    k1 = res$retlm$k1
    k2 = res$retlm$k2
    yalpha = k1*x^(-k2)
    #cat(sp, type, k1,paste0("x^(-",k2,")"),"\n")
    heap = data.frame(x=x,y=yalpha,species=rep(sp,N))
    heap[start:N,]
}

heapslaw <- function (x,y,p0) {
    objectFun <- function(p, x, y) {
        y.hat <- p[1] * x^(-p[2])
        J <- sqrt(sum((y - y.hat)^2))/length(x)
        return(J)
    }

    fit <- optim(p0, objectFun, gr = NULL, x, y, method = "L-BFGS-B", lower = c(0, 0), upper = c(p0[1]*3, 2))

    return(list(k1=fit$par[1], k2=fit$par[2]))
}

lmheapslaw <- function (x,y) {
    y[y==0] = 0.0001
    N = max(x)
    Y = matrix(y, ncol=(max(x)-min(x)+1), byrow=T)
    Y = colMeans(Y)
    #Y = matrixStats::colMedians(Y)
    X = min(x):max(x)
    model = lm(log(Y)~log(X))
    return(list(k1=exp(model$coefficients[1]),
                k2=abs(model$coefficients[2]), 
                rsquared=summary(model)$adj.r.squared))
}

heapslaw_gamma <- function (x,y,p0) {
    objectFun <- function(p, x, y) {
      y.hat <- p[1] * x^(p[2])
      J <- sqrt(sum((y - y.hat)^2))/length(x)
      return(J)
    }

    fit <- optim(p0, objectFun, gr = NULL, x, y, method = "L-BFGS-B", lower = c(0, 0), upper = c(p0[1]*3, 2))

    return(list(k1=fit$par[1], k2=fit$par[2]))
}

lmheapslaw_gamma <- function (x,y) {
    y[y==0] = 0.0001
    N = max(x)
    Y = matrix(y, ncol=(max(x)-min(x)+1), byrow=T)
    Y = colMeans(Y)
    #Y = matrixStats::colMedians(Y)
    X = min(x):max(x)
    model = lm(log(Y)~log(X))
    return(list(k1=exp(model$coefficients[1]),
                k2=abs(model$coefficients[2]), 
                rsquared=summary(model)$adj.r.squared))
}

get_gamma <- function (start, N, sp, type, start_fit=1, ...)  {
    x = 1:N
    res = get_fitting_gamma(sp, type, start_fit, ...)
    k1 = res$retlm$k1
    k2 = res$retlm$k2
    yalpha = k1*x^k2
    #cat(sp, type, k1,paste0("x^(",k2,")"),"\n")
    heap = data.frame(x=x,y=yalpha,species=rep(sp,N))
    heap[start:N,]
}


# Q: What happens if I normalized by genome size?
# A: Results of the alpha do not change only the values for k1 change
#lmheapslaw_norm <- function (x,y,sp) {
#    avg_len = unlist(species_df() %>% filter(species==sp) %>% select(avg_sum_len))
#    y[y==0] = 0.0001
#    N = max(x)
#    Y = matrix(y, ncol=(max(x)-min(x)+1), byrow=T)
#    Y = colMeans(Y)
#    Y = Y / avg_len
#    #Y = matrixStats::colMedians(Y)
#    X = min(x):max(x)
#    model = lm(log(Y)~log(X))
#    return(list(k1=exp(model$coefficients[1]),
#                k2=abs(model$coefficients[2]), 
#                rsquared=summary(model)$adj.r.squared))
#}

