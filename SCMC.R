library(abind)

# initial sampling on the hypercube
unrestricted = function(N, range) { # initial sample on the hypercube
  D = nrow(range)
  samp = NULL
  for (d in 1:D) samp = cbind(samp, runif(N, range[d,1], range[d,2]))
  return(samp)
}

# log-posterior / equivalent to the constraint if a uniform sample is drawn
logpost = function(sample, nu_t) { 
  term = constraint(sample)
  return(sum(pnorm(- term / nu_t, log = T)))
}

# sampling at each time step
Gibbs = function(x, q, d, a, nu_t, lpdent = lpdent) { 
  delta = rnorm(1,0,q[d])
  newx = x
  newx[d] = newx[d] + delta
  lpnum = logpost(newx, nu_t)
  ratio = lpnum - lpdent
  prob = min(1, exp(ratio))
  u = runif(1)
  if (u <= prob) {
    x = newx
    lpdent = lpnum
    a = a + 1
  }
  return(list(x = x, a = a, lpdent = lpdent))
}

# adaptive specification of the constraint parameter
adapt_seq = function(nu, nu0, t, N, sample, Wt) { 
  wt = c()
  for (i in 1:N) {
    term = constraint(sample[i,])
    cons1 = sum(pnorm(- term / nu, log = T))
    cons2 = sum(pnorm(- term / nu0, log = T))
    wt[i] = cons1 - cons2
  }
  Wt = Wt * exp(wt)
  Wt = Wt / sum(Wt)
  ESS = ifelse(sum(is.na(Wt)) == N, 0, 1 / sum(Wt ^ 2))
  return(ESS - (N / 2))
}

# main function: sampling from a constraint set
SCMC = function(N=1000, M=10, L=25, nuseq_T= 1e-3, lower=c(-1,-1), upper=c(1,1), qt = 1) {
  rge = cbind(lower,upper)
  D = nrow(rge)
  t = 1
  ESS = NULL
  nuseq = c(Inf)
  b = seq(1.5,.1,length = L)
  nuseq_0 = c(Inf, b^7)
  a = array(0,dim=c(L, D))
  samplet = array(dim=c(N, 1, D))
  lpdent = array(dim=c(N, 1))
  Wt = array(dim=c(N, 1))
  samplet[,1,] = unrestricted(N, range=rge)
  for (i in 1:N) lpdent[i, 1] = logpost(samplet[i, t,], nuseq[t])
  Wt[,1] = rep(1 / N, N)
  repeat {
    t = t+1
    newsample = samplet[, t-1,]
    newlpdent = lpdent[, t- 1]
    newWt = Wt[, t - 1]
    nuseq[t] = ifelse(adapt_seq(nu = nuseq_T, nu0 = nuseq[t-1], t = t, N = N, sample = newsample, Wt = newWt) > 0, nuseq_T, uniroot(adapt_seq, interval = c(ifelse(t < 3, min(.1, nuseq_0[t]), nuseq_T), ifelse(t != 2, nuseq[t-1], 1e5)), nu0 = nuseq[t-1], t = t, N = N, sample = newsample, Wt = newWt)$root)
    wt = c()
    for (i in 1:N) {
      term = constraint(newsample[i,])
      constraint1 = sum(pnorm(- term / nuseq[t], log = T)) 
      constraint2 = sum(pnorm(- term / nuseq[t-1], log = T))
      wt[i] = constraint1 - constraint2
    }
    newWt = newWt * exp(wt)
    newWt = newWt / sum(newWt)
    ESS[t] = 1 / sum(newWt ^ 2)
    index = sample(1:N, N, prob = newWt, replace = T)
    newsample = newsample[index,] 
    newlpdent = newlpdent[index]
    newWt = rep(1/N, N)
    q = apply(newsample, 2, sd) * qt
    a = abind(a, rep(0, D), along = 1)
    for(j in 1:M) {
      for (i in 1:N) {
        for (d in 1:D) {
          out = Gibbs(newsample[i,], q, d, a[t, d], nuseq[t], lpdent = newlpdent[i])
          newsample[i,] = out$x
          a[t,d] = out$a
          newlpdent[i] = out$lpdent
        }
      }
    }
    samplet = abind(samplet, newsample, along = 2)
    lpdent = abind(lpdent,newlpdent,along = 2)
    Wt = abind(Wt, newWt, along = 2)
    if (nuseq[t] <= nuseq_T) break
  }
  t_final = t
  sample = samplet[, t_final,]
  return(sample)
}

# constrained set example 1
crecsent = function(sample) { 
  out = c(sample[1] - sqrt(14 * sample[2] ^ 2 + 2), - (sample[1] - sqrt(33*sample[2] ^ 2 + 1)))
  return(out)
}

constraint = crecsent
l=c(-1,-1)
u=c(1,1)
sample0 = SCMC(N=1000, M=10, L=25, nuseq_T= 1e-3, lower=l, upper=u, qt = 1)
plot(sample0)

# constrained set example 2
mixture = function(sample, l = c(.4,.1,.1,.03), u = c(.6, .47, .47, .08)) {
  return(c(abs(sum(sample) - 1), l - sample, sample - u))
}

constraint=mixture
epsilon=1e-3
l=c(.4,.1,.1,.03) - 1e-3
u=c(.6,.47,.47,.08) + 1e-3
sample0 = SCMC( N=1000, M=10, L=25, nuseq_T= 1e-3, lower=l, upper=u, qt = 1)

