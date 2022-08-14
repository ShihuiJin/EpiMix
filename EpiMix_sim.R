#INLA code for simulation
#t1=Sys.time()
libray(INLA)
y=I[start:end] #case counts, estimation window: [start:end]
Infect=rep(0,end-start+1)
w=dlnorm(seq(1,50),1.132,0.742)
w=w/sum(w)
for(i in start:end)
{
  if(i>50){
    Infect[i-start+1]=sum(rev(w)*I[(i-50):(i-1)])
  }else{
    Infect[i-start+1]=sum(rev(w[1:(i-1)])*I[1:(i-1)])
  }
}

data=data.frame(cbind(y,Infect,X_sd[start:end-1,])) #X_sd: standardized X matrix for prediction
if(ncol(X_sd)==1) colnames(data)[3]=colnames(X_sd)='V1'
data$idx = 1:nrow(data)
f1 <- as.formula(paste0("y ~ 1 + f(idx, model='iid') + ", paste0(colnames(X_sd), collapse = " + ")))
ll=list()
for(i in 1:ncol(X_sd))
{
  ll[[colnames(X_sd)[i]]]=X_sd[start:end-1,i]
}
ll[["(Intercept)"]]=rep(1,end-start+1)
lc=inla.make.lincombs(ll)
model=inla(f1,data=data,family = "poisson",control.family=list(link='log'), E=Infect, lincomb =lc,control.compute=list(dic=TRUE, cpo=TRUE, return.marginals=TRUE, return.marginals.predictor=TRUE, config = TRUE),control.predictor = list(compute=TRUE))
#print(Sys.time()-t1)
