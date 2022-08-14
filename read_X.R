#read data for 2021 full-year version , including X matrix and cases 
if(region=='SG'){
  mobility_2020=read.csv(paste0('2020_',region,'_Region_Mobility_Report.csv'))
  mobility_2021=read.csv(paste0('2021_',region,'_Region_Mobility_Report.csv'))
}else if(region=='NZ'){
  mobility_2020=read.csv(paste0('2020_',region,'_Region_Mobility_Report.csv'))[1:321,]
  mobility_2021=read.csv(paste0('2021_',region,'_Region_Mobility_Report.csv'))[1:365,]
}

m=c(31,28,31,30,31,30,31,31,30,31,30,31)
m=cumsum(m)
#vaccination
vax=read.csv(paste0('vaccination_',region,'.csv'))
#2020 part
X_2020=mobility_2020[(ncol(mobility_2020)-5):ncol(mobility_2020)]
X_2020$vaccination=rep(0,nrow(X_2020))
X_2020$delta=rep(0,nrow(X_2020))
X_2021=mobility_2021[(ncol(mobility_2021)-5):ncol(mobility_2021)]
X_2021$vaccination=rep(0,nrow(X_2021))
if(region=='SG'){
  #p_vax=rev(taRifx::destring(vax$per.100.people))
  p_vax=rev(vax$per.100.people)
  X_2021$vaccination=p_vax[1:nrow(X_2021)]
}
if(region=="NZ"){
  p_vax=rev(taRifx::destring(vax$per.100.people)) #vaccination doses per 100 people
  X_2021$vaccination[(17+m[1]):nrow(X_2021)]=p_vax[1:(nrow(X_2021)-m[1]-16)] #for NZ, vaccination data from 17/02,2021
}

delta_count=as.matrix(read.csv(paste0('delta_count_',region,'.csv'))) #week data, starting from week 0
delta_day=rep(0,nrow(X_2021))
for(i in 1:nrow(X_2021))
{
  k=ceiling((i-2)/7)+1
  if(k<=nrow(delta_count)){
    delta_day[i]=delta_count[k,1]/delta_count[k,2]*100
    if(is.na(delta_day[i])==1)delta_day[i]=0
  }else{
    delta_day[i]=100
  }
}
X_2021$delta=delta_day

X=rbind(as.matrix(X_2020),as.matrix(X_2021[,(ncol(X_2021)-ncol(X_2020)+1):ncol(X_2021)]))

oxford_policy=as.matrix(read.csv(paste0('Oxford_new_',region,'.csv'))) #20 variables
for(i in 2:ncol(oxford_policy))
{
  X=cbind(X,taRifx::destring(oxford_policy[,i]))
}
colnames(X)[-c(1:8)]=colnames(oxford_policy)[-1]
X_sd=X

#normalizing for the 2021 part
for(i in 1:ncol(X))
{
  if(var(X[,i])>0){
    X_sd[,i]=(X[,i]-mean(X[nrow(X_2020)+1:nrow(X_2021),i]))/sqrt(var(X[nrow(X_2020)+1:nrow(X_2021),i]))
  }else{
    X_sd[,i]=0
  }
}

#variance
factor=rep(0,ncol(X))
for(i in 1:length(factor))
{
  factor[i]=sqrt(var(X[nrow(X_2020)+1:nrow(X_2021),i]))
}

X_sd=X_sd[,which(factor!=0)]

rm(delta_count,oxford_policy,vax,X_2020,X_2021)
rm(delta_day,i,k,m,p_vax)
#rm(mobility_2020,mobility_2021)

#read cases
if(region=='SG'){
  covid_2020=read.csv('covid_cases_breakdown_2020.csv')
  covid_2021=read.csv('covid_cases_breakdown_2021.csv')
  I=c(covid_2020$Total.Community[(nrow(covid_2020)-nrow(mobility_2020)+1):nrow(covid_2020)], taRifx::destring(covid_2021$Community))
  I_imp=c(covid_2020$Imported[(nrow(covid_2020)-nrow(mobility_2020)+1):nrow(covid_2020)], covid_2021$Imported)
  rm(covid_2020,covid_2021)
}else if(region=='NZ'){
  covid=read.csv(paste0('case_',region,'.csv'))
  case=rep(0,nrow(covid)+27+31)
  case[(28+31):(length(case))]=covid$local
  I=case[(366-nrow(mobility_2020)+1): (366+nrow(mobility_2021))]
  case[(28+31):(length(case))]=covid$imported
  I_imp=case[(366-nrow(mobility_2020)+1): (366+nrow(mobility_2021))]
  remove(case,covid)
}

