library(ncvreg)
#### importing the data ####
data=read.csv("riboflavin.csv",header=T,row.names = 1)
data=t(data)
dim(data)
y=data[,1]
X=data[,-1]



####data pre-processing####
var.col = apply(X,2,var)
var.col.sort = sort(var.col,decreasing = T)
var.col.sel = var.col.sort[1:500]
names.sel = names(var.col.sel)
X.subset = subset(X,select = names.sel)


cor=cor(X.subset,y)
names(cor)=names(var.col.sel)
cor.sort=sort(abs(cor),decreasing = T)
cor.sel=cor.sort[1:100]
names.col_sel=names(cor.sel)
X.subset = subset(X,select = names.col_sel)

Data = data.frame(y,X.subset)
Data.std = scale(Data,center=T,scale=T)
dim(Data.std)


####PMPL####
source('SCAD_MPL.R')

data=Data.std
y=data[,1]
X=data[,-1]

#### SCAD choose initial values SCAD.ebeta for PMPL
set.seed(57)
SCAD.cvfit = cv.ncvreg(X,y,penalty="SCAD",max.iter=100000)
SCAD.fit = SCAD.cvfit$fit
beta = SCAD.fit$beta[,SCAD.cvfit$min]
SCAD.ebeta = beta[-1]

#### PMPL_1:Select the bandwidth using the Rule of Thumb
SCAD.MPL.result = SCAD.MPL.BIC(y, X, SCAD.ebeta,0,1,0.01)
SCAD.MPL.ebeta = SCAD.MPL.result[-1]

#### PMPL_2
##The bandwidth output above is used as a priori.
##We use the maximum-estimated-likelihood cross-validation method 
##to select the bandwidth.
n=nrow(X)
bandwidth=c(0.147885727904182,seq(0.0933762829936264,0.434319059814676,length.out =49))#priori from PMPL_1
out=numeric(0)
for(band in 1:length(bandwidth)){
  h=bandwidth[band]
  SCAD.MPL.result = SCAD.MPL.BIC_h(y, X, SCAD.ebeta,0,1,0.01,h)
  SCAD.MPL.ebeta = SCAD.MPL.result[-1]
  residual0 = y-X%*%SCAD.MPL.ebeta
  a = numeric(n)
  for(i in 1:n){
    a[i] = sum( gaussionkernel((residual0-residual0[i])/h)/n/h )
  }
  out[band]=sum( log(a))
}
h=bandwidth[which.max(out)]
print(paste("h:",h))

SCAD.MPL.result = SCAD.MPL.BIC_h(y, X, SCAD.ebeta,0,1,0.01,h)
SCAD.MPL.ebeta = SCAD.MPL.result[-1]
names(SCAD.MPL.ebeta)=names(data[,-1])
p_sele=sum(SCAD.MPL.ebeta!=0)
print(paste("The number of selected variables:",p_sele))

###computing t value
sigma_hat=t(y-X%*%SCAD.MPL.ebeta)%*%(y-X%*%SCAD.MPL.ebeta)/(length(y)-p_sele)
sigma_hat=sqrt(sigma_hat)
X_del=X[,which(SCAD.MPL.ebeta!=0)]
c_jj=sqrt(diag(solve(t(X_del)%*%X_del)))
std=as.vector(sigma_hat)*c_jj
t_value=SCAD.MPL.ebeta[which(SCAD.MPL.ebeta!=0)]/std

result=cbind(SCAD.MPL.ebeta[SCAD.MPL.ebeta!=0],std,t_value)
print(result)
