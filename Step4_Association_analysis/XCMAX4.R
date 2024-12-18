XCMAX4 <- function(data){
  D=as.matrix(data[,1]) 
  G=as.matrix(data[,2])
  
  X=as.matrix(data[,-c(1,2)]) 
  gender <- as.matrix(data["gender"])
  n <- length(D) 
  Ind_AA <- rep(0,n)
  Ind_Aa <- rep(0,n) 
  Ind_AO <- rep(0,n) 
  Ind_AA[gender==1 & G==2] <- 1
  Ind_Aa[gender==1 & G==1] <- 1
  Ind_AO[gender==0 & G==1] <- 1
  
  esta <- summary(glm(D ~ X , family = binomial(link = "logit")))$coefficients[,1]
  
  X <- cbind(rep(1,n),X)
  
  estf=function(esta,X){
    re=1/(1+exp(-X%*%esta))
    return(re)
  }
  estpen <- estf(esta,X) 
  
  infora=function(estpen,X){
    l=dim(X)[2]
    Ia=matrix(0, nrow=l,ncol=l)
    for(i in 1:l){
      for(j in 1:l){
        Ia[i,j]=sum(X[,i]*X[,j]*(1-estpen)*estpen)
      }
    }
    return(Ia)
  }
  Ia=infora(estpen,X) 
  
  inforb=function(estpen,z1,z2,X,Ind_AO,Ind_Aa,Ind_AA){
    G = 2*Ind_AA + z1*Ind_Aa + z2*Ind_AO
    Ib=sum(G*G*(1-estpen)*estpen)
    return(Ib)
  }
  Ib11=inforb(estpen,1,1,X,Ind_AO,Ind_Aa,Ind_AA) 
  Ib02=inforb(estpen,0,2,X,Ind_AO,Ind_Aa,Ind_AA) 
  Ib12=inforb(estpen,1,2,X,Ind_AO,Ind_Aa,Ind_AA) 
  Ib22=inforb(estpen,2,2,X,Ind_AO,Ind_Aa,Ind_AA) 
  
  
  inforba=function(estpen,z1,z2,X,Ind_AO,Ind_Aa,Ind_AA){
    l=dim(X)[2]
    G = 2*Ind_AA + z1*Ind_Aa + z2*Ind_AO
    Iba=NULL
    for(i in 1:l){
      Iba[i]=sum(X[,i]*G*(1-estpen)*estpen)
    }
    return(Iba)
  }
  Iba11=inforba(estpen,1,1,X,Ind_AO,Ind_Aa,Ind_AA) 
  Iba02=inforba(estpen,0,2,X,Ind_AO,Ind_Aa,Ind_AA) 
  Iba12=inforba(estpen,1,2,X,Ind_AO,Ind_Aa,Ind_AA) 
  Iba22=inforba(estpen,2,2,X,Ind_AO,Ind_Aa,Ind_AA) 
  
  stest=function(estpen,D,z1,z2,Ind_AO,Ind_Aa,Ind_AA,Ib,Iba,Ia){
    G <- 2*Ind_AA + z1*Ind_Aa + z2*Ind_AO
    score <- sum(G*(D-estpen))
    variance <- Ib-t(as.matrix(Iba))%*%solve(Ia)%*%as.matrix(Iba)
    re=score/sqrt(variance) ## test statistic
    return(re)
  }
  s11=stest(estpen,D,1,1,Ind_AO,Ind_Aa,Ind_AA,Ib11,Iba11,Ia) 
  s02=stest(estpen,D,0,2,Ind_AO,Ind_Aa,Ind_AA,Ib02,Iba02,Ia) 
  s12=stest(estpen,D,1,2,Ind_AO,Ind_Aa,Ind_AA,Ib12,Iba12,Ia) 
  s22=stest(estpen,D,2,2,Ind_AO,Ind_Aa,Ind_AA,Ib22,Iba22,Ia)
  
  Iba11_02=rbind(Iba11,Iba02)
  Iba11_12=rbind(Iba11,Iba12)
  Iba11_22=rbind(Iba11,Iba22)
  Iba02_12=rbind(Iba02,Iba12)
  Iba02_22=rbind(Iba02,Iba22)
  Iba12_22=rbind(Iba12,Iba22)
  
  infor_dob=function(estpen,z11,z21,z12,z22,D,Ind_AO,Ind_Aa,Ind_AA){
    G1 = 2*Ind_AA + z11*Ind_Aa + z21*Ind_AO
    G2 = 2*Ind_AA + z12*Ind_Aa + z22*Ind_AO
    I1=sum(G1*G1*(1-estpen)*estpen)
    I12=sum(G1*G2*(1-estpen)*estpen)
    I2=sum(G2*G2*(1-estpen)*estpen)
    re=matrix(  c(I1,I12,I12,I2),nrow=2 )
    return(re)
  }
  inf11_02=infor_dob(estpen,1,1,0,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf11_12=infor_dob(estpen,1,1,1,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf11_22=infor_dob(estpen,1,1,2,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf02_12=infor_dob(estpen,0,2,1,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf02_22=infor_dob(estpen,0,2,2,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf12_22=infor_dob(estpen,1,2,2,2,D,Ind_AO,Ind_Aa,Ind_AA)
  
  C11_02=inf11_02 - as.matrix(Iba11_02)%*%solve(Ia)%*%t(as.matrix(Iba11_02)) 
  C11_12=inf11_12 - as.matrix(Iba11_12)%*%solve(Ia)%*%t(as.matrix(Iba11_12)) 
  C11_22=inf11_22 - as.matrix(Iba11_22)%*%solve(Ia)%*%t(as.matrix(Iba11_22)) 
  C02_12=inf02_12 - as.matrix(Iba02_12)%*%solve(Ia)%*%t(as.matrix(Iba02_12)) 
  C02_22=inf02_22 - as.matrix(Iba02_22)%*%solve(Ia)%*%t(as.matrix(Iba02_22)) 
  C12_22=inf12_22 - as.matrix(Iba12_22)%*%solve(Ia)%*%t(as.matrix(Iba12_22)) 
  
  RC11_02 <- C11_02[1,2] / sqrt( C11_02[1,1]* C11_02[2,2]) 
  RC11_12 <- C11_12[1,2] / sqrt( C11_12[1,1]* C11_12[2,2]) 
  RC11_22 <- C11_22[1,2] / sqrt( C11_22[1,1]* C11_22[2,2]) 
  RC02_12 <- C02_12[1,2] / sqrt( C02_12[1,1]* C02_12[2,2]) 
  RC02_22 <- C02_22[1,2] / sqrt( C02_22[1,1]* C02_22[2,2])
  RC12_22 <- C12_22[1,2] / sqrt( C12_22[1,1]* C12_22[2,2]) 
  
  vacov=matrix(c( 1,        RC11_02, RC11_12, RC11_22,
                  RC11_02,  1,       RC02_12, RC02_22,
                  RC11_12,  RC02_12, 1,       RC12_22,       
                  RC11_22,  RC02_22, RC12_22, 1       ),
               ncol=4,
               byrow = TRUE)
  
  zmax1 <- max(abs(s11),abs(s02),abs(s12),abs(s22))
  
  #rhombus formula
  p_rh <- function(zmax1,vacov,a,b,c,d){
    part1 <- 2*(pnorm(zmax1)-pnorm(-zmax1)-1)
    l12 <- acos(vacov[a,b]);l23 <- acos(vacov[b,c]);l34<-acos(vacov[c,d])
    part2 <- pnorm(l12*zmax1/2) + pnorm((pi-l12)*zmax1/2) - 1
    part3 <- pnorm(l23*zmax1/2) + pnorm((pi-l23)*zmax1/2) - 1
    part4 <- pnorm(l34*zmax1/2) + pnorm((pi-l34)*zmax1/2) - 1
    p_rh <- part1 + 4*dnorm(zmax1)/zmax1*(part2+part3+part4)
    return(p_rh)
  } 
  p_series <- c(p_rh(zmax1,vacov,1,2,3,4),
                p_rh(zmax1,vacov,1,2,4,3),
                p_rh(zmax1,vacov,1,3,2,4),
                p_rh(zmax1,vacov,1,3,4,2),
                p_rh(zmax1,vacov,1,4,2,3),
                p_rh(zmax1,vacov,1,4,3,2),
                p_rh(zmax1,vacov,2,1,3,4),
                p_rh(zmax1,vacov,2,1,4,3),
                p_rh(zmax1,vacov,2,3,1,4),
                p_rh(zmax1,vacov,2,4,1,3),
                p_rh(zmax1,vacov,3,1,2,4),
                p_rh(zmax1,vacov,3,2,1,4))
  p_rh <- min(p_series,1)
  
  return(list("statictic"=zmax1,
              "p-value"=p_rh))
 

}

