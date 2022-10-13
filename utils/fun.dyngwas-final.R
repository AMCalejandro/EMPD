
# This is the SCEBE algorithm that we downloaded from -> https://github.com/Myuan2019/SCEBE

####functions

# 计算时间差
timediff <- function(End_time, Start_time){
  return(difftime(End_time, Start_time, units = 'secs'))
}

##SPR for NEBE
slmm <- function(y1,y2, X){
  
  n1 <- length(y1)
  n<-n2<-n1
  
  y1bar <- mean(y1)
  y2bar <- mean(y2)
  
  y1c <- y1 - y1bar
  y2c <- y2 - y2bar
  s1 <- colSums(X)
  s2 <- colSums(X ^ 2)
  sxx <- (s2 - (s1 ^ 2) / n)
  
  b1 <- as.vector(colSums(y1c*X) / sxx)
  b2 <- as.vector(colSums(y2c*X) / sxx)
  y1_sum <- sum(y1)
  y2_sum <- sum(y2)
  y1X <- colSums(y1 * X)
  y2X <- colSums(y2 * X)
  
  #Residuals1 <- (y1 - X %*% diag(b1) - t(replicate(n1, a1))) ## still potential for improvement
  #b1_sd <- sqrt((colSums(Residuals1 ^ 2) / (n1 - 2)) / sxx)
  #Residuals2 <- (y2 - X %*% diag(b2) - t(replicate(n2, a2))) ## still potential for improvement
  #b2_sd <- sqrt((colSums(Residuals2 ^ 2) / (n2 - 2)) / sxx)
  
  sigma11 = sum(y1 ^ 2) - (y1_sum ^ 2 * s2 - 2 * s1 * y1_sum * y1X +
                             n * (y1X) ^ 2) / (n * sxx)
  var.gamma11 = sigma11 / (sxx * (n - 2))
  
  #sigma12 = sum(y1 * y2) - (y1_sum * y2_sum * s2 - y2_sum * s1 * y1X -
  #                           y1_sum * s1 * y2X + n * y1X * y2X) / (n * sxx)
  #var.gamma12 = sigma12 / (sxx * (n - 2))
  
  sigma22 = sum(y2 ^ 2) - (y2_sum ^ 2 * s2 - 2 * s1 * y2_sum * y2X +
                             n * (y2X) ^ 2) / (n * sxx)
  var.gamma22 = sigma22 / (sxx * (n - 2))
  
  b1_sd=sqrt(var.gamma11)
  b2_sd=sqrt(var.gamma22)
  
  
  t1_value <- b1 / b1_sd
  p1_value <- pt(q = -abs(t1_value), df =  n1 - 2) * 2 # Pr(>|t|)
  
  t2_value <- b2 / b2_sd
  p2_value <- pt(q = -abs(t2_value), df =  n2 - 2) * 2 # Pr(>|t|)
  
  
  result <- cbind(b1, b1_sd, p1_value,b2, b2_sd,  p2_value)
  colnames(result) <- c("b1_Est", "b1_Std",  "Pr1(>|t|)","b2_Est", "b2_Std",  "Pr2(>|t|)")
  rownames(result) <- colnames(X)
  return(result)
}
###########################


scebe_sim<-function(phenoData=myData,genoData=covs.N, fit0, Time="Time", pheno="Resp", method="scebe"){
  
  N = nrow(genoData)
  noCov = ncol(genoData)
  n_measure = phenoData %>% dplyr::group_by(ID) %>% dplyr::summarize(n=length(ID))%>%dplyr::select(n)%>%unlist()
  
  ##base model#### take out
  #fit0 <-try(lmer(Resp~Time+(Time|ID),data=phenoData),silent=TRUE)
  
  if(inherits(fit0,"try-error")){
    print("Base model didn't converge") #output <-  matrix(NA, nrow=1, ncol=9)
  } else {
    
    flme0<-summary(fit0)
    rf=ranef(fit0)$ID
    #r.name=rownames(rf)
    #rf1=data.frame(ID=as.factor(r.name),rf) %>% arrange(ID)
    #rf = rf1[, -which(names(rf1)=="ID")]
    # ind.parm <- merge(phenoData[!duplicated(phenoData$ID),c("ID",paste0("cov",1:noCov))],rf1,by="ID")
    # names(ind.parm)[c(2+noCov,3+noCov)]=c("Inter","Slope")
    names(rf)[c(1:2)]=c("Inter","Slope")
    
    
    if(method=="nebe"){
      
      etime.start<-Sys.time()###alletime
      {
        #allX=ind.parm[,2:(1+noCov)]
        #spr<-slmm(y1=ind.parm$Inter,y2=ind.parm$Slope,X=as.matrix(allX))
        spr<-slmm(y1=rf$Inter,y2=rf$Slope,X=as.matrix(genoData))
      }
      etime.end<-Sys.time()
      etime=timediff(etime.end,etime.start)
      etime=etime/noCov
      #petime = 0
      output=cbind(spr,etime)
      colnames(output)= c("est.ebe.b1", "se.ebe.b1", "p.ebe.b1",
                          "est.ebe.b2", "se.ebe.b2", "p.ebe.b2","etime")
    }#end of nebe
    
    
    if(method=="lme"){
      ##2.full model##
      
      
      
      
      allrst.lme=NULL
      
      for(i in 1:noCov){
        
        tmp = rep(genoData[,i],n_measure)### change 1: lognormal for WT 30% var
        phenoData1 = cbind(phenoData,tmp)
        names(phenoData1)[ncol(phenoData1)]=paste0("cov",i)
        
        form.lme = as.formula(paste(pheno, "~", Time, "|ID"))
        
        phenoData1 <- groupedData(form.lme,
                                  data = phenoData1,
                                  labels = list( x = "Time", y = "RESP"),
                                  order.groups = FALSE)
        
        ltime.start<-Sys.time()
        {
          fml <- as.formula(paste0(pheno, "~", Time, "*cov",i,"+(", Time, "|ID)"))
          fit1 <- try(lmer(fml, data = phenoData1),silent=TRUE)
          if(inherits(fit1,"try-error")){
            output<-  matrix(NA, nrow=1, ncol=7)
          } else {
            
            flme1<-summary(fit1)
            est.lme=coef(flme1)[3:4,]
            p.lme=2*pnorm(-abs(est.lme[,3]))
            rst.lme=c(est.lme[1,1:2],p.lme[1],est.lme[2,1:2],p.lme[2])
          }
        }#end of time start
        ltime.end<-Sys.time()
        ltime=timediff(ltime.end,ltime.start)
        #pltime = 0
        rst=c(rst.lme,ltime)
        rbind(allrst.lme,rst)->allrst.lme
      }#end of loop
      output=allrst.lme
      rownames(output)=paste0("cov",1:noCov)
      
      colnames(output)= c("est.lme.b1", "se.lme.b1", "p.lme.b1",
                          "est.lme.b2", "se.lme.b2", "p.lme.b2","ltime")
    }#end of lme
    
    if(method=="gallop"){
      ##GALLOP
      # Extract variance components and compute penalty matrix (P)
      pgtime.start<-Sys.time()
      {
        varcor = VarCorr(fit0)
        sig = attr(varcor, "sc") ##within subject error
        R = varcor$ID ##between subject error ##when nt==1 return error!
        P = ginv(R / sig ^ 2)
        
        
        # Put data in convenient arrays and vector; base model
        Zstar = cbind(1, phenoData[, Time])
        Xstar = Zstar
        Ystar = phenoData[,pheno]
        
        # Compute components of block-diagonal system with covariates  (without SNP)
        # Additionally compute and store object SS = Rot %*% Si, which is used later on
        A21 = matrix(NA, 2 * N , 2)
        q2  = rep(NA, 2 * N)
        SS = matrix(NA, 2 * N , 2)
        for (i in 1 : N ) {
          if (i == 1) {
            uk <- 1:n_measure[i]
          } else {
            uk = (sum(n_measure[1:(i - 1)]) + 1) : sum(n_measure[1:i])
          }
          u2 = (i - 1) * 2 + (1 : 2)
          Zstari = Zstar[uk, ]
          if(n_measure[i]==1) Zstari=matrix(Zstari,1)
          Si = crossprod(Zstari, Zstari)
          sv = svd(Si + P)
          Rot = sqrt(1 /sv$d) * sv$u ###trans matrix phi
          Q <-  tcrossprod(Rot, Zstari)
          
          SS[u2, ] = Rot %*% Si
          A21[u2, ] = Q %*% Xstar[uk,  ]
          q2[u2] = Q %*% Ystar[uk]
        }
        
        q1 = crossprod(Xstar, Ystar)
        A11 = crossprod(Xstar)
        # Solve the system (20)
        QQ = A11 - crossprod(A21)
        q = q1 - crossprod(A21, q2)
        sol = solve(QQ, q) #beta star
        blups = q2 - A21 %*% sol #theta
        
        # Compute sums of products per subject involved in the
        #crossprod(X, G), crossprod(G) and crossprod(G, y).
        #We use row-wise Kronecker product to avoid repeating SNP vector k times
        
        ex <- et <-  matrix(1, 1, 2)
        XTk = kronecker(et, Zstar) * kronecker(Xstar, ex)
        TTk = kronecker(et, Xstar) * kronecker(Xstar, et)
        Tyk = Ystar * Xstar
        XTs = matrix(0, N, ncol(XTk))
        TTs = matrix(0, N, ncol(TTk))
        Tys = matrix(0, N, 2)
        AtS = matrix(0, N, 2 * 2)
        for (i in 1:N) {
          if (i == 1) {
            uk = 1:n_measure[i]
          }
          else {
            uk = (sum(n_measure[1:(i - 1)]) + 1):sum(n_measure[1:i])
          }
          if(n_measure[i]==1) {
            XTs[i, ] <- matrix(XTk[uk, ],1)
            TTs[i, ] <- matrix(TTk[uk, ],1)
            Tys[i, ] <- matrix(Tyk[uk, ],1)
          }
          if(n_measure[i]!=1) {
            XTs[i, ] <- colSums(XTk[uk, ])
            TTs[i, ] <- colSums(TTk[uk, ])
            Tys[i, ] <- colSums(Tyk[uk, ])
          }
          u2 = (i - 1) * 2 + (1 : 2)
          AtS[i, ] = c(crossprod(A21[u2, ], SS[u2, ]))
        }
        
      }
      pgtime.end<-Sys.time()
      
      pgtime<-timediff(pgtime.end,pgtime.start)
      pgtime<-pgtime/noCov
      
      # Add WT
      allrst.gallop=NULL
      
      for(i in 1:noCov){
        gtime.start <- Sys.time()
        {
          wi = genoData[,i]
          lWT2 = rep(wi, each = 2)
          H1 = matrix(crossprod(wi, XTs), 2, 2) ##X'G
          H2 = lWT2 * SS ##phi Z'G
          AtH = matrix(crossprod(wi, AtS), 2, 2)##
          RR = H1 - AtH  ##H21 x'G-A21^tran' phi z'G
          Cfix = solve(QQ, RR)
          Cran = H2 - A21 %*% Cfix
          GtG = matrix(crossprod(wi ^ 2, TTs), 2, 2) #G'G
          Gty = matrix(crossprod(wi, Tys), 2, 1) #G'Y
          V = GtG - crossprod(H1, Cfix) - crossprod(H2, Cran)
          v = Gty - crossprod(H1, sol) - crossprod(H2, blups)
          Theta= solve(V, v)
          D = diag(solve(V))
          
          
          SE = sig * sqrt(D)
          Pval = 2 * pnorm(-abs(Theta / SE))
        }#end of gtime
        gtime.end <- Sys.time()
        gtime <- timediff(gtime.end, gtime.start)
        
        rst.gallop1<-cbind(Theta,SE,Pval)
        rst.gallop<-c(rst.gallop1[1,],rst.gallop1[2,],gtime,pgtime)
        
        
        allrst.gallop=rbind(allrst.gallop,rst.gallop)
      }# end of loop
      
      
      
      
      rownames(allrst.gallop)=paste0("cov",1:noCov)
      output=allrst.gallop
      colnames(output)= c("est.gal.b1", "se.gal.b1", "p.gal.b1",
                          "est.gal.b2", "se.gal.b2", "p.gal.b2","gtime","pgtime")
      
      
    }# end of gallop
    
    
    
    if(method=="scebe"){
      
      # allX=ind.parm[,2:(1+noCov)]
      spr<-slmm(y1=rf$Inter,y2=rf$Slope,X=as.matrix(genoData))
      
      
      Data=phenoData
      R=matrix(as.numeric(VarCorr(fit0)$ID),ncol=2)
      Ip=diag(1,2)
      
      pstime.start=Sys.time()
      {
        
        #allS=allW=list()
        #allM=allISMIS=list()
        alls11=alls12=alls21=alls22=NULL
        allw11=allw12=allw21=allw22=NULL
        allsms11=allsms12=allsms21=allsms22=NULL
        
        Time.lst = split(Data[[Time]], Data[["ID"]])
        
        for(i in 1:N){
          nt = n_measure[i]
          time <- Time.lst[[i]] #Data[Data$ID ==i, Time]##
          Zi <- cbind(1, time)
          Gi = diag(sigma(fit0) ^ 2, nt)
          ZinvGZ <- crossprod(Zi, solve(Gi, Zi)) ##when n measure==1, return error!
          Mi=R+ginv(ZinvGZ)
          
          SIGMAi=Zi%*%R%*%t(Zi)+Gi
          Wi <- crossprod(Zi, solve(SIGMAi, Zi))
          #allW[[i]]=Wi
          allw11=c(allw11,Wi[1,1])
          allw12=c(allw12,Wi[1,2])
          allw21=c(allw21,Wi[2,1])
          allw22=c(allw22,Wi[2,2])
          
          
          
          invSi <- R %*% ZinvGZ + Ip
          Si=solve(invSi)
          alls11<-c(alls11,Si[1,1])
          alls12<-c(alls12,Si[1,2])
          alls21<-c(alls21,Si[2,1])
          alls22<-c(alls22,Si[2,2])
          
          Ipi=diag(1,2)
          SMS=(Ipi-Si)%*%Mi%*%t(Ipi-Si)
          
          allsms11<-c(allsms11,SMS[1,1])
          allsms12<-c(allsms12,SMS[1,2])
          allsms21<-c(allsms21,SMS[2,1])
          allsms22<-c(allsms22,SMS[2,2])
          
          #allISMIS[[i]]<-SMS
          #allS[[i]]=Si
        }
        
        #C1 <- purrr::reduce(allW, `+`)
        C<-matrix(c(sum(allw11),sum(allw12),sum(allw21),sum(allw22)),2,byrow=T)
        
        inv_C<-solve(C)
        
        invc11<-inv_C[1,1]
        invc12<-inv_C[1,2]
        invc21<-inv_C[2,1]
        invc22<-inv_C[2,2]
        
      }
      pstime.end<-Sys.time()
      pstime<-timediff(pstime.end,pstime.start)
      
      
      ebe<-spr[,c(1,4)]
      
      stime.start<-Sys.time()
      {
        
        myX<-genoData #unique(phenoData[,5:(noCov+4)])
        Xc<-apply(myX,2,function(x){scale(x,scale=F)})
        
        XcX<-Xc*myX
        Xc2<-Xc^2
        Sxc2<-colSums(Xc2)
        
        xstar11<-colSums(allw11*Xc)
        xstar12<-colSums(allw12*Xc)
        xstar21<-colSums(allw21*Xc)
        xstar22<-colSums(allw22*Xc)
        
        A11<-colSums(Xc2*alls11)/Sxc2
        A12<-colSums(Xc2*alls12)/Sxc2
        A21<-colSums(Xc2*alls21)/Sxc2
        A22<-colSums(Xc2*alls22)/Sxc2
        
        B11<-colSums(Xc*alls11)/Sxc2
        B12<-colSums(Xc*alls12)/Sxc2
        B21<-colSums(Xc*alls21)/Sxc2
        B22<-colSums(Xc*alls22)/Sxc2
        
        Sxc22<-(Sxc2)^2
        
        D11=colSums(Xc2*allsms11)/Sxc22
        D12=colSums(Xc2*allsms12)/Sxc22
        D21=colSums(Xc2*allsms21)/Sxc22
        D22=colSums(Xc2*allsms22)/Sxc22
        
        
        
        ssh211<-B11*invc11*xstar11+B12*invc21*xstar11+B11*invc12*xstar21+B12*invc22*xstar21
        ssh212<-B11*invc11*xstar12+B12*invc21*xstar12+B11*invc12*xstar22+B12*invc22*xstar22
        ssh221<-B21*invc11*xstar11+B22*invc21*xstar11+B21*invc12*xstar21+B22*invc22*xstar21
        ssh222<-B21*invc11*xstar12+B22*invc21*xstar12+B21*invc12*xstar22+B22*invc22*xstar22
        
        ssh11=1-A11+ssh211
        ssh12=0-A12+ssh212
        ssh21=0-A21+ssh221
        ssh22=1-A22+ssh222
        
        detssh=ssh11*ssh22-ssh12*ssh21
        
        invssh11<-ssh22/detssh
        invssh12<-(-ssh12/detssh)
        invssh21<-(-ssh21/detssh)
        invssh22<-ssh11/detssh
        
        Vgamma211<-B11*invc11*B11+B12*invc21*B11+B11*invc12*B12+B12*invc22*B12
        #Vgamma211<-Vgamma211/SXc2
        Vgamma212<-B11*invc11*B21+B12*invc21*B21+B11*invc12*B22+B12*invc22*B22
        #Vgamma212<-Vgamma212/SXc2
        Vgamma221<-B21*invc11*B11+B22*invc21*B11+B21*invc12*B12+B22*invc22*B12
        #Vgamma221<-Vgamma221/SXc2
        Vgamma222<-B21*invc11*B21+B22*invc21*B21+B21*invc12*B22+B22*invc22*B22
        #Vgamma222<-Vgamma222/SXc2
        
        Vgamma11=D11-Vgamma211
        Vgamma12=D12-Vgamma212
        Vgamma21=D21-Vgamma221
        Vgamma22=D22-Vgamma222 ###VAR of NEBE
        
        Vssh11<-invssh11*Vgamma11*invssh11+invssh12*Vgamma21*invssh11+invssh11*Vgamma12*invssh12+invssh12*Vgamma22*invssh12
        Vssh12<-invssh11*Vgamma11*invssh21+invssh12*Vgamma21*invssh21+invssh11*Vgamma12*invssh22+invssh12*Vgamma22*invssh22
        Vssh21<-invssh21*Vgamma11*invssh11+invssh22*Vgamma21*invssh11+invssh21*Vgamma12*invssh12+invssh22*Vgamma22*invssh12
        Vssh22<-invssh21*Vgamma11*invssh21+invssh22*Vgamma21*invssh21+invssh21*Vgamma12*invssh22+invssh22*Vgamma22*invssh22
        
        
        est.ssh.inter<-invssh11*ebe[,1]+invssh12*ebe[,2]
        est.ssh.slope<-invssh21*ebe[,1]+invssh22*ebe[,2]
        
        t.ssh.inter<-est.ssh.inter/sqrt(Vssh11)
        t.ssh.slope<-est.ssh.slope/sqrt(Vssh22)
        
        p.ssh.inter<-2*pnorm(abs(t.ssh.inter),lower.tail=F)
        p.ssh.slope<-2*pnorm(abs(t.ssh.slope),lower.tail=F)
        
        
      }
      stime.end<-Sys.time()
      stime<-timediff(stime.end,stime.start)
      stime=stime/noCov
      pstime<-pstime/noCov
      output<-cbind(est.ssh.inter,sqrt(Vssh11),p.ssh.inter,est.ssh.slope,sqrt(Vssh22),p.ssh.slope,stime,pstime)
      
      colnames(output)= c("est.sc.b1", "se.sc.b1", "p.sc.b1",
                          "est.sc.b2", "se.sc.b2", "p.sc.b2","stime","pstime")
      
    }#end of scebe
    
    
    
  }#end of fit0
  return(output)
  
}#end of function