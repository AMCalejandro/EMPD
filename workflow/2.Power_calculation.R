library(faux)
library(lme4)

Ntests<-1000 #tolal N of SNPs tested  
NSims <-10000 #make Nsims at least 10 times more than N tests 

# We define a df to store the results for the intercept only model and the model with random slopes too
res_slopes<-as.data.frame(matrix(ncol=2, nrow=NSims))
names(res_slopes)<-c("B","P")
res_intercept<-as.data.frame(matrix(ncol=2, nrow=NSims))
names(res_intercept)<-c("B","P")


#Assumptions:
#age at baseline and UPDRS are independent; age ranged between 50 and 90
#UPDRS goes ((((UP and not down as Maryam initially said))))  withing 3 assessments with r1=-0.2 and is N(0,1)
#SNP is associated with the UPDRS rate of decline with r2=0.05


N<-3500
Rep<-3
Nt<-N*Rep
MAF_vector =seq(from = 0.001, to = 0.01, by = 0.001)
r1 <- 0.2 
r2 <- 0.1



dat<-as.data.frame(matrix(ncol=4, nrow=Nt))
names(dat)<-c("ID", "AGE", "UPDRS", "SNP")
snp<-as.data.frame(matrix(ncol=3, nrow=N))
names(snp)<-c("mUPDRS", "SNPrnd", "SNP")


# List to store different MAF powers
maf_list <- vector(mode='list', length= length(MAF_vector))
names(maf_list) = as.character(MAF_vector)


for (index in 1:length(MAF_vector)) {
  MAF = MAF_vector[index]
  for (sim in 1:NSims)
  {
    for (n in 1:N){
      dat$ID[((n-1)*Rep+1):(Rep*n)]<-n
    
      # Adding a visit number
      dat$Time[((n-1)*Rep+1):(Rep*n)]<- c(1:3)
      # Also, I have changed the age to make it more realistic
      dat$AGE[((n-1)*Rep+1)]<-runif(1,55,80)
      dat$AGE[((n-1)*Rep+2)]<-runif(1,dat$AGE[((n-1)*Rep+1)] + 1, dat$AGE[((n-1)*Rep+1)] + 2)
      dat$AGE[(n*Rep)]<-runif(1, dat$AGE[((n-1)*Rep+2)] + 1, dat$AGE[((n-1)*Rep+2)] + 2)
    }
    dat$UPDRS<-rnorm_pre(dat$AGE, mu = 15, sd = 3, r=r1)
    
    
    for (n in 1:N){snp$mUPDRS[n]<-mean(dat$UPDRS[((n-1)*Rep+1):(Rep*n)])}
    snp$SNPrnd<-rnorm_pre(snp$mUPDRS, mu = 0, sd = 1, r=r2)
    a<-sort(snp$SNPrnd); hom11<-a[N*MAF*MAF]; het12<-a[N*MAF*MAF+N*2*(MAF*(1-MAF))]
    b<-which(snp$SNPrnd<=hom11);snp$SNP[b]<-0
    b<-which(hom11<snp$SNPrnd & snp$SNPrnd<=het12);snp$SNP[b]<-1
    b<-which(het12<snp$SNPrnd);snp$SNP[b]<-2
    for (n in 1:N){dat$SNP[((n-1)*Rep+1):(Rep*n)]<-snp$SNP[n]}
    
    
    random_slopes <- lmer(scale(UPDRS) ~ 1+(1+Time|ID) + scale(AGE) + Time + SNP + SNP*Time, data=dat, REML=FALSE)
    coefs_slopes <- as.data.frame(coef(summary(random_slopes)))
    res_slopes$B[sim]<-coefs_slopes$Estimate[5]
    res_slopes$P[sim]<-2 * (1 - pnorm(abs(coefs_slopes$t[5])))
    
    random_intercept <- lmer(scale(UPDRS) ~ 1+(1+scale(AGE)|ID) + scale(AGE) + SNP, data=dat, REML=FALSE)
    coefs_intercept <- as.data.frame(coef(summary(random_intercept)))
    res_intercept$B[sim]<-coefs_intercept$Estimate[3]
    res_intercept$P[sim]<-2 * (1 - pnorm(abs(coefs_intercept$t[3])))
  }

  power_slopes <-length(which(res_slopes$P<=(0.05/Ntests)))/NSims
  power_intercept <-length(which(res_intercept$P<=(0.05/Ntests)))/NSims
  maf_list[[index]] = c(power_slopes, power_intercept)
}


