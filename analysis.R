rm(list=ls())
library(lmerTest);library(lme4);library(reshape);library(effects);library(matrixStats)
#code to plot error bars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
error.bar.horiz <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x+upper,y, x-lower, y, angle=90, code=3, length=length, ...)
}
#Code for VIFs
vif.mer <- function (fit) {
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

Data <- read.table("data.csv",sep=",",header=T)
head(Data)

#012 is now CIN, should become CNI (because the package can't deal with factors)
Data$Congruency <- Data$compatibility
with(Data,aggregate(ratings_raw,by=list(Congruency),mean))
Data$Congruency[Data$Congruency==2] <- 3
Data$Congruency[Data$Congruency==1] <- 2
Data$Congruency[Data$Congruency==3] <- 1
#Data$Congruencyz <- scale(Data$Congruency)

#Recode ratings so that 0=easy,359=hard
Data$ratings_raw <- 359-Data$ratings_raw

Data$ratings_raw_z <- NA;Data$resid_ratings <- NA
for(i in 1:length(table(Data$subject))){
  Data$ratings_raw_z[Data$subject==unique(Data$subject)[i]] <- scale(Data$ratings_raw[Data$subject==unique(Data$subject)[i]])  
}

subs <- unique(Data$sub)
N <- length(subs)

##################################
#1. Simple main effects analysis #
#1.1 RCZ
fit <- lmer(rcz~ratings_raw_z+(1|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
fit_a <- lmer(rcz~ratings_raw_z+(ratings_raw_z|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
anova(fit,fit_a)
anova(fit_a)

#To visualize this, aggregate RCZ and ratings in 20 bins, plot this and add the regression line
temp <- data.frame(effect('ratings_raw_z',fit_a,xlevels=list(ratings_raw_z=seq(-10,10,.5))));temp$rcz <- temp$fit;temp <- temp[,c(6,1,4,5)]
cuts <- 20
Data$rcz_bin <- NA;Data$ratings_raw_z_bin <- NA
for(i in 1:N){
  tempDat <- subset(Data,subject==subs[i])
  Data$rcz_bin[Data$subject==subs[i]] <- as.numeric(cut(tempDat$rcz, cuts))
  Data$ratings_raw_z_bin[Data$subject==subs[i]] <- as.numeric(cut(tempDat$ratings_raw_z, cuts))
}
avg_rcz_rating <- with(Data,aggregate(rcz,by=list(ratings_raw_z_bin=ratings_raw_z_bin,sub=subject),mean));avg_rcz_rating <- cast(avg_rcz_rating,sub~ratings_raw_z_bin)
avg_rating_rcz <- with(Data,aggregate(ratings_raw_z,by=list(ratings_raw_z_bin=ratings_raw_z_bin,sub=subject),mean));avg_rating_rcz <- cast(avg_rating_rcz,sub~ratings_raw_z_bin)

plot(colMeans(avg_rcz_rating,na.rm=T)~colMeans(avg_rating_rcz,na.rm=T),ylab="Rostro-Cingulate Zone (RCZ)",xlab="Difficulty ratings (z)",frame=F,ylim=c(-1,4),xlim=c(-3,5),pch=19)
polygon(x=c(temp$ratings_raw_z,rev(temp$ratings_raw_z)),y=c(temp$upper, rev(temp$lower)),border=F,col=rgb(.8,.8,.8,.5))
lines(temp$rcz~temp$ratings_raw_z,lty=3)
points(colMeans(avg_rcz_rating,na.rm=T)~colMeans(avg_rating_rcz,na.rm=T),pch=19)
error.bar(colMeans(avg_rating_rcz,na.rm=T),colMeans(avg_rcz_rating,na.rm=T),colSds(as.matrix(avg_rcz_rating),na.rm=T)/sqrt(N),length=0)
error.bar.horiz(colMeans(avg_rating_rcz,na.rm=T),colMeans(avg_rcz_rating,na.rm=T),colSds(as.matrix(avg_rating_rcz),na.rm=T)/sqrt(N),length=0)

#1.2. Insula
fit <- lmer(insula~ratings_raw_z+(1|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
fit_a <- lmer(insula~ratings_raw_z+(ratings_raw_z|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
anova(fit,fit_a)
anova(fit_a)

# Plot raw data and add the mixed model regression line
temp <- data.frame(effect('ratings_raw_z',fit_a,xlevels=list(ratings_raw_z=seq(-10,10,.5))));temp$insula <- temp$fit;temp <- temp[,c(6,1,4,5)]

#aggregate AI and ratings in 20 bins, plot this and add the regression line
cuts <- 20
Data$ai_bin <- NA;Data$ratings_raw_z_bin <- NA
for(i in 1:N){
  tempDat <- subset(Data,subject==subs[i])
  Data$ai_bin[Data$subject==subs[i]] <- as.numeric(cut(tempDat$insula, cuts))
  Data$ratings_raw_z_bin[Data$subject==subs[i]] <- as.numeric(cut(tempDat$ratings_raw_z, cuts))
}
avg_ai_rating <- with(Data,aggregate(insula,by=list(ratings_raw_z_bin=ratings_raw_z_bin,sub=subject),mean));avg_ai_rating <- cast(avg_ai_rating,sub~ratings_raw_z_bin)
avg_rating_ai <- with(Data,aggregate(ratings_raw_z,by=list(ratings_raw_z_bin=ratings_raw_z_bin,sub=subject),mean));avg_rating_ai <- cast(avg_rating_ai,sub~ratings_raw_z_bin)

plot(colMeans(avg_ai_rating,na.rm=T)~colMeans(avg_rating_ai,na.rm=T),ylab="Anterior Insula (AI)",xlab="Difficulty ratings (z)",frame=F,ylim=c(-1,4),xlim=c(-3,5),pch=19)
polygon(x=c(temp$ratings_raw_z,rev(temp$ratings_raw_z)),y=c(temp$upper, rev(temp$lower)),border=F,col=rgb(.8,.8,.8,.5))
lines(temp$insula~temp$ratings_raw_z,lty=3)
points(colMeans(avg_ai_rating,na.rm=T)~colMeans(avg_rating_ai,na.rm=T),pch=19)
error.bar(colMeans(avg_rating_ai,na.rm=T),colMeans(avg_ai_rating,na.rm=T),colSds(as.matrix(avg_ai_rating),na.rm=T)/sqrt(N),length=0)
error.bar.horiz(colMeans(avg_rating_ai,na.rm=T),colMeans(avg_ai_rating,na.rm=T),colSds(as.matrix(avg_rating_ai),na.rm=T)/sqrt(N),length=0)


##################################################
#2. More complex models controlling for CE and RT#
#2.1 RCZ
Data$Congruency <- as.factor(Data$Congruency)
fit <- lmer(rcz~rts+ratings_raw_z*Congruency+(1|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
fit_a <- lmer(rcz~rts+ratings_raw_z*Congruency+(rts|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
fit_b <- lmer(rcz~rts+ratings_raw_z*Congruency+(ratings_raw_z|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
fit_c <- lmer(rcz~rts+ratings_raw_z*Congruency+(Congruency|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead")) #doesn't run
anova(fit,fit_a)
anova(fit,fit_b)
fit2 <- lmer(rcz~rts+ratings_raw_z*Congruency+(rts+ratings_raw_z|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
anova(fit2)
library(multcomp)
contrast.matrix <- rbind("c vs n" = c(0,0,0,1,0,0,0),
                         "c vs i" = c(0,0,0,0,1,0,0),
                         "i vs n" = c(0,0,0,1,-1,0,0)
)
summary(glht(fit2, linfct=contrast.matrix), test=adjusted("none"))

plot(effect('ratings_raw_z:Congruency',fit2))

#2.2 AI
fit <- lmer(insula~rts+ratings_raw_z*Congruency+(1|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
fit_a <- lmer(insula~rts+ratings_raw_z*Congruency+(rts|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
fit_b <- lmer(insula~rts+ratings_raw_z*Congruency+(ratings_raw_z|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
fit_c <- lmer(insula~rts+ratings_raw_z*Congruency+(Congruency|subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead")) #doesn't work
anova(fit,fit_a)
anova(fit,fit_b) #yes
anova(fit_b)
library(multcomp)
contrast.matrix <- rbind("c vs n" = c(0,0,0,1,0,0,0),
                         "c vs i" = c(0,0,0,0,1,0,0),
                         "i vs n" = c(0,0,0,1,-1,0,0)
)
summary(glht(fit_b, linfct=contrast.matrix), test=adjusted("none"))
plot(effect('ratings_raw_z:Congruency',fit2))


#3. Causal Mediation Analysis ------------------------------------------------------
#Mediation sheat sheet:
#ACME: average causal mediation effects
#ADE: average direct effects
#"TotalEffect" tells how much influence there is from x->y, and this is then partitioned into
#acme(average): the mediated effect, and ade(average): the direct effect.
#Prop. Mediated (average) simply gives you ACME(average)/Total Effect

library('mediation');
detach("package:lmerTest", unload=TRUE) #detach lmertest or it gives problem when the mediation functions calls: getCall(medModel)[[1]] 
library(lme4)

#Data$Congruency <- as.factor(Data$Congruency)
Data$RCZ <- Data$rcz
Data$AI <- Data$insula
Data$rating_pure <- Data$ratings_raw
Data$Subject <- Data$subject
Data$rt <- Data$rts

table(Data$sub,Data$Congruency)

#to what extent are AI and RCZ correlated?
ai_rcz <- NA
for(i in unique(Data$subject)){
  ai_rcz[i] <- cor(Data$RCZ[Data$subject==i],Data$AI[Data$subject==i])
}
hist(ai_rcz)

#1. Is the influence of RCZ on rating mediated by AI
#construct the random effect structure depending on whether each random factor adds something
#fyi, there's a (small) bugg in the source code, random slopes cannot be added for factors so you need to code
#them as numeric (or not add random slopes, which in this case isn't needed!)
medModel0 <- lmer(AI~RCZ+Congruency+rt + (1 | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
medModel_a <- lmer(AI~RCZ+Congruency+rt + (RCZ | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
medModel_b <- lmer(AI~RCZ+Congruency+rt + (Congruency | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
medModel_c <- lmer(AI~RCZ+Congruency+rt + (rt | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
anova(medModel0,medModel_a) #p<.001
anova(medModel0,medModel_b) #ns
anova(medModel0,medModel_c) #ns
medModel <- lmer(AI~RCZ+Congruency+rt + (RCZ | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))

outModel0 <- lmer(rating_pure~AI+RCZ+Congruency+rt + (1 | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
outModel_a <- lmer(rating_pure~AI+RCZ+Congruency+rt + (AI | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
outModel_b <- lmer(rating_pure~AI+RCZ+Congruency+rt + (RCZ | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
outModel_d <- lmer(rating_pure~AI+RCZ+Congruency+rt + (Congruency | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
outModel_e <- lmer(rating_pure~AI+RCZ+Congruency+rt + (rt | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
anova(outModel0,outModel_a) #p<.001,BIC=61765
anova(outModel0,outModel_b) #p<.001,BIC=61720
anova(outModel0,outModel_d) #p<.001,BIC=61778
anova(outModel0,outModel_e) #p<.001,BIC=61565
outModel <- lmer(rating_pure~AI+RCZ+Congruency+rt + (rt | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))

#check diagnostics of both models
par(mfrow=c(1,1))
vif.mer(outModel)
plot(outModel)
qqnorm(resid(outModel))
summary(outModel)
anova(outModel)

vif.mer(medModel)
plot(medModel)
qqnorm(resid(medModel))
summary(medModel)
anova(medModel)

#Fit the mediation model (takes some time)
med <- mediate(model.m = medModel, model.y = outModel,
               treat="RCZ", mediator="AI",data=Data,sims=5000,
               covariates=c("Congruency","rt"))
summary(med)

#plot these results
par(mfrow=c(1,1))
par(pin=c(2, 1.915208))
plot(med,yaxt='n',xlab=c("Estimates"),main=c('AI-RCZ-Rating'),frame=F)
axis(2,at=1:3,labels=c("Total effect","Direct effect","Mediation effect"),las=1)


#2. Is the influence of AI on rating mediated by RCZ
#now predict rcz based on AI instead of reverse
medModel0 <- lmer(RCZ~AI+Congruency+rt + (1 | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
medModel_a <- lmer(RCZ~AI+Congruency+rt + (AI | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
medModel_b <- lmer(RCZ~AI+Congruency+rt + (Congruency | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
medModel_c <- lmer(RCZ~AI+Congruency+rt + (rt | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))
anova(medModel0,medModel_a) #p<.001,BIC=27329
anova(medModel0,medModel_b) #p=.50
anova(medModel0,medModel_c) #p<.001,BIC=27518
medModel <- lmer(RCZ~AI+Congruency+rt + (AI | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))

#outmodel is same as above
outModel <- lmer(rating_pure~AI+RCZ+Congruency+rt + (rt | Subject),data=Data,control = lmerControl(optimizer ="Nelder_Mead"))

#check diagnostics of both models
vif.mer(medModel)
plot(medModel)
qqnorm(resid(medModel))
anova(medModel)

#Fit the mediation model (takes some time)
med <- mediate(model.m = medModel, model.y = outModel,
               treat="AI", mediator="RCZ",data=medData,sims=5000,
               covariates=c("Congruency","rt"))
summary(med)

#plot these results
par(mfrow=c(1,1))
par(pin=c(2, 1.915208))
plot(med,yaxt='n',xlab=c("Estimates"),main=c('RCZ-AI-Rating'))
axis(2,at=1:3,labels=c("Total effect","Direct effect","Mediation effect"),las=1)
