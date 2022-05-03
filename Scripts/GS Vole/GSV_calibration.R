# Porsanger Analysis and Plots (abundance) -----
library(dplyr)
library(scales)
# Linear Models - only abundance - Porsanger
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Porsanger 
PorD <- readRDS("data/porsanger_mean_intervals.rds")
pormean <- dplyr::filter(PorD, species=="GRAASIDEMUS" & trapseason > 1 & !(station == "G19" & trapseason <= 6))
pormean$Abundance_HT  <- pormean$Abundance_HT/16 # divide by number of traps
pormean[,4:ncol(pormean)]<- log(pormean[,4:ncol(pormean)]+1)
r2.vec <- p.vec <- b.vec <- se.vec.b <- se.vec.int <- int.vec <- c()
mods <- list()
for(i in 1:13) 
{
  datacol <- pormean[, c(1:5, 5 + i)]
  varname <- names(datacol)[6]
  colnames(datacol)[6] = "var"

  m1 <- lm(var~Abundance_HT,data=datacol)
  s0 <- summary(m1)
  
  mods[[i]] <- m1
  r2.vec <- c(r2.vec,s0$r.squared)
  p.vec <- c(p.vec,lmp(m1))  
  b.vec <- c(b.vec,summary(m1)$coefficients[2,1] ) 
  se.vec.b <- c(se.vec.b,summary(m1)$coefficients[2,2] ) 
  se.vec.int <- c(se.vec.int,summary(m1)$coefficients[1,2] ) 
  int.vec <- c(int.vec,summary(m1)$coefficients[1,1] ) 
  #abline(b=0,a=0, lty=2)
  #abline(a=0,b=1, lty=2)
  print(i)
}

(porint <- data.frame(int=paste0("$",round(int.vec,2),"pm",round(se.vec.int,2),"$"),
                     coef=paste0("$",round(b.vec,2),"pm",round(se.vec.b,2),"$"),
                     R2=round(r2.vec,3),P=round(p.vec,3),
                     varnames=names(PorD)[6:ncol(PorD)]))

# Previous days ----- 
PorD2 <- readRDS("data/porsanger_mean_intervals_prewindow.rds")
porpre <- dplyr::filter(PorD2, species=="GRAASIDEMUS" & trapseason > 1 & !(station == "G19" & trapseason <= 6))
unique(porpre$station)
porpre$Abundance_HT  <- porpre$Abundance_HT/16 # divide by number of traps
porpre[,4:ncol(porpre)]<- log(porpre[,4:ncol(porpre)]+1)
r2.vecX <- p.vecX <- b.vecX <- se.vec.bX <- se.vec.intX <- int.vecX <- c()
modsX <- list()
for(i in 1:(ncol(porpre)-5)) 
{
  datacol <- porpre[, c(1:5, 5 + i)]
  varname <- names(datacol)[6]
  colnames(datacol)[6] = "var"
  
  m1 <- lm(var~Abundance_HT,data=datacol)
  s0 <- summary(m1)
  
  modsX[[i]] <- m1
  r2.vecX <- c(r2.vecX,s0$r.squared)
  p.vecX <- c(p.vecX,lmp(m1))  
  b.vecX <- c(b.vecX,summary(m1)$coefficients[2,1] ) 
  se.vec.bX <- c(se.vec.bX,summary(m1)$coefficients[2,2] ) 
  se.vec.intX <- c(se.vec.intX,summary(m1)$coefficients[1,2] ) 
  int.vecX <- c(int.vecX,summary(m1)$coefficients[1,1] ) 
  #abline(b=0,a=0, lty=2)
  #abline(a=0,b=1, lty=2)
  print(i)
}



(porpredf <- data.frame(int=paste0("$",round(int.vecX,2),"pm",round(se.vec.intX,2),"$"),
                        coef=paste0("$",round(b.vecX,2),"pm",round(se.vec.bX,2),"$"),
                     R2=round(r2.vecX,4),P=round(p.vecX,3),
                     varnames=names(PorD2)[6:ncol(PorD2)]))
# ACTUAL PLOT #####
pdf("Plots/R2_Reg_Porsanger.pdf", width=8*1.2,height=5*1.2)
par(mfrow=c(2,2),oma=c(0.5, 4, .5, .5))
plot(c(1,seq(3,23,2)),porint$R2[1:12], main="",
     xlab= "Time Window (days)",ylab=bquote(R^2), pch=19, col=alpha("gray70"), 
     type="b", ylim=c(0,.6), lty=2)
points(5,porint$R2[3], col=colorBlindBlack8[6], pch=19, cex=1.3)
mtext("Gray-sided Vole", cex=1.5, side=2, line = 5.5)

lines(1:12,porpredf$R2[1:12], xlab= "Time Window (days)",ylab=bquote(R^2), 
      pch=19, col="gray40", type="b", lty=2)
# points(4,porpredf$R2[4], col=alpha(4,1), pch=15, cex=1.5)

plot(pormean$int_4~pormean$Abundance_HT, pch=19, col=alpha(1,.8), 
     main=paste0(""), ylab="CT-index", xlab="CMR-Abundance")

mtext(bquote(R^2== .(round(r2.vec[4],2))), cex=1, line=0.3)

seqpred <- seq(min(pormean$Abundance_HT)-.1,max(pormean$Abundance_HT)+.1,.1)
pdf <- as.data.frame(predict(mods[[3]],newdata=data.frame(Abundance_HT=seqpred), interval = "confidence"))
polygon(c(seqpred,c(rev(seqpred))),c(pdf$lwr, rev(pdf$upr)),
        col = scales::alpha("gray50",0.5), border = NA)
abline(mods[[3]], lty=1, lwd=1.5)
dev.off()
# for the best model, do the predictions ----

cal_4apply <- function(y0,mod, mean.res=FALSE, interval="Wald")
{
  x0 <- investr::calibrate(mod,y0,interval=interval, mean.response = mean.res)
  return(x0[1:3])
}

## exponential -----
pdf("Plots/exp_porsanger.pdf", width=6,height=6)
{
  originaldata <- pormean
  # fit model
  modv <- mods[[3]]
  summary(modv)
  # predict grid
  seqx <- seq(min(originaldata$int_4)-.1, max(originaldata$int_4)+.1, 0.01)
  
  predicted_v <- matrix(unlist(t(sapply(seqx,cal_4apply, mod=modv))),ncol=3)
  #
  expredv <- exp(predicted_v)
  expseq <- exp(seqx)
  
  plot(expseq,expredv[,1], type="l", main="Gray-sided Vole",
       xlab="CT-counts (5-day interval)", ylab = "mean CMR-abundance",
       xlim=c(exp(min(originaldata$int_4)-.1),exp(max(originaldata$int_4)+.1)),
       ylim=c(0,12) # Pors
  )

  
  polygon(c(expseq,c(rev(expseq))),c(expredv[,2], rev(expredv[,3])),
          col = scales::alpha(colorBlindBlack8[6],0.7), border = NA)
  lines(expseq,exp(predicted_v[,1]), type="l", lty=1)
  lines(expseq,exp(predicted_v[,2]), type="l", lty=2)
  lines(expseq,exp(predicted_v[,3]), type="l", lty=2)
  points(exp(originaldata$int_4),exp(originaldata$Abundance_HT),pch=19, col=alpha("gray10",.8))
}
dev.off()
par(mfrow=c(1,2))
## log -----
pdf("Plots/log_preds.pdf", width=8.77,height=5) # needs also haakoya_calibration_final.R
par(mfrow=c(1,2))

# pdf("Plots/log_porsanger.pdf", width=6,height=6)
{
  originaldata <- pormean
  # fit model
  modv <- mods[[3]]
  summary(modv)
  # predict grid
  seqx <- seq(min(originaldata$int_4)-.1, max(originaldata$int_4)+.1, 0.01)
  
  predicted_v <- matrix(unlist(t(sapply(seqx,cal_4apply, mod=modv))),ncol=3)
  #
  predv <- predicted_v
  seq <- seqx
  
  plot(seq,predv[,1], type="l", 
       xlab="CT-index", ylab = "CMR-abundance",
       main="Gray-sided Vole"
       , ylim=c(-1,3), xlim=c(min(originaldata$int_4)-.1,max(originaldata$int_4)+.1)
       # xlim=c(0,1) # Pors
  )
  
  
  #mtext("log ()", side=2, line=2.5, font=1)
  #mtext("(log Growth rates)",side=2, line=2.5, cex=0.8)
  #mtext(varlab[v], side=1, line=2.5, font=1)
  #mtext("(log Growth rates)",side=1, line=3.5, cex=0.8)
  
  polygon(c(seq,c(rev(seq))),c(predv[,2], rev(predv[,3])),
          col = scales::alpha("gold",0.8), border = NA)
  lines(seq,predv[,1], type="l", lty=1)
  lines(seq,predv[,2], type="l", lty=2)
  lines(seq,predv[,3], type="l", lty=2)
  points(originaldata$int_4,originaldata$Abundance_HT,pch=19, col=alpha("gray10",.8))
}
dev.off()

# prediction ----
pred_general <-  function(dataset,ustation,mean.res=FALSE, gr=FALSE)
  # takes dataset, removes station and predicts for that station
{
  #start_time <- Sys.time()
  colnames(dataset)[4] <- "abundance"
  
  traindata <- filter(dataset, station != ustation)
  testdata <- filter(dataset, station == ustation)
  
  #fit model to traindata
  mod0 <- lm(var~abundance, data=traindata)
  predicted_x0 <- matrix(unlist(t(sapply(testdata$var,cal_4apply, mod=mod0))),ncol=3)
  #cbind(predicted_x0,testdata$abundance)
  
  coverage <- testdata$abundance<predicted_x0[,3] & testdata$abundance>predicted_x0[,2]
  if(gr==TRUE) cov_interval <- mapply(growthinterval,x0=predicted_x0[,1],y0=testdata$abundance)
  if(gr==FALSE) cov_interval <- mapply(densityinterval,x0=predicted_x0[,1],y0=testdata$abundance)
  
  Bias <- testdata$abundance-predicted_x0[,1]
  
  result <- data.frame(coverage=coverage*1, cov_interval=cov_interval,bias=Bias)
  
  return(result)
  
}



densityinterval <- function(x0,y0, min=.40, max=.75) 
{
  cond_a = abs(x0) < min & abs(y0) < min
  cond_b = x0 > max & y0 > max
  cond_c = (x0 >= min & y0 >= min) & (x0 <= max & y0 <= max)
  
  value <- ifelse(cond_a | cond_b | cond_c == TRUE, 1, 0)
  return(value)
}

library(doParallel)
unstations <- unique(originaldata$station)
originaldata <- porpre
imax=12
for(i in 1:imax)
{
  
  newGR <- originaldata[, c(1:5, 5 + i)]
  print(colnames(newGR)[6])
  
  colnames(newGR)[6] = "var"
  
  predDens <- foreach::foreach(s=1:length(unique(unstations)),.combine=rbind) %dopar% 
    pred_general(dataset=newGR,ustation=unstations[s],gr=FALSE)
  
  print(round(c(apply(predDens,2,mean),sqrt(mean(predDens$bias^2))),3))

}
