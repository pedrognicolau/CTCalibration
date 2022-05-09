# GS Vole Analysis and Plots (abundance) -----
library(dplyr)
library(scales)
# Linear Models - only abundance - Porsanger


# 1. IMPORT DATA -------
# Porsanger 
GSV1 <- readRDS("data/porsanger_cmrenc.rds")
# remove censored data as described in appendix
GSV2 <- dplyr::filter(GSV1, trapseason > 1 & !(station == "G19" & trapseason <= 6))
# divide by number of traps
GSV2$Abundance_HT  <- GSV2$Abundance_HT/16 
# compute log abundance +1
GSV2[,4:ncol(GSV2)]<- log(GSV2[,4:ncol(GSV2)]+1)

r2.vec <- p.vec <- b.vec <- se.vec.b <- se.vec.int <- int.vec <- c()
mods <- list()

# 2. RUN CALIBRATION MODELS FOR DIFFERENT TIME WINDOWS ------
## CMR-ENCOMPASSING ----
# compute p-value function
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
# run models for different windows and store values
for(i in 1:13) 
{
  datacol <- GSV2[, c(1:5, 5 + i)]
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
}

# -----
# data frame with relevant information from calibration regressions
# formated for LaTEX table
(CMRENC <- data.frame(varnames=names(GSV2)[6:ncol(GSV2)],int=paste0("$",round(int.vec,2),"pm",round(se.vec.int,2),"$"),
                      coef=paste0("$",round(b.vec,2),"pm",round(se.vec.b,2),"$"),
                      R2=round(r2.vec,3),P=round(p.vec,3)))

## REPEAT FOR CMR-PRECEEDING DAYS ----- 
GSV3 <- readRDS("data/porsanger_precmr.rds")
GSV4 <- dplyr::filter(GSV3, trapseason > 1 & !(station == "G19" & trapseason <= 6))
GSV4$Abundance_HT  <- GSV4$Abundance_HT/16 # divide by number of traps
GSV4[,4:ncol(GSV4)]<- log(GSV4[,4:ncol(GSV4)]+1)
r2.vecX <- p.vecX <- b.vecX <- se.vec.bX <- se.vec.intX <- int.vecX <- c()
modsX <- list()
# run models for different windows and store values
for(i in 1:(ncol(GSV4)-5)) 
{
  datacol <- GSV4[, c(1:5, 5 + i)]
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
}

# data frame with relevant information from calibration regressions
# formated for LaTEX table
(CMRPRE <- data.frame(int=paste0("$",round(int.vecX,2),"pm",round(se.vec.intX,2),"$"),
                        coef=paste0("$",round(b.vecX,2),"pm",round(se.vec.bX,2),"$"),
                     R2=round(r2.vecX,4),P=round(p.vecX,3),
                     varnames=names(GSV4)[6:ncol(GSV4)]))
# 3. PLOTS -----
# pdf("Plots/R2_Reg_Porsanger.pdf", width=8*1.2,height=5*1.2)

# plot R2 as a function of window
par(mfrow=c(2,2),oma=c(0.5, 4, .5, .5))
plot(c(1,seq(3,23,2)),CMRENC$R2[1:12], main="",
     xlab= "Time Window (days)",ylab=bquote(R^2), pch=19, col=alpha("gray70"), 
     type="b", ylim=c(0,.6), lty=2)
points(5,CMRENC$R2[3], col=4, pch=19, cex=1.3)
mtext("Gray-sided Vole", cex=1.5, side=2, line = 5.5)

lines(1:11,CMRPRE$R2, xlab= "Time Window (days)",ylab=bquote(R^2), 
      pch=19, col="gray40", type="b", lty=2)

# plot scatterplot + regression for optimized window
plot(GSV2$int_4~GSV2$Abundance_HT, pch=19, col=alpha(1,.8), 
     main=paste0(""), ylab="CT-index", xlab="CMR-Abundance")

mtext(bquote(R^2== .(round(r2.vec[4],2))), cex=1, line=0.3)

seqpred <- seq(min(GSV2$Abundance_HT)-.1,max(GSV2$Abundance_HT)+.1,.1)
pdf <- as.data.frame(predict(mods[[3]],newdata=data.frame(Abundance_HT=seqpred), interval = "confidence"))
polygon(c(seqpred,c(rev(seqpred))),c(pdf$lwr, rev(pdf$upr)),
        col = scales::alpha("gray50",0.5), border = NA)
abline(mods[[3]], lty=1, lwd=1.5)
# dev.off()
# 4. PREDICTIONS FOR OPTIMIZED WINDOW ----

cal_4apply <- function(y0,mod, mean.res=FALSE, interval="Wald")
  # Apply calibration function from INVESTR PACKAGE
{
  x0 <- investr::calibrate(mod,y0,interval=interval, mean.response = mean.res)
  return(x0[1:3])
}

# Obtain prediction intervals ----
#  pdf("Plots/log_preds.pdf", width=8.77,height=5) 
par(mfrow=c(1,2))

# pdf("Plots/log_porsanger.pdf", width=6,height=6)
{
  # optimized variable is from CMR-ENCOMPASSING data - GSV2
  originaldata <- GSV2
  # highest R^2 is model #3
  modv <- mods[[3]]
  # prediction grid
  seqx <- seq(min(originaldata$int_4)-.1, max(originaldata$int_4)+.1, 0.01)
  # predict values to obtain prediction intervals
  predicted_v <- matrix(unlist(t(sapply(seqx,cal_4apply, mod=modv))),ncol=3)
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
# dev.off()

# 5. PREDICTION METRICS (K-FOLD CROSS VALIDATION) ----

# prediction function
pred_general <-  function(dataset,ustation,mean.res=FALSE)
  # takes dataset, removes station and predicts for that station 'ustation'
{
  #start_time <- Sys.time()
  colnames(dataset)[4] <- "abundance"
  
  # separate train and test set
  traindata <- filter(dataset, station != ustation)
  testdata <- filter(dataset, station == ustation)
  
  # fit model to traindata
  mod0 <- lm(var~abundance, data=traindata)
  # predict
  predicted_x0 <- matrix(unlist(t(sapply(testdata$var,cal_4apply, mod=mod0))),ncol=3)
  #cbind(predicted_x0,testdata$abundance)
  # compute prediction metrics
  coverage <- testdata$abundance<predicted_x0[,3] & testdata$abundance>predicted_x0[,2]
  # compute ECM - predict correct abundance category
  cov_interval <- mapply(ECM,x0=predicted_x0[,1],y0=testdata$abundance)
  # bias
  Bias <- testdata$abundance-predicted_x0[,1]
  # lump results
  result <- data.frame(coverage=coverage*1, cov_interval=cov_interval,bias=Bias)
  
  return(result)
  
}


# ECM function
ECM <- function(x0,y0, min=0.12, 
                max=0.87)
  # quantiles min and max chosen according to quantiles of distribution
{
  cond_a = abs(x0) < min & abs(y0) < min
  cond_b = x0 > max & y0 > max
  cond_c = (x0 >= min & y0 >= min) & (x0 <= max & y0 <= max)
  
  value <- ifelse(cond_a | cond_b | cond_c == TRUE, 1, 0)
  return(value)
}

# run k-fold cross validation in parallel
library(doParallel)
# adjust dataset
originaldata <- GSV4
unstations <- unique(originaldata$station)
imax=12
for(i in 1:imax)
{
  
  newGR <- originaldata[, c(1:5, 5 + i)]
  print(colnames(newGR)[6])
  
  colnames(newGR)[6] = "var"
  
  predDens <- foreach::foreach(s=1:length(unique(unstations)),.combine=rbind) %dopar% 
    pred_general(dataset=newGR,ustation=unstations[s])
  
  print(round(c(apply(predDens,2,mean),sqrt(mean(predDens$bias^2))),3))

}
