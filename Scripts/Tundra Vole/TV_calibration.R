# Tundra vole Analysis and Plots (abundance) -----

# 1. IMPORT DATA -------
# read in data already scaled by number of traps
TV1 <- readRDS("data/haakoya_cmrenc.rds")
# compute log abundance +1
TV1[,3:ncol(TV1)]<- log(TV1[,3:ncol(TV1)]+1)
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

r2.vec <- p.vec <- b.vec <- se.vec.b <- se.vec.int <- int.vec <- c()
modsE <- list()
# run models for different windows and store values
for(i in 1:13) 
{
  datacol <- TV1[, c(1:5, 5 + i)]
  varname <- names(datacol)[6]
  colnames(datacol)[6] = "var"
  
  m1 <- lm(var~D50,data=datacol)
  s0 <- summary(m1)
  
  modsE[[i]] <- m1
  r2.vec <- c(r2.vec,s0$r.squared)
  p.vec <- c(p.vec,lmp(m1))  
  b.vec <- c(b.vec,summary(m1)$coefficients[2,1] ) 
  se.vec.b <- c(se.vec.b,summary(m1)$coefficients[2,2] ) 
  se.vec.int <- c(se.vec.int,summary(m1)$coefficients[1,2] ) 
  int.vec <- c(int.vec,summary(m1)$coefficients[1,1] ) 
  print(i)
}

# -----
# data frame with relevant information from calibration regressions
# formated for LaTEX table
(TV_ENC <- data.frame(varnames=names(PorD)[6:ncol(PorD)],int=paste0("$",round(int.vec,2),"pm",round(se.vec.int,2),"$"),
                      coef=paste0("$",round(b.vec,2),"pm",round(se.vec.b,2),"$"),
                      R2=round(r2.vec,3),P=round(p.vec,3)))

## REPEAT FOR CMR-PRECEEDING DAYS ----- 
TV2 <- readRDS("data/haakoya_precmr.rds")
TV2[,3:ncol(TV2)]<- log(TV2[,3:ncol(TV2)]+1)
r2.vecX <- p.vecX <- b.vecX <- se.vec.bX <- se.vec.intX <- int.vecX <- c()
modsX2 <- list()
i=1
# run models for different windows and store values
for(i in 1:(ncol(TV2)-5)) 
{
  datacol <- TV2[, c(1:5, 5 + i)]
  varname <- names(datacol)[6]
  colnames(datacol)[6] = "var"
  
  m1 <- lm(var~D50,data=datacol)
  s0 <- summary(m1)
  
  modsX2[[i]] <- m1
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

# data frame with relevant information from calibration regressions
# formated for LaTEX table
(TVPRE <- data.frame(int=paste0("$",round(int.vecX,2),"pm",round(se.vec.intX,2),"$"),
                        coef=paste0("$",round(b.vecX,2),"pm",round(se.vec.bX,2),"$"),
                        R2=round(r2.vecX,4),P=round(p.vecX,3),
                        varnames=names(TV2)[6:ncol(TV2)]))
# 3. PLOTS -----

# pdf("Plots/R2_Reg_Hakoya.pdf", width=8*1.2,height=5*1.2)
#par(mfrow=c(1,2),oma=c(0.5, 4, .5, .5))

# plot R2 as a function of window
plot(c(1,seq(3,23,2)),TV_ENC$R2[1:12], main="",
     xlab= "Time Window (days)",ylab=bquote(R^2), pch=19, col=alpha("gray70"), 
     type="b", ylim=c(0,.6), lty=2)

lines(1:12,TVPRE$R2[1:12], xlab= "Time Window (days)",ylab=bquote(R^2), 
      pch=19, col="gray40", type="b", lty=2)
points(1,TVPRE$R2[1], col=4, pch=19, cex=1.5)
mtext("Tundra Vole", cex=1.5, side=2, line = 5.5)

# plot scatterplot + regression for optimized window
plot(TV2$int_0~TV2$D50, pch=19, col=alpha(1,.8), 
     main=paste0(""), ylab="CT-index", xlab="CMR-Abundance")

mtext(bquote(R^2== .(round(r2.vecX[1],2))), cex=1, line=0.3)

seqpredx <- seq(min(TV1$D50)-.1,max(TV1$D50)+.3,.1)
pdf <- as.data.frame(predict(modsX2[[1]],newdata=data.frame(D50=seqpredx), interval = "confidence"))
polygon(c(seqpredx,c(rev(seqpredx))),c(pdf$lwr, rev(pdf$upr)),
        col = scales::alpha("gray50",0.5), border = NA)
abline(modsX2[[1]], lty=1, lwd=1.5)
# dev.off()

# 4. PREDICTIONS FOR OPTIMIZED WINDOW ----

cal_4apply <- function(y0,mod, mean.res=FALSE, interval="Wald")
  # Apply calibration function from INVESTR PACKAGE
{
  x0 <- investr::calibrate(mod,y0,interval=interval, mean.response = mean.res)
  return(x0[1:3])
}

#pdf("Plots/log_hakoya.pdf", width=8.77,height=5)
# pdf("Plots/log_hakoya.pdf", width=8.77,height=5)
{
  originaldata <- TV1
  # fit model
  modv <- modsX2[[1]]
  summary(modv)
  # predict grid
  seqx <- seq(min(originaldata$int_0)-.1, max(originaldata$int_0)+.1, 0.01)
  
  predicted_v <- matrix(unlist(t(sapply(seqx,cal_4apply, mod=modv))),ncol=3)
  #
  predv <- predicted_v
  seq <- seqx
  
  plot(seq,predv[,1], type="l",main="Tundra Vole",
       xlab="CT-index", ylab = "CMR-abundance",
       ylim=c(-1,3), xlim=c(min(originaldata$int_0)-.1,max(originaldata$int_0)+.1),
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
  points(originaldata$int_0,originaldata$D50,pch=19, col=alpha("gray10",.8))
}
#dev.off()

# 5. PREDICTION METRICS (K-FOLD CROSS VALIDATION) ----

# prediction function for a new station
pred_general <-  function(dataset,ustation,mean.res=FALSE)
  # takes dataset, removes station and predicts for that station
{
  #start_time <- Sys.time()
  colnames(dataset)[4] <- "abundance"
  
  # separate train and test set
  traindata <- filter(dataset, station != ustation)
  testdata <- filter(dataset, station == ustation)
  
  #fit model to traindata
  mod0 <- lm(var~abundance, data=traindata)
  # predict
  predicted_x0 <- matrix(unlist(t(sapply(testdata$var,cal_4apply, mod=mod0))),ncol=3)
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

ECM <- function(x0,y0, min=.28, max=.85) 
  # quantiles min and max chosen according to quantiles of distribution
{
  cond_a = abs(x0) < min & abs(y0) < min
  cond_b = x0 > max & y0 > max
  cond_c = (x0 >= min & y0 >= min) & (x0 <= max & y0 <= max)
  value <- 0
  if(cond_a == TRUE) value = 1
  
  value <- ifelse(cond_a | cond_b | cond_c == TRUE, 1, 0)
  return(value)
}

# run k-fold cross validation in parallel
library(doParallel)
# adjust dataset
originaldata <- TV1
unstations <- unique(originaldata$station)
imax=12
for(i in 1:imax)
{
  
  newGR <- originaldata[, c(1:5, 5 + i)]
  print(colnames(newGR)[6])
  
  colnames(newGR)[6] = "var"
  colnames(newGR)[4] = "abundance"
  
  predDens <- foreach::foreach(s=1:length(unique(unstations)),.combine=rbind) %dopar% 
    pred_general(dataset=newGR,ustation=unstations[s])
  
  print(round(c(apply(predDens,2,mean),sqrt(mean(predDens$bias^2))),3))
  # print(sqrt(mean(predDens$bias^2)))
  
}

summary(lm(TV2$int_1~TV2$D50))


# 6. Aggregated Analysis Tundra Vole & Plots #----
par(mfrow=c(1,1))
# pdf("Plots/hakoya_aggregated_lm.pdf", height=6,width=6)
originaldata <- TV2
# choose optimized window
# aggregate station data  by trap session
od2 <- aggregate(cbind(D50,int_0)~trapsession, data=originaldata, FUN=mean)

mx <- lm(int_0~D50, data=od2)
smx <- summary(mx)

## Plot ----
plot(od2$int_0~od2$D50, pch=19, col="gray30",cex=1.2, xlim=c(0.15,1.2), ylim=c(0.2,2.5),
     ylab="Aggregated CT-index",xlab="Aggregated CMR-Abundance", main="Tundra Vole")
mtext(bquote(R^2== .(round(smx$r.squared,2))), cex=1, line=0.3)
seqpred <- seq(min(od2$D50)-.1,max(od2$D50)+.3,.1)
pdf <- as.data.frame(predict(mx,newdata=data.frame(D50=seqpred), interval = "confidence"))
polygon(c(seqpred,c(rev(seqpred))),c(pdf$lwr, rev(pdf$upr)),
        col = scales::alpha("gray50",0.5), border = NA)
points(od2$int_0~od2$D50, pch=19, col="gray20",cex=1.2)
abline(mx, lty=1, lwd=1.5)
# dev.off()
