#####  REQUIREMENTS ####

library(vegan)
library(pROC)
library(rLims)

#Put path were output should be written here
prefix <- ""
while (prefix == ""){
  prefix <- readline("Where should the plots be stored? (e.g. /home/moi/put_images_here/)")
}
if (substr(prefix,nchar(prefix),nchar(prefix)) != "/"){
    prefix <- paste(prefix, "/", sep="")
}
prefix <- paste(prefix,"Coffee_Origin_IRvNMR_Images/", sep="")


ginv<-function(X, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}


Deriv1 <- function(x,y){
  y.prime <- diff(y) / diff(x)
  x.prime <- x[-length(x)] + diff(x)/2
  list(shift = x.prime,
       intensity = y.prime)
}

#function for computation of second derivative
Deriv2 <- function(x,y){
  h <- x[2] - x[1]
  Range <- 2:(length(x)-1)  # Drop first and last points
  list(shift = x[Range],
       intensity = (y[Range+1] - 2*y[Range] + y[Range-1]) / h^2)
}

msc <- function(spectrum, ...){
  UseMethod("msc", spectrum)    
}

msc.numeric <- function(spectrum, reference){
  Rm <- mean(reference)
  Sc <- spectrum - mean(spectrum);
  Rc <- reference - Rm;
  
  m = ginv(crossprod(Rc)) %*% crossprod(Rc, Sc);
  Sc / m + Rm
}

msc.matrix <- function(spectra, reference = apply(spectra, 2, median)){
  t(apply(spectra, 1, function(x) msc.numeric(x, reference = reference)));
}

integral.normalization<- function(spectra, value = 100){
  t(apply(spectra, 1, function(x) (x * value) / sum(abs(x))))
}

##Function for plotting Q2_average vs. Number of components with std. dev. error bars
plot.with.error.bars <- function(x, matrix, col = "black",
                                 ylim = range(c(apply(matrix, 2, mean) - apply(matrix, 2, sd),
                                                apply(matrix, 2, mean) + apply(matrix, 2, sd))), ...){
  avg <- apply(matrix,2,mean)
  sdev <- apply(matrix,2,sd)
  plot(x, avg, col = col, ylim = ylim, ...)
  arrows(x, avg - sdev, x, avg + sdev, length=0.05, angle=90, code=3, col = col)
} 

#Returns the indices of the top values of x that encompass the upper fraction of its total magnitude according to the given norm
upper <- function(x, fraction = 1, norm = function(x) sum(abs(x))){
  ordering <- order(x, decreasing = TRUE)
  n <- length(x)
  total <- norm(x)
  res <- c()
  partial <- 0
  partial.fraction <- 0
  for (index in ordering){
    if (partial.fraction < fraction){
      res <- c(res, index)
      partial = partial + x[index]
      partial.fraction <- partial / total
    }
    else{
      break;
    }
  }
  res
}

##Environment variables to determine how to compute certain things (normalization, centering and scaling, validation, number of components)
validating <- "kfold"  
valsize <- 7   
nComp <- 15
nModels <- 100
normalizing <- c(integral.normalization,msc,identity)
centering <- c(FALSE,FALSE,FALSE)
scaling <- c(TRUE,TRUE,TRUE) #c(FALSE,FALSE,FALSE)
centering1<-c(TRUE,TRUE,TRUE)
scaling1<- c(TRUE,TRUE,TRUE)
model <- c(lims.opls.model, lims.pls.model, lims.opls.model)
names <- c("IR", "NIR", "NMR")
color<-c("gray60","2","4")

nfigs <- 0

### READ DATA ###

##read metadata

metadata <- read.table("./data/metadata.csv",
                       header=TRUE,sep="\t")[1:317,]

#Change labels to -1 = Arabica, 1 = Robusta
metadata$species1<-metadata$species
levels(metadata$species1)[1]<-"Arabica"
levels(metadata$species1)<-c("1", "-1")

#set up 4 labels for scores plot: Robusta, Arabica Colombia, Arabica Brazil (neighbor), Arabica Peru (neighbor aussi),
#Other Arabicas
#groups<-factor(metadata$species1)
metadata$groups<-metadata$country
levels(metadata$groups)[levels(metadata$groups) %in% metadata$country[metadata$species == "Robusta"]] = "Robusta"
levels(metadata$groups)[!levels(metadata$groups) %in% c("Colombia", "Peru", "Brazil","Robusta")] = "Others"

## IR

#read data
irdata <- as.matrix(read.table("./data/mIRspectra.csv",
                               header=FALSE,sep=" ")[1:317,])
#nu scale for spectra
dx <- 0.964233
irscale <- seq(649.893311,3999.640625,dx)

###Data preprocessing
#filter samples with spectra in absorbance units
irdata <- irdata[!metadata$ID %in% c("AC1199", "AC1242", "AC1376", "AC1380", "AC1427"),]
metadata <- metadata[!metadata$ID %in% c("AC1199", "AC1242", "AC1376", "AC1380", "AC1427"),]

#remove duplicates
irdata <- irdata[!duplicated(metadata$ID),]
metadata <- metadata[!duplicated(metadata$ID),]

#sort
irdata <- irdata[order(metadata$ID),]
metadata <- metadata[order(metadata$ID),]

#spectra second derivatives
irdata.deriv2 <- t(apply(irdata,1,function(y) Deriv2(irscale,y)$intensity))
#irdata.deriv2<-irdata

#select nus and normalize
irdata.deriv2.wavelength.optimized1 <- irdata.deriv2[,(irscale > 800 & irscale < 1800)[2:3474]]#t(scale(t(irdata.deriv2[,(irscale>800 & irscale<1800)[2:3474]])))
irscale<- irscale[irscale>800 & irscale<1800]
irdata.deriv2.wavelength.optimized1 <- normalizing[[1]](irdata.deriv2.wavelength.optimized1)
#irdata.deriv2.wavelength.optimized1 <- t(lims.scaling(t(irdata.deriv2.wavelength.optimized1), center = FALSE,
#                                                        scale = function(x) max(abs(x))))
###   NIR    ###

#read data
id<-as.vector(as.matrix(read.csv("./data/NIRSamples.csv")))
dataNIR<-as.matrix(read.csv("./data/NIRspectra.csv", sep=""))
nirscale<-as.vector(as.matrix(read.csv("./data/NIRscale.csv", sep="", header=FALSE)))
#metadata <- read.table("~/CoffeeIR/FT-IR_ JCAMP-DX and .SPA files/sampleMetadata.csv",
# header=TRUE,sep="\t")[1:317,]

#remove duplicates
dataNIR <- dataNIR[!duplicated(id),]
id <- id[!duplicated(id)]

#filter common samples (with IR)
filter<-which(id %in% as.vector(metadata$ID),T)
dataNIR<-dataNIR[filter,]
id <- id[filter]
#nir.metadata<-metadata[which(id.selected%in%metadata$ID,T),]

#sort
dataNIR <- dataNIR[order(id),]
id <- id[order(id)]

#compute second derivative
nirdev2<-t(apply(dataNIR,1, function(y) Deriv2(nirscale,y)$intensity))

#select nus and normalize
nirdev2<-nirdev2[,1:900] ###scale
nirdev2 <- normalizing[[2]](nirdev2)
#nirdev2 <- t(lims.scaling(t(nirdev2),center = FALSE, scale = function(x) max(abs(x))))

#scale for second derivative
nirdev2scale<-as.numeric(nirscale[2:901])

##Read NMR

load("./data/dataNMR.rda")
nmr.spectra <- res$nmrData #vaya nombre...
nmr.id <- res$param$catalogID

#Select roasted samples
nmr.spectra <- nmr.spectra[res$param$presentation == "t",]
nmr.id <- nmr.id[res$param$presentation == "t"]

#Remove duplicates
nmr.spectra <- nmr.spectra[!duplicated(nmr.id),]
nmr.id <- nmr.id[!duplicated(nmr.id)]

#Select common data (should be all in IR)
filter <- which(nmr.id %in% as.vector(metadata$ID),T)
nmr.spectra <- nmr.spectra[filter,]
nmr.id <- nmr.id[filter]

#Back-filter: restrict metadata, IR, NIR to samples that were also in NMR
irdata.deriv2.wavelength.optimized1 <- irdata.deriv2.wavelength.optimized1[metadata$ID %in% nmr.id,]
nirdev2 <- nirdev2[metadata$ID %in% nmr.id,]
metadata <- metadata[metadata$ID %in% nmr.id,]

#Sort
nmr.spectra <- nmr.spectra[order(nmr.id),]
nmr.id <- nmr.id[order(nmr.id)]

#Normalize
nmr.spectra <- normalizing[[3]](nmr.spectra)
#nmr.spectra <- t(lims.scaling(t(nmr.spectra),center=FALSE,scale=function(x) max(abs(x))))

#Remove outliers that were detected on PCA
#AC1229 is outlier for IR, the rest were outliers for NMR
#1227 AC1290 "AC1339" nmr outliers
nmr.spectra <- nmr.spectra[!metadata$ID %in% c("AC1441", "AC1431", "AC1445", "AC1414", "AC1169", "AC1403",
                                               
                                               
                                               "AC1207", "AC1420", "AC1337", "AC1229",
                                               "AC1227", "AC1290", "AC1339"),]


nmr.id <- nmr.id[!metadata$ID %in% c("AC1441", "AC1431", "AC1445", "AC1414", "AC1169", "AC1403",
                                     
                                     
                                     "AC1207", "AC1420", "AC1337", "AC1229",
                                     "AC1227", "AC1290", "AC1339")]

# ### plot outliers #####
# 
# matplot(res$ppm, t(nmr.spectra[1:10,]), type="l", xlim=c(3,1), ylim=c(0,3e9))



irdata.deriv2.wavelength.optimized1 <-
  irdata.deriv2.wavelength.optimized1[!metadata$ID %in% c("AC1441", "AC1431", "AC1445", "AC1414", "AC1169", "AC1403",
                                                          "AC1207", "AC1420", "AC1337", "AC1229",
                                                          "AC1227", "AC1290", "AC1339"),]
nirdev2 <- nirdev2[!metadata$ID %in% c("AC1441", "AC1431", "AC1445", "AC1414", "AC1169", "AC1403",
                                       "AC1207", "AC1420", "AC1337", "AC1229",
                                       "AC1227", "AC1290", "AC1339"),]
metadata <- metadata[!metadata$ID %in% c("AC1441", "AC1431", "AC1445", "AC1414", "AC1169", "AC1403",
                                         "AC1207", "AC1420", "AC1337", "AC1229",
                                         
                                         "AC1227", "AC1290", "AC1339"),]
##ASSEMBLE DATA

#c<- which(metadata$species=="Arabica",T)[1:12]
#c1<-which(metadata$species=="100% ARABICA COLOMBIA",T)[1:12]
#d<-which(metadata$species=="Robusta",T)
#F<-c(c,c1,d)
#irFAR<-cbind(metadata$species[F],irdata.deriv2.wavelength.optimized1[F,])
#write.table(nmrFAR, file="nmrARF.csv")

#data <-  c(list(irdata.deriv2.wavelength.optimized1), list(nirdev2), list(nmr.spectra)) #list(nirdev2)
osc <- osc.lims(as.matrix((nirdev2)), as.matrix(as.numeric(metadata$species1)), nComp=1, training=1:101,
                xcenter = centering[[2]], xscale = scaling[[2]], ycenter = FALSE, 
                yscale = FALSE, accuracy = 1e-5, maxit = 100, 
                discriminant = FALSE)

data <- c(list(irdata.deriv2.wavelength.optimized1), list(osc$X), list(nmr.spectra))
scales <- c(irscale, nirdev2scale, res$ppm)
#y<-normalization(nmr.spectra, method = "pqn")
#data <- c(list(y$newXtrain))
#data <- c(list(osc$X))
#metadata<-metadata[F,]

#matplot(t(y$newXtrain), type="l")

### DATA PROCESSING ####
AR.Q2 <- c()
AR.ROC <- list(list("predicted" = c(), "observed" = c()),list("predicted" = c(), "observed" = c()),list("predicted" = c(), "observed" = c()))
CO.ROC <- list(list("predicted" = c(), "observed" = c()),list("predicted" = c(), "observed" = c()),list("predicted" = c(), "observed" = c()))
AR.Q2.densities <- c()
CO.Q2 <- c()
CO.Q2.densities <- c()
nComp <- 15

for (i in 1:3){
  predictors <- data[[i]]  
  
  ## ARABICA VS. ROBUSTA
  #tiff(paste(path,names[i],"1"))
  nfigs = nfigs + 1
  postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=10,width=8, horizontal=FALSE)
 # png(file = paste(prefix,nfigs,".png"), bg = "white")
  plot(c(), main = c(names[[i]], "Arabica vs. Robusta"), xlim=c(0,0), ylim=c(0,0))
  dev.off()
  #   ## PCA to detect outliers
  #   pca <- prcomp(predictors,scale=scaling[[i]],center=centering[[i]])
  #   plot(pca$x[,1],pca$x[,2],
  #        col=c("grey66", "black", "black", "black", "black")[as.numeric(metadata$groups)],
  #        pch=c(15,1,15,7,22)[as.numeric(metadata$groups)], #pch=c(16,15)[metadata$groups],
  #        xlab= paste("1st PC (", round((pca$sdev^2 / sum(pca$sdev^2))[1] * 100,1), "%)"), 
  #        ylab= paste("2nd PC (", round((pca$sdev^2 / sum(pca$sdev^2))[2] * 100,1), "%)"),
  #        cex=1.2, cex.axis=1.0,cex.lab=1.3)  ##png and pdf
  #   text(pca$x[,1:2], label=metadata$ID,cex=0.7, pos=1, col="red")
  #   abline(h=0);abline(v=0) 
  #   ordiellipse(cbind(pca$x[,1], pca$x[,2]),
  #               metadata$species1,
  #               conf=0.95, col="gray66", lwd=2)
  
  #saveDevice("pca")
  #}
  
  ##Model and validation
  
  #Next two lines: select just as many arabicas as there are robustas. The scores plots are expected to go
  ##gaga but everything else should be fine
  original.predictors <- predictors;
  selection <- sort(c(which(metadata$species1 == -1), sample(which(metadata$species1 == 1),length(which(metadata$species1 == -1)))))
  predictors <- predictors[selection,]
  ir.validation <-lims.validate(model = model[[i]], predictors = predictors, 
                                responses = as.matrix(as.numeric(as.character(metadata$species1)[selection])), 
                                nModels = nModels, nComp=nComp, ycenter = FALSE, yscale = FALSE,
                                codomain = c(-1,1), xscale = scaling[[i]], xcenter = centering[[i]], valsize = valsize,
                                method=validating)
  
  #   y1<-ir.validation$training.set[[1]]
  #   y2<-as.matrix(as.numeric(as.character(metadata$species1)))[y1,]
  #   X<-osc$X[y1,]
  #   X1<-irdata.deriv2.wavelength.optimized1[y1,]
  #   X2<-nmr.spectra[y1,]
  #   Xtrain = list(X,X1,X2);#list(X,X2,X3);
  #   Xtest<-osc$X[-y1,]
  #   X1test<-irdata.deriv2.wavelength.optimized1[-y1,]
  #   X2test<-nmr.spectra[-y1,]
  #   XTest = list(Xtest,X1test, X2test);
  #   m=c(dim(X)[2],dim(X1)[2], dim(X2)[2])
  #   multi<-multiblock(Xtrain,y2,m,3)
  #   
  #   ypred=predict1(Xtrain,multi);
  #   print(cbind(y2,ypred))
  #   R2 <- cov(y2,ypred)^2/(var(y2)*var(ypred))
  #   print(cbind("R2",R2))
  
  #   x <- Xtest
  #   Ytest<-as.matrix(as.numeric(as.character(metadata$species1)))[-y1,]
  #   ypred = predict1(x,multi) #Abdi verion 
  #   Q2 <- cov(Ytest,ypred)^2/(var(Ytest)*var(ypred))
  #   print(cbind("Q2: ",Q2))
  #   print("Confusion matrix test XWQ")
  # print(rbind(cbind("R/P",  "T", " F "),cbind("T  ",colSums(ypred>0&response>0),colSums(ypred<0&response>0)),cbind("F  ",colSums(ypred>0&response<0),colSums(ypred<0&response<0))))
  
#     ## Compute and save ROC curve (for all accumulated predictions) and AUCs for individual ROC curves
#     Pred1<-c(); W1<-c()
#     for (i in c(1:length(ir.validation$full.results))){ 
#       pred1<-(ir.validation$full.results[[i]]$predicted.values[,,nComp])
#       #  pred<-as.factor(as.numeric(ir.validation$full.results[[i]]$predicted.values[,,nComp]>0))
#       #levels(pred)<-c("-1","1")
#       #pred<-as.numeric(as.character(pred))
#       w1<-as.numeric(as.character(metadata$species1[-ir.validation$training.set[[i]]]))
#       
#       Pred1<-c(Pred1,pred1)
#       W1<-c(W1,w1)
#     }
#     
#     roc1 <- roc(W1,
#                 Pred1, percent=TRUE,
#                 # arguments for auc
#                 partial.auc=c(100, 90), partial.auc.correct=TRUE,
#                 partial.auc.focus="sens",
#                 # arguments for ci
#                 ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
#                 # arguments for plot
#                 plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
#                 print.auc=TRUE, show.thres=TRUE, legacy.axes = TRUE)
#     AR.ROC <- c(AR.ROC,roc1)
    
  
  ##Select optimum number of components
  
  #Plot Q2 vs. Nº of components to choose optimum
  plot.with.error.bars(1:nComp,ir.validation$Q2, pch=19, xlab = "Number of components", ylab = "Q2")
 # saveDevice(paste(names[i],"Q2" ))
  #Choose best number of components from the previous plot
  best.nComp <- 0
  while (best.nComp < 1){
    best.nComp <- readline("Choose number of components to proceed")
  }
  best.nComp <- as.numeric(best.nComp)
  
  ##Save data for ROC analysis
  j <- 0
  for (results in ir.validation$full.results){
    j = j + 1
   # AR.ROC[[i]]$predicted <- rbind(AR.ROC[[i]]$predicted, sign(results$predicted.values[,,best.nComp]))
    AR.ROC[[i]]$predicted <- rbind(AR.ROC[[i]]$predicted, (results$predicted.values[,,best.nComp][1:4]))
    AR.ROC[[i]]$observed <- rbind(AR.ROC[[i]]$observed, 
                                  as.numeric(as.character(metadata$species1[selection][-ir.validation$training.set[[j]]][1:4])))
  }
  
#   ARNIR<-c()
#   W<-c(); Pred<-c()
#   for (j in c(1:length(ir.validation$full.results))){ 
#     pred<-(ir.validation$full.results[[j]]$predicted.values[,,best.nComp])
#     #  pred<-as.factor(as.numeric(ir.validation$full.results[[i]]$predicted.values[,,nComp]>0))
#     #levels(pred)<-c("-1","1")
#     #pred<-as.numeric(as.character(pred))
#     w<-as.numeric(as.character(metadata$species1[-ir.validation$training.set[[j]]]))
#     #print(w)
#     #W<-c(W,w)
#     #Pred<-c(Pred,pred)
#     AR <- roc(w,
#               pred, percent=TRUE,
#               # arguments for auc
#               # partial.auc=c(100, 90), partial.auc.correct=TRUE,
#               # partial.auc.focus="sens",
#               # arguments for ci
#               ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
#               # arguments for plot
#               plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
#               print.auc=TRUE, show.thres=TRUE, legacy.axes = TRUE)
#     ARNIR<-c(ARNIR,list(AR))
#     
#   }
  ##Save Q2 for latter comparative plotting
  AR.Q2 <- c(AR.Q2, list(ir.validation$Q2))
  
  ##Compute Q2.density and select example model
  Q2.density <- density(ir.validation$Q2[,best.nComp])
  #  plot(Q2.density, xlab = "Q2 for optimum number of components") #Just for checking, eventually all Q2 densities will be in a single plot
  selected.model <- which.min(abs(ir.validation$Q2[,best.nComp] - Q2.density$x[which.max(Q2.density$y)])) #we draw
  best.model <- which.max(ir.validation$Q2[,best.nComp])
  
  ##Save Q2.density for latter comparative plotting
  AR.Q2.densities <- c(AR.Q2.densities, list(Q2.density))
  
  ##scores plot of model with Q2 in the maximum of the distribution
  par(mar=c(5.1,4.5,4,5.1))
  scoresX <- c()
  scoresY <- c()
  labelY <- ""
  if (identical(model[[i]], lims.pls.model)){
    scoresX <- ir.validation$full.results[[selected.model]]$xscores[,1]
    scoresY <- ir.validation$full.results[[selected.model]]$xscores[,2]
    labelY <- "Predictive Component"
  }
  else{
    if (identical(model[[i]], lims.opls.model)){
      scoresX <- ir.validation$full.results[[selected.model]]$xscores[,best.nComp]
      scoresY <- ir.validation$full.results[[selected.model]]$orthoscores[,1]
      labelY <- "Orthogonal Component"
    }
  }
  
  nfigs = nfigs + 1
  postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=10,width=8, horizontal=FALSE)
 # png(file = paste(prefix,nfigs,".png"), bg = "white")
  plot(scoresX,scoresY,
       xlab="Predictive Component",
       ylab=labelY, cex.lab=1.2,
       main = paste("Typical model: R2 = ",
                    substr(as.character(ir.validation$full.results[[selected.model]]$R2[best.nComp]),
                           start=1,stop=6), 
                    " Q2 = ",
                    substr(as.character(ir.validation$full.results[[selected.model]]$Q2[best.nComp]),
                           start=1,stop=6)),
       #squares = arabica, circles = robusta
       #black squares = colombia
       #grey squares = brazil
       #clear squares = peru
       #x squares = other arabicas
       col=c("grey66", "black", "black", "black",
             "black")[as.numeric(metadata$groups)[unlist(ir.validation$training.set[selected.model])]],
       pch=c(15,1,15,7,22)[as.numeric(metadata$groups)[unlist(ir.validation$training.set[selected.model])]])
  
  abline(h=0);abline(v=0) 
  ordiellipse(cbind(scoresX, scoresY), 
              metadata$species1[unlist(ir.validation$training.set[selected.model])],conf=0.95, col="gray66", lwd=2)
  
  dev.off()
  
  ##scores plot of model with best Q2
#   if (identical(model[[i]], lims.pls.model)){
#     scoresX <- ir.validation$full.results[[best.model]]$xscores[,1]
#     scoresY <- ir.validation$full.results[[best.model]]$xscores[,2]
#   }
#   else{
#     if (identical(model[[i]], lims.opls.model)){
#       scoresX <- ir.validation$full.results[[best.model]]$xscores[,best.nComp]
#       scoresY <- ir.validation$full.results[[best.model]]$orthoscores[,1]
#     }
#   }
  
  par(mar=c(5.1,4.5,4,5.1))
  
  if (identical(model[[i]], lims.pls.model)){
    scoresX <- ir.validation$full.results[[best.model]]$xscores[,1]
    scoresY <- ir.validation$full.results[[best.model]]$xscores[,2]
  }
  else{
    if (identical(model[[i]], lims.opls.model)){
      scoresX <- ir.validation$full.results[[best.model]]$xscores[,best.nComp]
      scoresY <- ir.validation$full.results[[best.model]]$orthoscores[,1]
    }
  }
  
  nfigs = nfigs + 1
  postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=10,width=8, horizontal=FALSE)
 # png(file = paste(prefix,nfigs,".png"), bg = "white")
  plot(scoresX, scoresY,
       xlab="Predictive Component",
       ylab=labelY, cex.lab=1.2,
       main = paste("Best model: R2 = ",
                    substr(as.character(ir.validation$full.results[[best.model]]$R2[best.nComp]),
                           start=1,stop=6), 
                    " Q2 = ",
                    substr(as.character(ir.validation$full.results[[best.model]]$Q2[best.nComp]),
                           start=1,stop=6)),
       #squares = arabica, circles = robusta
       #black squares = colombia
       #grey squares = brazil
       #clear squares = peru
       #x squares = other arabicas
       col=c("grey66", "black", "black", "black",
             "black")[as.numeric(metadata$groups)[unlist(ir.validation$training.set[best.model])]],
       pch=c(15,1,15,7,22)[as.numeric(metadata$groups)[unlist(ir.validation$training.set[best.model])]])
  
  abline(h=0);abline(v=0) 
  ordiellipse(cbind(scoresX, scoresY), 
              metadata$species1[unlist(ir.validation$training.set[best.model])],conf=0.95, col="gray66", lwd=2)
  
  dev.off()
  
  #saveDevice("irARopls")
  #   plot(ir.validation$full.results[[selected.model]]$xloadings[,1],
  #        ir.validation$full.results[[selected.model]]$xloadings[,2],
  #        type="p", main="")
  
  # plot(ir.validation$full.results[[best.model]]$xloadings[,1]/max(ir.validation$full.results[[best.model]]$xloadings[,1]),     
  #      ir.validation$full.results[[best.model]]$ortholoadings[,1]/(max(ir.validation$full.results[[best.model]]$ortholoadings[,1])),
  #        type="p", main="best")
  
  
  # 
  # plot(ir.validation$full.results[[selected.model]]$xloadings[,1],
  #      ir.validation$full.results[[selected.model]]$xloadings[,2],
  #      type="p", main="")
  # 
  # plot(ir.validation$full.results[[selected.model]]$ortholoadings[,2],
  #      ir.validation$full.results[[selected.model]]$xloadings[,1],
  #      type="p", main="")
  
  
  
  # loadings.quartiles <- matrix(c(1:(dim(predictors)[[2]]*5)),ncol=5);
  # 
  # current.loadings <- c(1:15);
  # for (i in c(1:dim(predictors)[[2]])){
  #   for (j in c(1:length(ir.validation$full.results))){
  #     current.loadings[[j]] <- ir.validation$full.results[[j]]$xloadings[[i]];
  #   }
  #   current.loadings <- sort(current.loadings);
  #   loadings.quartiles[i,1] <- current.loadings[[1]];
  #   loadings.quartiles[i,2] <- current.loadings[[25]];
  #   loadings.quartiles[i,3] <- current.loadings[[50]];
  #   loadings.quartiles[i,4] <- current.loadings[75];
  #   loadings.quartiles[i,5] <- current.loadings[100];
  # }
  # matplot(matrix(c(1:dim(predictors)[[2]]),ncol=1), loadings.quartiles, type="l");
  # 
  
  nirL<-t(do.call("cbind", lapply(ir.validation$full.results, function(x) x$xloadings[,1])))
  nirL2<-t(do.call("cbind", lapply(ir.validation$full.results, function(x) x$xloadings[,2])))
  nirL3<-t(do.call("cbind", lapply(ir.validation$full.results, function(x) x$xloadings[,3])))
  #t(nirdev2)
  
  #Checking loading stability: plot spectra, loadings, and relative error on loadings for first 2 loadings
  #Doing it sorted for easier comparison of high loading value vs. high variability
  #Filtering top 95% loadings
  ordering <- upper(apply(abs(nirL),2,median))
 # plot(apply(predictors,2,median)[ordering], type="l", main = "Spectrum ordered by loadings")

  nfigs = nfigs + 1
  postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=10,width=8, horizontal=FALSE)
 # png(file = paste(prefix,nfigs,".png"), bg = "white")
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), widths=c(2,1), heights=c(1,2))
  plot(apply(predictors,2,median)[ordering], type="l", main = "", 
       frame.plot = FALSE, ylab="", axes=FALSE, col=1, xlab="") #Spectrum ordered by loadings
  axis(2, cex=3.5, cex.lab=3.5, cex.axis=1.4)
  #grid(50,NA ,lwd = 2)
  #par(new=TRUE)
  plot(apply(abs(nirL),2,median)[ordering], type = "l", main = "",
       frame.plot = FALSE, axes=FALSE, ylab="", col=color[i], xlab="") #Median Loading 1 (Absolute value)
  
  #axis(3, cex=0.5, col=2)   # Draw the x-axis above the plot area
  axis(4, cex=3.5,cex.lab=3.5, col=color[i], cex.axis=1.4 )
  par(new=TRUE)
  # grid(50,NA ,lwd = 2)
  plot(apply(abs(nirL),2,function(x) sd(x) / mean(x))[ordering], type = "l",  
       main="", frame.plot = FALSE, ylab="", xlab="Order of decreasing loadings", 
       cex=3.5, cex.lab=3.5,cex.axis=1.4) #Loading 1 relative error (absolute values)
  #saveDevice(paste(names[i],"LoagOrder3" ))
  #  grid(50,NA ,lwd = 2)
  par(mfrow=c(1,1))
  dev.off()
  

#plot(apply(abs(nirL),2,function(x) sd(x) / mean(x))[ordering], type = "l",  
 #    main="Loading 1 relative error (absolute values)")

  plot(apply(nirL,2,function(x) sd(x) / mean(x))[ordering], type = "l",  
       main="Loading 1 relative error")
  plot(apply(nirL,2,function(x) sd(x) / mean(x))[ordering], type = "l",  ylim = c(-10,10),
       main="Loading 1 relative error (ZOOM)")
 
  #Another check for loading stability: correlations between first loadings for the different models, histogramize
  hist(apply(nirL, 1, function(x) cor(nirL[1,], x)), breaks = 20, xlab = "Correlation with model 1",
            main = "Loading 1 stability")
  
  ordering <- upper(apply(abs(nirL2),2,median))
#  plot(apply(predictors,2,median)[ordering], type="l", main = "Spectrum ordered by loadings")
  plot(apply(abs(nirL2),2,median)[ordering], type = "l", main = "Median Loading 2 (Absolute value)")
  plot(apply(nirL2,2,function(x) sd(x) / mean(x))[ordering], type = "l",  
       main="Loading 2 relative error")
  plot(apply(nirL2,2,function(x) sd(x) / mean(x))[ordering], type = "l",  ylim = c(-10,10),
       main="Loading 2 relative error (ZOOM)")
  plot(apply(abs(nirL2),2,function(x) sd(x) / mean(x))[ordering], type = "l",  
       main="Loading 2 relative error (absolute values)")
  #Another check for loading stability: correlations between first loadings for the different models, histogramize
  hist(apply(nirL2, 1, function(x) cor(nirL2[1,], x)), breaks = 20, xlab = "Correlation with model 1",
       main = "Loading 2 stability")
  
  ####COLOMBIA VS. OTHER ARABICAS
  predictors <- original.predictors;
  
  nfigs = nfigs + 1
  postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=10,width=8, horizontal=FALSE)
 # png(file = paste(prefix,nfigs,".png"), bg = "white")
  plot(c(), main = c(names[[i]], "Colombia vs. Other countries"), xlim=c(0,0), ylim=c(0,0))
  dev.off()
  
  ##Select Arabica Samples and build the corresponding predictors and responses matrices
  filter <- metadata$species1 == 1 & !metadata$ID %in% c("AC1436", "AC1437", "AC1432")
  predictors <- predictors[filter,]
  #   X<-osc$X[filter,]
  #   X1<-irdata.deriv2.wavelength.optimized1[filter,]
  #   X2<-nmr.spectra[filter,]
  #   
  responses <- metadata$groups[filter]
  levels(responses)[!levels(responses) == "Colombia"] = -1
  levels(responses)[levels(responses) == "Colombia"] = 1
  #filtered.metadata <- metadata[filter,]
  
  #   ## PCA to detect outliers
  #   pca <- prcomp(predictors,scale=scaling[[i]],center=centering[[i]])
  #   plot(pca$x[,1],pca$x[,2], 
  #        col=c("grey66", "black", "black", "black", "black")[as.numeric(metadata$groups[filter])], 
  #        pch=c(15,1,15,7,22)[as.numeric(metadata$groups[filter])], #pch=c(16,15)[metadata$groups],
  #        xlab= paste("1st PC (", round((pca$sdev^2 / sum(pca$sdev^2))[1] * 100,1), "%)"), 
  #        ylab= paste("2nd PC (", round((pca$sdev^2 / sum(pca$sdev^2))[2] * 100,1), "%)"),
  #        cex=1.2, cex.axis=1.0,cex.lab=1.3)  ##png and pdf
  #   text(pca$x[,1:2], label=metadata$ID[filter],cex=0.7, pos=1, col="red")
  #   abline(h=0);abline(v=0) 
  #   ordiellipse(cbind(pca$x[,1], pca$x[,2]),
  #               responses,
  #               conf=0.95, col="gray66", lwd=2)
  
  #saveDevice("pca")
  
  ## remove outliers ### 
  #sub outliers detected visually in the PCA scores plot
  #outliers<-which(ir.metadata$catalogID=="AC1427"|ir.metadata$catalogID=="AC1199"|ir.metadata$catalogID=="AC1380")
  #irdata.deriv2.wavelength.optimized1<-irdata.deriv2.wavelength.optimized1[-outliers,]
  #ir.metadata<-ir.metadata[-outliers,]
 
ir.validation <-lims.validate(model = model[[i]], predictors = predictors,
                              responses = as.matrix(as.numeric(as.character(responses))),
                              nModels = nModels, nComp=nComp,xscale = scaling1[[i]], ycenter = FALSE,
                              yscale = FALSE, xcenter = centering1[[i]], method=validating, valsize = valsize)
  
  ##Select optimum number of components
  
  #   y1<-ir.validation$training.set[[1]]
  #   y2<-as.matrix(as.numeric(as.character(responses)))[y1,]
  #   Xtest<-X[-y1,]
  #   X1test<-X1[-y1,]
  #   X2test<-X2[-y1,]
  #   XTest = list(Xtest,X1test, X2test);
  #   
  #   X<-X[y1,]
  #   X1<-X1[y1,]
  #   X2<-X2[y1,]
  #   Xtrain = list(X,X1,X2);#list(X,X2,X3);
  #   m=c(dim(X)[2],dim(X1)[2], dim(X2)[2])
  #   multi<-multiblock(Xtrain,y2,m,3)
  #   ypred=predict1(Xtrain,multi);
  #   print(cbind(y2,ypred))
  #   R2 <- cov(y2,ypred)^2/(var(y2)*var(ypred))
  #   print(cbind("R2",R2))
  #   
  #   x <- Xtest
  #   Ytest<-as.matrix(as.numeric(as.character(metadata$species1)))[-y1,]
  #   ypred = predict1(x,multi) #Abdi verion 
  #   Q2 <- cov(Ytest,ypred)^2/(var(Ytest)*var(ypred))
  #   print(cbind("Q2: ",Q2))
  #   print("Confusion matrix test XWQ")
  #   
  
#   W<-c(); Pred<-c()
#   for (j in c(1:length(ir.validation$full.results))){ 
#     pred<-(ir.validation$full.results[[j]]$predicted.values[,,best.nComp])
#     #  pred<-as.factor(as.numeric(ir.validation$full.results[[i]]$predicted.values[,,nComp]>0))
#     #levels(pred)<-c("-1","1")
#     #pred<-as.numeric(as.character(pred))
#     w<-as.numeric(as.character(responses[-ir.validation$training.set[[j]]]))
#     #print(w)
#     #W<-c(W,w)
#     #Pred<-c(Pred,pred)
#     roc1 <- roc(w,
#               pred, percent=TRUE,
#               # arguments for auc
#               #  partial.auc=c(100, 90), partial.auc.correct=TRUE,
#               # partial.auc.focus="sens",
#               # arguments for ci
#               ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
#               # arguments for plot
#               plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
#               print.auc=TRUE, show.thres=TRUE, legacy.axes = TRUE)
#   }
#   CO.ROC <- c(CO.ROC,roc1)
#   # AUC<-c(AUC,ROC$auc)  
#   rm(Pred) 
#   rm(conf)
  
  #Plot Q2 vs. Nº of components to choose optimum
  plot.with.error.bars(1:nComp,ir.validation$Q2, pch=19, xlab = "Number of components", 
                       ylab = "Q2")
  
  #Choose best number of components from the previous plot
  best.nComp <- 0
  while (best.nComp < 1){
    best.nComp <- readline("Choose number of components to proceed")
  }
  best.nComp <- as.numeric(best.nComp)
  
  ##Save data for ROC analysis
  j <- 0
  for (results in ir.validation$full.results){
    j = j + 1
 #   CO.ROC[[i]]$predicted <- rbind(CO.ROC[[i]]$predicted, sign(results$predicted.values[,,best.nComp]))
    CO.ROC[[i]]$predicted <- rbind(CO.ROC[[i]]$predicted, (results$predicted.values[,,best.nComp]))
 
    CO.ROC[[i]]$observed <- rbind(CO.ROC[[i]]$observed, 
                                  as.numeric(as.character(responses[-ir.validation$training.set[[j]]])))
  }
  
  ##Save Q2 for latter comparative plotting
  CO.Q2 <- c(CO.Q2, list(ir.validation$Q2))
  
  ##Compute Q2.density and select example model
  Q2.density <- density(ir.validation$Q2[,best.nComp])
  #  plot(Q2.density, xlab = "Q2 for optimum number of components") #Just for checking, eventually all Q2 densities will be in a single plot
  selected.model <- which.min(abs(ir.validation$Q2[,best.nComp] - Q2.density$x[which.max(Q2.density$y)])) #we draw
  best.model <- which.max(ir.validation$Q2[,best.nComp])
  
  ##Save Q2.density for latter comparative plotting
  CO.Q2.densities <- c(CO.Q2.densities, list(Q2.density))
  
  ##scores plot of model with Q2 in the maximum of the distribution #(for 2 components?)
  #Q2.density <- density(ir.validation$Q2[,2])
  #selected.model <- which.min(abs(ir.validation$Q2[,2] - Q2.density$x[which.max(Q2.density$y)])) #we draw
  par(mar=c(5.1,4.5,4,5.1))
  if (identical(model[[i]], lims.pls.model)){
    scoresX <- ir.validation$full.results[[selected.model]]$xscores[,1]
    scoresY <- ir.validation$full.results[[selected.model]]$xscores[,2]
    labelY <- "Predictive Component"
  }
  else{
    if (identical(model[[i]], lims.opls.model)){
      scoresX <- ir.validation$full.results[[selected.model]]$xscores[,best.nComp]
      scoresY <- ir.validation$full.results[[selected.model]]$orthoscores[,1]
      labelY <- "Orthogonal Component"
    }
  }
  
  nfigs = nfigs + 1
  postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=10,width=8, horizontal=FALSE)
 # png(file = paste(prefix,nfigs,".png"), bg = "white")
  plot(scoresX, scoresY,
       xlab="Predictive Component",
       ylab=labelY, cex.lab=1.2,
       main = paste("Typical model: R2 = ",
                    substr(as.character(ir.validation$full.results[[selected.model]]$R2[best.nComp]),
                           start=1,stop=6), 
                    " Q2 = ",
                    substr(as.character(ir.validation$full.results[[selected.model]]$Q2[best.nComp]),
                           start=1,stop=6)),
       #squares = arabica, circles = robusta
       #black squares = colombia
       #grey squares = brazil
       #clear squares = peru
       #x squares = other arabicas
       col=c("grey66", "black", "black", "black",
             "black")[as.numeric(metadata$groups)[filter][unlist(ir.validation$training.set[selected.model])]],
       pch=c(15,1,15,7,22)[as.numeric(metadata$groups)[filter][unlist(ir.validation$training.set[selected.model])]])
  
  
  abline(h=0);abline(v=0) 
  ordiellipse(cbind(scoresX, scoresY), 
              responses[unlist(ir.validation$training.set[selected.model])],conf=0.95, col="gray66", lwd=2)
  
  dev.off()
  
  #   plot(ir.validation$full.results[[selected.model]]$xloadings[,1],
  #        ir.validation$full.results[[selected.model]]$xloadings[,2],
  #        type="p", main="CO")
  #Scores plot of best model
  par(mar=c(5.1,4.5,4,5.1))
#   if (identical(model[[i]], lims.pls.model)){
#     scoresX <- ir.validation$full.results[[best.model]]$xscores[,1]
#     scoresY <- ir.validation$full.results[[selected.model]]$xscores[,2]
#   }
#   else{
#     if (identical(model[[i]], lims.opls.model)){
#       scoresX <- ir.validation$full.results[[selected.model]]$xscores[,best.nComp]
#       scoresY <- ir.validation$full.results[[selected.model]]$orthoscores[,1]
#     }
#   }
  
  if (identical(model[[i]], lims.pls.model)){
    scoresX <- ir.validation$full.results[[best.model]]$xscores[,1]
    scoresY <- ir.validation$full.results[[best.model]]$xscores[,2]
  }
  else{
    if (identical(model[[i]], lims.opls.model)){
      scoresX <- ir.validation$full.results[[best.model]]$xscores[,best.nComp]
      scoresY <- ir.validation$full.results[[best.model]]$orthoscores[,1]
    }
  }
  
  nfigs = nfigs + 1
  postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=10,width=8, horizontal=FALSE)
 # png(file = paste(prefix,nfigs,".png"), bg = "white")
  plot(scoresX, scoresY,
       xlab="Predictive Component",
       ylab=labelY, cex.lab=1.2,
       main = paste("Best model: R2 = ",
                    substr(as.character(ir.validation$full.results[[best.model]]$R2[best.nComp]),
                           start=1,stop=6), 
                    " Q2 = ",
                    substr(as.character(ir.validation$full.results[[best.model]]$Q2[best.nComp]),
                           start=1,stop=6)),
       #squares = arabica, circles = robusta
       #black squares = colombia
       #grey squares = brazil
       #clear squares = peru
       #x squares = other arabicas
       col=c("grey66", "black", "black", "black",
             "black")[as.numeric(metadata$groups)[filter][unlist(ir.validation$training.set[best.model])]],
       pch=c(15,1,15,7,22)[as.numeric(metadata$groups)[filter][unlist(ir.validation$training.set[best.model])]])
  
  
  abline(h=0);abline(v=0) 
  ordiellipse(cbind(scoresX, scoresY), 
              responses[unlist(ir.validation$training.set[best.model])],conf=0.95, col="gray66", lwd=2)
  
  dev.off()
  #saveDevice("irARopls")
  
  # loadings.quartiles <- matrix(c(1:(dim(predictors)[[2]]*5)),ncol=5);
  
  # current.loadings <- c(1:15);
  # for (i in c(1:dim(predictors)[[2]])){
  #   for (j in c(1:length(ir.validation$full.results))){
  #     current.loadings[[j]] <- ir.validation$full.results[[j]]$xloadings[[i]];
  #   }
  #   current.loadings <- sort(current.loadings);
  #   loadings.quartiles[i,1] <- current.loadings[[1]];
  #   loadings.quartiles[i,2] <- current.loadings[[25]];
  #   loadings.quartiles[i,3] <- current.loadings[[50]];
  #   loadings.quartiles[i,4] <- current.loadings[75];
  #   loadings.quartiles[i,5] <- current.loadings[100];
  # }
  #matplot(matrix(c(1:dim(predictors)[[2]]),ncol=1), loadings.quartiles, type="l");
  
  nirL<-t(do.call("cbind", lapply(ir.validation$full.results, function(x) x$xloadings[,1])))
  nirL2<-t(do.call("cbind", lapply(ir.validation$full.results, function(x) x$xloadings[,2])))
#  nirL3<-t(do.call("cbind", lapply(ir.validation$full.results, function(x) x$xloadings[,3])))
  #t(nirdev2)
  
  #Checking loading stability: plot spectra, loadings, and relative error on loadings for first 2 loadings
  #Doing it sorted for easier comparison of high loading value vs. high variability
  ordering <- upper(apply(abs(nirL),2,median))
  #par(mfrow=c(2,1))
  
  nfigs = nfigs + 1
  postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=10,width=8, horizontal=FALSE)
 # png(file = paste(prefix,nfigs,".png"), bg = "white")
  
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), widths=c(2,1), heights=c(1,2))
  #par(mar = c(4,4,0,0), oma = c(0,0, 0,0))
  plot(apply(predictors,2,median)[ordering], type="l", main = "", 
       frame.plot = FALSE, ylab="", axes=FALSE, col=1, xlab="") #Spectrum ordered by loadings
  axis(2, cex=3.5, cex=3.5, cex.axis=1.4)
  #grid(50,NA ,lwd = 2)
  #par(new=TRUE)
 # par(mar = c(4,4,0,0), oma = c(0,0, 0,0))
  plot(apply(abs(nirL),2,median)[ordering], type = "l", main = "",
       frame.plot = FALSE, axes=FALSE, ylab="", col=color[i], xlab="", cex=3.5, cex.lab=3.5) #Median Loading 1 (Absolute value)
  
  #axis(3, cex=0.5, col=2)   # Draw the x-axis above the plot area
  axis(4, cex.lab=3.5,cex=3.5,col=color[i],cex.axis=1.4)
  par(new=TRUE)
 # par(mar = c(4,4,0,0), oma = c(0,0, 0,0))
 # grid(50,NA ,lwd = 2)
  plot(apply(abs(nirL),2,function(x) sd(x) / mean(x))[ordering], type = "l",  
       main="", frame.plot = FALSE, ylab="", xlab="Order of decreasing loadings",cex.axis=1.4, cex.lab=3.5, cex=3.5) #Loading 1 relative error (absolute values)
  
  #  grid(50,NA ,lwd = 2)
  par(mfrow=c(1,1))
  
  dev.off()
    
  plot(apply(nirL,2,function(x) sd(x) / mean(x))[ordering], type = "l",  
       main="Loading 1 relative error", ylab="")
  
  par(mfrow=c(1,1))
  plot(apply(nirL,2,function(x) sd(x) / mean(x))[ordering], type = "l",  ylim = c(-10,10),
       main="Loading 1 relative error (ZOOM)")
  #Another check for loading stability: correlations between first loadings for the different models, histogramize
  hist(apply(nirL, 1, function(x) cor(nirL[1,], x)), breaks = 20, xlab = "Correlation with model 1",
       main = "Loading 1 stability")
  
  ordering <- upper(apply(abs(nirL2),2,median))
 # plot(apply(predictors,2,median)[ordering], type="l", main = "Spectrum ordered by loadings")
  plot(apply(abs(nirL2),2,median)[ordering], type = "l", main = "Median Loading 2 (Absolute value)")
  plot(apply(nirL2,2,function(x) sd(x) / mean(x))[ordering], type = "l",  
       main="Loading 2 relative error")
  plot(apply(nirL2,2,function(x) sd(x) / mean(x))[ordering], type = "l",  ylim = c(-10,10),
       main="Loading 2 relative error (ZOOM)")
  plot(apply(abs(nirL2),2,function(x) sd(x) / mean(x))[ordering], type = "l",  
       main="Loading 2 relative error (absolute values)")
  #Another check for loading stability: correlations between first loadings for the different models, histogramize
  hist(apply(nirL2, 1, function(x) cor(nirL[1,], x)), breaks = 20, xlab = "Correlation with model 1",
       main = "Loading 2 stability")

  
}

nfigs = nfigs + 1
postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=6,width=8, horizontal=FALSE)
#png(file = paste(prefix,nfigs,".png"), bg = "white")
plot.with.error.bars(c(1:15), AR.Q2[[1]], pch = 15, col = "black", xlab = "Number of components", ylab = "Q2",
                     ylim = c(-1,1),main="Arabica / Robusta", cex.axis=1.5, cex.lab=1.5, cex = 1.5)
par(new = T)
plot.with.error.bars(c(1:15), AR.Q2[[2]], pch = 15, col = "red", xlab = "Number of components", ylab = "Q2",
                     ylim = c(-1,1), cex.axis=1.5, cex.lab=1.5, cex=1.5)
par(new = T)
plot.with.error.bars(c(1:15), AR.Q2[[3]], pch = 15, col = "blue", xlab = "Number of components", ylab = "Q2",
                     ylim = c(-1,1), cex.axis=1.5, cex.lab=1.5, cex=1.5)
dev.off()

ymax <- max(as.numeric(lapply(AR.Q2.densities, function(x) max(x$y))))
xmin <- min(as.numeric(lapply(AR.Q2.densities, function(x) min(x$x))))
xmax <- max(as.numeric(lapply(AR.Q2.densities, function(x) max(x$x))))
nfigs = nfigs + 1
postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=6,width=8, horizontal=FALSE)
#png(file = paste(prefix,nfigs,".png"), bg = "white")
plot(AR.Q2.densities[[1]], xlab = "Q2", col="black", ylim = c(0,ymax), xlim = c(xmin,xmax), main="Arabica / Robusta")
lines(AR.Q2.densities[[2]], col="red")
lines(AR.Q2.densities[[3]], col="blue")
dev.off()

nfigs = nfigs + 1
postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=6,width=8, horizontal=FALSE)
#png(file = paste(prefix,nfigs,".png"), bg = "white")
plot.with.error.bars(c(1:15), CO.Q2[[1]], pch = 15, col = "black", xlab = "Number of components", ylab = "Q2",
                     ylim = c(-1,1), main = "Colombia / Others", cex.axis=1.5, cex.lab=1.5, cex=1.5)
par(new = T)
plot.with.error.bars(c(1:15), CO.Q2[[2]], pch = 15, col = "red", xlab = "Number of components", ylab = "Q2",
                     ylim = c(-1,1), cex.axis=1.5, cex.lab=1.5, cex=1.5)
par(new = T)
plot.with.error.bars(c(1:15), CO.Q2[[3]], pch = 15, col = "blue", xlab = "Number of components", ylab = "Q2",
                     ylim = c(-1,1), cex.axis=1.5, cex.lab=1.5, cex=1.5)
dev.off()

ymax <- max(as.numeric(lapply(CO.Q2.densities, function(x) max(x$y))))
xmin <- min(as.numeric(lapply(CO.Q2.densities, function(x) min(x$x))))
xmax <- max(as.numeric(lapply(CO.Q2.densities, function(x) max(x$x))))
nfigs = nfigs + 1
postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=6,width=8)
#png(file = paste(prefix,nfigs,".png"), bg = "white")
plot(CO.Q2.densities[[1]], xlab = "Q2", col="black", ylim = c(0,ymax), xlim = c(xmin,xmax), main="Colombia / Others")
lines(CO.Q2.densities[[2]], col="red")
lines(CO.Q2.densities[[3]], col="blue")
dev.off()

nfigs = nfigs + 1
postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=6,width=8, horizontal=FALSE)
#png(file = paste(prefix,nfigs,".png"), bg = "white")
   roc(as.numeric(AR.ROC[[1]]$observed), as.numeric(AR.ROC[[1]]$predicted),
       partial.auc=c(100, 90), partial.auc.correct=TRUE,percent=TRUE,
       partial.auc.focus="sens",
       # arguments for ci
       ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
       # arguments for plot
       plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
       print.auc=FALSE, show.thres=TRUE, legacy.axes = TRUE, col = "black", main = "Arabica / Robusta",
       cex.axis=1.5, cex.lab=1.5, cex=1.5)
  par(new = TRUE)
  roc(as.numeric(AR.ROC[[2]]$observed), as.numeric(AR.ROC[[2]]$predicted),
      partial.auc=c(100, 90), partial.auc.correct=TRUE,percent=TRUE,
      partial.auc.focus="sens",
      # arguments for ci
      ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
      # arguments for plot
      plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
      print.auc=FALSE, show.thres=TRUE, legacy.axes = TRUE, col = "red",
      cex.axis=1.5, cex.lab=1.5, cex=1.5)
  par(new = TRUE)
  roc(as.numeric(AR.ROC[[3]]$observed), as.numeric(AR.ROC[[3]]$predicted),
      partial.auc=c(100, 90), partial.auc.correct=TRUE,percent=TRUE,
      partial.auc.focus="sens",
      # arguments for ci
      ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
      # arguments for plot
      plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
      print.auc=FALSE, show.thres=TRUE, legacy.axes = TRUE, col = "blue",
      cex.axis=1.5, cex.lab=1.5, cex=1.5)

dev.off()

# roc(as.numeric(AR.ROC[[3]]$observed), as.numeric(AR.ROC[[3]]$predicted),
#     #partial.auc = c(100,90), partial.auc.correct = TRUE, partial.auc.focus = "sens", allow.invalid.partial.auc.correct = FALSE,
#     #   ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
#     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,percent=TRUE,
#     print.auc=FALSE, show.thres=TRUE, legacy.axes = TRUE, col = "blue")

  
#   roc(as.numeric(CO.ROC[[1]]$observed), as.numeric(CO.ROC[[1]]$predicted),
#       #partial.auc = c(100,90), partial.auc.correct = TRUE, partial.auc.focus = "sens", allow.invalid.partial.auc.correct = FALSE,
#     #  ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
#       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,percent=TRUE,
#       print.auc=FALSE, show.thres=TRUE, legacy.axes = TRUE, col = "black")

nfigs = nfigs + 1
postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=6,width=8, horizontal=FALSE)
#png(file = paste(prefix,nfigs,".png"), bg = "white")
roc(as.numeric(CO.ROC[[1]]$observed), as.numeric(CO.ROC[[1]]$predicted),
    partial.auc=c(100, 90), partial.auc.correct=TRUE,percent=TRUE,
    partial.auc.focus="sens",
    # arguments for ci
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    # arguments for plot
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=FALSE, show.thres=TRUE, legacy.axes = TRUE,col = "black", main = "Colombia / Others",
    cex.axis=1.5, cex.lab=1.5, cex=1.5)

  par(new = TRUE)
  roc(as.numeric(CO.ROC[[2]]$observed), as.numeric(CO.ROC[[2]]$predicted),
      partial.auc=c(100, 90), partial.auc.correct=TRUE,percent=TRUE,
      partial.auc.focus="sens",
      # arguments for ci
      ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
      # arguments for plot
      plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
      print.auc=FALSE, show.thres=TRUE, legacy.axes = TRUE,col = "red",
      cex.axis=1.5, cex.lab=1.5, cex=1.5)


  par(new = TRUE)
roc(as.numeric(CO.ROC[[3]]$observed), as.numeric(CO.ROC[[3]]$predicted),
    partial.auc=c(100, 90), partial.auc.correct=TRUE,percent=TRUE,
    partial.auc.focus="sens",
    # arguments for ci
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    # arguments for plot
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=FALSE, show.thres=TRUE, legacy.axes = TRUE,col = "blue",
    cex.axis=1.5, cex.lab=1.5, cex=1.5)

dev.off()

AUCs <- matrix(nrow = nrow(AR.ROC[[1]]$observed), ncol = 3)
k <- 0
for (tech in AR.ROC){
  k = k + 1
  for (i in 1:nrow(tech$observed)){
    AUCs[i,k] <- as.numeric(roc(tech$observed[i,], tech$predicted[i,], auc = TRUE, plot = FALSE)$auc)
  }
}

nfigs = nfigs + 1
postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=6,width=8, horizontal=FALSE)
#png(file = paste(prefix,nfigs,".png"), bg = "white")
boxplot(AUCs, names= c("mIR", "NIR", "NMR"), main = "Arabica / Robusta",
        cex.axis=1.5, cex.lab=1.5, cex=1.5)
dev.off()

# print(t.test(AUCs[,1], AUCs[,2]))
# print(t.test(AUCs[,1], AUCs[,3]))
# print(t.test(AUCs[,2], AUCs[,3]))

k <- 0
for (tech in CO.ROC){
  k = k + 1
  for (i in 1:nrow(tech$observed)){
    AUCs[i,k] <- as.numeric(roc(tech$observed[i,], tech$predicted[i,], auc = TRUE, plot = FALSE)$auc)
  }
}

nfigs = nfigs + 1
postscript(paste(prefix,nfigs,".eps", sep=""), paper="special",height=6,width=8, horizontal=FALSE)
#png(file = paste(prefix,nfigs,".png"), bg = "white")
boxplot(AUCs, names= c("mIR", "NIR", "NMR"), main = "Colombia / Others",
        cex.axis=1.5, cex.lab=1.5, cex=1.5)
dev.off()

print(t.test(AUCs[,1],AUCs[,2]))
print(t.test(AUCs[,1],AUCs[,3]))
print(t.test(AUCs[,2],AUCs[,3]))

# print(t.test(AUCs[,1], AUCs[,2]))
# print(t.test(AUCs[,1], AUCs[,3]))
# print(t.test(AUCs[,2], AUCs[,3]))



# plot(ROC[[1]])
# par(new = T)
# plot(ROC[[2]], col=2)
# par(new = T)
# plot(ROC[[3]], col=4)

# data<- 
#   c(list(osc$X),list(nmr.spectra),  list(irdata.deriv2.wavelength.optimized1))
# 
# y<-ir.validation$training.set[[1]]
# y1<-as.numeric(responses[y])
# Xtrain<-osc$X[y,]
# Xtrain1<-nmr.spectra[y,]
# Xtrain2<-irdata.deriv2.wavelength.optimized1[y,]
# data<- c(list(Xtrain),list(Xtrain1),list(Xtrain2))


#multi<-multiblock(data,y1,m,3)

