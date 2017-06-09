library(parallel)

.oneGenePerturbations=function(gene, network, refEXP, targetEXP, delta=0, measure="perturbation"){
  ##### an empty list to return in case of error
  l <- rep(0,length(rownames(targetEXP)))
  names(l) <- rownames(targetEXP)
  
  ##### return the list of 0s is the gene is not found in the network
  if(!(gene %in% targets(network))) return(l)
  
  ##### return the list of 0s is the gene is not found in the expression dataset
  if(!(gene %in% colnames(refEXP)) | !(gene %in% colnames(targetEXP)) ){
    print(paste("Couldn't compute perturbation for: ", gene, ". The gene was not found in the reference or the target expression data."))
    return(l)
  }
  
  print(paste("computing perturbations for gene ",gene))
  
  ##### compute the RMSE and get the best GRN in case multiple networks exist for the same gene (the case of some regulatory network inference methods like HLICORN)
  # listedgrn <- data.frame(t(network@GRN[which(network@GRN$Target==gene),1:3]),stringsAsFactors=FALSE)
  # results=lapply(listedgrn,.fitGRNcv,exp=refEXP)
  # # results=mclapply(listedgrn,.fitGRNcv,exp=refEXP, mc.cores = 2)
  # bestgrn <- listedgrn[which.min(sapply(results,f <- function(x){return(unlist(x)["RMSE"])}))]
  
  ##### or take it directly from networks infered with the hLicorn method using the CoRegNet R package (help saving computation time)
  listedgrn <- data.frame(network@GRN[which(network@GRN$Target==gene),],stringsAsFactors=FALSE)
  bestgrn <- t(listedgrn[which.min(listedgrn$RMSE),1:3])
  
  ##### estimate the gene perturbation score 
  perturbation <-  .estimatePerturbation(grn = unlist(bestgrn), refexp = refEXP, targetexp=targetEXP)
  
  ##### if an error occured (mainly the regulators of the network were not found in the expression data)
  if(is.null(perturbation)) return(l)
  
  ##### set to 0 all perturbations for which the confidence score is < then delta
  if (perturbation$Confidence < delta) return(l)
  
  ##### it is possible also to cancel just the perturbation scores for only the sample that are < than delta  
  # refs_to_0 <- which(perturbation$Confidence.perSample < delta)
  # perturbation$Confidence.perSample[refs_to_0] <- 0

  ##### we can provide values for perturbation (by default) and also for the percentage of perturbation, 
  ##### the expected expression for the gene according to its regulatory network, per sample ot the average 
  ##### confidence of the estimation, the RMSE and R2 of the model
  return(switch(measure,
                perturbation = perturbation$Target.perturbation,
                perturbationPercentage = perturbation$Target.perturbationPercentage,
                expectedExpression = perturbation$Target.pred,
                confidencePerSample = perturbation$Confidence.perSample,
                confidence = perturbation$Confidence,
                rmse = perturbation$RMSE,
                r2 = perturbation$R2))
}


.estimatePerturbation=function(grn,refexp,targetexp)
{
  
  ##### get all regulators (predictors) and gene (response) data and organise in a new dataframe
  #deal with several, one or no coregulators
  act = unique(unlist(strsplit(grn[2]," ")))
  act=act[which(act!="EMPTY" & !is.na(act))]
  rep = unique(unlist(strsplit(grn[3]," ")))
  rep=rep[which(rep!="EMPTY" & !is.na(rep))]
  
  if(!(all(c(act,rep) %in% colnames(refexp)) && all(c(act,rep) %in% colnames(targetexp)))){
    print(paste("Couldn't compute perturbation for: ", grn[1], ". One or multiple of its regulators were not found in the reference or the target expression data."))
    return()
  }
  
  X = refexp[,c(act,rep)]
  targetX=targetexp[,c(act,rep)]
  
  if(length(act)>1){
    coact = apply(refexp[,act],1,prod)
    X=data.frame(X,"coact"=as.numeric(coact))
    targetcoact = apply(targetexp[,act],1,prod)
    targetX = data.frame(targetX,"coact"=as.numeric(targetcoact))
  }
  if(length(rep)>1){
    corep = apply(refexp[,rep],1,prod)
    X=data.frame(X,"corep"=as.numeric(corep))
    targetcorep = apply(targetexp[,rep],1,prod)
    targetX = data.frame(targetX,"corep"=as.numeric(targetcorep))
  }
  
  # da contains all the necessary data with the first column y as the gene, response variable
  Y=refexp[,grn[1]]
  da = data.frame("y"= Y,X)
  # get fitted lm, coefficients and fitted values
  l=lm("y~.",da)
  fitted=l$fitted.values
  coefs=coef(l)
  coefs= coefs[2:length(coefs)]
  residues = Y-fitted
  numscores.RMSE = sqrt(sum(( (residues) ^2) )/length(fitted))
  # print(numscores.RMSE)
  numscores.R2 =  summary(l)$r.squared
  numscores.Confidence = 1-numscores.RMSE
  numscores.Confidence.perSample = 1-residues
  
  targetX <- data.frame(targetX)
  # rename the column to X when only the size of c(act,rep) is 1 since R transforms the dataframe into a vector and thus the name of the predictor (the TF) becomes the name of the learning variable ("X") (see above)
  # thus we rename the predictor as "X"
  if(!(length(c(act,rep))>1)){
    colnames(targetX) <- "X"
  }
  
  y=targetexp[,grn[1]]
  numscores.Target.pred = predict(l,targetX)
  numscores.Target.perturbation = numscores.Target.pred-y
  numscores.Target.perturbationPercentage = abs(numscores.Target.perturbation/numscores.Target.pred)
  
  # #test Remy's method on checking conserved gene networks
  # l=lm("y~.",data.frame(y, targetX))
  # ft=l$fitted.values
  # rs = y-fitted
  # RMSE = sqrt(sum(( (rs) ^2) )/length( ft))
  # print(summary(l)$r.squared)
  # if(RMSE>0.5) return()
  # # if(summary(l)$r.squared <= 0.5) return()
  
  if(length(act)==0){
    numscores.Coef.Acts=NA
    numscores.Coef.coActs=NA
  }else if(length(act)==1){
    numscores.Coef.Acts=coefs[1]
    numscores.Coef.coActs=NA
  }else{
    numscores.Coef.Acts=paste(coefs[1:length(act)],collapse= " ")
    numscores.Coef.coActs=coefs["coact"]
  }
  
  if(length(rep)==0){
    numscores.Coef.Reps=NA
    numscores.Coef.coReps=NA
  }else if(length(rep)==1){
    numscores.Coef.Reps=coefs[length(act)+1]
    numscores.Coef.coReps=NA
  }else{
    numscores.Coef.Reps=paste(coefs[(length(act)+1):length(rep)],collapse= " ")
    numscores.Coef.coReps=coefs["corep"]
  }
  numscores <- list(numscores.Coef.Acts, numscores.Coef.Reps, numscores.Coef.coActs, numscores.Coef.coReps, numscores.R2, numscores.RMSE, numscores.Confidence,
                    numscores.Confidence.perSample, numscores.Target.pred, numscores.Target.perturbation, numscores.Target.perturbationPercentage)
  names(numscores) <- c("Coef.Acts", "Coef.Reps", "Coef.coActs", "Coef.coReps", "R2", "RMSE", "Confidence",
                        "Confidence.perSample", "Target.pred", "Target.perturbation", "Target.perturbationPercentage")
  return(numscores)
}


.fitGRNcv=function(grn,exp)
{
  #get all regulators (predictors) and gene (response) data and organise in a new dataframe
  #deal with several, one or no coregulators
  act = unique(unlist(strsplit(grn[2]," ")))
  act=act[which(act!="EMPTY" & !is.na(act))]
  rep = unique(unlist(strsplit(grn[3]," ")))
  rep=rep[which(rep!="EMPTY" & !is.na(rep))]
  X = exp[,c(act,rep)]
  
  if(length(act)>1){
    coact = apply(exp[,act],1,prod)
    X=data.frame(X,"coact"=as.numeric(coact))
  }
  if(length(rep)>1){
    corep = apply(exp[,rep],1,prod)
    X=data.frame(X,"corep"=as.numeric(corep))
  }
  
  # da contains all the necessary data with the first column y as the gene, response variable
  Y=exp[,grn[1]]
  da = data.frame("y"= Y,X)
  
  # get 5-cross validated coefficients and fitted values
  Icol=sample(1:nrow(da))
  Ks=as.integer(cut(1:nrow(da),5))
  da=da[Icol,]
  cvlm=lapply(1:5,function(k){
    daL = da[which(Ks != k),]
    daP = da[which(Ks == k),]
    l=lm("y~.",daL)
    return(list("fitted"=predict(l,daP),"coefs"=coef(l)))
  })
  
  fitted=unlist(lapply(cvlm,function(x){return(x$fitted)}))[order(Icol)]
  coefs=apply(sapply(cvlm,function(x){return(x$coefs)}),1,mean)
  coefs= coefs[2:length(coefs)]
  
  R2=as.numeric(cor(Y,fitted)^2 )
  residues = Y-fitted
  RMSE = sqrt(sum(( (residues) ^2) )/length(fitted))
  
  
  
  numscores =c(rep.int(0,4),R2,RMSE)
  names(numscores)=c("Coef.Acts","Coef.Reps","Coef.coActs","Coef.coReps","R2","RMSE")
  
  if(length(act)==0){
    numscores["Coef.Acts"]=NA
    numscores["Coef.coActs"]=NA
  }else if(length(act)==1){
    numscores["Coef.Acts"]=coefs[1]
    numscores["Coef.coActs"]=NA
  }else{
    numscores["Coef.Acts"]=paste(coefs[1:length(act)],collapse= " ")
    numscores["Coef.coActs"]=coefs["coact"]
  }
  
  if(length(rep)==0){
    numscores["Coef.Reps"]=NA
    numscores["Coef.coReps"]=NA
  }else if(length(rep)==1){
    numscores["Coef.Reps"]=coefs[length(act)+1]
    numscores["Coef.coReps"]=NA
  }else{
    numscores["Coef.Reps"]=paste(coefs[(length(act)+1):length(rep)],collapse= " ")
    numscores["Coef.coReps"]=coefs["corep"]
  }
  return(numscores)
}


#function to compute regulator activity
.regulatorActivity=function (object, expData, minTarg = 5, withEvidences = FALSE, is.scaled = FALSE){
  adjlist = object@adjacencyList
  allgenes = unique(unlist(adjlist$bytf))  
  
  if(length(intersect(allgenes,rownames(expData))) < 0.1*length(allgenes)){
    warning(paste("More than 90% of the network genes are not in the expression data." ,
               "The influence of the regulators cannot be computed with these settings.\nNote" ,
               ": The expression data must be given with genes in line."))
  }
  
  if( sum(!allgenes %in% rownames(expData))>0){
    adjlist=.subsigrns(adjlist,rownames(expData))
  }    
  sampgenes = dimnames(expData)
  if(! is.scaled){
    expData = t(scale(t(expData),scale=FALSE))
    dimnames(expData)=sampgenes
  }
  
  tfs = names(adjlist$bytf)  
  
  tfScore = mclapply(tfs,function(tf){
    if(length(unlist(strsplit(tf," "))) > 1){
      cotf=unlist(strsplit(tf," "))
      acti = names(which(table(unlist(lapply(adjlist$bytf[cotf],function(x){return(x$act)})))== length(cotf)))
      repr =  names(which(table(unlist(lapply(adjlist$bytf[cotf],function(x){return(x$rep)})))== length(cotf)))
    }
    else{
      acti = adjlist$bytf[[tf]]$act
      repr =  adjlist$bytf[[tf]]$rep
    }
    
    if(length(repr)>=minTarg &(length(acti) >= minTarg )) {
      return(
        unlist( lapply(1:ncol(expData),function(tumor){
          g = (t.test( (expData[acti,tumor]) , (expData[repr,tumor]),alternative="two.sided"))
          return(g$statistic)
        }))
      )
    }
    else{ return(NULL)  }
  })
  
  names(tfScore) = unlist( tfs)
  tfscore = do.call(rbind,tfScore)
  colnames(tfscore) = colnames(expData)
  tfscore[tfscore<0] <- 0
  return(tfscore)
}

# used in the computation of TF activity
.subsigrns = function(sigrns,genes)
{
  
  subtg  =intersect(names(sigrns$bygene),genes)
  subtf = intersect(names(sigrns$bytf),genes)
  subsigrns = list()
  
  subsigrns$bytf = lapply(sigrns$bytf[subtf],function(subr){
    subn = list()
    subn$act = intersect(   subr$act  , subtg   )
    subn$rep = intersect(subr$rep,subtg)
    return(subn)
  })
  
  subsigrns$bygene = lapply(sigrns$bygene[subtg],function(subr){
    subn = list()
    subn$act = intersect(   subr$act  , subtf   )
    subn$rep = intersect(subr$rep,subtf)
    return(subn)
  })
  return(subsigrns)
}