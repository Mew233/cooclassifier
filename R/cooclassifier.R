cooclassifier<-function(data,n=50)
{
  data(staudt)
  data(label)
  
  scale.normalization<-function(matrix)
  {
    scale.norm<-t(apply(matrix,1,scale))#log(t(apply(matrix,1,scale))+10)
    return(scale.norm)
  }
  
  p.subgroup<-function(pred,ABC.pred,GCB.pred)
  {
    mean.ABC<-mean(ABC.pred)
    sd.ABC<-sd(ABC.pred)
    mean.GCB<-mean(GCB.pred)
    sd.GCB<-sd(GCB.pred)
    p.ABC<-dnorm(pred,mean.ABC,sd.ABC)/(dnorm(pred,mean.ABC,sd.ABC)+dnorm(pred,mean.GCB,sd.GCB))
    p.GCB<-1-p.ABC  
    p<-rbind(p.ABC,p.GCB)
    rownames(p)<-c("p.ABC","p.GCB")
    return(p)
  }
  
  staudt.scale<-scale.normalization(staudt)
  colnames(staudt.scale)<-colnames(staudt)
  
  ## select all labeled microarray data as train set
  ABC.train<-staudt.scale[,label=="ABC "]
  GCB.train<-staudt.scale[,label=="GCB "]
  
  ## select predictors
  p.value.scale<-vector()
  t.statistic.scale<-vector()
  for(i in 1:dim(ABC.train)[1])
  {
    t.test.result<-t.test(ABC.train[i,],GCB.train[i,])
    p.value.scale[i]<-t.test.result$p.value 
    t.statistic.scale[i]<-t.test.result$statistic 
  }
  names(p.value.scale)<-rownames(ABC.train)
  names(t.statistic.scale)<-rownames(ABC.train)
  FDR<-p.adjust(p.value.scale,"BH")
  t.statistic<-t.statistic.scale[rank(FDR)<=n]
  ABC.pred<-t.statistic %*% ABC.train[rank(FDR)<=n,]
  GCB.pred<-t.statistic %*% GCB.train[rank(FDR)<=n,]
  
  DLBCL<-as.matrix(data)
  rownames(DLBCL)<-toupper(rownames(DLBCL))
  
  ## prepare predictors
  pred<-names(FDR[rank(FDR)<=n])[names(FDR[rank(FDR)<=n]) %in% rownames(DLBCL)]
  DLBCL.pred<-scale.normalization(DLBCL[pred,])
  colnames(DLBCL.pred)<-colnames(DLBCL)
  DLBCL.LPS<-t.statistic[pred] %*% DLBCL.pred
  ABC.pred.test<-t.statistic[pred] %*% ABC.train[pred,]
  GCB.pred.test<-t.statistic[pred] %*% GCB.train[pred,]
  
  test.subgroup<-p.subgroup(DLBCL.LPS,ABC.pred.test,GCB.pred.test)
  sort.DLBCL<-sort(test.subgroup[1,],decreasing=T)
  
  list<-c(rep("unclassified",dim(DLBCL)[2]))
  list[test.subgroup[1,]>0.9]<-"ABC"
  list[test.subgroup[2,]>0.9]<-"GCB"
  names(list)<-colnames(test.subgroup)
  
  return(list)
}
