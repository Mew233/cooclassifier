library("ggplot2")
library("tidyr")
library("grid")
library("dplyr")
library(caret)
library(pheatmap)
library(forcats)
library("survival")
library("survminer")
library(matrixStats)

# Read in data
#mshap <- read.table("mshap_tcga_dlbc.csv",header=T,sep=",",row.names=1)
mshap <- read.table("mshap_lstaudt.csv",header=T,sep=",",row.names=1)
# scale.normalization<-function(matrix)
# {
#   scale.norm<-t(apply(matrix,1,scale))
#   return(scale.norm)
# }
# mshap.scale <- scale.normalization(mshap)
# mshap_3 <- dplyr::bind_cols(mshap,mshap_2)
#mshap <- mshap_3
# ht<-pheatmap::pheatmap(mshap_3, fontsize = 3,clustering_method = "complete")

# ind <- which(is.na(mshap), arr.ind=TRUE)
# mshap[ind] <- rowMeans(mshap,  na.rm = TRUE)[ind[,1]]
df <- t(mshap) %>% 
  as.data.frame() %>% 
  pivot_longer(everything(), names_to = "drug1") %>% 
  group_by(drug1) %>% 
  mutate(meanval = mean(value), stdev = sd(value))
df <- df %>% arrange(across("stdev", desc)) %>% distinct_(.,"drug1") %>% head(30)
#mshap <- mshap %>% filter(rowMeans(select(., where(is.numeric))) > 0.3)
# drugs <- c("Navitoclax","Saracatinib","Danusertib","Duvelisib","PONATINIB","Olaparib",
#            "TAK-733","Tivozanib","Dacomitinib","Tipifarnib","AZD6482","Veliparib",
#            "PF-3845","497839-62-0","REBASTINIB","Gedatolisib","Lapatinib","Pralatrexate")

mshap <- mshap%>% filter(row.names(mshap) %in% df$drug1)
ht<-pheatmap::pheatmap(mshap, fontsize = 3,clustering_method = "ward.D")
hc<-ht$tree_col
lbl <- cutree(hc, 2) # split gene dendrogram in 2 groups
effetive<-which(lbl==1) # grab genes of first group/effective ; second group/non-effective 
noneffetive<-which(lbl==2)
print(length(effetive))
print(length(noneffetive))
## 浅浅画一个pfs
metadata <- read.table("metadata_lstaudt.csv",header=T,sep=",")
lbl_df <- as.data.frame(lbl)
lbl_df$subject = rownames(lbl_df)
metadata <- left_join(metadata,lbl_df) %>% 
  mutate(type = case_when(
    lbl==1 ~ "predicted effective", 
    lbl==2 ~ "predicted noneffective"))

fit <- survfit(Surv(Progression_Free.Survival._PFS_.Time._yrs) ~ type, data = metadata)
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
)

# create train, test set
set.seed(42)
# Now Selecting 50% of data as sample from total 'n' samples of data 
split.data <- function(effective){
  sample <- sample.int(n = length(effective), size = floor(.5*length(effective)), replace = F)
  train <- effective[sample]
  test  <- effective[-sample]
  newList <- list("train" = train, "test" = test)
  return(newList)
}
newList_eff<-split.data(effetive)
newList_noneff<-split.data(noneffetive)
eff_train<-newList_eff$train
eff_test<-newList_eff$test
noneff_train<-newList_noneff$train
noneff_test<-newList_noneff$test
# create train label, and test label
create_label<-function(eff_train,noneff_train){
  temp1 <- as.data.frame(eff_train)
  temp1$Tumor = "Eff" 
  temp1 <- temp1%>%dplyr::select("Tumor")
  
  temp2 <- as.data.frame(noneff_train)
  temp2$Tumor = "Noneff" 
  temp2 <- temp2%>%dplyr::select("Tumor")
  return(bind_rows(temp1,temp2))
}

label_train <- create_label(eff_train,noneff_train)
label_test <- create_label(eff_test,noneff_test)

cooclassifier<-function(data=lstaudt_exp,n=50,cutoff=0.8)
{
  staudt.train<-data%>%dplyr::select(rownames(label_train))
  staudt.test<-data%>%dplyr::select(rownames(label_test))
  print(dim(staudt.train))
  scale.normalization<-function(matrix)
  {
    scale.norm<-t(apply(matrix,1,scale))##log(t(apply(matrix,1,scale))+10)
    return(scale.norm)
  }
  p.subgroup<-function(pred,Eff.pred,Noneff.pred)
  {
    mean.Eff<-mean(Eff.pred)
    sd.Eff<-sd(Eff.pred)
    mean.Noneff<-mean(Noneff.pred)
    sd.Noneff<-sd(Noneff.pred)
    p.Eff<-dnorm(pred,mean.Eff,sd.Eff)/(dnorm(pred,mean.Eff,sd.Eff)+dnorm(pred,mean.Noneff,sd.Noneff))
    p.Noneff<-1-p.Eff
    p<-rbind(p.Eff,p.Noneff)
    rownames(p)<-c("p.Eff","p.Noneff")
    return(p)
  }
  staudt.scale<-scale.normalization(staudt.train)
  #remove rows with na
  row.has.na <- apply(staudt.scale, 1, function(x){any(is.na(x))})
  print(sum(row.has.na))
  staudt.scale <- staudt.scale[!row.has.na,]
  colnames(staudt.scale)<-colnames(staudt.train)
  ## select all labeled microarray data as train set
  Eff.train<-staudt.scale[,label_train=="Eff"]
  Noneff.train<-staudt.scale[,label_train=="Noneff"]
  ## select predictors
  p.value.scale<-vector()
  t.statistic.scale<-vector()
  for(i in 1:dim(Eff.train)[1])
  {
    t.test.result<-t.test(Eff.train[i,],Noneff.train[i,])
    p.value.scale[i]<-t.test.result$p.value
    t.statistic.scale[i]<-t.test.result$statistic
  }
  names(p.value.scale)<-rownames(Eff.train)
  names(t.statistic.scale)<-rownames(Eff.train)
  FDR<-p.adjust(p.value.scale,"BH")
  t.statistic<-t.statistic.scale[rank(FDR)<=n]
  Eff.pred<-t.statistic %*% Eff.train[rank(FDR)<=n,]
  Noneff.pred<-t.statistic %*% Noneff.train[rank(FDR)<=n,]
  
  DLBCL<-as.matrix(staudt.test)
  rownames(DLBCL)<-toupper(rownames(DLBCL))
  ## prepare predictors
  pred<-names(FDR[rank(FDR)<=n])[names(FDR[rank(FDR)<=n]) %in% rownames(DLBCL)]
  print(pred)
  DLBCL.pred<-scale.normalization(DLBCL[pred,])
  colnames(DLBCL.pred)<-colnames(DLBCL)
  DLBCL.LPS<-t.statistic[pred] %*% DLBCL.pred
  Eff.pred.test<-t.statistic[pred] %*% Eff.train[pred,]
  Noneff.pred.test<-t.statistic[pred] %*% Noneff.train[pred,]
  test.subgroup<-p.subgroup(DLBCL.LPS,Eff.pred.test,Noneff.pred.test)
  sort.DLBCL<-sort(test.subgroup[1,],decreasing=T)
  list<-c(rep("unclassified",dim(DLBCL)[2]))
  list[test.subgroup[1,]>cutoff]<-"Eff"
  list[test.subgroup[2,]>cutoff]<-"Noneff"
  names(list)<-colnames(test.subgroup)
  
  totallist = vector('list', 2)
  totallist[[1]] <- list
  totallist[[2]] <- pred
  totallist[[3]] <- test.subgroup
  return(totallist)
}

lstaudt_exp <- read.csv("/Users/chengqi_xu/Documents/Elemento lab/synergyy_r/crc_organoid/orginal_lstaudt_DLBC_exp.csv")
#lstaudt_exp <- read.csv("/Users/chengqi_xu/Documents/Elemento lab/synergyy/data/cell_line_data/Customized/TCGA_PAN_RNA_TPM/tcga_DLBC_RNA_counts_csv_tpm.csv")
lstaudt_exp <- unique(lstaudt_exp)
rownames(lstaudt_exp) <- lstaudt_exp$V1
lstaudt_exp <- lstaudt_exp[, -1]
selected_genes <- read.table("/Users/chengqi_xu/Documents/Elemento lab/synergyy/results/selected_genes.txt",header=F)
lstaudt_exp<-lstaudt_exp %>% filter(row.names(lstaudt_exp) %in% selected_genes$V1)
rownames(lstaudt_exp) <- selected_genes$V1
result <- cooclassifier(data=lstaudt_exp,n=50,cutoff=0.7)

#Creates vectors having data points
expected_value <- factor(label_test$Tumor)
predicted_value <- factor(result[[1]])
#Creating confusion matrix
example <- confusionMatrix(data=expected_value, reference = predicted_value)
#Display results
print(example)
# pfs plot
predicted_value_df <- as.data.frame(predicted_value)
predicted_value_df$subject = rownames(predicted_value_df)
metadata_df <- left_join(metadata,predicted_value_df)%>%dplyr::filter(predicted_value %in% c("Eff","Noneff"))
fit <- survfit(Surv(Progression_Free.Survival._PFS_.Time._yrs) ~ predicted_value, data = metadata_df)

# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
          )

test.subgroup <- result[[3]]
a <- as.data.frame(test.subgroup[2,])
a$tumors <- rownames(a)
colnames(a) <- c("value","tumor")
a$type = "Noneff"
b <- as.data.frame(test.subgroup[1,])
b$tumors <- rownames(b)
colnames(b) <- c("value","tumor")
b$type = "Eff"
a %>%
  mutate(name = fct_reorder(tumor, value)) %>%
  ggplot( aes(x=name, y=value)) + geom_point(size = 0.5,color="#00008B")#AA336A

b %>%
  mutate(name = fct_reorder(tumor, value, .desc = TRUE)) %>%
  ggplot( aes(x=name, y=value)) + geom_point(size = 0.5,color="#AA336A")#

c<- rbind(a,b)
c %>%
  mutate(name = fct_reorder(tumor, value)) %>%
  ggplot( aes(x=name, y=value,fill=type)) + geom_point(size = 0.5)

#"


#exp
scale.normalization<-function(matrix)
{
  scale.norm<-log(t(apply(matrix,1,scale))+10)
  return(scale.norm)
}
lstaudt_exp.scale<-scale.normalization(lstaudt_exp)
colnames(lstaudt_exp.scale) <- colnames(lstaudt_exp) 
rownames(lstaudt_exp.scale) <- rownames(lstaudt_exp) 
pred <- result[[2]]
staudt.test <-as.data.frame(lstaudt_exp.scale)%>%dplyr::select(rownames(label_test))%>% dplyr::filter(row.names(.) %in% pred)
rownames(staudt.test) <- pred
pheatmap(staudt.test,fontsize = 3,cluster_cols = FALSE)

