library(reshape2); 
library(vegan); 
library(dplyr)
library(ade4); 
library(plotly)
library(compositions); 
library(pracma); 
library(DESeq2); 
library(fpc); 
library(tidyverse)
library(purrr)
library(cluster)
library(RColorBrewer)
library(ape)
setwd("/Users/taylorsehein/github/WHOI_repo")
raw <- read.delim("SP2018_iTags_2m_50.txt")
env_data <- read.delim("SP2018_SampleKey.txt")
head(raw)
colnames(raw)
str(raw)
summary(raw)
count.filtered <- raw
names(count.filtered)
seq_counts<-count.filtered[1:116]
tax_key<-count.filtered[c(1,117)]; head(tax_key[1:2,])
head(seq_counts[1:3])
row.names(seq_counts)<-seq_counts$OTU_ID
seq_counts$OTU_ID<-NULL; head(seq_counts[1:3])
covariance_matrix<-as.matrix(seq_counts)%*%t(seq_counts)
str(seq_counts)
cov_determinant<-det(covariance_matrix)
cov_determinant
log_rats<-data.frame(compositions::ilr(t(seq_counts)))
new_covdet<-det(as.matrix(log_rats)%*%t(log_rats))
cov_determinant #Original Count Data
new_covdet
lograt_pca<-prcomp(log_rats)
lograt_variances<-lograt_pca$sdev^2/sum(lograt_pca$sdev^2)
barplot(lograt_variances,
        main='Log-Ratio PCA Screeplot',
        xlab='Axis',
        ylab='% Variance',
        col=c(rep('black',1),rep('grey',18)))
legend('topright',fill=c('black','grey'),c('Should Present','??? Judgment Call'))
  
row.names(lograt_pca$x)
head(lograt_pca)

pca_lograt_frame<-data.frame(lograt_pca$x,
                             month=gsub('_.*','',rownames(lograt_pca$x)))
ggplot(pca_lograt_frame)+
  geom_point(aes(x=PC1,y=PC2,col=month))+
  ylab(paste0('PC2 ',round(lograt_variances[2]*100,2),'%'))+
  xlab(paste0('PC1 ',round(lograt_variances[1]*100,2),'%'))+
  scale_color_brewer(palette='Set1',name='Month')+
  ggtitle('Log-Ratio PCA Ordination')+
  coord_fixed(ratio=lograt_variances[2]/lograt_variances[1])+
  theme_bw()

jac_dmat<-vegdist(t(seq_counts),method="jaccard")
pcoa_jac<-ape::pcoa(jac_dmat)
samp_no<-dim(seq_counts)[2]
jac_variances<-pcoa_jac$values$Relative_eig
barplot(jac_variances,
        main='Jaccard PCoA Screeplot',
        xlab='Axis',
        ylab='% Variance',
        col=c(rep('black',3),rep('grey',16)),
        cex.lab=1.5,cex.main=1.5)
legend('topright',fill=c('black','grey'),c('Should Present','Unnecessary to Present'),cex=1.5)

euc_dmat<-dist(log_rats) 
pcoa_euc<-ape::pcoa(euc_dmat)
euc_variances<-pcoa_euc$values$Relative_eig
barplot(euc_variances,
        main='Euclidean PCoA Screeplot',
        xlab='Axis',
        ylab='% Variance',
        col=c(rep('black',2),rep('darkgrey',2),rep('lightgrey',17)),
        cex.main=1.5,cex.lab=1.5)
legend('topright',fill=c('black','darkgrey','lightgrey'),c('Should Present','Questionable','Unnecessary'),cex=1.5)

##Error, suggests scatter3d
pcoa_jac_frame<-data.frame(pcoa_jac$vectors,month=gsub('_.*','',rownames(pcoa_jac$vectors)))
eigenvalues<-round(jac_variances,4)*100
plot_ly(pcoa_jac_frame,x=~Axis.1,y=~Axis.2,z=~Axis.3,colors=~brewer.pal(6,'Set1'),color=~month)%>%
  layout(title='PCoA Jaccard Distance',
         scene=list(xaxis=list(title=paste0('Axis1 ',eigenvalues[1],'%'),
                               scale=eigenvalues[1]/100),
                    yaxis=list(title=paste0('Axis2 ',eigenvalues[2],'%'),
                               scale=eigenvalues[2]/100),
                    zaxis=list(title=paste0('Axis3 ',eigenvalues[3],'%'),
                               scale=eigenvalues[3]/100),
                    aspectratio = list(x=3, y=3*eigenvalues[2]/eigenvalues[1], z=3*eigenvalues[3]/eigenvalues[1])))

##Error, suggests scatter3d
pcoa_euc_frame<-data.frame(pcoa_euc$vectors,month=gsub('_.*','',rownames(pcoa_euc$vectors)))
euc_eigenvalues<-round(euc_variances,4)*100
plot_ly(pcoa_euc_frame,x=~Axis.1,y=~Axis.2,z=~Axis.3,colors=~brewer.pal(6,'Set1'),color=~month)%>%
  layout(title='PCoA Euclidean Distance',
         scene=list(xaxis=list(title=paste0('Axis1 ',euc_eigenvalues[1],'%'),
                               scale=euc_eigenvalues[1]/100),
                    yaxis=list(title=paste0('Axis2 ',euc_eigenvalues[2],'%'),
                               scale=euc_eigenvalues[2]/100),
                    zaxis=list(title=paste0('Axis3 ',euc_eigenvalues[3],'%'),
                               scale=euc_eigenvalues[3]/100),
                    aspectratio = list(x=3, 
                                       y=3*euc_eigenvalues[2]/euc_eigenvalues[1], 
                                       z=3*euc_eigenvalues[3]/euc_eigenvalues[1])))

cluster_ex<-hclust(vegdist(t(seq_counts),method='jaccard'),method="average")
plot(cluster_ex,main='Jaccard Hierarchical Clustering',xlab='',sub='')
set.seed(071510)
euc_nmds<-metaMDS(euc_dmat,k=2,autotransform=FALSE)
jac_nmds<-metaMDS(jac_dmat,k=2,autotransform=FALSE)
euc_nmds$stress
jac_nmds$stress
euc_frame<-data.frame(euc_nmds$points,
                      month=gsub('_.*','',rownames(log_rats)))
jac_frame<-data.frame(jac_nmds$points,
                      month=gsub('_.*','',rownames(log_rats)))
ggplot(euc_frame,aes(x=MDS1,y=MDS2,col=month))+
  geom_point(size=2)+
  scale_color_brewer(palette='Set1',name='Month')+
  theme_bw()+ggtitle('Euclidean Distance NMDS')
ggplot(jac_frame,aes(x=MDS1,y=MDS2,col=month))+
  geom_point(size=2)+
  scale_color_brewer(palette='Set1',name='Month')+
  theme_bw()+ggtitle('Jaccard Distance NMDS')

within_seq_means<-apply(t(seq_counts),2,mean)
within_seq_vars<-apply(t(seq_counts),2,var)
plot(within_seq_means,within_seq_vars,log='xy',
     main='Heteroskedastic Data',
     xlab='Mean # Counts',
     ylab='Var # Counts')
