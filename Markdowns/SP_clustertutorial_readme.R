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
transformed_counts<-DESeq2::varianceStabilizingTransformation(as.matrix(seq_counts_duplicate))

seq_counts_duplicate[seq_counts_duplicate == 0] <- 1
View(seq_counts_duplicate)

within_trans_means<-apply(transformed_counts,1,mean)
within_trans_vars<-apply(transformed_counts,1,var)
# Plot:
plot(within_trans_means,within_trans_vars)

transformed_detrended<-apply(transformed_counts,1,pracma::detrend)
trans_dt_scaled<-apply(transformed_detrended,2,scale)
trans_dt_scaled<-apply(transformed_detrended,2,scale)
rownames(trans_dt_scaled)<-colnames(transformed_counts)
hist(transformed_counts,
     main='VST Sequence Count Observations',
     xlab='# Total Observations',
     ylab='Frequency')
hist(transformed_detrended,
     main='VST+Detrended Data',
     xlab='Total Observations Accounting for Linear Trends',
     ylab='Frequency')
hist(trans_dt_scaled,
     main='VST+Detrended+Scaled Data (Z-scores)',
     xlab='Expression Level Relative to per-OTU Avg',
     ylab='Frequency')
temporal_dmat<-dist(t(trans_dt_scaled))
n_clusts<-2:20
hc_full_cluster<-hclust(temporal_dmat)
hc_clusts<-lapply(n_clusts,function(x) cutree(hc_full_cluster,x))
kmed_clusts<-lapply(n_clusts, function(x) cluster::pam(temporal_dmat,k=x))
hc_stats<-lapply(hc_clusts,function(x) fpc::cluster.stats(temporal_dmat,
                                                          clustering=x))
kmed_stats<-lapply(kmed_clusts, function(x) fpc::cluster.stats(temporal_dmat,
                                                               clustering=x$clustering))
ripping_stats<-function(list,func){
  ## Essentially all this function does is implements a function (func) on a list
  ## and coerces the output to a column vector (this will be handy when we want to make a data frame)
  output<-do.call(rbind,lapply(list,func))
  return(output)
}
func_list<-rep(list(function(x) x$cluster.number,
                    function(x) x$within.cluster.ss,
                    function(x) x$avg.silwidth,
                    function(x) x$ch),2)
stats_list<-rep(list(hc_stats,kmed_stats),each=4)
collected_stats<-purrr::map2(stats_list,func_list,ripping_stats)
nclusts<-rep(n_clusts,length(collected_stats))
method<-rep(c('hc','kmed'),each=length(n_clusts)*length(collected_stats)/2)
ind_name<-rep(c('n','ss','sil','ch'),each=length(n_clusts)) %>%
  rep(2)
index_frame<-data.frame(index=do.call(rbind,collected_stats),
                        nc=nclusts,
                        Method=method,
                        ind=ind_name)
index_frame %>%
  filter(ind=='ss') %>%
  ggplot(aes(x=nc,y=index,col=Method)) +
  geom_point() +
  geom_line(aes(group=Method)) +
  ylab('Within Cluster Sum Square Error') +
  xlab('Number Clusters')+
  theme_bw()
index_frame %>%
  filter(ind=='ch') %>%
  ggplot(aes(x=nc,y=index,col=Method)) +
  geom_point() +
  geom_line(aes(group=Method)) +
  ylab('C-H Stat') +
  xlab('Number Clusters')+
  theme_bw()
silhouette_profile_kmed8<-cluster::silhouette(kmed_clusts[[3]])
silhouette_frame<-data.frame(cluster=silhouette_profile_kmed8[,1],
                             neighbor=silhouette_profile_kmed8[,2],
                             dist=silhouette_profile_kmed8[,3],
                             otu_id=rownames(silhouette_profile_kmed8))
new_silframe<-silhouette_frame %>%
  arrange(dist) %>%
  group_by(cluster) %>%
  mutate(idno=1:n(),
         tot_num=n())
ggplot(new_silframe,aes(x=idno,y=dist))+
  geom_bar(stat='identity') +
  coord_flip() +
  ggtitle('Silhouette Profiles for Individual Clusters') +
  facet_wrap(~cluster)
medoid_otus<-kmed_clusts[[3]]$medoids
medoid_dynamics<-trans_dt_scaled[,medoid_otus]
head(silhouette_frame)
colnames(tax_key)[1]<-"otu_id"
tax_bycluster<-left_join(silhouette_frame, tax_key, by="otu_id") 
str(tax_key)
tax_key$otu_id <- as.factor(tax_key$otu_id)
str(tax_key)
head(tax_bycluster[1:2,])
tax_bycluster_sum<-tax_bycluster %>%
  group_by(cluster, level_0) %>%
  summarise(richness =  n()) %>% #Get richness information
  group_by(level_0) %>%
  mutate(n_per_tax=sum(richness)) %>% #Figure out total number of OTUs assigned to each taxon
  as.data.frame
tax_bycluster_sum$level_0[which(tax_bycluster_sum$level_0=='Opisthokonts')]<-'Opisthokont' # Cleaning this up
head(tax_bycluster_sum)
tax_color<-c("#67000d","#e31a1c","#dd3497","#fcbba1","#fed976","#fc8d59","#a63603","#addd8e","#7f2704","#238b45","#a1d99b","#081d58","#1f78b4","#a6cee3","#8c6bb1","#9e9ac8","#984ea3","#081d58","#662506","#ffffff","#969696","#525252","#000000")
ggplot(tax_bycluster_sum, aes(x=cluster, y=richness, fill=level_0))+
  geom_bar(color="black", stat="identity")+
  coord_flip()+
  scale_fill_brewer(palette='Set3')+
  theme_minimal()+
  scale_x_continuous(breaks=1:4,labels=as.character(1:4),name='Cluster #')+
  scale_y_continuous(name='# Species')+
  theme(legend.position='bottom',text=element_text(size=12))
medoid_dyn_long<-as.data.frame(t(medoid_dynamics)) %>%
  mutate(otu_id=colnames(medoid_dynamics),
         clust_num=1:4) %>%
  gather(timepoint,z_score,'Mar_19_rep1':'Oct_10_rep2') %>%
  mutate(time_numeric=as.numeric(gsub('_.*','',timepoint)),
         plotting_timepoint=paste('Month',ceiling(time_numeric/6),gsub('_.*','',timepoint)))
ggplot()+
  geom_point(size=4,
             shape=21, 
             color="white", 
             aes(fill=factor(clust_num),
                 x=time_numeric,
                 y=z_score),
             data=medoid_dyn_long)+
  geom_line(data=medoid_dyn_long,
            aes(x=time_numeric,
                y=z_score,
                group=otu_id,
                col=factor(clust_num)),size=1.25,linetype=5)+
  facet_wrap(~clust_num,ncol=2)+
  scale_color_brewer(palette='Set2',name='Cluster #')+
  scale_fill_brewer(palette='Set2',guide=FALSE)+
  scale_x_continuous(breaks=1:8,labels=unique(medoid_dyn_long$plotting_timepoint))+
  scale_y_continuous(name='Z-score Expression',limits=c(-2,4.5))+
  theme(axis.text.x=element_text(angle=90,hjust=0.5),
        panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color='gray'),
        axis.title.x=element_blank())+
  geom_rect(data=data.frame(ymin=rep(-2,4),
                            ymax=rep(4.5,4),
                            xmin=seq(1,8,by=1),
                            xmax=c(seq(1,8,by=1),8)),
            aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax),
            col='gray',
            alpha=0.25)
