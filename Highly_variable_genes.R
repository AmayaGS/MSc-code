#====================================================================================

# Dissertation
# Amaya Syed - ID: 190805496

#====================================================================================

# cleaning any previous objects
rm(list = ls())
ls()
rm()

install.packages("tidyverse")

install.packages("doMPI")

# install necessary packages
install.packages("doParallel")

# loading necessary packages
library(rstudioapi)
library(matrixStats)
library(igraph)
library(cluster)
library(Rtsne)
library(foreach)
library(doParallel)
library(parallel)
library(tidyverse)
library("tidyr")
library("ggplot2")
library("data.table")
library("Rtsne")
library("matrixStats")
library("dbscan")
library("class")
library("pals")
library("foreach")
library("doParallel")
library("parallel")
library("doMPI")
numCores<-detectCores()

registerDoParallel(numCores)


# Setting work directory
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))


# loading the whole dataset to identify the of genes with highest variance
mat_all <- read.csv("rna_seq_peac1.csv", header= TRUE, sep='\t')

expr <- as.data.frame(mat_all)
rownames(expr) <- expr[,1]
expr[,1] <- NULL

# loading the subsetted data to identify and keep the list of the genes with the highest variance
mat = read.csv("rna_seq_peac1.txt", header= TRUE, sep=',')
mat_t = t(mat)
df <- as.data.frame(mat_t)
names(df) <- df %>% slice(1) %>% unlist()
df <- df %>% slice(-1)

expr <- data.frame(lapply(df, as.numeric), check.names=F, row.names = rownames(df))
expr <- expr[rowMeans(as.matrix(expr)) >= 1,]
expr <-na.omit(expr)

#############################

# startn with expr from whole data or subsetted data. This code is taken from Nikolai Oskolkov's github on scRNAseq -

mean_expr <- as.numeric(rowMeans(expr,na.rm=TRUE))
log10_mean_expr<-log10(mean_expr)
names(log10_mean_expr)<-rownames(expr)
sd_expr<-rowSds(as.matrix(expr),na.rm=TRUE)
cv_squared_expr<-(sd_expr/mean_expr)^2
log10_cv_squared_expr<-log10(cv_squared_expr)
names(log10_cv_squared_expr)<-rownames(expr)
plot(log10_cv_squared_expr~log10_mean_expr,pch=20,cex=0.5,xlab="LOG10 ( MEAN EXPRESSION )", ylab="LOG10 ( CV2 )")

# plot log(mean) vs log(cv)
N_bin<-100
bin_size<-(max(log10_mean_expr)-min(log10_mean_expr))/N_bin
mean_cv_squared_points<-vector(length=N_bin)
sd_cv_squared_points<-vector(length=N_bin)
mean_expr_points<-vector(length=N_bin)
a<-min(log10_mean_expr)
for(i in 1:N_bin)
{
        b<-a+bin_size
        mean_cv_squared_points[i]<-mean(log10_cv_squared_expr[log10_mean_expr>=a & log10_mean_expr<=b],na.rm=TRUE)
        sd_cv_squared_points[i]<-sd(log10_cv_squared_expr[log10_mean_expr>=a & log10_mean_expr<=b],na.rm=TRUE)
        mean_expr_points[i]<-mean(log10_mean_expr[log10_mean_expr>=a & log10_mean_expr<=b],na.rm=TRUE)
        a<-b
}
points(mean_cv_squared_points ~ mean_expr_points,col="red",pch=20)
boundary_cv_squared_points<-mean_cv_squared_points + 2*sd_cv_squared_points

# fit loess curve to the data
fit_mean_cv_squared_points <- loess(mean_cv_squared_points[is.finite(mean_cv_squared_points)&
                                                                   is.finite(mean_expr_points)]~
                                            mean_expr_points[is.finite(mean_cv_squared_points)&
                                                                     is.finite(mean_expr_points)],span=1)

# find boundaries over which genes are considered highly variable
fit_boundary_cv_squared_points <- loess(boundary_cv_squared_points[is.finite(boundary_cv_squared_points)&
                                                                           is.finite(mean_expr_points)]~
                                                mean_expr_points[is.finite(boundary_cv_squared_points)
                                                                 &is.finite(mean_expr_points)],span=1)
# subset those genes
pred_cv_squared_expr <- predict(fit_mean_cv_squared_points, log10_mean_expr)
pred_boundary_cv_squared_expr <- predict(fit_boundary_cv_squared_points, log10_mean_expr)

expr_pruned <- expr[log10_cv_squared_expr >= pred_boundary_cv_squared_expr,]
expr_pruned <- expr_pruned[grepl("NA",rownames(expr_pruned))==FALSE,]
text(log10_mean_expr[rownames(expr_pruned)],log10_cv_squared_expr[rownames(expr_pruned)],
     labels=rownames(expr_pruned),cex=0.5)
lines(fit_mean_cv_squared_points$fitted~
              mean_expr_points[is.finite(mean_expr_points)],col="red",lwd=3)

expr_hv <- t(expr_pruned)

# keep highly variable genes
hv_genes <- rownames(expr_pruned)

# write them to a csv file
write.csv(expr_hv, 'dataset_expr.csv')
write.csv(hv_genes, 'dataset_hv_genes.csv', row.names=FALSE)



# heatmap


#exp.data <- read.csv(file = "data/GSE68849-expression.csv",
#                     stringsAsFactors = FALSE)

exp.long <- pivot_longer(data = hv_genes,
                         cols = -c(subject, treatment),
                         names_to = "gene",
                         values_to = "expression")

exp.long$log.expression <- log(exp.long$expression)

exp.heatmap <- ggplot(data = exp.long, mapping = aes(x = subject,
                                                     y = gene,
                                                     fill = log.expression)) +
        geom_tile() +
        xlab(label = "Subject") +
        facet_grid(~ treatment, switch = "x", scales = "free_x", space = "free_x") +
        theme(axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(filename = "output/expression-heatmap.pdf", plot = exp.heatmap)
