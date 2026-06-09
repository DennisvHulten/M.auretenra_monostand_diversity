# de novo DAPC
# Alejandra Hernandez

## Dependencies ==============================================================
suppressMessages(library("adegenet"))
suppressMessages(library("vcfR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("svglite"))
suppressMessages(library("gridExtra"))


## Main code =================================================================
## Directory
setwd("/Users/alejhernandez/Dropbox/Rare species/scolymia/dart/F - Overall genetic structure whole/F3 - PCA_DAPC")


## Main code =================================================================
## Import with vcfR and popfile (as medatada, no prior)
coral.genlight <- vcfR2genlight(read.vcfR("sco_c3_80.vcf"))          #643 loci with more than 2 alelles, ignored

# Import metadata
popfile <- read.table("sco_popfile_c3_80.txt", header=FALSE, sep="\t")
colnames(popfile)[1] <- "indv"
colnames(popfile)[2] <- "pop"


####################################################################################################
#K=2
grp2 <- find.clusters(coral.genlight, stat = "BIC", 
                      n.iter = 100000, n.pca = 120, 
                      choose.n.clust = 2,
                      n.start = 1000, set.seed(2))

#Choose the number of clusters (>=2 - assess BIC): 
#2

pop2 <- grp2$grp
max_PCAs2 <- as.integer(length(pop2) / 3) # as <= N/3 advised
dapc2 <- dapc(coral.genlight, pop2, n.pca = max_PCAs2, n.da = 10)
optimum_score2 <- optim.a.score(dapc2)              #1 PCs

dapc.pop2 <- dapc(coral.genlight, pop2, n.pca = optimum_score2$best, n.da = 10)

#Extract coordinates and add metadata
dapc_coord2 <- as.data.frame(dapc.pop2$ind.coord)
dapc_assig2 <- as.data.frame(dapc.pop2$assign)                    #Posterior assignment
dapc_coord2$indv <- rownames(dapc_coord2)
dapc_coord2 <- bind_cols(dapc_coord2, dapc_assig2)
colnames(dapc_coord2)[3] <- "assignment"                    #Hardcoded
dapc_coord2$pop <- with(popfile,pop[match(dapc_coord2$indv,indv)])

#Rename assignment
dapc_coord2$assignment <- as.character(dapc_coord2$assignment)
dapc_coord2$assignment[dapc_coord2$assignment == "1"] <- "Cluster 1"
dapc_coord2$assignment[dapc_coord2$assignment == "2"] <- "Cluster 2"

#Setting site order 
dapc_coord2$pop <- as.factor(dapc_coord2$pop)
levels(dapc_coord2$pop)
dapc_coord2$pop <- factor(dapc_coord2$pop, levels = c("KAL","HUL","EST","SNA","SEA","MAR","DIR"))

# #Setting site shapes
colors <- c("#2e294e","#f46036","#d7263d","#1b998b","#0496ff","#c5d86d","#00f5d4",
             "#001219", "#005f73","#0a9396", "#94d2bd", "#e9d8a6","#ee9b00", "#ca6702", "#bb3e03",
             "#ae2012", "#9b2226", "#6a040f", "#03071e","#e3d5ca")

#Extract variance explained by eigenvectors
percent_dapc.pop2 <- as.data.frame(dapc.pop2$eig/sum(dapc.pop2$eig)*100)



#Common across graphs
clean_theme <-  theme(axis.title = element_text(size=16),
                      axis.line = element_line(colour = "black"),
                      axis.text = element_text(size=16),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.grid = element_blank(),
                      panel.background = element_blank())

# clean_theme2 <-  theme(axis.title = element_text(size=16),
#                        axis.line = element_line(colour = "black"),
#                        axis.text = element_text(size=16),
#                        panel.grid.major = element_blank(),
#                        panel.grid.minor = element_blank(),
#                        panel.border = element_blank(),
#                        panel.grid = element_blank(),
#                        panel.background = element_blank(),
#                        legend.position ="none")
# 
# theme2 <- theme(legend.position = "none",
#                 axis.text = element_text(size = 24),
#                 axis.title = element_text(size = 24),
#                 axis.ticks = element_line(),
#                 panel.grid = element_blank(),
#                 panel.border = element_rect(size = 1, fill = NA),
#                 panel.background = element_blank(),
#                 axis.text.y = element_text(size=24),
#                 axis.text.x = element_text(size=24))


#Graph 
DAPC_K2 <- ggplot(data = dapc_coord2,aes(x=LD1, fill=assignment)) + 
  geom_area(stat="bin", alpha=0.6,binwidth = 1)+ 
  scale_fill_manual("K = 2", values=colors)+
  clean_theme+
  labs(x = "Discriminant function 1", y = "Density")+
  theme(legend.position=c(0.9,0.9))#+
#  scale_x_continuous(limits = c(-8,4), breaks= c(-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4))

DAPC_K2


DAPC_K2_scatter <- scatter(dapc.pop2, col=colors)

#ggsave(plot = DAPC_K2, filename = "denovo_DAPC_K2_sco_whole.png",dpi = 300, limitsize = FALSE, width = 12, height = 8)



########################################################################################
#K=3
grp3 <- find.clusters(coral.genlight, stat = "BIC", 
                      n.iter = 100000, n.pca = 120, 
                      choose.n.clust = 3,
                      n.start = 1000, set.seed(2))

#Choose the number of clusters (>=2 - assess BIC): 
#3

pop3 <- grp3$grp
max_PCAs3 <- as.integer(length(pop3) / 3) # as <= N/3 advised
dapc3 <- dapc(coral.genlight, pop3, n.pca = max_PCAs3, n.da = 10)
optimum_score3 <- optim.a.score(dapc3)              #1 PCs

dapc.pop3 <- dapc(coral.genlight, pop3, n.pca = optimum_score3$best, n.da = 10)

#Extract coordinates and add metadata
dapc_coord3 <- as.data.frame(dapc.pop3$ind.coord)
dapc_assig3 <- as.data.frame(dapc.pop3$assign)                    #Posterior assignment
dapc_coord3$indv <- rownames(dapc_coord3)
dapc_coord3 <- bind_cols(dapc_coord3, dapc_assig3)
colnames(dapc_coord3)[3] <- "assignment"                    #Hardcoded
dapc_coord3$pop <- with(popfile,pop[match(dapc_coord3$indv,indv)])

#Rename assignment
dapc_coord3$assignment <- as.character(dapc_coord3$assignment)
dapc_coord3$assignment[dapc_coord3$assignment == "1"] <- "Cluster 1"
dapc_coord3$assignment[dapc_coord3$assignment == "2"] <- "Cluster 2"
dapc_coord3$assignment[dapc_coord3$assignment == "3"] <- "Cluster 3"

#Setting site order 
dapc_coord3$pop <- as.factor(dapc_coord3$pop)
levels(dapc_coord3$pop)
dapc_coord3$pop <- factor(dapc_coord3$pop, levels = c("KAL","HUL","EST","SNA","SEA","MAR","DIR"))

#Graph 
DAPC_K3 <- ggplot(data = dapc_coord3,aes(x=LD1, fill=assignment)) + 
  geom_area(stat="bin", alpha=0.6,binwidth = 1)+ 
  scale_fill_manual("K = 2", values=colors)+
  clean_theme+
  labs(x = "Discriminant function 1", y = "Density")+
  theme(legend.position=c(0.9,0.9))#+
#  scale_x_continuous(limits = c(-8,4), breaks= c(-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4))

DAPC_K3


DAPC_K3_scatter <- scatter(dapc.pop3, col=colors)

#ggsave(plot = DAPC_K3, filename = "denovo_DAPC_K3_sco_whole.png",dpi = 300, limitsize = FALSE, width = 12, height = 8)



########################################################################################
#K=4
grp4 <- find.clusters(coral.genlight, stat = "BIC", 
                      n.iter = 100000, n.pca = 120, 
                      choose.n.clust = 4,
                      n.start = 1000, set.seed(2))

#Choose the number of clusters (>=2 - assess BIC): 
#4

pop4 <- grp4$grp
max_PCAs4 <- as.integer(length(pop4) / 3) # as <= N/3 advised
dapc4 <- dapc(coral.genlight, pop4, n.pca = max_PCAs4, n.da = 10)
optimum_score4 <- optim.a.score(dapc4)              #22

dapc.pop4 <- dapc(coral.genlight, pop4, n.pca = optimum_score4$best, n.da = 10)

#Extract coordinates and add metadata
dapc_coord4 <- as.data.frame(dapc.pop4$ind.coord)
dapc_assig4 <- as.data.frame(dapc.pop4$assign)                    #Posterior assignment
dapc_coord4$indv <- rownames(dapc_coord4)
dapc_coord4 <- bind_cols(dapc_coord4, dapc_assig4)
colnames(dapc_coord4)[5] <- "assignment"                    #Hardcoded
dapc_coord4$pop <- with(popfile,pop[match(dapc_coord4$indv,indv)])

#Rename assignment
dapc_coord4$assignment <- as.character(dapc_coord4$assignment)
dapc_coord4$assignment[dapc_coord4$assignment == "1"] <- "Cluster 1"
dapc_coord4$assignment[dapc_coord4$assignment == "2"] <- "Cluster 2"
dapc_coord4$assignment[dapc_coord4$assignment == "3"] <- "Cluster 3"
dapc_coord4$assignment[dapc_coord4$assignment == "4"] <- "Cluster 4"

#Setting site order 
dapc_coord4$pop <- as.factor(dapc_coord4$pop)
levels(dapc_coord4$pop)
dapc_coord4$pop <- factor(dapc_coord4$pop, levels = c("KAL","HUL","EST","SNA","SEA","MAR","DIR"))

#Graph 
DAPC_K4 <- ggplot() + geom_point(data = dapc_coord4, 
                                 aes(x=LD1, y=LD2, color=as.character(assignment),
                                     alpha =0.99),size = 3.5) + 
  scale_color_manual("K = 4", values=colors)+
  #  scale_fill_manual("pop", values=def_shapes)+
  clean_theme+
  labs(x = "LD1", y = "LD2")+#theme(legend.position=c(0.1,0.8))+
  ggtitle("K = 4")

DAPC_K4


DAPC_K4_scatter <- scatter(dapc.pop4, col=colors)


ggsave(plot = DAPC_K4, filename = "denovo_DAPC_K4_sco_whole.png",dpi = 300, limitsize = FALSE, width = 12, height = 8)



########################################################################################
#K=5
grp5 <- find.clusters(coral.genlight, stat = "BIC", 
                      n.iter = 100000, n.pca = 120, 
                      choose.n.clust = 5,
                      n.start = 1000, set.seed(2))

#Choose the number of clusters (>=2 - assess BIC): 
#5


pop5 <- grp5$grp
max_PCAs5 <- as.integer(length(pop5) / 3) # as <= N/3 advised
dapc5 <- dapc(coral.genlight, pop5, n.pca = max_PCAs5, n.da = 10)
optimum_score5 <- optim.a.score(dapc5)              #1

dapc.pop5 <- dapc(coral.genlight, pop5, n.pca = optimum_score5$best, n.da = 10)

#Extract coordinates and add metadata
dapc_coord5 <- as.data.frame(dapc.pop5$ind.coord)
dapc_assig5 <- as.data.frame(dapc.pop5$assign)                    #Posterior assignment
dapc_coord5$indv <- rownames(dapc_coord5)
dapc_coord5 <- bind_cols(dapc_coord5, dapc_assig5)
colnames(dapc_coord5)[3] <- "assignment"
dapc_coord5$pop <- with(popfile,pop[match(dapc_coord5$indv,indv)])

#Rename assignment
dapc_coord5$assignment <- as.character(dapc_coord5$assignment)
dapc_coord5$assignment[dapc_coord5$assignment == "1"] <- "Cluster 1"
dapc_coord5$assignment[dapc_coord5$assignment == "2"] <- "Cluster 2"
dapc_coord5$assignment[dapc_coord5$assignment == "3"] <- "Cluster 3"
dapc_coord5$assignment[dapc_coord5$assignment == "4"] <- "Cluster 4"
dapc_coord5$assignment[dapc_coord5$assignment == "5"] <- "Cluster 5"

#Setting site order 
dapc_coord5$pop <- as.factor(dapc_coord5$pop)
levels(dapc_coord5$pop)
dapc_coord5$pop <- factor(dapc_coord5$pop, levels = c("KAL","HUL","EST","SNA","SEA","MAR","DIR"))

#Graph 
DAPC_K5 <- ggplot(data = dapc_coord5,aes(x=LD1, fill=assignment)) + 
  geom_area(stat="bin", alpha=0.6,binwidth = 1)+ 
  scale_fill_manual("K = 5", values=colors)+
  clean_theme+
  labs(x = "Discriminant function 1", y = "Density")+
  theme(legend.position=c(0.9,0.9))#+
#  scale_x_continuous(limits = c(-8,4), breaks= c(-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4))

DAPC_K5


DAPC_K5_scatter <- scatter(dapc.pop5, col=colors)

#ggsave(plot = DAPC_K5, filename = "denovo_DAPC_K5_sco_whole.png",dpi = 300, limitsize = FALSE, width = 12, height = 8)



########################################################################################
#K=6
grp6 <- find.clusters(coral.genlight, stat = "BIC", 
                      n.iter = 100000, n.pca = 120, 
                      choose.n.clust = 6,
                      n.start = 1000, set.seed(2))

#Choose the number of clusters (>=2 - assess BIC): 
#6

pop6 <- grp6$grp
max_PCAs6 <- as.integer(length(pop6) / 3) # as <= N/3 advised
dapc6 <- dapc(coral.genlight, pop6, n.pca = max_PCAs6, n.da = 10)
optimum_score6 <- optim.a.score(dapc6)              #21

dapc.pop6 <- dapc(coral.genlight, pop6, n.pca = optimum_score6$best, n.da = 10)

#Extract coordinates and add metadata
dapc_coord6 <- as.data.frame(dapc.pop6$ind.coord)
dapc_assig6 <- as.data.frame(dapc.pop6$assign)                    #Posterior assignment
dapc_coord6$indv <- rownames(dapc_coord6)
dapc_coord6 <- bind_cols(dapc_coord6, dapc_assig6)
colnames(dapc_coord6)[7] <- "assignment"
dapc_coord6$pop <- with(popfile,pop[match(dapc_coord6$indv,indv)])

#Rename assignment
dapc_coord6$assignment <- as.character(dapc_coord6$assignment)
dapc_coord6$assignment[dapc_coord6$assignment == "1"] <- "Cluster 1"
dapc_coord6$assignment[dapc_coord6$assignment == "2"] <- "Cluster 2"
dapc_coord6$assignment[dapc_coord6$assignment == "3"] <- "Cluster 3"
dapc_coord6$assignment[dapc_coord6$assignment == "4"] <- "Cluster 4"
dapc_coord6$assignment[dapc_coord6$assignment == "5"] <- "Cluster 5"
dapc_coord6$assignment[dapc_coord6$assignment == "6"] <- "Cluster 6"

#Setting site order 
dapc_coord6$pop <- as.factor(dapc_coord6$pop)
levels(dapc_coord6$pop)
dapc_coord6$pop <- factor(dapc_coord6$pop, levels = c("KAL","HUL","EST","SNA","SEA","MAR","DIR"))

#Graph 
DAPC_K6 <- ggplot() + geom_point(data = dapc_coord6, 
                                 aes(x=LD1, y=LD2, color=as.character(assignment),
                                     alpha =0.99),size = 3.5) + 
  scale_color_manual("K = 6", values=colors)+
  #  scale_fill_manual("pop", values=def_shapes)+
  clean_theme+
  labs(x = "LD1", y = "LD2")+#theme(legend.position=c(0.1,0.8))+
  ggtitle("K = 6")

DAPC_K6


DAPC_K6_scatter <- scatter(dapc.pop6, col=colors)

#ggsave(plot = DAPC_K6, filename = "denovo_DAPC_K6_sco_whole.png",dpi = 300, limitsize = FALSE, width = 12, height = 8)



########################################################################################
#K=7
grp7 <- find.clusters(coral.genlight, stat = "BIC", 
                      n.iter = 100000, n.pca = 120, 
                      choose.n.clust = 7,
                      n.start = 1000, set.seed(2))

#Choose the number of clusters (>=2 - assess BIC): 
#7

pop7 <- grp7$grp
max_PCAs7 <- as.integer(length(pop7) / 3) # as <= N/3 advised
dapc7 <- dapc(coral.genlight, pop7, n.pca = max_PCAs7, n.da = 10)
optimum_score7 <- optim.a.score(dapc7)              #18

dapc.pop7 <- dapc(coral.genlight, pop7, n.pca = optimum_score7$best, n.da = 10)

#Extract coordinates and add metadata
dapc_coord7 <- as.data.frame(dapc.pop7$ind.coord)
dapc_assig7 <- as.data.frame(dapc.pop7$assign)                    #Posterior assignment
dapc_coord7$indv <- rownames(dapc_coord7)
dapc_coord7 <- bind_cols(dapc_coord7, dapc_assig7)
colnames(dapc_coord7)[8] <- "assignment"
dapc_coord7$pop <- with(popfile,pop[match(dapc_coord7$indv,indv)])

#Rename assignment
dapc_coord7$assignment <- as.character(dapc_coord7$assignment)
dapc_coord7$assignment[dapc_coord7$assignment == "1"] <- "Cluster 1"
dapc_coord7$assignment[dapc_coord7$assignment == "2"] <- "Cluster 2"
dapc_coord7$assignment[dapc_coord7$assignment == "3"] <- "Cluster 3"
dapc_coord7$assignment[dapc_coord7$assignment == "4"] <- "Cluster 4"
dapc_coord7$assignment[dapc_coord7$assignment == "5"] <- "Cluster 5"
dapc_coord7$assignment[dapc_coord7$assignment == "6"] <- "Cluster 6"
dapc_coord7$assignment[dapc_coord7$assignment == "7"] <- "Cluster 7"


#Setting site order 
dapc_coord7$pop <- as.factor(dapc_coord7$pop)
levels(dapc_coord7$pop)
dapc_coord7$pop <- factor(dapc_coord7$pop, levels = c("KAL","HUL","EST","SNA","SEA","MAR","DIR"))

#Graph 
DAPC_K7 <- ggplot() + geom_point(data = dapc_coord7, 
                                 aes(x=LD1, y=LD2, color=as.character(assignment),
                                     alpha =0.99),size = 3.5) + 
  scale_color_manual("K = 7", values=colors)+
  #  scale_fill_manual("pop", values=def_shapes)+
  clean_theme+
  labs(x = "LD1", y = "LD2")+#theme(legend.position=c(0.1,0.8))+
  ggtitle("K = 7")

DAPC_K7


DAPC_K7_scatter <- scatter(dapc.pop7, col=colors)

ggsave(plot = DAPC_K7, filename = "denovo_DAPC_K7_sco_whole.png",dpi = 300, limitsize = FALSE, width = 12, height = 8)
