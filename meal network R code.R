################################################################

# PREG HEI T1 BREAKFAST NETWORK

################################################################

# import data
mydata <- read.sas7bdat("C:/Users/schwedhelmramc2/Documents/PEAS data_M1/HEI/PREG/SAS/preg_hei_t1_occ1.sas7bdat")


# skeptic transformation for non-normally distributed data
nphuge<-huge.npn(mydata, npn.func = "skeptic")
# if data is normally distributed, use correlation or covariance matrix
cor_mydata <- cor(mydata)


# see what happens when we don't regularize the correlation matrix
nphugeL = huge(nphuge,lambda=0.0001 ,method = "glasso")
nphugeL
plot(nphugeL)


# see how regularization/penalty parameter works
nphugeL = huge(nphuge,lambda=NULL ,method = "glasso")
nphugeL
plot(nphugeL)


# choosing optimal lambda from cross-validated glasso
res<-screen_cv.glasso(nphuge, fold = 5, use.package = "glasso",verbose = TRUE)


# visualizing with optimal lambda from cross-validation
nphugeL = huge(nphuge,lambda=0.2285652   ,method = "glasso")
nphugeL
plot (nphugeL)


# generating graph
res.lasso.m<-glasso(nphuge,rho=0.2285652) #put here optimal lambda obtained from previous line
AM <- res.lasso.m$wi != 0
diag(AM) <- F
g.lasso.m.ug <- as(AM, "graphNEL")
nodes(g.lasso.m.ug)<-names(mydata)  
glmat.m<-as(g.lasso.m.ug,"matrix")
iglmat.m<-as(glmat.m,"igraph")

# preparing edge weights
# this is our regularized correlation matrix
corrmat.m<-(glmat.m*nphuge)

corrmat_tri <- corrmat.m
corrmat_tri[lower.tri(corrmat_tri)] <- NA
corrm <- melt(corrmat_tri)
corrm1 <- na.omit(corrm)
corrm2 <- corrm1[!(corrm1$value ==0),]
# order edges in the same order as in igraph
my_weights <- corrm2[order(corrm2$Var1),]

# set (absolute) weights
my_weights2 <- my_weights
my_weights2$abs_value <- abs(my_weights2$value)
# export this to add edge attributes in cytoscape
my_weights2$Negative <- ifelse(my_weights2$value == my_weights2$abs_value, FALSE, TRUE)
write.csv(my_weights2,file="edgeW_PREG_HEI_T1_OCC1.csv")

iglmat.m0 <- set_edge_attr(iglmat.m, "weight", value=my_weights2$abs_value)
E(iglmat.m0)$weight <- my_weights2$abs_value

# plot graph
plot(iglmat.m, edge.width=E(iglmat.m0)$weight*4, vertex.size=4,vertex.label.dist=0.3, vertex.color="green", 
     edge.color="black")


# export igraph to sif (to open with cytoscape)
igraphToSif(iglmat.m, "PREG_HEI_T1_occ1.sif", "and")


# LOUVAIN - obtaining communities with edge weights. 
louv <- cluster_louvain(iglmat.m, weights = E(iglmat.m0)$weight)
plot(louv,iglmat.m,vertex.size=2)


# checking and comparing properties of communities
modularity(louv)
length(louv)
# vs
louvain(abs(corrmat.m))


# printing each community
print(louv[1])
print(louv[2])
print(louv[3])
print(louv[4])
print(louv[5])
print(louv[6])
print(louv[7])
print(louv[8])
print(louv[9])
print(louv[10])
print(louv[11])    



# participation coeff
participation(corrmat.m, comm =("louvain"))
# clustering coeff
clustcoeff(corrmat.m)


## LINKCOMM COMMUNITIES ##

#you get linked communities 
lc <- getLinkCommunities(my_weights) 
# you get comminuities
print(getNodesIn(lc, clusterids = 1, type = "names"))
print(getNodesIn(lc, clusterids = 2, type = "names"))
print(getNodesIn(lc, clusterids = 3, type = "names")) 
print(getNodesIn(lc, clusterids = 4, type = "names")) 
print(getNodesIn(lc, clusterids = 5, type = "names"))

# get community centrality
getCommunityCentrality(lc)


################################################################

# PREG HEI T3 BREAKFAST NETWORK

################################################################

# import data
mydata <- read.sas7bdat("C:/Users/schwedhelmramc2/Documents/PEAS data_M1/HEI/PREG/SAS/preg_hei_t3_occ1.sas7bdat")


# skeptic transformation for non-normally distributed data
nphuge<-huge.npn(mydata, npn.func = "skeptic")

# choosing optimal lambda from cross-validated glasso
res<-screen_cv.glasso(nphuge, fold = 5, plot.it = TRUE, se = TRUE, use.package = "glasso",verbose = TRUE)

# generating graph
res.lasso.m<-glasso(nphuge,rho=0.1811154) #put here optimal lambda obtained from previous line
AM <- res.lasso.m$wi != 0
diag(AM) <- F
g.lasso.m.ug <- as(AM, "graphNEL")
nodes(g.lasso.m.ug)<-names(mydata)  
glmat.m<-as(g.lasso.m.ug,"matrix")
iglmat.m<-as(glmat.m,"igraph")

# preparing edge weights
# this is our regularized correlation matrix
corrmat.m<-(glmat.m*nphuge)
corrmat_tri <- corrmat.m
corrmat_tri[lower.tri(corrmat_tri)] <- NA
corrm <- melt(corrmat_tri)
corrm1 <- na.omit(corrm)
corrm2 <- corrm1[!(corrm1$value ==0),]
# order edges in the same order as in igraph
my_weights <- corrm2[order(corrm2$Var1),]

# set (absolute) weights
my_weights2 <- my_weights
my_weights2$abs_value <- abs(my_weights2$value)
# export this to add edge attributes in cytoscape
my_weights2$Negative <- ifelse(my_weights2$value == my_weights2$abs_value, FALSE, TRUE)
write.csv(my_weights2,file="edgeW_PREG_HEI_T3_OCC1.csv")
iglmat.m0 <- set_edge_attr(iglmat.m, "weight", value=my_weights2$abs_value)
E(iglmat.m0)$weight <- my_weights2$abs_value

# plot graph
plot(iglmat.m, edge.width=E(iglmat.m0)$weight*4, vertex.size=4,vertex.label.dist=0.3, vertex.color="green", 
     edge.color="black")


# export igraph to sif (to open with cytoscape)
igraphToSif(iglmat.m, "PREG_HEI_T13_occ1.sif", "and")


# LOUVAIN - obtaining communities with edge weights. 
louv <- cluster_louvain(iglmat.m, weights = E(iglmat.m0)$weight)
plot(louv,iglmat.m,vertex.size=2)


# checking and comparing properties of communities
modularity(louv)
length(louv)
# vs
louvain(abs(corrmat.m))


# printing each community
print(louv[1])
print(louv[2])
print(louv[3])
print(louv[4])
print(louv[5])
print(louv[6])




# participation coeff
# participation(corrmat.m, comm =("louvain"))
# clustering coeff
# clustcoeff(corrmat.m)