library(testthat)
library(WGPAN)
library('dplyr')

test_check("WGPAN")


ID <- c("Amanita","Polyporus","Coprinus", "Daldinia", "Entoloma", "Flaviporus")
myData<- matrix(1:36, nrow = 6, ncol = 6,  dimnames = list(ID, ID))
diag(myData)<- 0 #set diagonal to 0
clades = c("Agaricales", "Polyporales", "Agaricales", "Xylariales", "Agaricales", "Polyporales")
cladeList = data.frame(ID, clades)
x<- myData
idist = na.omit(reshape::melt(x))
colnames(idist) <- c("a", "b", "dist")
allSelfDists = data.frame()
clade_loop <- unique(cladeList$clades)
ID_names <- cladeList$ID
cladeAssignments <-cladeList$clades
WithinDistanceTable <- withinGroupDistances(distanceMatrix = x, cladeAssignments = cladeAssignments, ID_names = ID_names)
boxplot(distance ~ clade, data = WithinDistanceTable, las=3 , main = "within-clade distances 1-11-17", ylab = "SNPs/MB" )
withinSDtable<- WithinDistanceTable %>% group_by(clade)%>%summarize(sd = sd(distance), mean = mean(distance))
bstrees<- generateBootstrapTrees(distanceMatrix = x, numberOfTrees = 100, withinSDtable= withinSDtable, cladeAssignments = cladeAssignments, ID_names = ID_names)
