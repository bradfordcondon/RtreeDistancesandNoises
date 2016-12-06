#All R functions related to noising tree.
require(ape)
require(phytools)
require(distory)

##
#written 7-24-16
#goal- to take a distance matrix, SD, and number of trees, 
#and return the bootstrap trees for this SD via normal noising
##
bootStrapWithNoiseFromNorm <- function(originalMatrix, standardDev, 
                                       numTrees = 1000, rootTree = NULL) {
  print("warning: rooting is not currently supported")
bstrees = vector("list", 0)  #create blank object to put noised trees 
class(bstrees) <- "multiPhylo"
loop = c(1:numTrees)  #number of trees and intervals to make
for (i in loop) {
  x <- originalMatrix #reset matrix
  rowsneeded = row.names(originalMatrix)#reset taxa list
  for (j in rowsneeded) {
    rowsneeded=  rowsneeded[rowsneeded != j]    #remove self from list
    for (k in rowsneeded) {
      oldval = x[j,k]
      newval = oldval + rnorm(1, mean= 0, sd = standardDev)
      x[j, k ] <-newval
      x[k,j] <- newval
    }
  }
  x[x <0]< - 0  #set negative to zero
  treex<-  nj(as.dist(x))
  treex$edge.length[treex$edge.length <0] <- abs(treex$edge.length[treex$edge.length<0])
  if (!is.null(rootTree)){
  treex= multi2di(treex)
  }
  bstrees[[i]] <- treex 
  # message = paste("finished with tree number" , i)
  # print(message)
}
return(bstrees)
}
####
####

##
#written 7-24-16
#Goal: To take a range of SDs and a distance matrix as an input, noisea given number of rep trees, 
#and to return the bootstrap values for each node in list form as an output
##
bootStrapWithNoiseOverRange <- function(originalMatrix,rangeMin = 1, rangeMax = 1000, 
                                       numTrees = 1000, rootTree = NULL) {
  sdRange = c(rangeMin:rangeMax) #range of SD to check
  loop = c(1:numTrees)  #number of trees and intervals to make
  bootstrapList = vector("list", 0)#blank list of bootstrap values
  mainTree = nj(as.dist(originalMatrix))
  mainTree$edge.length[mainTree$edge.length <0] <- abs(mainTree$edge.length[mainTree$edge.length<0])
  for (a in sdRange) {
    bstrees = vector("list", 0)  #create blank object to put noised trees 
    class(bstrees) <- "multiPhylo"
    for (i in loop) {
      x <- originalMatrix #reset matrix
      for (j in rowsneeded) {
        others =  rowsneeded[rowsneeded != j]    #remove self from list
        for (k in others) {
          oldval = x[j,k]
          newval = oldval + rnorm(1, mean= 0, sd = a)
          x[j, k ] <-newval
          x[k,j] <- newval
        }
      }
    x[x <0]< - 0  #set negative to zero
    treex<-  nj(as.dist(x))
    treex$edge.length[treex$edge.length <0] <- abs(treex$edge.length[treex$edge.length<0])
      if (!is.null(rootTree)){
        treex= root(treex, rootTree)
      }
    bstrees[[i]] <- treex 
    }
    #perform bootstrap analysis
      if (!is.null(rootTree)){ 
    clade_support <- prop.clades(mainTree, bstrees, rooted = TRUE)  #generate BS values based on trees
      } else{ 
      clade_support <- prop.clades(mainTree, bstrees, rooted = FALSE)  #generate BS values based on trees
      }
    clade_support = clade_support/numTrees*100  #convert to percent
    clade_support[is.na(clade_support)]<- 0  #set NA to zero
    bootstrapList[[a]]<-clade_support
     message = paste("finished with BS trees from SD " , a)
    print(message)
  }
  return(bootstrapList)
}


##readAndCombineTrees
#imports directory of tree files, trims and roots as necessary, returns multiphylo object
#Last modified 7-21-16

readAndCombineTrees<- function(pathToTrees, taxaToExclude = NULL, taxonToRoot = NULL){
file_list <-  list.files(path = pathToTrees) # List files to import.
allMyTrees = vector("list", 2)
class(allMyTrees) <- "multiPhylo" #create blank multiPhylo Object
loop = as.list(c(1:length(file_list))) #make a loop as long as the file list
for (i in (loop)) {
  filename = file_list[i]
  treeFile = paste(pathToTrees,"/",filename, sep="")
  itree <- read.tree(treeFile, tree.names = filename)
  itree$tip.label <-   gsub("NEW_GENOMES/FASTA/", "", perl=TRUE, itree$tip) #trim node labels
  itree$tip.label <-   gsub("_.*", "", perl=TRUE, itree$tip)
  if (is.null(taxaToExclude) == FALSE) {
    for (j in taxaToExclude){ #drop tips matching taxaToExclude
      itree= drop.tip(itree, c(j, trim.internal = TRUE )) #delete tips as needed
    }
  }
  if (is.null(taxonToRoot) == FALSE) {
    itree= root(itree, taxonToRoot, resolve.root=TRUE)
  }
  #convert to absolute value
  itree$edge.length[itree$edge.length <0] <- abs(itree$edge.length[itree$edge.length<0]) 
  allMyTrees[[i]] <- itree
}
names(allMyTrees) = file_list
return(allMyTrees)
}

# #####determining distances
#Function: treesToCladeCompare.  
#input: Multiphylo object, list of clade assignments
#output: A list of two dataframes.  The first documenting each clade's average in and out distances for each tree.  The second is a simplified, average in/out score for each tree.
#Last modified- 7-15-16 
treesToCladeCompare<- function(multiPhyloTrees, cladeListFile = "list.txt", referenceTree = NULL, refCompareMethod = "edgeset") {
clade_list <- read.table(file = cladeListFile, header=FALSE)
colnames(clade_list)<- c("ID", "host", "cladenum") #assumes 3 column clade list.  clade num isnt necesary
clade_list$ID <- gsub(".*._v_", "", perl=TRUE, clade_list$ID)#trim the ID names, assuming names are FARMAN style
clade_list$ID <- gsub("_.*", "", perl=TRUE, clade_list$ID)
clade_loop <- unique(clade_list$host)
#remove taxa in the list that aren't in the tree
 taxaToKeep = multiPhyloTrees[[1]]$tip.label
clade_list = clade_list[which(clade_list$ID %in% taxaToKeep),]
allTreesData<- data.frame()
simpleTreesData<- data.frame()
if (is.null(referenceTree) == FALSE) {
  refDistTracker = data.frame(stringsAsFactors=FALSE)
  compare = vector("list", 2)
  class(compare) <- "multiPhylo" #create blank multiPhylo Object
  compare[[1]]<- referenceTree #store ref tree in slot1 of a multiphylo object.
}  
loop = as.list(c(1:length(multiPhyloTrees))) #Set up tree loop
 for (b in loop) {
   itree = multiPhyloTrees[[b]]  #retrieve distance matrix for each tree
   idist= cophenetic.phylo(itree)
   in_tracker = character() #reset trackers for this individual tree
   out_tracker = character()
   clade_means = data.frame()
   if (is.null(referenceTree) == FALSE) {
       compare[[2]] <- itree #put this tree in our comparison slot
       distance = c(distance= as.character(dist.multiPhylo(compare, method= refCompareMethod)), treename = as.character(names(multiPhyloTrees[b]))) 
       refDistTracker=  rbind(refDistTracker, distance, stringsAsFactors=FALSE)
     }
     for (i in clade_loop){
       clade_in = character()
       clade_out = character()
       #retrieve all genomes matching i as character
       itaxa <- clade_list[which(clade_list$host == i),1]
       nottaxa <- clade_list[which(clade_list$host != i),1]
       for (a in itaxa) {
         in_distancesa= idist[a,itaxa]
         out_distancesa = idist[a,nottaxa]
         in_tracker = c(in_tracker, in_distancesa)
         out_tracker = c(out_tracker, out_distancesa)
         clade_in = c(clade_in, in_distancesa)
         clade_out = c(clade_out, out_distancesa)
       }
       #now, store the average in and out distances for this particular clade
       clade_data= data.frame("avg_in" = mean(as.numeric(clade_in)), "avg_out"= mean(as.numeric(clade_out)),  "clade"=i, check.names=FALSE)
       clade_means =  rbind(clade_means, clade_data)
     }
     #Finally, assign the mean in, mean out, min out, and max in distances for the tree
     toadd = data.frame(clade_means, tree=names(multiPhyloTrees[b]))
     simpleadd = data.frame("in_avg" = mean(clade_means$avg_in, na.rm=TRUE), "out_avg" = mean(clade_means$avg_out, na.rm=TRUE), "tree" = names(multiPhyloTrees[b]))
     allTreesData = rbind(allTreesData, toadd)
     simpleTreesData = rbind(simpleTreesData, simpleadd)
 }#finish looping through trees

simpleTreesData= cbind(simpleTreesData, inPRank = rank(simpleTreesData$in_avg)/length(simpleTreesData$in_avg)*100, outPRank= rank(simpleTreesData$out_avg)/length(simpleTreesData$out_avg)*100)
if (is.null(referenceTree) == FALSE) {
    colnames(refDistTracker) = c("distance", "treeName")
    myOutput <- list("cladeByCladeData"=allTreesData, "treeByTreeSummary"= simpleTreesData, "distanceToReferenceTree" = refDistTracker)
  }
  else {
myOutput <- list("cladeByCladeData"=allTreesData, "treeByTreeSummary"= simpleTreesData, "distanceToReferenceTree" = NULL)
  }
  return(myOutput)
}

###
#Trim tree A based on tips in tree B.
#Intended purpose is for cutting back the reference tree to remove exlcuded taxa
#input: Two tree files- a tree to trim, and a tree to serve as teh template for trimming tips.
###
trimTreeToTree<- function(treeToTrim, limitingTree) {
  unTrimTips = treeToTrim$tip
  limTips = limitingTree$tip
  tipsToDrop <- unTrimTips[!unTrimTips %in% limTips]
  for (i in tipsToDrop){
    treeToTrim =  drop.tip(treeToTrim, c(i, trim.internal = TRUE )) #delete tips as needed
  }
  finalTree <-treeToTrim
  return(finalTree) 
}
###Trim tree based on list of tips


trimTreeFromList<- function(treeToTrim, listToDrop) {
  unTrimTips = treeToTrim$tip
  limTips <- listToDrop
  tipsToDrop <- unTrimTips[!unTrimTips %in% limTips]
  for (i in tipsToDrop){
    treeToTrim =  drop.tip(treeToTrim, c(i, trim.internal = TRUE )) #delete tips as needed
  }
  finalTree <-treeToTrim
  return(finalTree) 
}
#####Generate report (histogram, quantile cutoffs) for in/out clade distances
#Not sure if this is working as intended
#7-28-16

computeQuantiles <- function (treeDistanceData, cutOffNum, greaterThan = TRUE){
    totalNum = length(dataToConvert[,1])
    if (is.null(greaterThan) == FALSE) {
      x = ceiling(percentCutoff*totalNum/100)
      sort(dataToConvert[1:x,1])    
    }
    if (is.null(greaterThan) == TRUE) {
      x = floor(percentCutoff*totalNum/100)
    }
    return(x)
}



#calculate the average distance between groups given a distance matrix and a clade list
#12-5-16

withinGroupDistancesTwo <- function(distanceMatrix, clade_list){
  require('dplyr')
  #remove items from clade_list not in tree
  clade_list<- filter(clade_list, Idcorrected %in% row.names(distanceMatrix))

  idist = na.omit(melt(distanceMatrix))
  colnames(idist) <- c("a", "b", "dist")
  
  allSelfDists = data.frame()
  clade_loop <- unique(clade_list$hostGenus)
  #loop through all clade self comparisons, add to DF
  for (i in clade_loop){
    #retrieve all genomes matching i as character
    itaxa <- as.character(clade_list[which(clade_list$hostGenus == i),2])
    if (length(itaxa > 1)){
      for (a in itaxa) {     
        intaxa =  itaxa[itaxa != a]    #remove variable from list
        #retrieve values matching a
        ataxa= idist[idist$a %in% a,]
        #calculate in-distances
        in_distancesa= ataxa[ataxa$b %in% intaxa, 3]
        if(length(in_distancesa > 0)){
        toAdd = cbind(as.numeric(in_distancesa), as.character(i), a)
        }
        
        allSelfDists = rbind(allSelfDists, toAdd)
      }
      }
    }
  colnames(allSelfDists) <- c("distance", "clade")
  allSelfDists$distance<- as.numeric(as.character(allSelfDists$distance))
  return(allSelfDists)
}







#calculate the average distance between groups given a distance matrix and a clade list
#8-15-16

withinGroupDistances <- function(distanceMatrix, clade_list){
idist = na.omit(melt(distanceMatrix))
colnames(idist) <- c("a", "b", "dist")
allSelfDists = data.frame()
clade_loop <- unique(clade_list$hostGenus)
#loop through all clade self comparisons, add to DF
for (i in clade_loop){
  #retrieve all genomes matching i as character
  itaxa <- as.character(clade_list[which(clade_list$hostGenus == i),2])
  if (length(itaxa < 2)){ #Skip this clade if theres 1 or 0 taxa.
  }else{
  for (a in itaxa) {     
    intaxa =  itaxa[itaxa != a]    #remove variable from list
    #retrieve values matching a
    ataxa= idist[idist$a %in% a, ]
    #calculate in-distances
    in_distancesa= ataxa[ataxa$b %in% intaxa, 3] 
    toAdd = cbind(as.numeric(in_distancesa), as.character(i))
    allSelfDists = rbind(allSelfDists, toAdd)
  }
  }
}
colnames(allSelfDists) <- c("distance", "clade")
allSelfDists$distance<- as.numeric(as.character(allSelfDists$distance))
return(allSelfDists)

}


###
#12-2-16
#Function to write 1,000 bootstrap trees, where within-clade distances (only!) are noised by the normal distribution with SD= that clade's within SD.
#returns bootstrap trees
###

bootstrapFromNormInCladeOnly<- function(startingTree, clade_list, withinSDtable, ntrees = 1000  ){
master<-startingTree
bstrees = vector("list", 0)  #create blank object to put noised trees in
class(bstrees) <- "multiPhylo"
loop = c(1:ntrees)  #number of trees and intervals to make
for (i in loop) {  
  x <- master #reset matrix
  for (a in as.character(clade_list$Idcorrected)) {
    #lookup this clade
    thisClade <-  as.character(unlist(clade_list %>% filter(Idcorrected == as.character(a)) %>% select (hostGenus )))
    #determine SD to noise from
    thisSD <- withinSDtable %>% filter(clade == thisClade) %>% select (sd)
    #get all taxa in this clade excluding itself
    matchingTaxa <- as.list(clade_list %>% filter(hostGenus == thisClade) %>% select (Idcorrected) %>% filter(Idcorrected != as.character(a)))
    for (b in matchingTaxa$Idcorrected){  
      sublista <- idist[which(idist$a == a  ) ,]
      sublist <- sublista[which((sublista$b) == b ),]  #retrieve the distance between these two strains
      thisDist = as.numeric(sublist[,3]) 
      newVal = (thisDist + rnorm(1, mean = 0, sd = as.numeric(thisSD)))   #randomly sample 1 value from rnorm of this SD
      x[b,a] <-  newVal
      x[a,b] <- newVal
    }
  }   
  x <- abs(x)  #Convert negative to positive
  treex<-  nj(as.dist(x))
  treex$edge.length[treex$edge.length <0] <- abs(treex$edge.length[treex$edge.length<0]) 
  bstrees[[i]] <- treex 
  message = paste("finished with tree number" , i)
  print(message)
}
return(bstrees)
}
