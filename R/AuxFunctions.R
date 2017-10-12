library(ape)
library(phytools)
library(distory)

## readAndCombineTrees imports directory of tree files, trims and roots as necessary, returns
## multiphylo object
#' Collates a folder of tree files into a single multiPhylo tree object.
#'
#' @param pathToTrees The path to your trees folder.
#' @param taxaToExclude Any taxa you would like excluded from the tree.  Useful to enforce conformity.
#' @param taxonToRoot Taxon to root all trees at.
#' @examples
#' readAndCombineTrees(c(/usr/local/files/trees), c(unicorns), c(hamsters) )
#'
#' @return Output is a multiPhylo object.
#' 
#'

readAndCombineTrees <- function(pathToTrees, taxaToExclude = NULL, taxonToRoot = NULL) {
    file_list <- list.files(path = pathToTrees)  # List files to import.
    allMyTrees = vector("list", 2)
    class(allMyTrees) <- "multiPhylo"  #create blank multiPhylo Object
    loop = as.list(c(1:length(file_list)))  #make a loop as long as the file list
    for (i in (loop)) {
        filename = file_list[i]
        treeFile = paste(pathToTrees, "/", filename, sep = "")
        itree <- read.tree(treeFile, tree.names = filename)
        itree$tip.label <- gsub("NEW_GENOMES/FASTA/", "", perl = TRUE, itree$tip)  #trim node labels
        itree$tip.label <- gsub("_.*", "", perl = TRUE, itree$tip)
        if (is.null(taxaToExclude) == FALSE) {
            for (j in taxaToExclude) {
                # drop tips matching taxaToExclude
                itree = drop.tip(itree, c(j, trim.internal = TRUE))  #delete tips as needed
            }
        }
        if (is.null(taxonToRoot) == FALSE) {
            itree = root(itree, taxonToRoot, resolve.root = TRUE)
        }
        # convert to absolute value
        itree$edge.length[itree$edge.length < 0] <- abs(itree$edge.length[itree$edge.length < 0])
        allMyTrees[[i]] <- itree
    }
    names(allMyTrees) = file_list
    return(allMyTrees)
}

# #####determining distances Function: treesToCladeCompare.  input: Multiphylo object, list of
# clade assignments output: A list of two (or three) dataframes.  The first documenting each
# clade's average in and out distances for each tree.  The second is a simplified, average
# in/out score for each tree.
#' @param multiPhyloTrees Your tree object from readAndCombineTrees.
#' @param cladeListFile A list of clades that your taxa belong to.
#' @param referenceTree A reference Tree.
#' @param refCompareMethod A comparison method for the reference tree.
#'
#' 
#'

treesToCladeCompare <- function(multiPhyloTrees, cladeListFile = "list.txt", referenceTree = NULL, 
    refCompareMethod = "edgeset") {
    clade_list <- read.table(file = cladeListFile, header = FALSE)
    colnames(clade_list) <- c("ID", "host", "cladenum")  #assumes 3 column clade list.  clade num isnt necesary
    clade_list$ID <- gsub(".*._v_", "", perl = TRUE, clade_list$ID)  #trim the ID names, assuming names are FARMAN style
    clade_list$ID <- gsub("_.*", "", perl = TRUE, clade_list$ID)
    clade_loop <- unique(clade_list$host)
    # remove taxa in the list that aren't in the tree
    taxaToKeep = multiPhyloTrees[[1]]$tip.label
    clade_list = clade_list[which(clade_list$ID %in% taxaToKeep), ]
    allTreesData <- data.frame()
    simpleTreesData <- data.frame()
    if (is.null(referenceTree) == FALSE) {
        refDistTracker = data.frame(stringsAsFactors = FALSE)
        compare = vector("list", 2)
        class(compare) <- "multiPhylo"  #create blank multiPhylo Object
        compare[[1]] <- referenceTree  #store ref tree in slot1 of a multiphylo object.
    }
    loop = as.list(c(1:length(multiPhyloTrees)))  #Set up tree loop
    for (b in loop) {
        itree = multiPhyloTrees[[b]]  #retrieve distance matrix for each tree
        idist = cophenetic.phylo(itree)
        in_tracker = character()  #reset trackers for this individual tree
        out_tracker = character()
        clade_means = data.frame()
        if (is.null(referenceTree) == FALSE) {
            compare[[2]] <- itree  #put this tree in our comparison slot
            distance = c(distance = as.character(dist.multiPhylo(compare, method = refCompareMethod)), 
                treename = as.character(names(multiPhyloTrees[b])))
            refDistTracker = rbind(refDistTracker, distance, stringsAsFactors = FALSE)
        }
        for (i in clade_loop) {
            clade_in = character()
            clade_out = character()
            # retrieve all genomes matching i as character
            itaxa <- clade_list[which(clade_list$host == i), 1]
            nottaxa <- clade_list[which(clade_list$host != i), 1]
            for (a in itaxa) {
                in_distancesa = idist[a, itaxa]
                out_distancesa = idist[a, nottaxa]
                in_tracker = c(in_tracker, in_distancesa)
                out_tracker = c(out_tracker, out_distancesa)
                clade_in = c(clade_in, in_distancesa)
                clade_out = c(clade_out, out_distancesa)
            }
            # now, store the average in and out distances for this particular clade
            clade_data = data.frame(avg_in = mean(as.numeric(clade_in)), avg_out = mean(as.numeric(clade_out)), 
                clade = i, check.names = FALSE)
            clade_means = rbind(clade_means, clade_data)
        }
        # Finally, assign the mean in, mean out, min out, and max in distances for the tree
        if (nrow(clade_means) > 0) {
            toadd = data.frame(clade_means, tree = names(multiPhyloTrees[b]))
        }
        simpleadd = data.frame(in_avg = mean(clade_means$avg_in, na.rm = TRUE), out_avg = mean(clade_means$avg_out, 
            na.rm = TRUE), tree = names(multiPhyloTrees[b]))
        allTreesData = rbind(allTreesData, toadd)
        simpleTreesData = rbind(simpleTreesData, simpleadd)
        
    }  #finish looping through trees
    
    simpleTreesData = cbind(simpleTreesData, inPRank = rank(simpleTreesData$in_avg)/length(simpleTreesData$in_avg) * 
        100, outPRank = rank(simpleTreesData$out_avg)/length(simpleTreesData$out_avg) * 100)
    
    
    if (is.null(referenceTree) == FALSE) {
        colnames(refDistTracker) = c("distance", "treeName")
        myOutput <- list(cladeByCladeData = allTreesData, treeByTreeSummary = simpleTreesData, distanceToReferenceTree = refDistTracker)
    } else {
        myOutput <- list(cladeByCladeData = allTreesData, treeByTreeSummary = simpleTreesData, distanceToReferenceTree = NULL)
    }
    
    
    
    return(myOutput)
}

### Trim tree A based on tips in tree B. Intended purpose is for cutting back the reference tree
### to remove exlcuded taxa
trimTreeToTree <- function(treeToTrim, limitingTree) {
    unTrimTips = treeToTrim$tip
    limTips = limitingTree$tip
    tipsToDrop <- unTrimTips[!unTrimTips %in% limTips]
    for (i in tipsToDrop) {
        treeToTrim = drop.tip(treeToTrim, c(i, trim.internal = TRUE))  #delete tips as needed
    }
    finalTree <- treeToTrim
    return(finalTree)
}
##### Generate report (histogram, quantile cutoffs) for in/out clade distances

computeQuantiles <- function(treeDistanceData, cutOffNum, greaterThan = TRUE) {
    totalNum = length(dataToConvert[, 1])
    
    if (is.null(greaterThan) == FALSE) {
        x = ceiling(percentCutoff * totalNum/100)
        sort(dataToConvert[1:x, 1])
        
    }
    if (is.null(greaterThan) == TRUE) {
        x = floor(percentCutoff * totalNum/100)
    }
    
    # meltedDataframe = return (meltedDataFrame)
}

## AnalyzeTips goal: determine which taxa are commonly missing from tips determine distribution
## of tip #'s so trees can be filtered.
analyzeTips <- function(multiPhyloObject) {
    tipSumCounter <- data_frame()
    tipLabelCounter <- data_frame(tip = character(0), count = numeric(0))
    loop <- c(1:length(multiPhyloObject))
    outputObject <- NULL
    for (i in loop) {
        thisTree = multiPhyloObject[[i]]
        thisTreeTips <- thisTree$tip.label
        toAdd <- cbind(names(multiPhyloObject)[i], as.numeric(length(thisTreeTips)))
        tipSumCounter <- rbind(tipSumCounter, toAdd)
        for (thisTip in thisTreeTips) {
            if (thisTip %in% tipLabelCounter[, 1]) {
                # if the column does exist, add +1.
                tipLabelCounter[tipLabelCounter$tip == thisTip, 2] <- as.numeric(tipLabelCounter[tipLabelCounter$tip == 
                  thisTip, 2]) + 1
                
            } else {
                # if it doesnt, create it, set it to 1
                toAddtwo <- cbind(tip = thisTip, count = 1, stringsAsFactors = FALSE)
                tipLabelCounter <- rbind(tipLabelCounter, toAddtwo, stringsAsFactors = FALSE)
            }
        }
    }
    colnames(tipSumCounter) <- c("treeName", "tipNumber")
    outputObject <- list(tipSumCounter, tipLabelCounter)
    return(outputObject)
}



# Function: treesToCladeCompareSpecific.  input: Multiphylo object, list of clade assignments,
# list of clades to focus on output: A list of two (or three) dataframes.  The first documenting
# each clade's average in and out distances for each tree.  The second is a simplified, average
# in/out score for each tree.  too convoluted- abandoned for now.
treesToCladeCompareKey <- function(multiPhyloTrees, cladeListFile = "list.txt", referenceTree = NULL, 
    refCompareMethod = "edgeset", keyClades = NULL) {
    clade_list <- read.table(file = cladeListFile, header = FALSE)
    colnames(clade_list) <- c("ID", "host", "cladenum")  #assumes 3 column clade list.  clade num isnt necesary
    clade_list$ID <- gsub(".*._v_", "", perl = TRUE, clade_list$ID)  #trim the ID names, assuming names are FARMAN style
    clade_list$ID <- gsub("_.*", "", perl = TRUE, clade_list$ID)
    # remove NAs.
    clade_list <- clade_list[!is.na(clade_list$ID), ]
    clade_loop <- unique(clade_list$host)
    # remove taxa in the list that aren't in the tree
    taxaToKeep = multiPhyloTrees[[1]]$tip.label
    clade_list = clade_list[which(clade_list$ID %in% taxaToKeep), ]
    allTreesData <- data.frame()
    simpleTreesData <- data.frame()
    keyCladesOutput <- data.frame()
    if (is.null(referenceTree) == FALSE) {
        refDistTracker = data.frame(stringsAsFactors = FALSE)
        compare = vector("list", 2)
        class(compare) <- "multiPhylo"  #create blank multiPhylo Object
        compare[[1]] <- referenceTree  #store ref tree in slot1 of a multiphylo object.
    }
    loop = as.list(c(1:length(multiPhyloTrees)))  #Set up tree loop
    for (b in loop) {
        itree = multiPhyloTrees[[b]]  #retrieve distance matrix for each tree
        idist = cophenetic.phylo(itree)
        in_tracker = character()  #reset trackers for this individual tree
        out_tracker = character()
        clade_means = data.frame()
        treeKeys = data.frame()
        
        if (is.null(referenceTree) == FALSE) {
            compare[[2]] <- itree  #put this tree in our comparison slot
            distance = c(distance = as.character(dist.multiPhylo(compare, method = refCompareMethod)), 
                treename = as.character(names(multiPhyloTrees[b])))
            refDistTracker = rbind(refDistTracker, distance, stringsAsFactors = FALSE)
        }
        for (i in clade_loop) {
            # retrieve all genomes matching i as character
            itaxa <- clade_list[which(clade_list$host == i), 1]
            nottaxa <- clade_list[which(clade_list$host != i), 1]
            for (a in itaxa) {
                in_distancesa = idist[a, itaxa]
                out_distancesa = idist[a, nottaxa]
                in_tracker = c(in_tracker, in_distancesa)
                out_tracker = c(out_tracker, out_distancesa)
                clade_in = c(clade_in, in_distancesa)
                clade_out = c(clade_out, out_distancesa)
            }
            if (is.null(keyClades == FALSE)) {
                # if running keyclades routine... if this clade is a key clade
                if (i %in% keyClades) {
                  otherKeys <- clade_list[which(clade_list$host != j & clade_list$host %in% keyClades), 
                    1]
                  for (a in itaxa) {
                    inKeyDists = mean(idist[a, itaxa])
                    outKeyDists = mean(idist[a, otherKeys])
                    treeKeys <- rbind(cbind(inKeyDists, outKeyDists))
                  }
                  
                }
            }
            
            # now, store the average in and out distances for this particular clade
            clade_data = data.frame(avg_in = mean(as.numeric(clade_in)), avg_out = mean(as.numeric(clade_out)), 
                clade = i, check.names = FALSE)
            clade_means = rbind(clade_means, clade_data)
        }
        # Finally, assign the mean in, mean out, min out, and max in distances for the tree
        if (nrow(clade_means) > 0) {
            toadd = data.frame(clade_means, tree = names(multiPhyloTrees[b]))
        }
        simpleadd = data.frame(in_avg = mean(clade_means$avg_in, na.rm = TRUE), out_avg = mean(clade_means$avg_out, 
            na.rm = TRUE), tree = names(multiPhyloTrees[b]))
        allTreesData = rbind(allTreesData, toadd)
        simpleTreesData = rbind(simpleTreesData, simpleadd)
        # add this tree's key info
        if (is.null(keyClades == FALSE)) {
            addKey = data.frame(out = mean(treeKeys[, 2]), `in` = mean(treeKeys[, 1]), tree = names(multiPhyloTrees[b]))
            keyCladesOutput = rbind(keyCladesOutput, addKey)
        }
    }  #finish looping through trees
    
    simpleTreesData = cbind(simpleTreesData, inPRank = rank(simpleTreesData$in_avg)/length(simpleTreesData$in_avg) * 
        100, outPRank = rank(simpleTreesData$out_avg)/length(simpleTreesData$out_avg) * 100)
    
    if (is.null(referenceTree) == FALSE) {
        colnames(refDistTracker) = c("distance", "treeName")
        myOutput <- list(cladeByCladeData = allTreesData, treeByTreeSummary = simpleTreesData, distanceToReferenceTree = refDistTracker)
    } else {
        myOutput <- list(cladeByCladeData = allTreesData, treeByTreeSummary = simpleTreesData, distanceToReferenceTree = NULL)
    }
    if (is.null(keyClades) == FALSE) {
        # add key clades info to list
        myOutput[4] <- keyCladesOutput
    }
    return(myOutput)
}



multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel ncol: Number of columns of plots nrow: Number of rows needed, calculated from #
        # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots == 1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
        }
    }
}
