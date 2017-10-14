#' Collates a folder of tree files into a single multiPhylo tree object.
#'
#' @param pathToTrees The path to your trees folder.
#' @param taxaToExclude Any taxa you would like excluded from the tree.  Useful to enforce conformity.
#' @param taxonToRoot Taxon to root all trees at.
#'
#' @return Output is a multiPhylo object.
#'
#' @export
#'

readAndCombineTrees <- function(pathToTrees, taxaToExclude = NULL, taxonToRoot = NULL) {
    fileList <- list.files(path = pathToTrees)  # List files to import.
    allMyTrees = vector("list", 2)
    class(allMyTrees) <- "multiPhylo"  #create blank multiPhylo Object
    loop = as.list(c(1:length(fileList)))  #make a loop as long as the file list
    for (i in (loop)) {
        fileName = fileList[i]
        treeFile = paste(pathToTrees, "/", fileName, sep = "")
        itree <- ape::read.tree(treeFile, tree.names = fileName)
        itree$tip.label <- gsub("NEW_GENOMES/FASTA/", "", perl = TRUE, itree$tip)  #trim node labels
        itree$tip.label <- gsub("_.*", "", perl = TRUE, itree$tip)
        if (is.null(taxaToExclude) == FALSE) {
            for (j in taxaToExclude) {
                # drop tips matching taxaToExclude
                itree = ape::drop.tip(itree, c(j, trim.internal = TRUE))  #delete tips as needed
            }
        }
        if (is.null(taxonToRoot) == FALSE) {
            itree = ape::root(itree, taxonToRoot, resolve.root = TRUE)
        }
        # convert to absolute value
        itree$edge.length[itree$edge.length < 0] <- abs(itree$edge.length[itree$edge.length < 0])
        allMyTrees[[i]] <- itree
    }
    names(allMyTrees) = fileList
    return(allMyTrees)
}

#' Trims the tips of a target tree A to match the tips of source tree B.
#'
#' @param treeToTrim The target tree that you would like to drop tips on.
#' @param limitingTree The source tree that will determine which tips are dropped.
#'
#' @return Output is the trimmed tree of treeToTrim.
#' @export
#'

trimTreeToTree <- function(treeToTrim, limitingTree) {
    unTrimTips = treeToTrim$tip
    limTips = limitingTree$tip
    tipsToDrop <- unTrimTips[!unTrimTips %in% limTips]
    for (i in tipsToDrop) {
        treeToTrim = ape::drop.tip(treeToTrim, c(i, trim.internal = TRUE))  #delete tips as needed
    }
    finalTree <- treeToTrim
    return(finalTree)
}
### Trim tree based on list of tips
#' Drops a set of tips from a tree.
#'
#' @param treeToTrim The target tree that you would like to drop tips on.
#' @param listToDrop A list of tips to drop.
#'
#' @return Output is the trimmed tree of treeToTrim.
#' @export
#'

trimTreeFromList <- function(treeToTrim, listToDrop) {
    unTrimTips = treeToTrim$tip
    limTips <- listToDrop
    tipsToDrop <- unTrimTips[!unTrimTips %in% limTips]
    for (i in tipsToDrop) {
        treeToTrim = ape::drop.tip(treeToTrim, c(i, trim.internal = TRUE))  #delete tips as needed
    }
    finalTree <- treeToTrim
    return(finalTree)
}


# #####determining distances Function: treesToCladeCompare.  input: Multiphylo object, list of clade assignments output: A list of two dataframes.  The first
# documenting each clade's average in and out distances for each tree.  The second is a simplified, average in/out score for each tree.

#' Calculate in and out distances between and within clades for a set of trees.
#' @param multiPhyloTrees The tree object.
#' @param cladeListFile The file.
#' @param referenceTree A reference tree.
#' @param refCompareMethod The method to use for dist.multiPhylo.
#' @export

treesToCladeCompare <- function(multiPhyloTrees, cladeListFile = "list.txt", referenceTree = NULL, refCompareMethod = "edgeset") {
    cladeList <- read.table(file = cladeListFile, header = FALSE)
    colnames(cladeList) <- c("ID", "host", "cladenum")  #assumes 3 column clade list.  clade num isnt necesary
    cladeList$ID <- gsub(".*._v_", "", perl = TRUE, cladeList$ID)  #trim the ID names, assuming names are FARMAN style
    cladeList$ID <- gsub("_.*", "", perl = TRUE, cladeList$ID)
    cladeLoop <- unique(cladeList$host)
    # remove taxa in the list that aren't in the tree
    taxaToKeep = multiPhyloTrees[[1]]$tip.label
    cladeList = cladeList[which(cladeList$ID %in% taxaToKeep), ]
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
        inTracker = character()  #reset trackers for this individual tree
        outTracker = character()
        cladeMeans = data.frame()
        if (is.null(referenceTree) == FALSE) {
            compare[[2]] <- itree  #put this tree in our comparison slot
            distance = c(distance = as.character(dist.multiPhylo(compare, method = refCompareMethod)), treename = as.character(names(multiPhyloTrees[b])))
            refDistTracker = rbind(refDistTracker, distance, stringsAsFactors = FALSE)
        }
        for (i in cladeLoop) {
            cladeIn = character()
            cladeOut = character()
            # retrieve all genomes matching i as character
            itaxa <- cladeList[which(cladeList$host == i), 1]
            nottaxa <- cladeList[which(cladeList$host != i), 1]
            for (a in itaxa) {
                inDistancesA = idist[a, itaxa]
                outDistancesA = idist[a, nottaxa]
                inTracker = c(inTracker, inDistancesA)
                outTracker = c(outTracker, outDistancesA)
                cladeIn = c(cladeIn, inDistancesA)
                cladeOut = c(cladeOut, outDistancesA)
            }
            # now, store the average in and out distances for this particular clade
            cladeData = data.frame(averageIn = mean(as.numeric(cladeIn)), averageOut = mean(as.numeric(cladeOut)), clade = i, check.names = FALSE)
            cladeMeans = rbind(cladeMeans, cladeData)
        }
        # Finally, assign the mean in, mean out, min out, and max in distances for the tree
        toadd = data.frame(cladeMeans, tree = names(multiPhyloTrees[b]))
        simpleadd = data.frame(inAverage = mean(cladeMeans$averageIn, na.rm = TRUE), outAverage = mean(cladeMeans$averageOut, na.rm = TRUE), tree = names(multiPhyloTrees[b]))
        allTreesData = rbind(allTreesData, toadd)
        simpleTreesData = rbind(simpleTreesData, simpleadd)
    }  #finish looping through trees

    simpleTreesData = cbind(simpleTreesData, inPRank = rank(simpleTreesData$inAverage)/length(simpleTreesData$inAverage) * 100, outPRank = rank(simpleTreesData$outAverage)/length(simpleTreesData$outAverage) *
        100)
    if (is.null(referenceTree) == FALSE) {
        colnames(refDistTracker) = c("distance", "treeName")
        myOutput <- list(cladeByCladeData = allTreesData, treeByTreeSummary = simpleTreesData, distanceToReferenceTree = refDistTracker)
    } else {
        myOutput <- list(cladeByCladeData = allTreesData, treeByTreeSummary = simpleTreesData, distanceToReferenceTree = NULL)
    }
    return(myOutput)
}

#' Determines which tips are missing from trees, and how often.
#' @param multiPhyloObject The tree object.
#' @export

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
                tipLabelCounter[tipLabelCounter$tip == thisTip, 2] <- as.numeric(tipLabelCounter[tipLabelCounter$tip == thisTip, 2]) + 1

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
