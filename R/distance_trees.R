####### calculate the average distance between groups given a distance matrix and a clade list
#' Calculate the average distance between groups given a distance matrix and a clade
#' @param distanceMatrix The distance matrix
#' @param cladeAssignments The clade assignments
#' @param ID_names  List of IDs
#' @export


withinGroupDistances <- function(distanceMatrix, cladeAssignments, ID_names) {
    cladeList = data.frame(ID_names, cladeAssignments)
    idist = na.omit(reshape::melt(distanceMatrix))
    colnames(idist) <- c("a", "b", "dist")  #break matrix down into sets of pairwise distances between taxon a & b where dist = dist.
    allSelfDists = data.frame()
    clade_loop <- unique(cladeAssignments)
    # loop through all clade self comparisons, add to DF
    for (i in clade_loop) {
        # retrieve all genomes matching i as character
        itaxa <- as.character(cladeList[which(cladeList$cladeAssignments == i), 1])  # this is highly problematic
        if (length(itaxa > 1)) {
            for (a in itaxa) {
                intaxa = itaxa[itaxa != a]  #remove variable from list
                # retrieve values matching a
                ataxa = idist[idist$a %in% a, ]
                # calculate in-distances
                inDistancesA = ataxa[ataxa$b %in% intaxa, 3]
                if (length(inDistancesA > 0)) {
                  toAdd = cbind(as.numeric(inDistancesA), as.character(i), a)
                }
                allSelfDists = rbind(allSelfDists, toAdd)
            }
        }
    }
    colnames(allSelfDists) <- c("distance", "clade")
    allSelfDists$distance <- as.numeric(as.character(allSelfDists$distance))
    return(allSelfDists)
}

########
#' Bootstraps
#' @param distanceMatrix hk
#' @param numberOfTrees hkj
#' @param withinSDtable hkj
#' @param cladeAssignments hkj
#' @param ID_names hkj
#'
generateBootstrapTrees <- function(distanceMatrix, numberOfTrees, withinSDtable, cladeAssignments, ID_names) {

    cladeList = data.frame(ID_names, cladeAssignments)
    bstrees = vector("list", 0)  #create blank object to put noised trees in
    class(bstrees) <- "multiPhylo"
    loop = c(1:numberOfTrees)  #number of trees and intervals to make
    for (i in loop) {
        x <- distanceMatrix  #reset matrix
        for (a in as.character(cladeList$ID_names)) {
            # lookup this clade
            thisClade <- as.character(unlist(cladeList %>% filter(ID_names == as.character(a)) %>% select(cladeAssignments)))
            # determine SD to noise from
            thisSD <- withinSDtable %>% filter(clade == thisClade) %>% select(sd)
            # get all taxa in this clade excluding itself
            matchingTaxa <- as.list(cladeList %>% filter(cladeAssignments == thisClade) %>% select(ID_names) %>% filter(ID_names != as.character(a)))

            for (b in matchingTaxa$ID) {
                sublista <- idist[which(idist$a == a), ]
                sublist <- sublista[which((sublista$b) == b), ]  #retrieve the distance between these two strains
                thisDist = as.numeric(sublist[, 3])
                newVal = (thisDist + rnorm(1, mean = 0, sd = as.numeric(thisSD)))  #randomly sample 1 value from rnorm of this SD
                x[b, a] <- newVal
                x[a, b] <- newVal
            }
        }
        x <- abs(x)  #Convert negative to positive
        treex <- ape::nj(as.dist(x))
        treex$edge.length[treex$edge.length < 0] <- abs(treex$edge.length[treex$edge.length < 0])
        bstrees[[i]] <- treex
        message = paste("finished with tree number", i)
        print(message)
    }
    return(bstrees)
}



#takes a distance matrix, SD, and number of trees, and return the bootstrap trees for this SD via normal noising
########
#' noise
#' @param originalMatrix hkjh
#' @param standardDev kjhjkh
#' @param numTrees hkjh
#' @param rootTree kjhkj
#'
bootStrapWithNoiseFromNorm <- function(originalMatrix, standardDev, numTrees = 1000, rootTree = NULL) {
    bstrees = vector("list", 0)  #create blank object to put noised trees
    class(bstrees) <- "multiPhylo"
    loop = c(1:numTrees)  #number of trees and intervals to make
    for (i in loop) {
        x <- originalMatrix  #reset matrix
        rowsneeded = row.names(originalMatrix)  #reset taxa list
        for (j in rowsneeded) {
            rowsneeded = rowsneeded[rowsneeded != j]  #remove self from list
            for (k in rowsneeded) {
                oldval = x[j, k]
                newval = oldval + rnorm(1, mean = 0, sd = standardDev)
                x[j, k] <- newval
                x[k, j] <- newval
            }
        }
        x[x < 0] < -0  #set negative to zero
        treex <- nj(as.dist(x))
        treex$edge.length[treex$edge.length < 0] <- abs(treex$edge.length[treex$edge.length < 0])
        if (!is.null(rootTree)) {
            treex = ape::multi2di(treex)
        }
        bstrees[[i]] <- treex
        message = paste("finished with tree number", i)
        print(message)
    }
    return(bstrees)
}

#Goal: To take a range of SDs and a distance matrix as an input, noisea given number of rep trees, and to return the bootstrap values for each node in
## list form as an output
########
#' Bootstraps
#' @param originalMatrix hkjh
#' @param rangeMin hkj
#' @param rangeMax khj
#' @param numTrees khj
#' @param rootTree hkj
#'
bootStrapWithNoiseOverRange <- function(originalMatrix, rangeMin = 1, rangeMax = 1000, numTrees = 1000, rootTree = NULL) {
    sdRange = c(rangeMin:rangeMax)  #range of SD to check
    loop = c(1:numTrees)  #number of trees and intervals to make
    bootstrapList = vector("list", 0)  #blank list of bootstrap values
    mainTree = ape::nj(as.dist(originalMatrix))
    mainTree$edge.length[mainTree$edge.length < 0] <- abs(mainTree$edge.length[mainTree$edge.length < 0])
    for (a in sdRange) {
        bstrees = vector("list", 0)  #create blank object to put noised trees
        class(bstrees) <- "multiPhylo"
        for (i in loop) {
            x <- originalMatrix  #reset matrix
            for (j in rowsneeded) {
                others = rowsneeded[rowsneeded != j]  #remove self from list
                for (k in others) {
                  oldval = x[j, k]
                  newval = oldval + rnorm(1, mean = 0, sd = a)
                  x[j, k] <- newval
                  x[k, j] <- newval
                }
            }
            x[x < 0] < -0  #set negative to zero
            treex <- nj(as.dist(x))
            treex$edge.length[treex$edge.length < 0] <- abs(treex$edge.length[treex$edge.length < 0])
            if (!is.null(rootTree)) {
                treex = ape::root(treex, rootTree)
            }
            bstrees[[i]] <- treex
        }
        # perform bootstrap analysis
        if (!is.null(rootTree)) {
            clade_support <- prop.clades(mainTree, bstrees, rooted = TRUE)  #generate BS values based on trees
        } else {
            clade_support <- prop.clades(mainTree, bstrees, rooted = FALSE)  #generate BS values based on trees
        }
        clade_support = clade_support/numTrees * 100  #convert to percent
        clade_support[is.na(clade_support)] <- 0  #set NA to zero
        bootstrapList[[a]] <- clade_support
        message = paste("finished with BS trees from SD ", a)
        print(message)
    }
    return(bootstrapList)
}

