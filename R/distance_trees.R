#' Calculate the average distance between groups given a distance matrix and a clade
#' @param distanceMatrix The distance matrix
#' @param cladeAssignments The clade assignments
#' @param ID_names  List of IDs
#' @export


withinGroupDistances <- function(distanceMatrix, cladeAssignments, ID_names) {
    cladeList = data.frame(ID_names, cladeAssignments)
    idist = stats::na.omit(reshape::melt(distanceMatrix))
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
    colnames(allSelfDists) <- c("distance", "clade", "ID")
    allSelfDists$distance <- as.numeric(as.character(allSelfDists$distance))
    return(allSelfDists)
}

########
#' Generates bootstraps trees based on within-group standard deviation.
#' @param distanceMatrix A pairwise distance matrix.
#' @param numberOfTrees The number of bootstrap trees to generate. Default: 1000
#' @param withinSDtable The within clade SD table.
#' @param cladeAssignments Ordered list of clade assignments for taxa.
#' @param ID_names List of taxa names corresponding to the clade assignments.
#' @export
#'
generateBootstrapTrees <- function(distanceMatrix, numberOfTrees = 1000, withinSDtable, cladeAssignments, ID_names) {

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
                newVal = (thisDist + stats::rnorm(1, mean = 0, sd = as.numeric(thisSD)))  #randomly sample 1 value from rnorm of this SD
                x[b, a] <- newVal
                x[a, b] <- newVal
            }
        }
        x <- abs(x)  #Convert negative to positive
        treex <- ape::nj(stats::as.dist(x))
        treex$edge.length[treex$edge.length < 0] <- abs(treex$edge.length[treex$edge.length < 0])
        bstrees[[i]] <- treex
        message = paste("finished with tree number", i)
        print(message)
    }
    return(bstrees)
}
