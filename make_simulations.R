library(NELSI)
library(phangorn)


# number of heterotachy events
n_events <- c(0, 1, 2, 5)
number_of_taxa <- c(10, 100)


rec_bind <- function(lmat){
    if(length(lmat) == 1) return(lmat[[1]])
    if(length(lmat) == 2) return(cbind(lmat[[1]], lmat[[2]]))
    if(length(lmat) > 2) return(cbind(lmat[[1]], rec_bind(lmat[-1])))
    #lmat <- lapply(letters[1:5], function(x) matrix(x, 5, 5))
}




for(nevent in n_events){
    for(ntax in number_of_taxa){
        for(i in 1:10){
            tr <- rtree(ntax)
            n <- nevent
         #sample_branches <- function(tr, n)
            nodes_nested <- vector()
            sampled_nodes <- vector()
            attempts <- 0
            repeat{
                if(length(sampled_nodes) >= n) break
                internal_nodes <- unique(tr$edge[, 1])
                root <- tr$edge[!tr$edge[, 1] %in% tr$edge[, 2], 1][1]
                internal_nodes <- internal_nodes[!internal_nodes %in% root]
                sampled <- sample(internal_nodes, 1)
                descendants <- get.descending.nodes.branches(tr, sampled)$descending.nodes
                if(!any(descendants %in% nodes_nested)){
                    nodes_nested <- c(nodes_nested, descendants, sampled)
                    sampled_nodes <- c(sampled_nodes, sampled)
                }else{
                    attempts <- attempts + 1
                }
                if(attempts > 20) break
            }
            write.tree(tr, file = paste0('simulations_gamma/ntax', ntax, '_event', nevent,
                                         '_true_tree_sim', i, '_gamma.tree'))
            gamma_rates <- phangorn:::discrete.gamma(0.1, 4)
            part1_500 <- rec_bind(lapply(gamma_rates,
                                         function(rate) as.DNAbin(simSeq(tr, l = 125, rate = rate))))
            part1_250 <- rec_bind(lapply(gamma_rates,
                                         function(rate) as.DNAbin(simSeq(tr, l = 62, rate = rate))))
            if(nevent > 0){
                tr$edge.length[tr$edge[, 2] %in% sampled_nodes] <- tr$edge.length[tr$edge[, 2] %in% sampled_nodes] * 10
                tr$node.label <- rep('original', tr$Nnode)
                tr$node.label[tr$edge[tr$edge[, 2] %in% sampled_nodes, 2]-length(tr$tip.label)] <- 'heterotachy'
                write.tree(tr, file = paste0('simulations_gamma/ntax', ntax, '_event', nevent,
                                         '_hetero_tree_sim', i, '_gamma.tree'))
            }
            part2_500 <- rec_bind(lapply(gamma_rates,
                                         function(rate) as.DNAbin(simSeq(tr, l = 500, rate = rate))))
            part2_250 <- rec_bind(lapply(gamma_rates,
                                         function(rate) as.DNAbin(simSeq(tr, l = 250, rate = rate))))
            aln_1000 <- cbind(part1_500, part2_500)
            aln_500 <- cbind(part1_250, part2_250)
            write.dna(aln_1000, file = paste0('simulations_gamma/1000nt_ntax', ntax, '_event_', nevent,
                      '_hetero_sim', i, '_gamma.fasta'), format = 'fasta', nbcol = -1, colsep = '')
            write.dna(aln_500, file = paste0('simulations_gamma/500nt_ntax', ntax, '_event_', nevent,
                      '_hetero_sim', i, '_gamma.fasta'), format = 'fasta', nbcol = -1, colsep = '')
        }
    }
}

