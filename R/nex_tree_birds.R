

tree_bird <- ape::read.nexus("./trees/bird_species.nex")
plot(tree_bird)

birds <- ape::consensus(tree_bird, p = 1, check.labels = TRUE) # Looks shit

birds <- ape::consensus(tree_bird, p = 0.5, check.labels = TRUE) 
plot(birds)

http://blog.phytools.org/2016/04/consensus-methods-and-computing-average.html

birds <- phytools::averageTree(tree_bird, method="symmetric.difference") 

tree_bird[[1]]$tip.label