# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

source("tsdist.R")

NUM_CORES = max(1, detectCores() - 1)

# ================================================================================
# Time series clustering using community detection.
# Parameters:
#       dist: Distance matrix between every pair of time series. Many distance
#           Functions are implemented in tsClust.R
#       kOrEps: Parameter for the network construction method. 
#       method: Network construction method. Knn connects the k (kOrEps) most 
#           similar nodes extracted from the dist matrix. Eps connects the nodes 
#           that have distance lower than eps (kOrEps).
#       communityDetectionFunc: Community detection function. The igraph package
#           has some community detection algorithms implemented. 
#       cutat: If the community detection algorhtms is hierarchical, this value
#           cuts the dendrogram to achieve the desired number of clusters. Default
#           value is -1 (no cut).
# ================================================================================
ts.community.detection <- function(dist, kOrEps, method=c("knn", "eps"), 
                                   communityDetectionFunc=cluster_louvain, cutat=-1){
    if (method[[1]]=="knn") {
        net=net.knn.create(dist, kOrEps)
    } else if (method[[1]]=="eps") {
        net=net.epsilon.create(dist, kOrEps)   
    } 
    communities=communityDetectionFunc(net)
    if (cutat > 0) {
        communities=cutat(communities, cutat)
    } else {
        communities=membership(communities)
    }
    return(communities)
}

net.knn.create <- function(dist, k) {
    d = as.matrix(dist) + diag(Inf, nrow(dist), ncol(dist));
    net = graph.empty(nrow(dist), directed=F);
    knnRows = apply(d, 1, function(x) which(x %in% sort(x, method="quick")[1:k])[1:k]);
    knnRows = matrix(knnRows, nrow=k);
    for (j in 1:ncol(knnRows)) {
        for (i in 1:nrow(knnRows)) {
            if (!are.connected(net, j, knnRows[i,j])) {
                net = add.edges(net, c(j, knnRows[i,j]));
            }            
        }
    }
    net
}

net.epsilon.create <- function(dist, epsilon) {
    n = matrix(0, ncol(dist), nrow(dist))
    n[dist < epsilon] = 1;
    return(graph.adjacency(n, mode="undirected", diag=F));
}

# ================================================================================
# Traditional clustering algorithms
# ================================================================================
tsclust.traditional.clust <- function(distMatrix, k, alg=c("agnes", "diana", "pam", "complete_linkage", "single_linkage", 
                                         "average_linkage", "median_linkage", "centroid_linkage", "dbscan")) {
    rand = NULL
    clustering = NULL
    if (alg=="agnes") {
        clustering = agnes(as.dist(distMatrix))
    }
    if (alg=="diana") {
        clustering = diana(as.dist(distMatrix))
    }
    if (alg == "pam") {
        clustering = as.integer(pam(distMatrix, k, diss=T, cluster.only=T, do.swap=F))
    }
    if (alg == "complete_linkage") {
        clustering = hclust(as.dist(distMatrix), method="complete")
    }
    if (alg == "single_linkage") {
        clustering = hclust(as.dist(distMatrix), method="single")
    }
    if (alg == "average_linkage") {
        clustering = hclust(as.dist(distMatrix), method="average")
    }
    if (alg == "median_linkage") {
        clustering = hclust(as.dist(distMatrix), method="median")
    }
    if (alg == "centroid_linkage") {
        clustering = hclust(as.dist(distMatrix), method="centroid")
    }
    if (alg == "dbscan") {
        clustering = dbscan(as.dist(distMatrix), mean(distMatrix), method="dist")$cluster + 1
    }
    return(clustering)
}