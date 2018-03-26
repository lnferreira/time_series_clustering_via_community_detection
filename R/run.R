source("tsDist.R")
source("tsClustAlgs.R")

# read dataset
ts_dataset = read.csv("../dataset/dataset.csv")
# get time series labels
labels = ts_dataset$class
# remove time series labels
ts_dataset$class = NULL
# transform it into a list 
ts_list = lapply(1:nrow(ts_dataset), function(i) ts_dataset[i,])
# for every plot, the user should press 'enter' to see the next plot
par(ask=T)
# plot the dataset
matplot(t(ts_dataset), pch=1, lty=1, t="l", xlab="time", ylab="value", main="Time series dataset", col=c(rep(2,10),rep(4,10)))
# calculate the distances for every pair of time series using DTW. 
# You can choose anyone distance function in tsDist.R or write your own one.
ts_dist = dist.parallel(tsList = ts_list, distFunc = tsdiss.dtw, cores = 1)
# Dist matrix normalization
ts_dist = dist.normalize(ts_dist)
# plot the distance matrix
heatmap(ts_dist, main = "Distance matrix")
# create the network connecting the k = 5 nearest neighbors (shortest distances)
net = net.knn.create(dist = ts_dist, k = 5)
# get the net layout, just to plot the network nodes in the same position
net_layout = layout_components(net)
# plot the network
plot(net, layout=net_layout)
# community detection using 
communities = cluster_louvain(net)
# plot communities
plot(communities, net, layout=net_layout)
# This is the main function:
# The output is similar to the kmeans function from the stats package.
# Given a distance matrix, it constructs the network, applies a community detection 
# algorithm, and returns the membership of each node.
clustering = ts.community.detection(dist = ts_dist, kOrEps = 5, method = "knn")
print(clustering)