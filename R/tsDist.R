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

library("TSclust")
library("TSdist") 
library("entropy")
library("DescTools")
library("compiler")
library("scales")
library("parallel")
library("igraph")

NUM_CORES = max(1, detectCores() - 1)

# ==============================================================================
# Parallel Dist with mcmapply
# ==============================================================================
dist.parallel <- function(tsList, distFunc=tsdiss.euclidean, cores=NUM_CORES, outputFilePath) {
    distFuncCompiled <- cmpfun(distFunc)
    tsListLength = length(tsList)
    combs = combn(tsListLength, 2)
    d = mcmapply(dist.parallel2.compute, x=combs[1,], y=combs[2,], 
                 MoreArgs=list(tsList=tsList, distFunc=distFuncCompiled), mc.cores=cores)
    dist = matrix(0, tsListLength, tsListLength)
    dist[lower.tri(dist)] = d
    dist = as.matrix(as.dist(dist))
    if(!missing(outputFilePath)) {
        cat("\nsaving file:", outputFilePath, "\n")
        saveRDS(dist, outputFilePath)
    }
    dist
}

dist.parallel2.compute <- function(x, y, tsList, distFunc) { 
    distFunc(tsList[[x]], tsList[[y]])
}

# ==============================================================================
# Distance Matrix normalization
# ==============================================================================
dist.normalize <- function(dist) {
    distNorm = matrix(0, nrow(dist), ncol(dist))
    d = dist[upper.tri(dist)] 
    d = rescale(d)
    distNorm[upper.tri(distNorm)] = d
    distNorm = distNorm + t(distNorm)
    distNorm
}

# ==============================================================================
# [X] Person Correlation
# ==============================================================================
tsdiss.correlation <- function(ts1, ts2) {
    CorDistance(ts1, ts2)
} 

# ==============================================================================
# [ ] Mutual Information
# ==============================================================================
tsdiss.mutualInformation <- function(ts1, ts2) {
    MutInf(ts1,ts2)
}

# ==============================================================================
# [X] DTW
# ==============================================================================
tsdiss.dtw <- function(ts1, ts2) {
    diss.DTWARP(ts1,ts2)
} 

# ==============================================================================
# [X] Euclidean L_2
# ==============================================================================
tsdiss.euclidean <- function(ts1, ts2) {
    diss.EUCL(ts1, ts2)
} 

# ==============================================================================
# [X] Manhattan L_1 (City Block)
# ==============================================================================
tsdiss.manhattan <- function(ts1, ts2) {
    sum(abs(ts1-ts2))
} 

# ==============================================================================
# [X] Infinite Norm L_{\inf} (chebyshev)
# ==============================================================================
tsdiss.infiniteNorm <- function(ts1, ts2) {
    max(abs(ts1-ts2))
} 

# ==============================================================================
# [X] Sorensen
# ==============================================================================
tsdiss.sorensen <- function(ts1, ts2) {
    sum(abs(ts1-ts2))/sum(abs(ts1+ts2))
} 

# ==============================================================================
# [X] Gower
# ==============================================================================
tsdiss.gower <- function(ts1, ts2) {
    sum(abs(ts1-ts2)) / max(length(ts1), length(ts2))
}

# ==============================================================================
# [-] Soergel
#     obs: todos valores iguais a 2 pois series estao normalizadas (z-score)
# ==============================================================================
tsdiss.soergel <- function(ts1, ts2) {
    sum(abs(ts1-ts2))/sum(mapply(max, ts1, ts2))
}

# ==============================================================================
# [-] Kulczynski Dissimilarity
#     obs: todos valores iguais a -2 pois series estao normalizadas (z-score)
# ==============================================================================
tsdiss.kulczynski <- function(ts1, ts2) {
    sum(abs(ts1-ts2))/sum(mapply(min, ts1, ts2))
}

# ==============================================================================
# [X] Canberra
# ==============================================================================
tsdiss.canberra <- function(ts1, ts2) {
    sum(abs(ts1-ts2)/(ts1+ts2))
}

# ==============================================================================
# [X] Lorentzian
# ==============================================================================
tsdiss.lorentzian <- function(ts1, ts2) {
    sum(log(1+abs(ts1-ts2)))
}

# ==============================================================================
# [X] Intersection Similarity
# ==============================================================================
tssim.intersection <- function(ts1, ts2) {
    sum(mapply(min, ts1, ts2))
}

# ==============================================================================
# [X] Intersection Dissimilarity == manhattan / 2
# ==============================================================================
tsdiss.intersection <- function(ts1, ts2) {
    sum(abs(ts1-ts2))/2
}

# ==============================================================================
# [X] Wave Hedges
# ==============================================================================
tsdiss.waveHedges <- function(ts1, ts2) {
    sum(log(1+abs(ts1-ts2)))
}

# ==============================================================================
# [X] Czekanowski similarity
# ==============================================================================
tssim.czekanowski <- function(ts1, ts2) {
    2 * sum(mapply(min, ts1, ts2)) / sum(ts1+ts2)
}

# ==============================================================================
# [X] Czekanowski dissimilarity
# ==============================================================================
tsdiss.czekanowski <- function(ts1, ts2) {
    sum(abs(ts1-ts2)) / sum(ts1+ts2)
}

# ==============================================================================
# [X] Motyka similarity
# ==============================================================================
tssim.motyka <- function(ts1, ts2) {
    sum(mapply(min, ts1, ts2)) / sum(ts1+ts2)
}

# ==============================================================================
# [X] Motyka dissimilarity
# ==============================================================================
tsdiss.motyka <- function(ts1, ts2) {
    sum(mapply(max, ts1, ts2)) / sum(ts1+ts2)
}

# ==============================================================================
# [-] Kulczynski Similarity == 1 / Kulczynski Dissimilarity
#     obs: todos valores iguais pois series estao normalizadas (z-score)
# ==============================================================================
tssim.kulczynski <- function(ts1, ts2) {
    sum(mapply(min, ts1, ts2)) / sum(abs(ts1-ts2))
}

# ==============================================================================
# [-] Ruzicka
#     obs: todos valores iguais pois series estao normalizadas (z-score)
# ==============================================================================
tssim.ruzicka <- function(ts1, ts2) {
    sum(mapply(min, ts1, ts2)) / sum(mapply(max, ts1, ts2))
}

# ==============================================================================
# [-] Tanimoto
#     obs: todos valores iguais pois series estao normalizadas (z-score)
# ==============================================================================
tsdiss.tanimoto <- function(ts1, ts2) {
    sumTs1 = sum(ts1)
    sumTs2 = sum(ts2)
    sumMinTs1Ts2 = sum(mapply(min, ts1, ts2))
    (sumTs1 + sumTs2 - (2 * sumMinTs1Ts2)) / (sumTs1 + sumTs2 - sumMinTs1Ts2)
}

# ==============================================================================
# [ ] Inner product similarity
# ==============================================================================
tssim.innerProduct <- function(ts1, ts2) {
    sum(ts1*ts2)
}

# ==============================================================================
# [ ] Harmonic mean similarity
# ==============================================================================
tssim.harmonicMean <- function(ts1, ts2) {
    2*sum((ts1*ts2)/(ts1+ts2))
}

# ==============================================================================
# [ ] Cosine similarity
# ==============================================================================
tssim.cosine <- function(ts1, ts2) {
    sum(ts1*ts2)/(sqrt(sum(ts1^2))*sqrt(sum(ts2^2)))
}

# ==============================================================================
# [ ] Kumar-Hassebrook (Jaccard) Revisar se é igual ao jaccardSim mesmo
# ==============================================================================
tssim.kumarHasserbrook <- function(ts1, ts2) {
    sumTs1TimesTs2 = sum(ts1*ts2)
    sumTs1TimesTs2 / ((sum(ts1^2)) + (sum(ts2^2)) - sumTs1TimesTs2)
}

# ==============================================================================
# [ ] Jaccard Dissimilarity
# ==============================================================================
tsdiss.jaccard <- function(ts1, ts2) {
    sum((ts1-ts2)^2) / ((sum(ts1^2)) + (sum(ts2^2)) - sum(ts1*ts2))
}

# ==============================================================================
# [ ] Dice similarity
# ==============================================================================
tssim.dice <- function(ts1, ts2) {
    2*sum(ts1*ts2) / (sum(ts1^2) + sum(ts2^2))
}

# ==============================================================================
# [ ] Dice dissimilarity
# ==============================================================================
tsdiss.dice <- function(ts1, ts2) {
    sum((ts1-ts2)^2) / (sum(ts1^2) + sum(ts2^2))
}

# ==============================================================================
# [ ] Fidelity similarity
# ==============================================================================
tssim.fidelity <- function(ts1, ts2) {
    sum(sqrt(ts1*ts2))
}

# ==============================================================================
# [ ] Bhattacharyya dissimilarity
# ==============================================================================
tsdiss.bhattacharyya <- function(ts1, ts2) {
    -log(sum(sqrt(ts1*ts2)))
}

# ==============================================================================
# [ ] Hellinger dissimilarity
# ==============================================================================
tsdiss.hellinger <- function(ts1, ts2) {
    2*sqrt(1-sum(sqrt(ts1*ts2)))
}

# ==============================================================================
# [ ] Matusita dissimilarity
# ==============================================================================
tsdiss.matusita <- function(ts1, ts2) {
    sqrt(1-sum(sqrt(ts1*ts2)))
}

# ==============================================================================
# [ ] Squared-chord dissimilarity
# ==============================================================================
tsdiss.squaredChord <- function(ts1, ts2) {
    sum((sqrt(ts1) - sqrt(ts2))^2)
}

# ==============================================================================
# [ ] Squared-chord similarity
# ==============================================================================
tssim.squaredChord <- function(ts1, ts2) {
    2*sum(sqrt(ts1*ts2))-1
}

# ==============================================================================
# [ ] Squared Euclidean dissimilarity
# ==============================================================================
tsdiss.squaredEuclidean <- function(ts1, ts2) {
    sum((ts1-ts2)^2)
}

# ==============================================================================
# [ ] Pearson Chi Squared dissimilarity
# ==============================================================================
tsdiss.pearsonChiSquared <- function(ts1, ts2) {
    sum(((ts1-ts2)^2)/ts2)
}

# ==============================================================================
# [ ] Neyman Chi Squared dissimilarity
# ==============================================================================
tsdiss.neymanChiSquared <- function(ts1, ts2) {
    sum(((ts1-ts2)^2)/ts1)
}

# ==============================================================================
# [ ] Squared Chi Squared dissimilarity
# ==============================================================================
tsdiss.squaredChiSquared <- function(ts1, ts2) {
    sum(((ts1-ts2)^2)/(ts1+ts2))    
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.probabilisticSymmetricChiSquared <- function(ts1, ts2) {
    2*sum(((ts1-ts2)^2)/(ts1+ts2))    
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.divergence <- function(ts1, ts2) {
    2*sum((ts1-ts2)^2 / (ts1+ts2)^2)
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.clark <- function(ts1, ts2) {
    sqrt(sum((abs(ts1-ts2)/(ts1+ts2))^2))
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.additiveSymmetricChiSquared <- function(ts1, ts2) {
    sum((((ts1-ts2)^2)*(ts1+ts2))/(ts1*ts2))
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.kullbackLeibler <- function(ts1, ts2) {
    sum(ts1*log(ts1/ts2))
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.jeffreys <- function(ts1, ts2) {
    sum((ts1-ts2) * log(ts1/ts2))
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.kDivergence <- function(ts1, ts2) {
    sum(ts1*log((2*ts1)/(ts1+ts2)))
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.topsoe <- function(ts1, ts2) {
    sum(ts1*log((2*ts1)/(ts1+ts2)) + ts2*log((2*ts2)/(ts1+ts2)))
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.jensenShannon <- function(ts1, ts2) {
    (sum(ts1*log((2*ts1)/(ts1+ts2))) + sum(ts2*log((2*ts2)/(ts1+ts2)))) / 2
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.jensenDifference <- function(ts1, ts2) {
    sum(((ts1*log(ts1) + ts2*log(ts2))/2) - (((ts1+ts2)/2)*log((ts1+ts2)/2)))
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.taneja <- function(ts1, ts2) {
    sum(((ts1+ts2)/2)*log((ts1+ts2)/(2*sqrt(ts1*ts2))))
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.kumarJohnson <- function(ts1, ts2) {
    sum(((ts1^2 - ts2^2)^2) / (2*(ts1*ts2)^(3/2)))
}

# ==============================================================================
# [ ] 
# ==============================================================================
tsdiss.avgL1LInf <- function(ts1, ts2) {
    (sum(abs(ts1-ts2)) + max(abs(ts1-ts2))) / 2
}

# ==============================================================================
# [ ] NEW DISTANCES:
# ==============================================================================

tsdiss.acfd <- function(ts1, ts2) {
    ACFDistance(ts1,ts2)
}

tsdiss.arlpcCeps <- function(ts1, ts2) {
    ARLPCCepsDistance(ts1, ts2)
}

tsdiss.armah <- function(ts1, ts2) {
    as.numeric(ARMahDistance(ts1,ts2)$p_value)
}

tsdiss.arpic <- function(ts1, ts2) {
    ARPicDistance(ts1, ts2)
}

# lag.max=(min(length(x), length(y))-1)
tsdiss.ccor <- function(ts1, ts2) {
    CCorDistance(ts1, ts2)
}

tsdiss.cdmd <- function(ts1, ts2) {
    CDMDistance(ts1, ts2)
}

tsdiss.cid <- function(ts1, ts2) {
    CIDDistance(ts1, ts2)
}

tsdiss.cort <- function(ts1, ts2) {
    CortDistance(ts1, ts2)
}

tsdiss.dissim <- function(ts1, ts2) {
    DissimDistance(ts1,ts2)
}

tsdiss.fourierDist <- function(ts1, ts2) {
    FourierDistance(ts1, ts2)
}

tsdiss.intper <- function(ts1, ts2) {
    IntPerDistance(ts1, ts2)
}

tsdiss.ncd <- function(ts1, ts2) {
    NCDDistance(ts1, ts2)
}

tsdiss.pacf <- function(ts1, ts2) {
    PACFDistance(ts1, ts2)
}

tsdiss.pdcd <- function(ts1, ts2) {
    PDCDistance(ts1, ts2)
}

tsdiss.sts <- function(ts1, ts2) {
    STSDistance(ts1, ts2)
}

tsdiss.per <- function(ts1, ts2) {
    PerDistance(ts1, ts2)
}

tsdiss.wav <- function(ts1, ts2) {
    WavDistance(ts1, ts2)
}