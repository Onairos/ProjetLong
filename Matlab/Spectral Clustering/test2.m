load data;

sigma = 5;
nbclusters = 3;

[clusters, evalues, evectors] = spcl(x', nbclusters, sigma, 'kmean', [2 2]);