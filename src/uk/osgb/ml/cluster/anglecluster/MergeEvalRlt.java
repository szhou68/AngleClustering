package uk.osgb.ml.cluster.anglecluster;

public interface MergeEvalRlt extends Comparable{
	AngCluster getCluster(); // from_cluster of the proposed merge
	double getVariance(); // increase of within-cluster variation if merged 
}
