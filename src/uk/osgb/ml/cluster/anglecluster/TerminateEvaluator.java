package uk.osgb.ml.cluster.anglecluster;

import java.util.Collection;

public interface TerminateEvaluator<T> {
	boolean terminate(Collection<AngCluster<T>> clus, MergeEvalRlt rlt);
}
