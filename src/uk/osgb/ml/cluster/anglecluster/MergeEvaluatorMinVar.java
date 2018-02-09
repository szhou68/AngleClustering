package uk.osgb.ml.cluster.anglecluster;

import java.util.Collection;
import java.util.Vector;

public class MergeEvaluatorMinVar implements MergeEvaluator {

	/* (non-Javadoc)
	 * @see uk.osgb.ml.cluster.anglecluster.MergeEvaluator#evaluate(uk.osgb.ml.cluster.anglecluster.AngCluster, uk.osgb.ml.cluster.anglecluster.AngCluster)
	 */
	@Override
	public MergeEvalRlt evaluate(AngCluster cFrom, AngCluster cTo) {
		AngRef angF = cFrom.getAngF();
		AngRef angT = cTo.getAngT();
		Collection<AngRef> refs = new Vector<AngRef>(cFrom.members.size() + cTo.members.size());
		refs.addAll(cFrom.getMembers());
		refs.addAll(cTo.getMembers());
		//
		double varFrom = AngClusterGenerator.compRefColAngleVariation(cFrom.getMembers(), cFrom.getAngF(), cFrom.getAngT());
		double varTo = AngClusterGenerator.compRefColAngleVariation(cTo.getMembers(), cTo.getAngF(), cTo.getAngT());
		//
		double variance = AngClusterGenerator.compRefColAngleVariation(refs, angF, angT);
		// return the increase of circular variance if this merge takes place
		return new MergeEvalRltMinVar(variance - varFrom - varTo, cFrom.getWeight() + cTo.getWeight(), cFrom);
	}

}
