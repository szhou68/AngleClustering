package uk.osgb.ml.cluster.anglecluster;

import java.util.Collection;

public class TerminateEvaluatorVar implements TerminateEvaluator {
	double totalVar = 0.0;
	double varTol = 0.0;
	//
	public TerminateEvaluatorVar(){}
	public TerminateEvaluatorVar(double tv, double vt){
		if(tv>= 0.0){
			totalVar = tv;
		}
		if(vt > 0.0 && vt <= 1.0){
			varTol = vt;
		}
	}
	public TerminateEvaluatorVar(Collection<? extends AngObject> refs, double vt){
		totalVar = AngClusterGenerator.compRefColAngleVariation(refs, null, null);
		if(vt > 0.0 && vt <= 1.0){
			varTol = vt;
		}
	}
	@Override
	public boolean terminate(Collection clus, MergeEvalRlt rlt) {
		if(varTol == 1.0){
			return false;
		}
		double withinVar = AngClusterGenerator.getWithinCluserVariation(clus);
		double varInc = ((MergeEvalRltMinVar)rlt).getVariance();
		if(totalVar == 0.0){
			return true;
		}
		if((withinVar+varInc) / totalVar > varTol){ 
			return true;
		}else{
			return false;
		}
	}

}
