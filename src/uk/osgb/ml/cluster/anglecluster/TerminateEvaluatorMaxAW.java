/*  
 * 
 * 
 * 
 * 
 */
package uk.osgb.ml.cluster.anglecluster;

import java.util.Collection;

import uk.osgb.utilities.AngleUtility;

public class TerminateEvaluatorMaxAW implements TerminateEvaluator {
	//
	double totalVar = 0.0;
	double varTol = 1.0;
	double maxAngWidth = 0.0;
	
	public TerminateEvaluatorMaxAW(double maxWidthDeg){
		maxAngWidth = Math.toRadians(maxWidthDeg);
	}
	public TerminateEvaluatorMaxAW(double totalVar, double varTol, double maxWidthDeg){
		this.totalVar = totalVar;
		this.varTol = varTol;
		maxAngWidth = Math.toRadians(maxWidthDeg);
	}
	@Override
	public boolean terminate(Collection clus, MergeEvalRlt rlt) {
		AngCluster cluFrom = rlt.getCluster();
		AngCluster cluTo = cluFrom.next;
		if(cluFrom == cluTo){
			return true;
		}
		double angF = cluFrom.getAngF().getAng();
		double angT = cluTo.getAngT().getAng();
		double width = AngleUtility.angFromTo(angF, angT);
		if(varTol == 1.0) {// consider angular width only
			if(width > maxAngWidth){
				return true;
			}else{
				return false;
			}
		}else{// even if within max width, 
			if(width <= maxAngWidth) {
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
			}else {
				return true;
			}

		}
	}

}
