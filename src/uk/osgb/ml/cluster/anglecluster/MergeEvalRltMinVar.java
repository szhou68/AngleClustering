package uk.osgb.ml.cluster.anglecluster;

public class MergeEvalRltMinVar implements MergeEvalRlt {
	double variance; // variance if merged
	double weightSum;
	AngCluster cluFrom;
	
	public MergeEvalRltMinVar(double v, double w, AngCluster cf){
		variance = v;
		weightSum = w;
		cluFrom = cf;
	}
	@Override
	public int compareTo(Object o) {
		MergeEvalRltMinVar other = (MergeEvalRltMinVar)o;
		if(variance < other.variance){
			return -1;
		}else if(variance > other.variance){
			return 1;
		}else{
			if(weightSum < other.weightSum){
				return -1;
			}else if(weightSum > other.weightSum){
				return 1;
			}else{
				return cluFrom.compareTo(other.cluFrom);
			}
		}
	}
	//
	public double getVariance() {
		return variance;
	}
	public void setVariance(double variance) {
		this.variance = variance;
	}
	public double getWeightSum() {
		return weightSum;
	}
	public void setWeightSum(double weightSum) {
		this.weightSum = weightSum;
	}
	public AngCluster getCluFrom() {
		return cluFrom;
	}
	public void setCluFrom(AngCluster cluFrom) {
		this.cluFrom = cluFrom;
	}
	@Override
	public AngCluster getCluster() {
		// TODO Auto-generated method stub
		return cluFrom;
	}

}
