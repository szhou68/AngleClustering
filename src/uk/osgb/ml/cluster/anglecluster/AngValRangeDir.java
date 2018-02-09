/** directed angle range, angF to angT in CCW direction
 */
package uk.osgb.ml.cluster.anglecluster;

import uk.osgb.utilities.AngleUtility;

public class AngValRangeDir<T> implements Comparable{
	double angF, angT, angDiff;
	double maxWeight, totalWeight;
	AngRef<T> refF, refT;

	public AngValRangeDir(double angF, AngRef<T> refF, double angT, AngRef<T> refT) {
		super();
		this.angF = angF;
		this.refF = refF;
		this.angT = angT;
		this.refT = refT;
		angDiff = AngleUtility.angFromTo(angF, angT);
		maxWeight = refF.weight > refT.weight?refF.weight:refT.weight;
		totalWeight = refF.weight + refT.weight;
	}
	@Override
	public int compareTo(Object o) {
		AngValRangeDir other = (AngValRangeDir)o;
		if(angDiff < other.angDiff){
			return -1;
		}else if(angDiff > other.angDiff){
			return 1;
		}else{// same angDiff
			
			if(totalWeight < other.totalWeight){
				return -1;
			}else if(totalWeight > other.totalWeight){
				return 1;
			}
			if(maxWeight < other.maxWeight){
				return -1;
			}else if(maxWeight > other.maxWeight){
				return 1;
			}
			if(angF < other.angF){
				return -1;
			}else if(angF > other.angF){
				return 1;
			}else{
				if(angT < other.angT){
					return -1;
				}else if(angT > other.angT){
					return 1;
				}else{
					if(refF.compareTo(other.refF) < 0){
						return -1;
					}else if(refF.compareTo(other.refF) > 0){
						return 1;
					}else if(refT.compareTo(other.refT) < 0){
						return -1;
					}else if(refT.compareTo(other.refT) > 0){
						return 1;
					}
				}
			}
		}
		return 0;
	}
	//
	public double getAngF() {
		return angF;
	}
	//
	public void setAngF(double angF) {
		this.angF = angF;
	}
	public double getAngT() {
		return angT;
	}
	public void setAngT(double angT) {
		this.angT = angT;
	}
	public double getAngDiff() {
		return angDiff;
	}
	public void setAngDiff(double angDiff) {
		this.angDiff = angDiff;
	}
	public AngRef<T> getRefF() {
		return refF;
	}
	public void setRefF(AngRef<T> refF) {
		this.refF = refF;
	}
	public AngRef<T> getRefT() {
		return refT;
	}
	public void setRefT(AngRef<T> refT) {
		this.refT = refT;
	}
}
