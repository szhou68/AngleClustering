/** an AngValue represents an numeric angle value in (-PI, PI].
 *  There may be multiple angle references to this value
 * 
 * 
 */
package uk.osgb.ml.cluster.anglecluster;

import java.util.TreeSet;

import uk.osgb.utilities.AngleUtility;

public class AngValue<T> extends AngObject<T> {
	double ang = -Math.PI;
	double totalWeight = 0.0;
	TreeSet<AngRef<T>> angRefs;
	//
	public AngValue(double angVal){
		ang = angVal;
	}
	public AngValue(AngRef<T> ref){
		ang = ref.getAng();
		angRefs = new TreeSet<AngRef<T>>();
		angRefs.add(ref);
	}
	public boolean addRef(AngRef<T> ref){
		if(angRefs == null){
			ang = ref.getAng();
			angRefs = new TreeSet<AngRef<T>>();
			angRefs.add(ref);
			return true;
		}else{
			if(ang == ref.getAng()){
				angRefs.add(ref);
				return true;
			}else{
				return false;
			}
		}
	}
	public double compTotalWeight(){
		if(angRefs!=null){
			for(AngRef<T> ref:angRefs){
				totalWeight+=ref.getWeight();
			}
		}
		return totalWeight;
	}
	public int getNumRefs(){
		if(angRefs!=null){
			return angRefs.size();
		}else{
			return 0;
		}
	}
	@Override
	public double angDistance(AngObject<T> other) {
		return AngleUtility.minAngDiff(ang, other.getAng());
	}

	@Override
	public long getId() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setId(long id) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double getAng() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setAng(double ang) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double getWeight() {
		// TODO Auto-generated method stub
		return totalWeight;
	}

	@Override
	public void setWeight(double weight) {
		// TODO Auto-generated method stub
		
	}

}
