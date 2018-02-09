/** Angle Reference, with angle value in (-PI, PI], ID, weight, and an optional object ref 
 * 
 */
package uk.osgb.ml.cluster.anglecluster;

import uk.osgb.utilities.AngleUtility;

public class AngRef<T> extends AngObject<T>{
	static long idCnt = 0;
	long id;
	double ang;
	double weight;
	T objRef = null;
	AngCluster<T> cluster = null;

	/**
	 * @param ang
	 * @param weight
	 * @param objRef
	 */
	public AngRef(double ang, double weight, T objRef) {
		super();
		this.id = idCnt++;
		this.ang = ang;
		this.weight = weight;
		this.objRef = objRef;
	}
	
	/**
	 * @param ang
	 * @param weight
	 */
	public AngRef(double ang, double weight) {
		super();
		this.id = idCnt++;
		this.ang = ang;
		this.weight = weight;
	}
	
	/**
	 * @param ang
	 * @param weight
	 * @param _id
	 */
	public AngRef(double ang, double weight, long _id){
		this.id = _id;
		this.ang = ang;
		this.weight = weight;
	}
	//
	/* (non-Javadoc)
	 * @see com.osgb.data.utilities.AngObject#angDistance(com.osgb.data.utilities.AngRef)
	 */
	@Override
	public double angDistance(AngObject other){
		return AngleUtility.minAngDiff(ang, other.getAng());
	}
	
	// getter setter
	/* (non-Javadoc)
	 * @see uk.osgb.ml.cluster.anglecluster.AngObject#getId()
	 */
	@Override
	public long getId() {
		return id;
	}
	/* (non-Javadoc)
	 * @see com.osgb.data.utilities.AngObject#setId(long)
	 */
	@Override
	public void setId(long id) {
		this.id = id;
	}
	/* (non-Javadoc)
	 * @see com.osgb.data.utilities.AngObject#getAng()
	 */
	@Override
	public double getAng() {
		return ang;
	}
	/* (non-Javadoc)
	 * @see com.osgb.data.utilities.AngObject#setAng(double)
	 */
	@Override
	public void setAng(double ang) {
		this.ang = ang;
	}
	public void turnAng(double ang){
		this.ang = AngleUtility.angleTurn(this.ang, ang);
	}
	/* (non-Javadoc)
	 * @see com.osgb.data.utilities.AngObject#getWeight()
	 */
	@Override
	public double getWeight() {
		return weight;
	}
	/* (non-Javadoc)
	 * @see com.osgb.data.utilities.AngObject#setWeight(double)
	 */
	@Override
	public void setWeight(double weight) {
		this.weight = weight;
	}
	/**
	 * @return
	 */
	public T getObjRef() {
		return objRef;
	}
	/**
	 * @param objRef
	 */
	public void setObjRef(T objRef) {
		this.objRef = objRef;
	}
	/**
	 * @return
	 */
	public AngCluster<T> getCluster() {
		return cluster;
	}
	/**
	 * @param cluster
	 */
	public void setCluster(AngCluster<T> cluster) {
		this.cluster = cluster;
	}
	//
	/**
	 * @param refQ
	 * @param refF
	 * @param refT
	 * @return
	 */
	public static boolean refInRangeClosed(AngRef refQ, AngRef refF, AngRef refT){
		if(refQ.compareTo(refF) == 0 || refQ.compareTo(refT) == 0){
			return true;
		}
		return refInRangeOpen(refQ, refF, refT);
	}
	//
	/**
	 * @param refQ
	 * @param refF
	 * @param refT
	 * @return
	 */
	public static boolean refInRangeOpen(AngRef refQ, AngRef refF, AngRef refT){
		if(AngleUtility.containsAng(refQ.ang, refF.ang, refT.ang)){
			return true;
		}
		if(refQ.ang == refF.ang || refQ.ang == refT.ang){
			if(refT.compareTo(refF) > 0){
				if(refQ.compareTo(refF) > 0 && refT.compareTo(refQ) >0){
					return true;
				}
			}else if(refT.compareTo(refF) < 0){
				if(refQ.compareTo(refF) > 0 && refT.compareTo(refQ) > 0){
					return true;
				}
			}
		}
		return false;
	}
	//
	/**
	 * @param mod
	 * @return
	 */
	public AngRef<T> modulo(int mod){
		AngRef<T> newRef = null;
		double newAng = ang;
		if(mod == 2){// PI
			if(ang < 0.0){
				newAng = ang + Math.PI;
			}
			newRef = new AngRef(newAng, weight, id);
			newRef.setObjRef(objRef);
		}else if (mod == 4){
			double halfPI = Math.PI*0.5;
			if(ang >= halfPI){
				newAng = ang - halfPI;
			}else if(ang <= - halfPI){
				newAng = ang + Math.PI;
			}else if( ang < 0.0){
				newAng = ang + halfPI;
			}
			newRef = new AngRef(newAng, weight, id);
			newRef.setObjRef(objRef);
		}
		return newRef;
	}
	public static void main(String[] args){
		AngRef ref = new AngRef((-165.0/180)*Math.PI, 1.0, null);
		AngRef newRef = ref.modulo(2);
		System.out.println("done");
		
	}
}
