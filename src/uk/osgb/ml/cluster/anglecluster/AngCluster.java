/** cluster of angle references
 * 
 */
package uk.osgb.ml.cluster.anglecluster;

import java.util.Iterator;
import java.util.TreeSet;
import uk.osgb.utilities.AngleUtility;


/**
 * @author Sheng.Zhou@os.uk
 *
 * @param <T>
 */
public class AngCluster<T> extends AngObject<T>{
	static long cIdCnt = 0;
	long cId;
	double angMean;
	double weightSum = 0.0;
	double lenSum = 0.0; // weights projected to mean angle of the cluster
	double variance = 0.0;
	AngCluster<T> next, prev;
	TreeSet<AngRef<T>> members = new TreeSet<AngRef<T>>(); // ordered by ref ang, weight and then id
	//
	AngRef<T> angF, angT; // range of angles covered by this cluster
	//
	public AngCluster(AngRef<T> ref){
		angMean = ref.getAng();
		this.cId = cIdCnt++;
		members.add(ref);
	}
	public AngCluster(double initAng){
		angMean = initAng;
		this.cId = cIdCnt++;
	}
	public AngCluster(){
		this.cId = cIdCnt++;
	}
	public boolean addAngleRef(AngRef<T> ang){
		ang.setCluster(this);
		return members.add(ang);
	}
	//
	public boolean removeAngleRef(AngRef<T> ang){
		return members.remove(ang);
	}
	//
	public AngCluster<T> getNext() {
		return next;
	}
	public void setNext(AngCluster<T> next) {
		this.next = next;
	}
	public AngCluster<T> getPrev() {
		return prev;
	}
	public void setPrev(AngCluster<T> prev) {
		this.prev = prev;
	}
	
	/** after all angles are added, this method should be called to find the cluster boundary (maximum "external" angle between the these two angles)
	 *  this should work even if the angle coverage is over PI 
	 * @return
	 */
	public boolean compFromToAngles(){
		if(members.isEmpty()){
			return false;
		}
		AngRef<T> minRef = members.first();
		if(members.size() == 1){
			angF = minRef;
			angT = minRef;
			return true;
		}
		AngRef<T> maxRef = members.last();
		//
		double maxAng = AngleUtility.angFromTo(maxRef.getAng(), minRef.getAng());
		angF = minRef;
		angT = maxRef;
		
		AngRef<T> ref0 = members.first();
		Iterator<AngRef<T>> iter = members.iterator();
		iter.next(); // skip first on (ref0)
		while(iter.hasNext()){
			AngRef<T> ref1 = iter.next();
			double angDiff = AngleUtility.angFromTo(ref0.getAng(), ref1.getAng());
			if(angDiff > maxAng){
				maxAng = angDiff;
				angF = ref1;
				angT = ref0;
			}
			ref0 = ref1;
		}
		return true;
 	}
	//
	public void merge(AngCluster<T> nextClu){
		angT = nextClu.angT;
		members.addAll(nextClu.getMembers());
		next = nextClu.next;
		next.prev = this;
		if(prev == nextClu){
			prev = this;
		}
		compWeightSum();
		compProjectedWeightSum();
		compClusterMean();
		compVariation();
		
	}
	//
	public double getAngWidth(){
		return AngleUtility.angFromTo(angF.ang, angT.ang);
	}
	//
	public double compWeightSum(){
		double sum = 0.0;
		if(members.isEmpty()){
			return sum;
		}
		for(AngRef<T> ref: members){
			sum+= ref.getWeight();
		}
		return weightSum = sum;
	}
	//
	public double compProjectedWeightSum(){
		double sum = 0.0;
		if(members.isEmpty()){
			return sum;
		}
		double angFtoM = AngleUtility.angFromTo(angF.getAng(), angMean);
		for(AngRef<T> ref: members){
			double weight = ref.getWeight();
			double ang = ref.getAng();
			double angDiffM = AngleUtility.angFromTo(angMean, ang);
			double projW =weight*Math.cos(angDiffM); // projection on angMean, maybe negative
			sum+=projW;
		}
		return lenSum = sum;
		
	}
	//
	public double compClusterMean(){
		if(members.isEmpty()){
			return Double.MAX_VALUE;// not initilaised			
		}
		if(members.size() ==1){
			weightSum = angF.getWeight();
			return angMean = angF.getAng();
		}else{
			return angMean = AngClusterGenerator.compRefColAngleMean(members, angF, angT);
		}
	}
	
	public double compVariation(){
		if(members.isEmpty()){
			return 0.0;// not initilaised			
		}
		if(members.size() ==1){
			return 0.0;
		}
		//
		double varSum = 0.0;
		varSum = AngClusterGenerator.compRefColAngleVariation(members, angF, angT);
		return variance = varSum;
	}
	//
	public long getId() {
		return cId;
	}

	public void setId(long cId) {
		this.cId = cId;
	}

	public double getAng() {
		return angMean;
	}

	public void setAng(double ang) {
		this.angMean = ang;
	}

	public double getWeight() {
		return weightSum;
	}

	public void setWeight(double weightSum) {
		this.weightSum = weightSum;
	}

	public double getVariance() {
		return variance;
	}

	public void setVariance(double variance) {
		this.variance = variance;
	}

	public TreeSet<AngRef<T>> getMembers() {
		return members;
	}

	public void setMembers(TreeSet<AngRef<T>> members) {
		this.members = members;
	}

	public AngRef<T> getAngF() {
		return angF;
	}
	public void setAngF(AngRef<T> angF) {
		this.angF = angF;
	}
	public AngRef<T> getAngT() {
		return angT;
	}
	public void setAngT(AngRef<T> angT) {
		this.angT = angT;
	}
	@Override
	public double angDistance(AngObject<T> other) {
		return AngleUtility.minAngDiff(angMean, other.getAng());
	}
	public void clearAngRefs(){
		members.clear();
	}
	public double getLenSum() {
		return lenSum;
	}
	public void setLenSum(double lenSum) {
		this.lenSum = lenSum;
	}
	public void report(){
		System.out.println(members.size() + " mean: "+angMean*180.0/Math.PI+ " weight: "+weightSum + " projected weight: "+ lenSum + " variance: "+variance);
	}
	public String report(int precision){
		String fmt = "%."+precision+"f";
		String rtn = members.size()+", " + String.format(fmt, angMean*180.0/Math.PI) + ", " + String.format(fmt, weightSum) + ", " + String.format(fmt, lenSum) +", " + String.format(fmt,  variance);
		return rtn;
	}

}
