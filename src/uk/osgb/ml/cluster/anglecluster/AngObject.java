/** Abstract class for object with a SINGLE angle value, an ID and a weight
 * 
 */
package uk.osgb.ml.cluster.anglecluster;

public abstract class AngObject<T> implements Comparable {

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Object o){
		AngObject<T> other = (AngObject<T>)o;
		long id = this.getId();
		long id2 = other.getId();
		double weight = this.getWeight();
		double weight2 = other.getWeight();
		double ang = this.getAng();
		double ang2 = other.getAng();
		if(ang < ang2){
			return -1;
		}else if(ang > ang2){
			return 1;
		}else{
			if(weight < weight2){
				return -1;
			}else if(weight > weight2){
				return 1;
			}else{
				if(id < id2){
					return -1;
				}else if(id > id2){
					return 1;
				}
			}
		}
		return 0;
	}
	//
	public abstract double angDistance(AngObject<T> other);

	public abstract long getId();

	public abstract void setId(long id);

	public abstract double getAng();

	public abstract void setAng(double ang);

	public abstract double getWeight();

	public abstract void setWeight(double weight);
}
