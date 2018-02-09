/**
 *  compare angle and then ID only
 * 
 */
package uk.osgb.ml.cluster.anglecluster;

import java.util.Comparator;

public class AngRefComparatorID implements Comparator<AngRef> {

	@Override
	public int compare(AngRef o1, AngRef o2) {
		long id = o1.id;
		long id2 = o2.id;
		double ang = o1.ang;
		double ang2 = o2.ang;
		if(ang < ang2){
			return -1;
		}else if(ang > ang2){
			return 1;
		}else{
			if(id < id2){
				return -1;
			}else if(id > id2){
				return 1;
			}
		}
		return 0;
	}

}
