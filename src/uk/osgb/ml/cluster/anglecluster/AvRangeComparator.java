/*
 *  compare two angle ranges, assuming they are not overlapping except on boundary
 * 
 */
package uk.osgb.ml.cluster.anglecluster;

import java.util.Comparator;

public class AvRangeComparator implements Comparator<AngValRangeDir> {

	public int compare(AngValRangeDir o1, AngValRangeDir o2) {
		double angF1 = o1.getAngF();
		double angT1 = o1.getAngT();
		double angF2 = o2.getAngF();
		double angT2 = o2.getAngT();
		AngRef refF1 = o1.getRefF();
		AngRef refT1 = o1.getRefT();
		AngRef refF2 = o2.getRefF();
		AngRef refT2 = o2.getRefT();
		if(angT1 < angT2){
			return -1;
		}else if(angT1 > angT2){
			return 1;
		}else{
			if(refT1.compareTo(refT2) > 0){
				return 1;
			}else if(refT1.compareTo(refT2) < 0){
				return -1;
			}
		}
		return 0;

	}

}
