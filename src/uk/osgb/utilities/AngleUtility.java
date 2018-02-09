package uk.osgb.utilities;

public class AngleUtility {
	public static double rad60 = Math.PI / 3.0; // PI == 3.141592653589793
	public static double rad90 = Math.PI * 0.5;
	//public static double rad120 = 2*rad60 - 0.000000000000001;
	public static double rad120 = 2*rad60;
	public static double rad22_5 = Math.PI*0.125;
	public static double rad45 = Math.PI*0.25;
	public static double rad15 = Math.PI/12.0;
	public static double rad10 = Math.PI/18.0;
	public static double rad05 = Math.PI/36.0;
//
	//
	//  ang1/ang2 are in [-PI, PI]
	// return value in degree
	//
	/**
	 * @param ang1 angle in (-PI, PI]
	 * @param ang2 angle in (-PI, PI]
	 * @return the minimum angle between ang1 and ang2 (converted to degree) in [0, 180]
	 */
	public static double minAngDiffDeg(double ang1, double ang2){
		double diff = Math.max(ang1, ang2) - Math.min(ang1, ang2);
		if(diff > Math.PI)
			diff = 2*Math.PI - diff;
		if(diff < 0.0)
			diff = 0.0;
		return 180*diff/Math.PI;
	}
	// [0, PI]
	/**
	 * @param ang1
	 * @param ang2
	 * @return
	 */
	public static double minAngDiff(double ang1, double ang2){
		double diff = Math.max(ang1, ang2) - Math.min(ang1, ang2);
		if(diff > Math.PI)
			diff = 2*Math.PI - diff;
		if(diff < 0.0)
			diff = 0.0;
		return diff;
	}
	//
	//  angF/angT are in [-PI, PI], returns the angle (always > 0) turned from angF to angT
	//
	/**
	 * @param angF
	 * @param angT
	 * @return FromTo angle in [0.0, 2*PI)
	 */
	public static double angFromTo(double angF, double angT){
		double ang = angT - angF;
		if(ang < 0.0){
			ang = 2*Math.PI + ang;
		}
		return ang;
	}
	// turning turnAng from base, return angle value in [-PI, PI]
	
	/**
	 * @param base the base angle, in [-PI, PI]
	 * @param turnAng angle to turn, in [-2PI, 2PI], postive for CCW and negative for CW turn
	 * @return angle after turning, in [-PI, PI]
	 */
	public static double angleTurn(double base, double turnAng){
		double ang = base + turnAng;
		if(ang >= 2*Math.PI){
			ang = ang - 2*Math.PI; 
		}else if (ang > Math.PI){
			ang = ang - 2*Math.PI;
		}else if(ang <=-2*Math.PI){
			ang = ang + 2*Math.PI;
		}else if(ang < -Math.PI){
			ang = 2*Math.PI + ang;
		}
		return ang;
	}
	//
	//
	/**
	 * @param angF
	 * @param angT
	 * @return
	 */
	public static double angFromToDeg(double angF, double angT){
		return 180*angFromTo(angF, angT)/Math.PI;
	}
	// turn ang by 90 degree ccw or cw
	/**
	 * @param ang
	 * @param ccw
	 * @return
	 */
	public static double angTurn90Deg(double ang, boolean ccw){
		double newAng = ang;
		if(ccw){
			newAng = ang + 0.5*Math.PI;
			if(newAng > Math.PI){
				newAng = newAng - 2*Math.PI;
			}
		}else{
			newAng = ang - 0.5*Math.PI;
			if(newAng < -Math.PI){
				newAng = 2*Math.PI + newAng;
			}
		}
		return 180*newAng/Math.PI;
	}
	//
	/**
	 * @param ang
	 * @param ccw
	 * @return
	 */
	public static double angTurnHalfPI(double ang, boolean ccw){
		double newAng = ang;
		if(ccw){
			newAng = ang + 0.5*Math.PI;
			if(newAng > Math.PI){
				newAng = newAng - 2*Math.PI;
			}
		}else{
			newAng = ang - 0.5*Math.PI;
			if(newAng < -Math.PI){
				newAng = 2*Math.PI + newAng;
			}
		}
		return newAng;
		
	}
	//
	public static double angSum(double ang1, double ang2){
		double newAng = ang1 + ang2;
		if(newAng > Math.PI){
			newAng = newAng - 2*Math.PI;
		}else if(newAng < -Math.PI){
			newAng = 2*Math.PI + newAng;
		}
		return newAng;
	}
	/**
	 * @param ang
	 * @return
	 */
	public static double angReverse(double ang){
		if(ang > 0){
			return ang - Math.PI;
		}else{
			return Math.PI + ang;
		}
	}
	// the angle bisector of angle(angF, angT)
	/**
	 * @param angF
	 * @param angT
	 * @return
	 */
	public static double angBisector(double angF, double angT){
		double ang = angFromTo(angF, angT) * 0.5;
		ang +=angF;
		if(ang > Math.PI){
			ang = ang - 2*Math.PI;
		}
		return ang;
		
	}
	/** bisector the minimum angle between ang1 and ang2, if angle is PI, use from 1 to 2
	 * @param ang1
	 * @param ang2
	 * @return
	 */
	public static double angBisectorMin(double ang1, double ang2){
		double diff12 = angFromTo(ang1, ang2);
		if(diff12 <=Math.PI){
			return angBisector(ang1, ang2);
		}else{
			return angBisector(ang2, ang1);
		}
	}
	//
	// if angQ is between angF and angT (turning from angF to angT) - note that angF/angT is not included
	//
	/** assume angF != angT
	 * @param angQ query angle
	 * @param angF FROM angle
	 * @param angT TO angle
	 * @return
	 */
	public static boolean containsAng(double angQ, double angF, double angT){
		if(angQ == angF || angQ == angT){
			return false;
		}
		//if(Math.abs(angFromTo(angF, angT) - angFromTo(angF, angQ) - angFromTo(angQ, angT)) <= Double.MIN_VALUE){
		if(angFromTo(angF, angQ ) < angFromTo(angF, angT)){
			return true;
		}else{
			return false;
		}
	}
	/**
	 * @param theta
	 * @return one of eight directions: E, NE, N, NW, W, SW, S, SE
	 */
	public static String ang2Dir8(double theta){
		if(theta >=rad22_5){
			if(theta >= 3*rad22_5){
				if(theta > 5*rad22_5){
					if(theta > 7*rad22_5){
						return "W"; //[157.5, 180]
					}else{
						return "NW"; // [112.5, 157.5)
					}
				}else{
					return "N"; // [67.5, 112.5)
				}
			}else{
				return "NE"; // [22.5, 67.5)
			}
		}else if(theta >= -rad22_5){
			return "E"; //[-22.5, 22.5)
		}else if(theta>=-3*rad22_5){
			return "SE"; //[-67.5, -22.5)
		}else if(theta>= -5*rad22_5){
			return "S"; // [-112.5, -67.5)
		}else if(theta >= -7*rad22_5){
			return "SW"; 
		}else{
			return "W";
		}
	}
	
	//
	public static void main(String[] agrs){
		double a[] = {0, 30};
		double w[] = {2, 1};
		//double a[] = {-180, 180};
		//double w[] = {1, 1};
		double wSum = 0.0;
		for(double i:w){
			wSum+=i;
		}
		double x = 0.0;
		double y = 0.0;
		for(int i = 0; i < a.length; ++i){
			y+=w[i]*Math.sin(Math.toRadians(a[i]));
			x+=w[i]*Math.cos(Math.toRadians(a[i]));
		}
		y/=wSum;
		x/=wSum;
		double g = (Math.PI/wSum) / Math.sin(Math.PI/wSum);
		double r = Math.sqrt(x*x + y*y);
		double cosA = x / r;
		double sinA = y / r;
		double theta = Math.atan2(y, x);
		double gr = g*r;
		System.out.println("x: "+x+" y: "+ y + " r: "+r+ " mean: "+ Math.toDegrees(theta) + " 1-r:" + (1-r));
		System.out.println("g: "+g + " corrected r: " + gr);
		System.out.println("Done");
	}

}
