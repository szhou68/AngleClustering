// various generic utility methods, using built-in data types only
package uk.osgb.utilities;

import java.awt.geom.Point2D;

import com.vividsolutions.jts.geom.util.AffineTransformation;

public class GeomUtility {
	// sphericity of 3D object
	/**
	 * @param volume
	 * @param surfArea
	 * @return
	 */
	public static double sphericity(double volume, double surfArea){
		return Math.cbrt(36*Math.PI*volume*volume) / surfArea;
	}
	// roundness of 2D object
	/**
	 * @param area
	 * @param perimeter
	 * @return
	 */
	public static double roundness(double area, double perimeter){
		return 2*Math.sqrt(Math.PI*area)/perimeter;
	}
	/********************************************	
	 *
	 * affine transformation 
	 *
	 *********************************************/	
	// rotate around (x0, y0) so that (x1, y1) is on x-axis (x', 0)
	/**
	 * @param x0
	 * @param y0
	 * @param x1
	 * @param y1
	 * @return
	 */
	public static AffineTransformation compAT(double x0, double y0, double x1, double y1){
		AffineTransformation at = new AffineTransformation();
		//at.rotate(x1-x0, y1-y0, x0, y0); // this be OK?
		at.translate(-x0, -y0);
		double len = Math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
		double xdiff = x1-x0;
		double ydiff = y1-y0;
		double cosTheta = xdiff / len;
		double sinTheta = - ydiff / len;
		at.rotate(sinTheta, cosTheta);
		return at;
	}

	/**
	 * @param theta
	 * @return
	 */
	public static AffineTransformation compAT(double theta){
		AffineTransformation at = new AffineTransformation();
		at.rotate(theta);
		return at;
	}
	// input: connected line segments v1{(x1,y1), (xi,yi)} and v2 {(xi,yi), (x2,y2)}
	// output: the angle turning from vector{xi-x1, yi-y1} to (x2-xi, y2-yi) 
	// the absolute angle value is in the range of [0, PI]
	// if turning CCW or collinear, positive value is returned 
	// if turning CW, NEGATIVE value is returned
	//
	public static double getVectorJoiningAngleDIR(double x1, double y1, double xi, double yi, double x2, double y2){
		double v1x = xi-x1;
		double v1y = yi-y1;
		double v2x = x2-xi;
		double v2y = y2-yi;
		// get the absolute angle
		double ang = Math.acos((v1x*v2x+v1y*v2y)/Math.sqrt((v1x*v1x+v1y*v1y)*(v2x*v2x+v2y*v2y)));
		// check direction using cross product
		return (v1x*v2y - v1y*v2x)<0.0 ? -ang:ang; 
	}
	// 
	// input: connected line segments v1{(x1,y1), (xi,yi)} and v2 {(xi,yi), (x2,y2)}
	// output: the angle between vectors{xi-x1, yi-y1} and (x2-xi, y2-yi) in the range of [0, PI]
	//
	public static double getVectorJoiningAngle(double x1, double y1, double xi, double yi, double x2, double y2){
		double v1x = xi-x1;
		double v1y = yi-y1;
		double v2x = x2-xi;
		double v2y = y2-yi;
		return Math.acos((v1x*v2x+v1y*v2y)/Math.sqrt((v1x*v1x+v1y*v1y)*(v2x*v2x+v2y*v2y)));
	}
	//
	// angle between (x1s-y1s x1e-y1e) and (x2s-y2s x2e-y2e)
	// move x2s, y2s to x1e y1e and then call 
	//
	public static double getVectorJoiningAngle(double x1s, double y1s, double x1e, double y1e, double x2s, double y2s, double x2e, double y2e){
		return getVectorJoiningAngle(x1s, y1s, x1e, y1e, x1e+x2e-x2s, y1e+y2e-y2s);
	}
	//
	// compute the length of a simple polyline represented by a double array [x0][y0][x1][y1][...]
	//
	public static double getLength(double[] simpln){
		double len = 0.0;
		if(simpln!=null){
			int coordnum = simpln.length;
			if(coordnum > 3){
				int cnt = 0;
				double x0 = simpln[cnt++];
				double y0 = simpln[cnt++];
				while(cnt < coordnum){
					double x1 = simpln[cnt++];
					double y1 = simpln[cnt++];
					double xdiff = x1-x0;
					double ydiff = y1-y0;
					len += Math.sqrt(xdiff*xdiff+ydiff*ydiff);
					x0 = x1;
					y0 = y1;
				}
			}
		}
		return len;
	}
	//
	// find the position on the polyline where the ratio between the distance from the starting vertex to the position
	// and the length of the polyline is 'ratio'
	//
	// return the index of p0(the vertex before the point)
	//
	// position of the middle point is returned via midpt
	//
	// the segment containing midpt is returned via p0 and p1
	//
	public static int getPointAlongLine(double[] simpln, double ratio, Point2D.Double midpt, Point2D.Double p0, Point2D.Double p1){
		int rtn = -1;
		if(ratio < 0.0 || ratio > 1.0){
			return rtn;
		}
		if(simpln!=null){
			double len = getLength(simpln);
			if(len>0.0){
				double intlen = len*ratio;
				int coordnum = simpln.length;
				double lentmp = 0.0;
				int cnt = 0;
				double x0 = simpln[cnt++];
				double y0 = simpln[cnt++];
				while(cnt < coordnum){
					rtn++;
					double x1 = simpln[cnt++];
					double y1 = simpln[cnt++];
					double xdiff = x1-x0;
					double ydiff = y1-y0;
					double seglen = Math.sqrt(xdiff*xdiff+ydiff*ydiff);
					if(lentmp+ seglen >= intlen){
						double k = (intlen - lentmp) / seglen;
						if(midpt!=null){
							double x = x0 + k*(x1-x0);
							double y = y0 + k*(y1-y0);
							midpt.setLocation(x, y);
							if(p0!=null){
								p0.setLocation(x0, y0);
							}
							if(p1!=null){
								p1.setLocation(x1, y1);
							}
						}
						break;
					}else{
						x0 = x1;
						y0 = y1;
					}
				}
			}
		}
		return rtn;
	}
	//
	// input vector base -> ve
	// return a new vector by turning the vector around base for 90 degree CCW or CW
	// if retunit is true,  returned vector will be normalised (i.e. length = 1.0) 
	//
	public static Point2D.Double getPerpendicularVector(double xb, double yb, double xe, double ye, boolean ccw, boolean retunit){
		double x = xe - xb;
		double y = ye - yb;
		if(retunit){
			double len = Math.sqrt(x*x + y*y);
			x /= len;
			y /= len;
		}
		if(ccw){
			return new Point2D.Double(xb-y, yb+x);
		}else{
			return new Point2D.Double(xb+y, yb-x);
		}
	}
	//
	//positive if they are CCW, negative if CW, 0.0 if collinear, 5 minus and 2 times operation
	//
	public static double triAreaDouble2( double x0, double y0, double x1, double y1, double x2, double y2){
		return (x1 - x0)*(y2 - y0) - (y1 - y0)*(x2 - x0);
	}
	// same as above, but with inverted sign, 5 plus/minus, 3 times
	public static double triAreaDouble3(double ax, double ay, double bx, double by, double cx, double cy){
		return ax*(cy - by) + bx*(ay - cy) + cx * (by - ay);
	}
}
