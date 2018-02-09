/** 
 * Various utility functions operating on JTS geometries
 * 
 */
package uk.osgb.utilities;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;


import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.util.AffineTransformation;
import com.vividsolutions.jts.geom.util.NoninvertibleTransformationException;
import com.vividsolutions.jts.math.Vector2D;
import com.vividsolutions.jts.operation.distance.DistanceOp;
import com.vividsolutions.jts.simplify.DouglasPeuckerSimplifier;

import uk.osgb.utilities.GeomUtility;
import uk.osgb.utilities.AngleUtility;

public class JTSUtility {

	public static double OSMINX = 0.0;
	public static double OSMINY = 0.0;
	public static double OSMAXX = 700000.0;
	public static double OSMAXY = 1300000.0;

	public static double rad60 = Math.PI / 3.0; // PI == 3.141592653589793
	public static double rad90 = Math.PI * 0.5;
	//public static double rad120 = 2*rad60 - 0.000000000000001;
	public static double rad120 = 2*rad60;
	public static double rad22_5 = Math.PI*0.125;
	public static double rad45 = Math.PI*0.25;
	public static double rad15 = Math.PI/12.0;
	public static double rad10 = Math.PI/18.0;
	public static double rad05 = Math.PI/36.0;

	//public static double CDTTolerance = 0.001;


	/********************************************	
	 *
	 * polygon 
	 *
	 *********************************************/	
	//
	/** longest diagonal in a polygon, return as Coordinate[2] 
	 * @param geom
	 * @return
	 */
	public static Coordinate[] compLongestDiagonal(Geometry geom){
		Geometry hull = geom.convexHull();
		String gtype = hull.getGeometryType();
		Coordinate[] diagonal = new Coordinate[2];
		if(gtype.compareToIgnoreCase("POLYGON") == 0){
			Coordinate[] coords = hull.getCoordinates(); // coords[0] equals to coords[coords.length-1]
			double clbase = -1.0;
			for(int i = 0; i < coords.length-2; ++i){
				Coordinate ci = coords[i];
				for(int j = i+1; j < coords.length-1; ++j){
					Coordinate cj = coords[j];
					double cl = ci.distance(cj);
					if(cl > clbase){
						clbase = cl;
						diagonal[0] = ci;
						diagonal[1] = cj;
					}
				}
			}
			return diagonal;
		}else if(gtype.equalsIgnoreCase("LINESTRING") || gtype.equalsIgnoreCase("LinearRing")){
			return hull.getCoordinates(); //?
		}else if(gtype.compareToIgnoreCase("POINT") == 0){
			diagonal = new Coordinate[2];
			diagonal[0] = diagonal[1] = hull.getCoordinate();
			return diagonal;
		}
		return null;
	}
	/** longest diagonal in a polygon, coordinates set to sp and ep, length is returned
	 * @param geom
	 * @param sp
	 * @param ep
	 * @return
	 */
	public static double compLongestDiagonal(Geometry geom, Coordinate sp, Coordinate ep){
		Geometry hull = geom.convexHull();
		String gtype = hull.getGeometryType();
		if(gtype.compareToIgnoreCase("POLYGON") == 0){
			Coordinate[] coords = hull.getCoordinates(); // coords[0] equals to coords[coords.length-1]
			double clbase = -1.0;
			for(int i = 0; i < coords.length-2; ++i){
				Coordinate ci = coords[i];
				for(int j = i+1; j < coords.length-1; ++j){
					Coordinate cj = coords[j];
					double cl = ci.distance(cj);
					if(cl > clbase){
						clbase = cl;
						sp.setCoordinate(ci);
						ep.setCoordinate(cj);
					}
				}
			}
			return sp.distance(ep);
		}else if(gtype.equalsIgnoreCase("LINESTRING") || gtype.equalsIgnoreCase("LinearRing")){
			Coordinate[] dig = hull.getCoordinates();
			sp.setCoordinate(dig[0]);
			ep.setCoordinate(dig[1]);
			return sp.distance(ep);
		}else if(gtype.compareToIgnoreCase("POINT") == 0){
			sp.setCoordinate(hull.getCoordinate());
			ep.setCoordinate(hull.getCoordinate());
			return 0.0;
		}
		return 0.0;
	}
	/**
	 * @param geom
	 * @param sp
	 * @param ep
	 * @return
	 */
	public static double getLongestSegment(Geometry geom, Coordinate sp, Coordinate ep){
		String gtype = geom.getGeometryType();
		if(gtype.equalsIgnoreCase("Polygon")){
			return getLongestSegment((Polygon)geom, sp, ep);
		}else if(gtype.equalsIgnoreCase("LineString") || gtype.equalsIgnoreCase("LinearRing")){
			return getLongestSegment((LineString)geom, sp, ep);
		}else{
			return 0.0;
		}
	}
	//
	/**
	 * @param plg
	 * @param sp
	 * @param ep
	 * @return
	 */
	public static double getLongestSegment(Polygon plg, Coordinate sp, Coordinate ep){
		double len = 0.0;
		// exterior ring
		LineString extRing = plg.getExteriorRing();
		len = getLongestSegment(extRing, sp, ep);
		if(plg.getNumInteriorRing() > 0){
			Coordinate sp2 = new Coordinate();
			Coordinate ep2 = new Coordinate();
			for(int i = 0; i < plg.getNumInteriorRing(); ++i){
				LineString intRing = plg.getInteriorRingN(i);
				double len2 = getLongestSegment(intRing, sp2, ep2);
				if(len2 > len){
					len = len2;
					sp.setCoordinate(sp2);
					ep.setCoordinate(ep2);
				}
			}
		}
		return len;
	}

	/** to complete: get the subset of vertices (from sidTH vertex to eidTH vertex inclusive)
	 *  in linestring or exterior ring of polgyon
	 *    
	 * @param geom LineString or Polygon
	 * @param sid
	 * @param eid
	 * @param isCCW
	 * @param seq
	 * @return
	 */
	public static int getVertexSeq(Geometry geom, int sid, int eid, boolean isCCW, Vector<Coordinate> seq){
		// to be completed
		return 0;
	}
	//
	//
	/**
	 * @param ls
	 * @param sp
	 * @param ep
	 * @return
	 */
	public static double getLongestSegment(LineString ls, Coordinate sp, Coordinate ep){
		double len = 0.0;
		Coordinate[] coords = ls.getCoordinates();
		Coordinate sc = coords[0];
		for(int i = 1; i < coords.length; ++i){
			Coordinate ec = coords[i];
			double segLen = sc.distance(ec);
			if(segLen > len){
				len = segLen;
				sp.setCoordinate(sc);
				ep.setCoordinate(ec);
			}
			sc = ec;
		}
		return len;
	}
	/**
	 * @param geom
	 * @param tol
	 * @param cor0
	 * @param cor1
	 * @param cor2
	 * @return
	 */
	public static double getLargestVerticalCornerDistSum(Geometry geom, double tol, Coordinate cor0, Coordinate cor1, Coordinate cor2){
		double distSum = 0.0;
		String gtype = geom.getGeometryType();
		if(gtype.equalsIgnoreCase("Polygon")){
			return getLargestVerticalCornerDistSum((Polygon)geom, tol, cor0, cor1, cor2);
		}else if(gtype.equalsIgnoreCase("LINESTRING") || gtype.equalsIgnoreCase("LinearRing")){
			return getLargestVerticalCornerDistSum((LineString)geom, tol, cor0, cor1, cor2);
		}else{
			return distSum;
		}
	}
	/**
	 * @param plg
	 * @param tol
	 * @param cor0
	 * @param cor1
	 * @param cor2
	 * @return
	 */
	public static double getLargestVerticalCornerDistSum(Polygon plg, double tol, Coordinate cor0, Coordinate cor1, Coordinate cor2){
		double distSum = 0.0;
		LineString extRing = plg.getExteriorRing();
		distSum = getLargestVerticalCornerDistSum(extRing, tol, cor0, cor1, cor2);
		if(plg.getNumInteriorRing() > 0){
			Coordinate cor00 = new Coordinate();
			Coordinate cor11 = new Coordinate();
			Coordinate cor22 = new Coordinate();
			for(int i = 0; i < plg.getNumInteriorRing(); ++i){
				LineString intRing = plg.getInteriorRingN(i);
				double distSum2 = getLargestVerticalCornerDistSum(intRing, tol, cor00, cor11, cor22);
				if(distSum2 > distSum){
					distSum = distSum2;
					cor0.setCoordinate(cor00);
					cor1.setCoordinate(cor11);
					cor2.setCoordinate(cor22);
				}
			}
		}
		return distSum;
	}
	//
	/**
	 * @param ls
	 * @param tol
	 * @param cor0
	 * @param cor1
	 * @param cor2
	 * @return
	 */
	public static double getLargestVerticalCornerDistSum(LineString ls, double tol, Coordinate cor0, Coordinate cor1, Coordinate cor2){
		double distSum = 0.0;
		Coordinate[] coords = ls.getCoordinates();
		if(coords.length < 3){
			return distSum;
		}
		if(coords[0].equals2D(coords[coords.length-1])){// closed
			Coordinate[] tmp = new Coordinate[coords.length+1];
			for(int i = 0; i < coords.length; ++i){
				tmp[i] = coords[i];
			}
			tmp[tmp.length-1] = coords[1];
			coords = tmp;
		}
		Coordinate c0 = coords[0];
		Coordinate c1 = coords[1];
		for(int i = 2; i < coords.length; ++i){
			Coordinate c2 = coords[i];
			double dist1 = c0.distance(c1);
			double dist2 = c1.distance(c2);
			if(dist1 > 0.0 && dist2 > 0.0){
				double ang0 = Math.atan2(c1.y - c0.y, c1.x - c0.x);
				double ang1 = Math.atan2(c2.y - c1.y, c2.x - c1.x);
				double angDiff = AngleUtility.minAngDiff(AngleUtility.minAngDiff(ang0, ang1), Math.PI*0.5);
				if(angDiff < tol){// vertical corner

					double distTotal = dist1+dist2;

					if(distTotal > distSum){
						distSum = distTotal;
						cor0.setCoordinate(c0);
						cor1.setCoordinate(c1);
						cor2.setCoordinate(c2);
					}
				}
			}
			c0 = c1;
			c1 = c2;
		}
		return distSum;
	}
	/**
	 * @param plg
	 * @param tol
	 * @param cor0
	 * @param cor1
	 * @param cor2
	 * @return
	 */
	public static double getLargestVerticalCornerDistSumWeighted(Polygon plg, double tol, Coordinate cor0, Coordinate cor1, Coordinate cor2){
		double distSum = 0.0;
		LineString extRing = plg.getExteriorRing();
		distSum = getLargestVerticalCornerDistSumWeighted(extRing, tol, cor0, cor1, cor2);
		if(plg.getNumInteriorRing() > 0){
			Coordinate cor00 = new Coordinate();
			Coordinate cor11 = new Coordinate();
			Coordinate cor22 = new Coordinate();
			for(int i = 0; i < plg.getNumInteriorRing(); ++i){
				LineString intRing = plg.getInteriorRingN(i);
				double distSum2 = getLargestVerticalCornerDistSumWeighted(intRing, tol, cor00, cor11, cor22);
				if(distSum2 > distSum){
					distSum = distSum2;
					cor0.setCoordinate(cor00);
					cor1.setCoordinate(cor11);
					cor2.setCoordinate(cor22);
				}
			}
		}
		return distSum;
	}
	/**
	 * @param geom
	 * @param tol
	 * @param cor0
	 * @param cor1
	 * @param cor2
	 * @return
	 */
	public static double getLargestVerticalCornerDistSumWeighted(Geometry geom, double tol, Coordinate cor0, Coordinate cor1, Coordinate cor2){
		double distSum = 0.0;
		String gtype = geom.getGeometryType();
		if(gtype.equalsIgnoreCase("Polygon")){
			return getLargestVerticalCornerDistSumWeighted((Polygon)geom, tol, cor0, cor1, cor2);
		}else if(gtype.equalsIgnoreCase("LINESTRING") || gtype.equalsIgnoreCase("LinearRing")){
			return getLargestVerticalCornerDistSumWeighted((LineString)geom, tol, cor0, cor1, cor2);
		}else{
			return distSum;
		}
	}
	// this one favours corner where the lengths of two edges forming the corner are more similar
	/** this one favours corner where the lengths of two edges forming the corner are more similar by comparing: dist1+dist2+dist1*dist2
	 * @param ls
	 * @param tol
	 * @param cor0
	 * @param cor1
	 * @param cor2
	 * @return
	 */
	public static double getLargestVerticalCornerDistSumWeighted(LineString ls, double tol, Coordinate cor0, Coordinate cor1, Coordinate cor2){
		double distSum = 0.0;
		Coordinate[] coords = ls.getCoordinates();
		if(coords.length < 3){
			return distSum;
		}
		if(coords[0].equals2D(coords[coords.length-1])){// closed
			Coordinate[] tmp = new Coordinate[coords.length+1];
			for(int i = 0; i < coords.length; ++i){
				tmp[i] = coords[i];
			}
			tmp[tmp.length-1] = coords[1];
			coords = tmp;
		}
		Coordinate c0 = coords[0];
		Coordinate c1 = coords[1];
		for(int i = 2; i < coords.length; ++i){
			Coordinate c2 = coords[i];
			double dist1 = c0.distance(c1);
			double dist2 = c1.distance(c2);
			if(dist1 > 0.0 && dist2 > 0.0){
				double ang0 = Math.atan2(c1.y - c0.y, c1.x - c0.x);
				double ang1 = Math.atan2(c2.y - c1.y, c2.x - c1.x);
				double angDiff = AngleUtility.minAngDiff(AngleUtility.minAngDiff(ang0, ang1), Math.PI*0.5);
				if(angDiff < tol){
					double distTotal = dist1+dist2+dist1*dist2;

					if(distTotal > distSum){
						distSum = distTotal;
						cor0.setCoordinate(c0);
						cor1.setCoordinate(c1);
						cor2.setCoordinate(c2);
					}
				}
			}
			c0 = c1;
			c1 = c2;
		}
		return distSum;
	}
	/******************************************
	 * 
	 *  Oriented Bounding Box
	 * 
	 ******************************************/
	//
	// Bounding Rectangle on longest edge
	//
	/** compute the bounding rectangle aligned to the longest segment on the geometry 
	 * @param geom 
	 * @param br
	 * @return
	 */
	public static Coordinate compBRonLE(Geometry geom, Coordinate[] br){
		Coordinate sp = new Coordinate();
		Coordinate ep = new Coordinate();
		getLongestSegment(geom, sp, ep);
		return compBR(geom, sp, ep, br);
	}
	//
	// BR on edges of the largest vertical corner
	/** computing the bounding rectangle aligned to the longer edge of the largest vertical corner of the geometry
	 * @param geom 
	 * @param tol tolerance (in radian) to decide if a corner is vertical
	 * @param br four corners of the resulting bounding rectangle
	 * @return width (x) and height (y) of the bounding rectangle are returned as a coordinate object
	 */
	public static Coordinate compBRonVC(Geometry geom, double tol, Coordinate[] br){
		Coordinate cor0 = new Coordinate();
		Coordinate cor1 = new Coordinate();
		Coordinate cor2 = new Coordinate();
		double sum = getLargestVerticalCornerDistSum(geom, tol, cor0, cor1, cor2);
		//
		if(sum > 0.0){
			double dist1 = cor0.distance(cor1);
			double dist2 = cor1.distance(cor2);
			if(dist1>=dist2){
				return compBR(geom, cor0, cor1, br);
			}else{
				return compBR(geom, cor1, cor2, br);
			}
		}else{ // no vertical corner
			return compBRonLE(geom, br);
		}
	}
	//
	/** BR on edges of the largest vertical corner (weighted)
	 * @param geom
	 * @param tol
	 * @param br
	 * @return
	 */
	public static Coordinate compBRonVCWeighted(Geometry geom, double tol, Coordinate[] br){
		Coordinate cor0 = new Coordinate();
		Coordinate cor1 = new Coordinate();
		Coordinate cor2 = new Coordinate();
		double sum = getLargestVerticalCornerDistSumWeighted(geom, tol, cor0, cor1, cor2);
		if(sum > 0.0){
			//
			double dist1 = cor0.distance(cor1);
			double dist2 = cor1.distance(cor2);
			if(dist1>=dist2){
				return compBR(geom, cor0, cor1, br);
			}else{
				return compBR(geom, cor1, cor2, br);
			}
		}else{// no vertical corner
			return compBRonLE(geom, br);
		}
	}
	//
	// BR along the lonest diagonal
	//
	/** compute the bounding rectangle aligned to the longest diagonal of the convex hull of the geometry
	 * @param geom
	 * @param br
	 * @return
	 */
	public static Coordinate compBRonLD(Geometry geom, Coordinate[] br){
		Coordinate[] diagonal = compLongestDiagonal(geom);
		return compBR(geom, diagonal[0], diagonal[1], br);

	}
	//
	// bounding box along a vector sp-ep
	// use transformation defined by sp-ep to compute MBR 
	//
	/** compute bounding box along a vector sp--ep
	 * @param geom
	 * @param sp
	 * @param ep
	 * @param br
	 * @return
	 */
	public static Coordinate compBR(Geometry geom, Coordinate sp, Coordinate ep, Coordinate[] br){
		return compBR(geom, sp.x, sp.y, ep.x, ep.y, br);
	}
	//
	public static Coordinate compBR(Geometry geom, double axisOrient, Coordinate[] br){
		Coordinate[] coords = geom.getCoordinates();
		if(coords==null || coords.length == 0){
			return null;
		}
		Coordinate sp = coords[0];
		double ex = sp.x + Math.cos(axisOrient);
		double ey = sp.y + Math.sin(axisOrient);
		return compBR(geom, sp.x, sp.y, ex, ey, br);
	}
	//
	/** compute bounding box along a vector (sx, sy) -- (ex, ey)
	 * @param geom
	 * @param sx
	 * @param sy
	 * @param ex
	 * @param ey
	 * @param br FOUR corners of the bounding box
	 * @return
	 */
	public static Coordinate compBR(Geometry geom, double sx, double sy, double ex, double ey, Coordinate[] br){
		AffineTransformation at = GeomUtility.compAT(sx, sy, ex, ey);
		Envelope e = null;
		Coordinate dest = new Coordinate(0.0, 0.0);
		for(Coordinate ver:geom.getCoordinates()){
			at.transform(new Coordinate(ver.x, ver.y), dest);
			if(e!=null){
				e.expandToInclude(dest);
			}else{
				e = new Envelope(dest.x, dest.x, dest.y, dest.y);
			}
		}
		Coordinate extent = null;
		if(e!=null){
			extent = new Coordinate(e.getWidth(), e.getHeight());
			if(br[0]!=null){
				br[0].x = e.getMinX();
				br[0].y = e.getMinY();
			}else{
				br[0] = new Coordinate(e.getMinX(), e.getMinY());
			}
			if(br[1]!=null){
				br[1].x = e.getMaxX();
				br[1].y = e.getMinY();
			}else{
				br[1] = new Coordinate(e.getMaxX(), e.getMinY());
			}
			if(br[2]!=null){
				br[2].x = e.getMaxX();
				br[2].y = e.getMaxY();
			}else{
				br[2] = new Coordinate(e.getMaxX(), e.getMaxY());
			}
			if(br[3]!=null){
				br[3].x = e.getMinX();
				br[3].y = e.getMaxY();
			}else{
				br[3] = new Coordinate(e.getMinX(), e.getMaxY());
			}
			try {
				AffineTransformation atInv = at.getInverse();
				atInv.transform(br[0], br[0]);
				atInv.transform(br[1], br[1]);
				atInv.transform(br[2], br[2]);
				atInv.transform(br[3], br[3]);
			} catch (NoninvertibleTransformationException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
		return extent;
	}
	// to be completed
	// BR align to the pair of furtherest parallel edges (with length restriction)
	private static Coordinate compBRParallel(Geometry geom, double minEdgeLen, Coordinate[] br){
		return null;
	}
	//
	// to be completed
	/** compute the largest interiror rectangle aligned to X-Y axis 
	 * @param geom
	 * @param cor1
	 * @param cor2
	 * @param cor3
	 * @return
	 */
	private static Coordinate compMaxIntRectXY(Geometry geom, Coordinate cor1, Coordinate cor2,Coordinate cor3){
		return null;
	}
	//
	// to be completed
	/** compute the largest interior rectangle aligned to vector sp-ep
	 * @param geom
	 * @param sp
	 * @param ep
	 * @param cor1
	 * @param cor2
	 * @param cor3
	 * @return
	 */
	private static Coordinate compMaxIntRectAligned(Geometry geom, Coordinate sp, Coordinate ep, Coordinate cor1, Coordinate cor2,Coordinate cor3){
		return null;
	}
	//
	// transform centroid of geometry to coordinate origin and the longer side the x-axis
	//
	/**
	 * @param geom
	 * @param mmbr
	 * @return
	 */
	public static Geometry normaliseGeometryByMMBR(Geometry geom, Geometry mmbr){
		Coordinate[] coords = mmbr.getCoordinates();
		Coordinate centroid = new Coordinate(0.5*(coords[0].x + coords[2].x), 0.5*(coords[0].y + coords[2].y));
		Coordinate endPt;
		double dist01 = coords[0].distance(coords[1]);
		double dist12 = coords[1].distance(coords[2]);
		if(dist01>=dist12){
			endPt = new Coordinate(0.5*(coords[1].x + coords[2].x), 0.5*(coords[1].y + coords[2].y));
		}else{
			endPt = new Coordinate(0.5*(coords[2].x + coords[3].x), 0.5*(coords[2].y + coords[3].y));
		}
		AffineTransformation at = GeomUtility.compAT(centroid.x, centroid.y, endPt.x, endPt.y);
		return at.transform(geom);
	}
	//
	/** move and rotate a geometry to the centerid of its aligned minimum bounding rectangle while the longer edges parallel to x-axis 
	 * @param geom
	 * @param minEdgeLen
	 * @param useCH
	 * @return
	 */
	public static Geometry normaliseGeometryByMMBR(Geometry geom, double minEdgeLen, boolean useCH){
		String gtype;
		Geometry hull = null;
		if(useCH){
			hull = geom.convexHull();
			gtype = hull.getGeometryType();
		}else{
			gtype = geom.getGeometryType();
		}
		// if the convex hull contains 3 or more points, a Polygon; 2 points, a LineString; 1 point, a Point; 0 points, an empty GeometryCollection.

		if(gtype.compareToIgnoreCase("POLYGON") == 0){
			AffineTransformation at = null;
			double area = Double.MAX_VALUE;
			Envelope minE = null;
			AffineTransformation minAT = null;
			Coordinate[] vertices = useCH?hull.getCoordinates():geom.getBoundary().getCoordinates();
			Coordinate v0, v1;
			for(int i = 0; i < vertices.length-1;++i){
				v0 = vertices[i];
				v1 = vertices[i+1];
				double len = v0.distance(v1);
				if(len==0.0 || len < minEdgeLen){
					continue;
				}
				// move origin to v0 and rotate so that v0-v1 is on X-axis
				at = GeomUtility.compAT(v0.x, v0.y, v1.x, v1.y);
				Envelope e = new Envelope(0.0,0.0,0.0,0.0);
				Coordinate dest = new Coordinate(0.0, 0.0);
				for(Coordinate ver:vertices){
					at.transform(new Coordinate(ver.x, ver.y), dest);
					e.expandToInclude(dest);
				}
				double ea = e.getArea();
				if(ea<area){
					area = ea;
					minE = e;
					minAT = at;
				}
			}
			if(minE!=null){
				try {
					Coordinate c0 = new Coordinate((minE.getMinX()+minE.getMaxX())*0.5, (minE.getMinY()+minE.getMaxY())*0.5);
					Coordinate c1;
					if(minE.getWidth() >= minE.getHeight()){
						c1 = new Coordinate(minE.getMaxX(), (minE.getMinY()+minE.getMaxY())*0.5);
					}else{
						c1 = new Coordinate((minE.getMinX()+minE.getMaxX())*0.5, minE.getMaxY());
					}
					AffineTransformation minATInv = minAT.getInverse();
					//
					minATInv.transform(c0, c0);
					minATInv.transform(c1, c1);
					//
					at = GeomUtility.compAT(c0.x, c0.y, c1.x, c1.y);
					return at.transform(geom);
				} catch (NoninvertibleTransformationException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				return null;
			}else{
				return null;
			}
		}else if(gtype.compareToIgnoreCase("POINT")==0){
			return geom;
		}
		return null;	
	}
	/***********************************************
	 * 
	 * 
	 *  Minimum Minimum Bounding Rectangle
	 * 
	 * 
	 ***********************************************/
	//
	// compute the mmbr on a give orientation defined by sinTheta and cosTheta
	//
	/**
	 * @param geom
	 * @param sp
	 * @param ep
	 * @param cor0
	 * @param cor1
	 * @param cor2
	 * @param cor3
	 * @return
	 */
	public static Coordinate compMinMBR(Geometry geom, Coordinate sp, Coordinate ep, Coordinate cor0, Coordinate cor1, Coordinate cor2,Coordinate cor3){
		AffineTransformation at = GeomUtility.compAT(sp.x, sp.y, ep.x, ep.y);
		Envelope env = null;
		Coordinate[] coords = geom.getCoordinates();
		Coordinate dest = new Coordinate(0.0, 0.0);
		for(Coordinate coord:coords){
			at.transform(coord, dest);
			if(env !=null){
				env.expandToInclude(dest);
			}else{
				env = new Envelope(dest, dest);
			}
		}
		if(env!=null){
			try {
				Coordinate c0 = new Coordinate(env.getMinX(), env.getMinY());
				Coordinate c1 = new Coordinate(env.getMaxX(), env.getMinY());
				Coordinate c2 = new Coordinate(env.getMaxX(), env.getMaxY());
				Coordinate c3 = new Coordinate(env.getMinX(), env.getMaxY());

				AffineTransformation atInv = at.getInverse();

				atInv.transform(c0, c0);
				atInv.transform(c1, c1);
				atInv.transform(c2, c2);
				atInv.transform(c3, c3);

				if(cor0!=null){
					cor0.setCoordinate(c0);
				}
				if(cor1!=null){
					cor1.setCoordinate(c1);
				}
				if(cor2!=null){
					cor2.setCoordinate(c2);
				}
				if(cor3!=null){
					cor3.setCoordinate(c3);
				}
			} catch (NoninvertibleTransformationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return new Coordinate(env.getWidth(), env.getHeight());
		}
		return null;
	}
	//
	// set four corners of mmbr to cor0-cor3 and return width and height of the mmbr as a coordinate
	//
	/**
	 * @param geom
	 * @param minEdgeLen
	 * @param useCH
	 * @param cor0
	 * @param cor1
	 * @param cor2
	 * @param cor3
	 * @return
	 */
	public static Coordinate compMinMBR(Geometry geom, double minEdgeLen, boolean useCH, Coordinate cor0, Coordinate cor1, Coordinate cor2,Coordinate cor3){
		//
		String gtype;
		Geometry hull = null;
		if(useCH){
			hull = geom.convexHull();
			gtype = hull.getGeometryType();
		}else{
			gtype = geom.getGeometryType();
		}
		// if the convex hull contains 3 or more points, a Polygon; 2 points, a LineString; 1 point, a Point; 0 points, an empty GeometryCollection.

		if(gtype.compareToIgnoreCase("POLYGON") == 0){
			AffineTransformation at = null;
			double area = Double.MAX_VALUE;
			Envelope minE = null;
			AffineTransformation minAT = null;
			Coordinate[] vertices = useCH?hull.getCoordinates():geom.getBoundary().getCoordinates();
			Coordinate v0, v1;
			for(int i = 0; i < vertices.length-1;++i){
				v0 = vertices[i];
				v1 = vertices[i+1];
				double len = v0.distance(v1);
				if(len==0.0 || len < minEdgeLen){
					continue;
				}
				// move origin to v0 and rotate so that v0-v1 is on X-axis
				at = GeomUtility.compAT(v0.x, v0.y, v1.x, v1.y);
				Envelope e = new Envelope(0.0,0.0,0.0,0.0);
				Coordinate dest = new Coordinate(0.0, 0.0);
				for(Coordinate ver:vertices){
					at.transform(new Coordinate(ver.x, ver.y), dest);
					e.expandToInclude(dest);
				}
				double ea = e.getArea();
				if(ea<area){
					area = ea;
					minE = e;
					minAT = at;
				}
			}
			if(minE!=null){
				try {
					Coordinate c0 = new Coordinate(minE.getMinX(), minE.getMinY());
					Coordinate c1 = new Coordinate(minE.getMaxX(), minE.getMinY());
					Coordinate c2 = new Coordinate(minE.getMaxX(), minE.getMaxY());
					Coordinate c3 = new Coordinate(minE.getMinX(), minE.getMaxY());

					AffineTransformation minATInv = minAT.getInverse();

					minATInv.transform(c0, c0);
					minATInv.transform(c1, c1);
					minATInv.transform(c2, c2);
					minATInv.transform(c3, c3);

					if(cor0!=null){
						cor0.setCoordinate(c0);
					}
					if(cor1!=null){
						cor1.setCoordinate(c1);
					}
					if(cor2!=null){
						cor2.setCoordinate(c2);
					}
					if(cor3!=null){
						cor3.setCoordinate(c3);
					}
				} catch (NoninvertibleTransformationException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				return new Coordinate(minE.getWidth(), minE.getHeight());
			}else{
				return null;
			}
		}else if(gtype.equalsIgnoreCase("LINESTRING")  ||gtype.equalsIgnoreCase("LinearRing")){// only set cor0 and cor1
			Coordinate[] coords = geom.getCoordinates();
			if(coords[0].equals2D(coords[coords.length-1])){
				Polygon plg = geom.getFactory().createPolygon(coords);
				return compMinMBR(plg, minEdgeLen, useCH, cor0, cor1, cor2, cor3);
			}
			Envelope env = geom.getEnvelopeInternal();
			if(cor0!=null){
				cor0.x = env.getMinX();
				cor0.y = env.getMinX();
			}
			if(cor1!=null){
				cor1.x = env.getMaxX();
				cor1.y = env.getMaxY();
			}
			return new Coordinate(cor0.distance(cor1), 0.0);
		}else if(gtype.compareToIgnoreCase("POINT")==0){
			cor0.x = geom.getCoordinate().x;
			cor0.y = geom.getCoordinate().y;
			return new Coordinate(0,0);
		}
		return null;
	}
	//
	// 
	//
	// compute the minimum MBR of JTS Geom and return an array of 4 coordinates in CCW order
	//
	/**
	 * @param geom
	 * @param useCH
	 * @return
	 */
	public static Coordinate[] compMinMBRCoord(Geometry geom, boolean useCH){
		Coordinate[] coords = new Coordinate[4];
		for(int i = 0; i < coords.length;++i){
			coords[i] = new Coordinate();
		}
		Coordinate dim = compMinMBR(geom, useCH, coords[0], coords[1], coords[2], coords[3]);
		if(dim!=null){
			if(dim.x == 0.0){// point
				coords[1].setCoordinate(coords[0]);
				coords[2].setCoordinate(coords[0]);
				coords[3].setCoordinate(coords[0]);
			}else if(dim.y == 0.0){
				coords[2].setCoordinate(coords[1]);
				coords[3].setCoordinate(coords[0]);
			}
			return coords;
		}else{
			return null;
		}
	}
	//
	// find the minimum mbr and return the width and height of the mmbr
	// this implementation is NOT the fastest. A possible improvement may be to compute all the affine transformation first and remove the redundant ones (prior to compute AFT, normanlise to: v0x < v1x or v0x == v1x and v0y < v1y)
	// if useCH is false, all segments on the boundary will be used for testing - more desirable for building?
	//
	/**
	 * @param geom
	 * @param useCH
	 * @param cor0
	 * @param cor1
	 * @param cor2
	 * @param cor3
	 * @return
	 */
	public static Coordinate compMinMBR(Geometry geom, boolean useCH, Coordinate cor0, Coordinate cor1, Coordinate cor2,Coordinate cor3){
		return compMinMBR(geom, 0.0, useCH,  cor0, cor1, cor2, cor3);
	}

	/**
	 * @param geom
	 * @param useCH
	 * @return
	 */
	public static Geometry compMinMBRGeom(Geometry geom, boolean useCH){
		Coordinate[] coords = new Coordinate[5];
		for(int i = 0; i < coords.length;++i){
			coords[i] = new Coordinate();
		}
		Coordinate dim = compMinMBR(geom, useCH, coords[0], coords[1], coords[2], coords[3]);
		GeometryFactory gf = geom.getFactory();
		if(dim!=null){
			if(dim.x == 0.0){
				coords[1] = coords[2] = coords[3]  = coords[0];
			}else if(dim.y == 0){// linestring
				coords[2] = coords[1];
				coords[3] = coords[0];
			}
			coords[4].setCoordinate(coords[0]);
			return gf.createPolygon(coords);//?
		}else{
			return null;
		}
	}		
	/*******************
	 * 	shape indicator
	 ******************/
	//
	// orthogonality of 2D polygons: Area(CV)/2*Area(MMBR(geom) in [0, 1] - 0 for circles and 1 for rectangles???
	//
	//
	private static double orthogonality(Geometry geom){
		//
		String gtype = geom.getGeometryType();
		Coordinate[] coords = null;
		if(gtype.compareToIgnoreCase("POLYGON") == 0){
			coords = geom.getBoundary().getCoordinates();
		}else if(gtype.equalsIgnoreCase("LINESTRING")  ||gtype.equalsIgnoreCase("LinearRing")){
			coords = geom.getBoundary().getCoordinates();
			if(!coords[0].equals2D(coords[coords.length-1])){
				return 0.0;
			}
		}else{
			return 0.0;
		}
		//
		Coordinate sp = coords[0];
		Coordinate[] varray = new Coordinate[coords.length-1];
		for(int i = 1; i < coords.length-2; ++i){
			Coordinate ep = coords[i];
			double x = ep.x - sp.x;
			double y = ep.y - sp.y;


		}
		//
		Coordinate[] corners = new Coordinate[4];
		for(int i = 0; i < 4; ++i){
			corners[i] = new Coordinate();
		}
		Coordinate mmbrDim = compMinMBR(geom, false, corners[0], corners[1], corners[2], corners[3]);
		double areaMMBR = mmbrDim.x * mmbrDim.y;

		return 1.0;
	}
	/*************************************
	 * 
	 *  orientation relating to MMBR-based sides
	 * 
	 * 
	 * 	
	 **************************************/
	/**
	 * check which "quadrant" pt falls into regarding coords
	 *  
	 * The quadrants are defined by the medial axis of the rectangle
	 * 
	 * coords is FOUR coordinates c0/c1/c2/c3 forming a rectangle in CCW order
	 * return 0 for c0-c1 side; 1 for c1-c2 side, 2 for c2-c3 side and 3 for c3-c0 side
	 *
	 * @param coords Coordinate[4] which 
	 * @param pt
	 * @return
	 */
	public static int computeOrientQuadrant(Coordinate[] coords, Coordinate pt){
		double dist01 = coords[0].distance(coords[1]);
		double dist12 = coords[2].distance(coords[1]);
		Coordinate p0, p1;
		if(dist01 >=dist12){// c0-c1 is the longer side
			double x01 = coords[0].x + (coords[1].x - coords[0].x)*dist12/dist01;
			double y01 = coords[0].y + (coords[1].y - coords[0].y)*dist12/dist01;
			// first internal pt
			double x0 = (x01 + coords[3].x)*0.5;
			double y0 = (y01 + coords[3].y)*0.5;
			p0 = new Coordinate(x0, y0);
			//
			double x10 = coords[1].x - (coords[1].x - coords[0].x)*dist12/dist01;
			double y10 = coords[1].y - (coords[1].y - coords[0].y)*dist12/dist01;
			//second interpt:
			double x1 = (x10 + coords[2].x)*0.5;
			double y1 = (y10 + coords[2].y)*0.5;
			p1 = new Coordinate(x1, y1);
			if(areaDoubleTriangle(pt, coords[0], p0) >=0){
				if(areaDoubleTriangle(pt, coords[1], p1) < 0){
					return 0;
				}else{
					if(areaDoubleTriangle(pt, coords[2], p1) < 0){
						return 1;
					}else{
						return 2;
					}
				}
			}else{
				if(areaDoubleTriangle(pt, coords[3], p0) >=0){
					return 3;
				}else{
					return 2;
				}
			}
		}else{// dist01<=dist12
			double x12 = coords[1].x + (coords[2].x - coords[1].x) * dist01/dist12;
			double y12 = coords[1].y + (coords[2].y - coords[1].y) * dist01/dist12;
			double x0 = (coords[0].x + x12)*0.5;
			double y0 = (coords[0].y + y12)*0.5;
			p0 = new Coordinate(x0, y0);
			double x21 = coords[2].x - (coords[2].x - coords[1].x)*dist01/dist12;
			double y21 = coords[2].y - (coords[2].y - coords[1].y)*dist01/dist12;
			double x1 = (coords[3].x + x21)*0.5;
			double y1 = (coords[3].y + y21)*0.5;
			p1 = new Coordinate(x1, y1);
			if(areaDoubleTriangle(pt, coords[0], p0) >=0){
				if(areaDoubleTriangle(pt, coords[1], p0) < 0){
					return 0;
				}else{
					if(areaDoubleTriangle(pt, coords[2], p1) < 0){
						return 1;
					}else{
						return 2;
					}
				}
			}else{
				if(areaDoubleTriangle(pt, coords[3], p1) >= 0){
					return 3;
				}else{
					return 2;
				}
			}
		}
	}


	//
	// find the shortest vector between geom1 and a collection of other geometries
	//
	/**
	 * @param geom1
	 * @param geometries
	 * @return
	 */
	public static Coordinate[] compShortestVector(Geometry geom1, Collection geometries){
		double minD = Double.MAX_VALUE;
		Iterator iter = geometries.iterator();
		Coordinate[] rlt = null;
		while(iter.hasNext()){
			Geometry geom2 = (Geometry)iter.next();
			Coordinate[] vec = DistanceOp.nearestPoints(geom1, geom2);
			double dist = vec[0].distance(vec[1]);
			if(dist < minD){
				minD = dist;
				rlt = vec;
			}
		}
		return rlt;
	}
	//
	// assume pt is on interior of ls (i.e. does not overlapping two endpoints 
	//

	public static LineString[] splitLineString(LineString ls, Coordinate coord, double tol){
		Coordinate[] coords = ls.getCoordinates();
		Point pt = ls.getFactory().createPoint(coord);
		Coordinate[] seg = new Coordinate[2];
		seg[0] = coords[0];
		int i = 1;
		boolean onVertex = false;
		for(; i < coords.length;++i){
			seg[1] = coords[i];
			if(coord.equals2D(seg[1])){
				onVertex = true;
				break;
			}
			LineString lnSeg = ls.getFactory().createLineString(seg);
			double dist = lnSeg.distance(pt);
			if(dist < tol){
				break;
			}
			seg[0] = seg[1];
		}
		if(i==coords.length){
			return null;
		}
		int numPt1 = i+1;
		int numPt2 = onVertex?coords.length - i:coords.length-i+1;
		Coordinate[] ls1 = new Coordinate[numPt1];
		Coordinate[] ls2 = new Coordinate[numPt2];
		if(onVertex){
			for(int j = 0; j < numPt1; ++j){
				ls1[j] = new Coordinate(coords[j]);
			}
			for(int j = 0; j < numPt2; ++j){
				ls2[j] = new Coordinate(coords[i++]);
			}
		}else{
			for(int j = 0; j < numPt1-1; ++j){
				ls1[j] = new Coordinate(coords[j]);
			}
			ls1[numPt1-1] = new Coordinate(coord);
			ls2[0] = new Coordinate(coord);
			for(int j = 1; j < numPt2; ++j){
				ls2[j] = new Coordinate(coords[i++]);
			}
		}
		LineString[] rlt = new LineString[2];
		rlt[0] = ls.getFactory().createLineString(ls1);
		rlt[1] = ls.getFactory().createLineString(ls2);
		return rlt;
	}		
	/*************************
	 *  construction
	 ************************* */
	//
	// building a linestring or polygon, COPYing coords
	//
	public static Geometry coords2Geom(Coordinate[] coords, GeometryFactory gf, boolean isPlg){
		Coordinate[] coordArray = null;
		if(isPlg){
			int sz = coords.length;
			boolean isClosed = false;
			if(coords[0].equals2D(coords[sz-1])){
				isClosed = true;
			}
			if(isClosed){
				coordArray = new Coordinate[sz];
			}else{
				coordArray = new Coordinate[sz+1];
			}
			for(int i = 0; i < sz;++i){
				coordArray[i] = new Coordinate(coords[i]);
			}
			if(!isClosed){
				coordArray[sz] = new Coordinate(coordArray[0]);
			}
			return gf.createPolygon(coordArray);
		}else{
			int sz = coords.length;
			coordArray = new Coordinate[sz];
			for(int i = 0; i < sz;++i){
				coordArray[i] = new Coordinate(coords[i]);
			}
			return gf.createLineString(coordArray);
		}
	}
	//
	/**
	 * @param ls
	 * @param maxSegLen
	 * @param preserveVer
	 * @return
	 */
	public static Collection densify(LineString ls, double maxSegLen, boolean preserveVer){
		List rlt = new ArrayList();
		Coordinate[] coords = ls.getCoordinates();
		Coordinate sp, ep;
		int numVer = ls.getNumPoints();
		if(preserveVer){
			for(int i = 0; i < numVer -1; ++i){
				sp = coords[i];
				ep = coords[i+1];
				rlt.add(sp);
				compSegmentation(sp, ep, 0.0, maxSegLen, true, rlt);
			}
			//if(coords[0].compareTo(coords[numVer-1])!=0){
			rlt.add(coords[numVer-1]);
			//}
		}else{
			rlt.add(ls.getCoordinateN(0));
			double slen = 0.0;
			for(int i = 0; i < numVer-1; ++i){
				sp = coords[i];
				ep = coords[i+1];
				slen = compSegmentation(sp, ep, slen, maxSegLen, false, rlt);
			}
			//if(coords[0].compareTo(coords[numVer-1])!=0){
			rlt.add(coords[numVer-1]);
			//}
		}
		return rlt;
	}
	//
	// vertices are preserved by default
	//
	/**
	 * @param ls
	 * @param maxSegLen
	 * @param minSegNum
	 * @return
	 */
	public static Collection densify(LineString ls, double maxSegLen, int minSegNum){
		List rlt = new ArrayList();
		Coordinate[] coords = ls.getCoordinates();
		Coordinate sp, ep;
		int numVer = ls.getNumPoints();
		for(int i = 0; i < numVer -1; ++i){
			sp = coords[i];
			ep = coords[i+1];
			//rlt.add(sp);
			//compSegmentation(sp, ep, 0.0, maxSegLen, true, rlt);
			densifySegment(sp, ep, maxSegLen, minSegNum, false, rlt);
		}
		rlt.add(coords[numVer-1]);
		return rlt;
	}
	//
	// newly created segmentation points are returned via coords
	//
	/**
	 * @param sp
	 * @param ep
	 * @param startLen
	 * @param segLen
	 * @param preserveVer
	 * @param coords
	 * @return
	 */
	public static double compSegmentation(Coordinate sp, Coordinate ep, double startLen, double segLen, boolean preserveVer, List coords){
		double len = sp.distance(ep);
		double dx = ep.x - sp.x;
		double dy = ep.y - sp.y;
		if(preserveVer){// ignore startLen, segLen is the maximum segment length, only generated internal vertices are added to coords list
			double slen = len;
			if(slen <= segLen)
				return 0.0;
			int numSeg = 2;
			for(; ;){
				slen  = len / numSeg;
				if(slen <= segLen)
					break;
				else
					++numSeg;
			}
			double sdx = dx / numSeg;
			double sdy = dy / numSeg;
			double x = sp.x;
			double y = sp.y;
			for(int i = 1; i < numSeg; ++i){
				x += sdx;
				y += sdy;
				//System.out.println(x+", "+y);
				coords.add(new Coordinate(x, y));
			}
		}else{ // if startLen == 0.0, sp will be added to coords but ep will never be added
			if(segLen > 0.0){
				double locLen = startLen;
				if(locLen == 0.0){
					locLen+=segLen;
				}
				while(locLen < len){
					double ratio = locLen / len;
					double x = sp.x + dx*ratio;
					double y = sp.y + dy*ratio;
					//System.out.println(x+", "+y);
					coords.add(new Coordinate(x, y));
					locLen+=segLen;
				}
				return locLen - len;
			}else{
				return startLen;
			}
		}
		return 0.0;
	}
	//newly created segmentation points are returned via coords
	/**
	 * @param sp
	 * @param ep
	 * @param maxSegLen
	 * @param addEP
	 * @param coords
	 * @return
	 */
	public static double densifySegment(Coordinate sp, Coordinate ep, double maxSegLen, boolean addEP, List coords){
		double len = sp.distance(ep);
		if(len <= maxSegLen){
			return len;
		}
		int segNum = (int)(len / maxSegLen)+1;
		double segLen = len / (double)segNum;
		double xstep = (ep.x - sp.x)/(double)segNum;
		double ystep = (ep.y - sp.y)/(double)segNum;
		if(addEP){
			coords.add(sp);
		}
		for(int i = 1; i < segNum;++i){
			double x = sp.x + xstep*(double)i;
			double y = sp.y + ystep*(double)i;
			coords.add(new Coordinate(x, y));
		}
		if(addEP){
			coords.add(ep);
		}
		//
		return segLen;
	}
	////newly created segmentation points are returned via coords
	/**
	 * @param sp
	 * @param ep
	 * @param segNum
	 * @param addEP
	 * @param coords
	 * @return
	 */
	public static double densifySegment(Coordinate sp, Coordinate ep, int segNum, boolean addEP, List coords){
		double len = sp.distance(ep);
		if(segNum <=1){
			return len;
		}
		double segLen = len / (double)segNum;
		double xstep = (ep.x - sp.x)/(double)segNum;
		double ystep = (ep.y - sp.y)/(double)segNum;
		if(addEP){
			coords.add(sp);
		}
		for(int i = 1; i < segNum;++i){
			double x = sp.x + xstep*(double)i;
			double y = sp.y + ystep*(double)i;
			coords.add(new Coordinate(x, y));
		}
		if(addEP){
			coords.add(ep);
		}
		//
		return segLen;
	}
	//newly created segmentation points are returned via coords
	/**
	 * @param sp
	 * @param ep
	 * @param maxSegLen
	 * @param minSegNum
	 * @param addEP
	 * @param coords
	 * @return
	 */
	public static double densifySegment(Coordinate sp, Coordinate ep, double maxSegLen, int minSegNum, boolean addEP, List coords){
		double len = sp.distance(ep);
		if(len <= maxSegLen){//
			return densifySegment(sp, ep, minSegNum, addEP, coords);
		}
		int segNum = (int)(len / maxSegLen)+1;
		if(segNum < minSegNum){
			return densifySegment(sp, ep, minSegNum, addEP, coords);
		}
		double segLen = len / (double)segNum;
		double xstep = (ep.x - sp.x)/(double)segNum;
		double ystep = (ep.y - sp.y)/(double)segNum;
		if(addEP){
			coords.add(sp);
		}
		for(int i = 1; i < segNum;++i){
			double x = sp.x + xstep*(double)i;
			double y = sp.y + ystep*(double)i;
			coords.add(new Coordinate(x, y));
		}
		if(addEP){
			coords.add(ep);
		}
		//
		return segLen;
	}

	//
	/**
	 * @param geom
	 * @return
	 */
	public static Geometry rmvCollinearVertices(Geometry geom){
		return DouglasPeuckerSimplifier.simplify(geom, 0.0);
	}
	// get a segment sp-ep on exterior boundary of the polygon where polygon interior is on the LEFT of sp-ep 
	/**
	 * @param plg
	 * @param sp
	 * @param ep
	 * @return
	 */
	public static int getExteriorSeg(Polygon plg, Coordinate sp, Coordinate ep){
		Coordinate[] coords = plg.getExteriorRing().getCoordinates();
		Coordinate p1 = coords[0], p2 = coords[1];
		for(int i = 1; i < coords.length-1; ++i){
			if(p1.equals2D(p2)){
				p1 = coords[i];
				p2 = coords[i+1];
			}else{
				break;
			}
		}
		int ori = JTSUtility.orientation(coords); 
		if( ori == -1){//ccw
			sp.setCoordinate(p1);
			ep.setCoordinate(p2);
		}else{// cw
			sp.setCoordinate(p2);
			ep.setCoordinate(p1);
		}
		return ori;
	}
	//
	/**
	 * @param gc
	 * @param point
	 * @param excludePt
	 * @param excludeLs
	 * @param excludePlg
	 * @return
	 */
	public static Geometry getNearestPart(GeometryCollection gc, Point point, boolean excludePt, boolean excludeLs, boolean excludePlg){
		double distN = Double.MAX_VALUE;
		Geometry rlt = null;
		for(int i = 0; i < gc.getNumGeometries(); ++i){
			Geometry geom = gc.getGeometryN(i);
			String gtype = geom.getGeometryType();
			if(excludePt && gtype.equalsIgnoreCase("Point")){
				continue;
			}else if(excludeLs && gtype.equalsIgnoreCase("LineString")){
				continue;
			}else if(excludePlg && gtype.equalsIgnoreCase("Polygon")){
				continue;
			}else if(gtype.equals("GeometryCollection")){
				geom = getNearestPart((GeometryCollection)geom, point, excludePt, excludeLs, excludePlg);
			}
			double dist = geom.distance(point);
			if(dist < distN){
				distN = dist;
				rlt = geom;
			}
		}
		return rlt;
	}
	//
	public static Geometry rect2Geometry(double minx, double miny, double maxx, double maxy, GeometryFactory gf){
		Coordinate[] coords = new Coordinate[5];
		coords[0] = new Coordinate(minx, miny);
		coords[1] = new Coordinate(maxx, miny);
		coords[2] = new Coordinate(maxx, maxy);
		coords[3] = new Coordinate(minx, maxy);
		coords[4] = new Coordinate(minx, miny);
		return gf.createPolygon(coords);
	}
	//
	public static Geometry rect2GeometryDim(double width, double height, double cenx, double ceny, GeometryFactory gf) {
		double halfW = width*0.5;
		double halfH = height*0.5;
		return rect2Geometry(cenx - halfW, ceny - halfH, cenx + halfW, ceny + halfH, gf);
	}

	/*****************************
	 * geometric computation
	 ******************************/
	// angle of line seg in rad
	public static double compAngle(Coordinate sp, Coordinate ep){
		double xdiff = ep.x - sp.x;
		double ydiff = ep.y - sp.y;
		return Math.atan2(ydiff, xdiff);
	}

	// knowing the two line segments are intersecting
	/**
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @param x3
	 * @param y3
	 * @param x4
	 * @param y4
	 * @return
	 */
	public static Coordinate compLineIntersectionDir(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
		double xa = x2 - x1;
		double ya = y2 - y1;
		double xb = x3 - x4;
		double yb = y3 - y4;
		double xc = x1 - x3;
		double yc = y1 - y3;
		double den = ya*xb - xa*yb;
		if(den == 0.0){// collinear
			return null;
		}
		double num1 = yb*xc - xb*yc;
		double t1 = num1 / den;
		return new Coordinate(x1 + xa*t1, y1 + ya*t1);
	}
	// polgyon vertex sequence CW or CCW?
	//// assume coords[0] equals coords[coords.length-1]
	// -1: ccw; 0: undefined (line segment); 1: cw
	//
	/**
	 * @param coords
	 * @return
	 */
	public static int orientation(Coordinate[] coords){
		double area2 = areaDouble(coords);
		if(area2 < 0.0){
			return -1;
		}else if(area2 > 0.0){
			return 1;
		}else{
			return 0;
		}
	}
	// negative if CCW
	/**
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	public static double areaDoubleTriangle(Coordinate a, Coordinate b, Coordinate c){
		return a.x*(c.y-b.y) + b.x*(a.y - c.y) + c.x * (b.y - a.y);
	}
	// assume coords[0] equals coords[coords.length-1]
	/**
	 * @param coords
	 * @return 2*area of coords, negative if ccw
	 */
	public static double areaDouble(Coordinate[] coords){
		double rlt = 0.0;
		if(coords.length > 3){
			Coordinate sp = coords[0];
			for(int i = 1; i < coords.length; ++i){
				Coordinate ep = coords[i];
				rlt+= (sp.y + ep.y)*(ep.x - sp.x);
				sp = ep;
			}
		}
		return rlt;
	}
	/** orientation of a line segment sc-ec in radian (-PI, PI]
	 * @param sc
	 * @param ec
	 * @return
	 */
	public static double segOrientationRad(Coordinate sc, Coordinate ec){
		return segOrientationRad(sc.x, sc.y, ec.x, ec.y);
	}
	//
	public static double segOrientationRad(double xs, double ys, double xe, double ye){
		double xdiff = xe - xs;
		double ydiff = ye - ys;
		return Math.atan2(ydiff, xdiff);
	}
	/** orientation of a line segment sc-ec in degree (-180, 180)
	 * @param sc
	 * @param ec
	 * @return
	 */
	public static double segOrientationDeg(Coordinate sc, Coordinate ec){
		return 180*segOrientationRad(sc, ec) / Math.PI;
	}
	//
	// if c0-c1 and c2-c3 are near vertical 
	//
	public static boolean isNearVertical(Coordinate c0, Coordinate c1,Coordinate c2,Coordinate c3, double errDeg ){
		double ang0 = compAngle(c0, c1);
		double ang1 = compAngle(c2, c3);
		double diff = AngleUtility.minAngDiffDeg(ang0, ang1);
		if(Math.abs(90 - diff) <= errDeg){
			return true;
		}else{
			return false;
		}
	}
	public static boolean isNearParalell(Coordinate c0, Coordinate c1,Coordinate c2,Coordinate c3, double errDeg ){
		double ang0 = JTSUtility.compAngle(c0, c1);
		double ang1 = JTSUtility.compAngle(c2, c3);
		double diff = AngleUtility.minAngDiffDeg(ang0, ang1);
		if(diff <= errDeg){
			return true;
		}else{
			return false;
		}
	}
	/** angle in radian for line segment from sc to ec
	 * @param sc
	 * @param ec
	 * @return
	 */
	public static double compEdgeAngleRad(Coordinate sc, Coordinate ec){
		double xdff = ec.x - sc.x;
		double ydiff = ec.y - sc.y;
		return Math.atan2(ydiff, xdff);
	}

	/** normalised angle for line c0-c1 while the line goes from left to right or from bottom up
	 * @param c0
	 * @param c1
	 * @return
	 */
	public static double compEdgeAngleNormRad(Coordinate c0, Coordinate c1){
		if(c0.x < c1.x){
			return compEdgeAngleRad(c0, c1);
		}else if(c0.x > c1.x){
			return compEdgeAngleRad(c1, c0);
		}else{
			if(c0.y <= c1.y){
				return compEdgeAngleRad(c0, c1);
			}else{
				return compEdgeAngleRad(c1, c0);
			}
		}
	}	
	/** given a front orientation frontOrient, compute the front side of rectangle mmbr
	 *  the centre of mmbr is used as the base of orientation
	 * @param frontOrient
	 * @param mmbr Coordinate[4], an aligned MBR
	 * @return front side index (0 for side mmbr[0]-mmbr[1]
	 */
	public static int compFrontIdxCorner(double frontOrient, Coordinate[] mmbr){
		if(frontOrient > Math.PI || frontOrient < - Math.PI){
			return -1; // error
		}
		Coordinate cen = new Coordinate((mmbr[0].x + mmbr[2].x)*0.5, (mmbr[0].y + mmbr[2].y)*0.5);
		double ang0 = JTSUtility.compEdgeAngleRad(cen, mmbr[0]);
		double ang1 = JTSUtility.compEdgeAngleRad(cen, mmbr[1]);
		if(AngleUtility.containsAng(frontOrient, ang0, ang1)){
			return 0;
		}
		double ang2 = JTSUtility.compEdgeAngleRad(cen, mmbr[2]);
		if(AngleUtility.containsAng(frontOrient, ang1, ang2)){
			return 1;
		}
		double ang3 = JTSUtility.compEdgeAngleRad(cen, mmbr[3]);
		if(AngleUtility.containsAng(frontOrient, ang2, ang3)){
			return 2;
		}
		if(AngleUtility.containsAng(frontOrient, ang3, ang0)){
			return 3;
		}
		int shortSide = mmbr[0].distance(mmbr[1]) <= mmbr[1].distance(mmbr[2])?0:1;
		if(frontOrient == ang0){
			return shortSide == 0?0:3;
		}
		if(frontOrient == ang1){
			return shortSide == 0?0:1;
		}
		if(frontOrient == ang2){
			return shortSide == 0?2:1;
		}
		if(frontOrient == ang3){
			return shortSide == 0?2:3;
		}
		return -1; // error
	}
	//
	/**Given a front orientation frontOrient, use edges orientations to align frontOrient, 
	 * the angle of the longest edge with smallest angle difference to frontOrient is returned
	 * 
	 * @param frontOrient
	 * @param mmbr
	 * @return
	 */
	public static double compFrontAngleByEdge(double frontOrient, Coordinate[] mmbr){
		boolean longE01 = mmbr[0].distance(mmbr[1]) > mmbr[1].distance(mmbr[2]);
		double ang01 = compEdgeAngleRad(mmbr[0], mmbr[1]);
		double ang12 = compEdgeAngleRad(mmbr[1], mmbr[2]);
		double ang23 = AngleUtility.angReverse(ang01);
		double ang30 = AngleUtility.angReverse(ang12);
		if(AngleUtility.containsAng(frontOrient, ang30, ang01)){
			double ang = AngleUtility.minAngDiff(frontOrient, ang01);
			if(ang < rad45 || (ang == rad45 && longE01)){
				return ang01;
			}else{
				return ang30;
			}
		}else if(AngleUtility.containsAng(frontOrient, ang01, ang12)){
			double ang = AngleUtility.minAngDiff(frontOrient, ang01);
			if(ang < rad45 || (ang == rad45 && longE01)){
				return ang01;
			}else{
				return ang12;
			}
		}else if(AngleUtility.containsAng(frontOrient, ang12, ang23)){
			double ang = AngleUtility.minAngDiff(frontOrient, ang23);
			if(ang < rad45 || (ang == rad45 && longE01)){
				return ang23;
			}else{
				return ang12;
			}
		}else{
			double ang = AngleUtility.minAngDiff(frontOrient, ang23);
			if(ang < rad45 || (ang == rad45 && longE01)){
				return ang23;
			}else{
				return ang30;
			}
		}
	}

	/** compute the front orientation (coming out of the front side of mmbr at right angle) based on a front direction reference (not necessarily at a right angle to any side)
	 * @param frontRef: angle value in radian
	 * @param mmbr four coordinates in CCW, representing the mmbr of a geometry
	 * @return
	 */

	public static double compFrontAngleByCorner(double frontRef, Coordinate[] mmbr){

		int sideIdx = compFrontIdxCorner(frontRef, mmbr);

		if(sideIdx == 0){// side 0-1
			return JTSUtility.compEdgeAngleRad(mmbr[3], mmbr[0]);
		}
		if(sideIdx == 1){// side 1-2
			return JTSUtility.compEdgeAngleRad(mmbr[0], mmbr[1]);
		}
		if(sideIdx == 2){// side 2-3
			return JTSUtility.compEdgeAngleRad(mmbr[1], mmbr[2]);
		}
		if(sideIdx == 3){// side 3-0
			return JTSUtility.compEdgeAngleRad(mmbr[2], mmbr[3]);
		}
		return Double.MAX_VALUE; // error computing front side idx
	}
	/******************************************
	 * 
	 *  Triangle
	 * 
	 ******************************************/

	// assuming corA-corB-corC are in ccw order
	//
	/**
	 * @param corA
	 * @param corB
	 * @param corC
	 * @return
	 */
	public static Coordinate compFermatPoint(Coordinate corA, Coordinate corB, Coordinate corC){
		double a = corC.distance(corB);
		double aa = a*a;
		double b = corC.distance(corA);
		double bb = b*b;
		double c = corA.distance(corB);
		double cc = c*c;
		double maxDist = Math.max(aa, Math.max(bb, cc)); 
		double maxAng;
		if(aa == maxDist){
			maxAng = 0.5*Math.acos((bb+cc-aa)/(b*c));
		}else if (bb == maxDist){
			maxAng = 0.5*Math.acos((aa+cc-bb)/(a*c));
		}else{
			maxAng = 0.5*Math.acos((aa+bb-cc)/(a*b));
		}
		if(maxAng < rad120){
			Vector2D vca = new Vector2D(corC.x-corA.x, corC.y-corA.y);
			vca = vca.rotate(rad60);
			Vector2D vba = new Vector2D(corB.x - corA.x, corB.y - corA.y);
			vba = vba.rotate(-rad60);
			// find intersection point of corB-vca and corC-vba
			double xa = vca.getX() + corA.x - corB.x;
			double ya = vca.getY() + corA.y - corB.y;
			double xb = corC.x - vba.getX() - corA.x;
			double yb = corC.y - vba.getY() - corA.y;
			double xc = corB.x - corC.x;
			double yc = corB.y - corC.y;
			double den = ya*xb - xa*yb;
			double num1 = yb*xc - xb*yc;
			double t1 = num1 / den;
			return new Coordinate(corB.x + xa*t1, corB.y + ya*t1);
		}else{
			return null;
		}
	}
	/**
	 * @param coords
	 * @return
	 */
	public static Coordinate compFermatPoint(Coordinate[] coords){
		if(coords.length!=3){
			return null;
		}
		return compFermatPoint(coords[0], coords[1], coords[2]);
	}
	/**
	 * @param theta
	 * @return one of eight directions: E, NE, N, NW, W, SW, S, SE
	 */
	public static String ang2Dir(double theta){
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
	public static int checkPlineConnectivity(Coordinate[] carray0, Coordinate[] carray1){
		int rtn = 0;
		Coordinate v0s = carray0[0];
		Coordinate v0e = carray0[carray0.length-1];
		Coordinate v1s = carray1[0];
		Coordinate v1e = carray1[carray1.length-1];
		if(v0e.equals2D(v1s)){
			rtn = 0; // both original order
			if(v0s.equals2D(v1e)){
				rtn+=4; // ring
			}
		}else if(v0e.equals2D(v1e)){
			rtn = 2; // a0 original, a1 reversed
			if(v0s.equals2D(v1s)){
				rtn+=4;
			}
		}else if(v0s.equals2D(v1s)){
			rtn = 1; // a0 reverse, v1 original
			if(v0e.equals2D(v1e)){
				rtn+=4;
			}
		}else if(v0s.equals2D(v1e)){
			rtn = 3; // both reverse
			if(v0e.equals2D(v1s)){
				rtn+=4;
			}
		}else{
			return -1;
		}
		return rtn;
	}
	//
	// note: returned array maybe a ring, just one common coord is omitted, the other closes the ring
	//
	public static Coordinate[] mergeCoords(Coordinate[] carray0, Coordinate[] carray1, int code){
		if(code< 0)
			return null;
		boolean c0rev = (code & 1) != 0; // 0b001
		boolean c1rev = (code & 2) != 0; // 0b010
		boolean isRing = (code & 4) != 0; // 0b100

		int numVer = carray0.length+carray1.length-1;
		
		Coordinate[] newArray = new Coordinate[numVer];
		int cnt = 0;
		if(c0rev){
			for(int i = carray0.length-1; i >=0; --i){
				newArray[cnt] = carray0[i];
				++cnt;
			}
		}else{
			for(int i = 0; i < carray0.length; ++i){
				newArray[cnt] = carray0[i];
				++cnt;
			}
		}
		if(c1rev){
			for(int i = carray1.length-2; i >= 0; --i){
				newArray[cnt] = carray1[i];
				++cnt;
			}
		}else{
			for(int i = 1; i < carray1.length; ++i){
				newArray[cnt] = carray1[i];
				++cnt;
			}
		}
		return newArray;
	}
	//
	public static boolean checkTouch(Geometry ln, Geometry boundary){
		if(ln.getGeometryType().compareToIgnoreCase("LINESTRING")==0){
			return checkTouchLS((LineString)ln, boundary);
		}else if(ln.getGeometryType().compareToIgnoreCase("MMULTILINESTRING")==0){
			int numParts = ln.getNumGeometries();
			for(int i = 0; i < numParts; ++i){
				LineString ls = (LineString)ln.getGeometryN(i);
				if(checkTouchLS(ls, boundary)){
					return true;
				}
			}
			return false;
		}
		return false;
	}
	
	public static boolean checkTouchLS(LineString ln, Geometry boundary){
		int numVerLn = ln.getNumPoints();
		Coordinate sp = ln.getCoordinateN(0);
		Coordinate ep = ln.getCoordinateN(numVerLn-1);
		Coordinate[] coords = boundary.getCoordinates();
		for(Coordinate coord:coords){
			if(coord.equals2D(sp))
				return true;
			else if(coord.equals2D(ep))
				return true;
		}
		return false;
	}
	//	
}

