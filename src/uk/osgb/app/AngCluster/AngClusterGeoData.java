package uk.osgb.app.AngCluster;
import java.util.Collection;
import java.util.List;
import java.util.Vector;

import uk.osgb.utilities.JTSUtility;
import uk.osgb.utilities.GeomUtility;
import uk.osgb.utilities.AngleUtility;
import uk.osgb.utilities.JUMPUtility;

import uk.osgb.ml.cluster.anglecluster.AngCluster;
import uk.osgb.ml.cluster.anglecluster.AngClusterGenerator;

import com.vividsolutions.jts.algorithm.ConvexHull;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jump.feature.AttributeType;
import com.vividsolutions.jump.feature.BasicFeature;
import com.vividsolutions.jump.feature.Feature;
import com.vividsolutions.jump.feature.FeatureCollection;
import com.vividsolutions.jump.feature.FeatureDataset;
import com.vividsolutions.jump.feature.FeatureSchema;

/** this class provides various wrapper to cluster building edges in geo dataset
 * 
 * @author SZhou
 *
 */
public class AngClusterGeoData {
	public static double ANGREVTOL = Math.PI/60.0;
	//
	/**
	 * @param feats a collection of JUMP features
	 * @param fn_out output shapefile pathname
	 * @param fid_nm FID field name of the input features
	 * @param varTol maximum within-cluster variance ratio over total angular variance. 1.0 if not used for control. if < 1.0, will be checked even if the maximum angular width is not reached
	 * @param maxAW maximum angular width (in DEGREE) for generated clusters. >=360 if not used for control
	 * @param modulo modulo operation to be applied. 0 for none, 2 for PI and 4 for 0.5PI
	 * @param pertAngDeg perturbance angle  
	 */
	public static void clusteringJUMPFeatures(Collection<Feature> feats, String fn_out, String fid_nm, double varTol, double maxAW, int modulo, double pertAngDeg){
		//
		boolean nullFid = false;
		if(fid_nm == null){
			nullFid = true;
		}
		// original geometry copied over to new feature with clustering informaiton
		FeatureSchema sch = new FeatureSchema();
		sch.addAttribute("GEOMETRY",  AttributeType.GEOMETRY);
		sch.addAttribute("FID", AttributeType.STRING);
		sch.addAttribute("NUMSIDES", AttributeType.INTEGER);
		sch.addAttribute("ROUNDNESS", AttributeType.DOUBLE);
		sch.addAttribute("SQUARENESS", AttributeType.DOUBLE);
		// vectors for cluster members
		FeatureSchema schClu = new FeatureSchema();
		schClu.addAttribute("GEOMETRY",  AttributeType.GEOMETRY);
		schClu.addAttribute("FID", AttributeType.STRING);
		schClu.addAttribute("EDGECOUNT", AttributeType.INTEGER);
		schClu.addAttribute("DIRECTION", AttributeType.DOUBLE);
		schClu.addAttribute("BOUNDARY_R", AttributeType.DOUBLE);
		schClu.addAttribute("BOUNDARY_L", AttributeType.DOUBLE);
		schClu.addAttribute("ANGLEWIDTH", AttributeType.DOUBLE);
		schClu.addAttribute("VARIANCE", AttributeType.DOUBLE);
		schClu.addAttribute("WEIGHT", AttributeType.DOUBLE);
		schClu.addAttribute("LENGTH_PROJ", AttributeType.DOUBLE);

		FeatureCollection fc = new FeatureDataset(sch);//
		FeatureCollection fcClu = new FeatureDataset(schClu);
		int i = 0;
		for(Feature feat:feats){
			String fid = nullFid?Integer.toString(i++):feat.getString(fid_nm);
			Geometry geom = ((Polygon)feat.getGeometry()).getExteriorRing();
			Coordinate centroid = geom.getCentroid().getCoordinate();
			double area = Math.abs(JTSUtility.areaDouble(geom.getCoordinates())*0.5);
			double perimeter = geom.getLength();
			Geometry mmbr = JTSUtility.compMinMBRGeom(geom, false);
			double areaMmbr = mmbr.getArea();
			double roundness = GeomUtility.roundness(area, perimeter);
			double squareness = area/Math.abs(JTSUtility.areaDouble(mmbr.getCoordinates())*0.5);;
			
			Collection<AngCluster> clus2 = clusteringJTSGeometry(geom, false, varTol, maxAW, 0.001, modulo, pertAngDeg);
			
			AngClusterGenerator.report(clus2);
			String clutxt = "";
			for(AngCluster clu:clus2){
				Feature featClu = new BasicFeature(schClu);
				Coordinate[] coord = new Coordinate[2];
				featClu.setAttribute("FID", fid);
				featClu.setAttribute("EDGECOUNT", clu.getMembers().size());
				double ang = clu.getAng();
				double angR = clu.getAngF().getAng();
				double angL = clu.getAngT().getAng();
				double angW = clu.getAngWidth();
				double projLen = clu.getLenSum();
				double weight = clu.getWeight();
				double variance = clu.getVariance();
				double xe = centroid.x + projLen*Math.cos(ang);
				double ye = centroid.y + projLen*Math.sin(ang);
				coord[0] = new Coordinate(centroid);
				coord[1] = new Coordinate(xe, ye);
				Geometry vec = geom.getFactory().createLineString(coord);
				featClu.setGeometry(vec);
				featClu.setAttribute("DIRECTION", ang*180.0/Math.PI); // in degree
				featClu.setAttribute("BOUNDARY_R", angR*180/Math.PI);
				featClu.setAttribute("BOUNDARY_L", angL*180/Math.PI);
				featClu.setAttribute("ANGLEWIDTH", angW*180/Math.PI);
				featClu.setAttribute("VARIANCE", variance);
				featClu.setAttribute("WEIGHT", weight);
				featClu.setAttribute("LENGTH_PROJ", projLen);
				fcClu.add(featClu);
			}
			//Geometry hull = clusterCH(clus2);
			//double areaHull = Math.abs(JTSUtility.areaDouble(hull.getCoordinates()));
			Feature feature = new BasicFeature(sch);
			feature.setAttribute("GEOMETRY", feat.getGeometry());
			feature.setAttribute("NUMSIDES", clus2.size());
			feature.setAttribute("ROUNDNESS", roundness);
			feature.setAttribute("SQUARENESS", squareness);
			feature.setAttribute("FID", fid);
			fc.add(feature);

			if(i%100 == 0){
				System.out.println(i + " of " + feats.size());
			}
		}
		System.out.println("exporting clustering result to ShapeFile "+fn_out +"......");
		JUMPUtility.exportShapeFile(fcClu, fn_out);
		// export polygons with additional information (number of sides, roundness, squarness) if required
		String fn_out_main = fn_out.substring(0, fn_out.length()-4);
		JUMPUtility.exportShapeFile(fc, fn_out_main+"_plg.shp");
		System.out.println("Done...");		
	}
	/*********************
	 * 
	 *  
	 * 
	 ********************/
	/****************
	 * convert a JTS Polygon edge segments to a vector of orientation and a vector of length 
	 * @param plg
	 * @param oriVec
	 * @param lenVec
	 * @param useHole interior rings representing holes are also processed
	 * @param reverse
	 * @return
	 */
	private static boolean JTSPolygon2Vector(Polygon plg, Vector<Double> oriVec, Vector<Double> lenVec, boolean useHole, boolean reverse){
		int numHoles = plg.getNumInteriorRing();
		JTSRing2Vector(plg.getExteriorRing().getCoordinates(), oriVec, lenVec, reverse);
		if(useHole && numHoles> 0){
			for(int i = 0; i < numHoles; ++i){
				JTSRing2Vector(plg.getInteriorRingN(i).getCoordinates(), oriVec, lenVec, !reverse);
			}
		}
		return true;
	}
	
	private static boolean JTSMultiPolygon2Vector(MultiPolygon mulPlg, Vector<Double> oriVec, Vector<Double> lenVec, boolean useHole, boolean reverse){
		int numPlg = mulPlg.getNumGeometries();
		for(int i = 0; i < numPlg; ++i){
			JTSPolygon2Vector((Polygon)mulPlg.getGeometryN(i), oriVec, lenVec, useHole, reverse);
		}
		return true;
	}
	
	private static boolean JTSLineString2Vector(LineString ls, Vector<Double> oriVec, Vector<Double> lenVec, boolean reverse){
		return JTSRing2Vector(ls.getCoordinates(), oriVec, lenVec, reverse);
	}
	private static boolean JTSMultiLineString2Vector(MultiLineString multiLs, Vector<Double> oriVec, Vector<Double> lenVec, boolean reverse){
		int numLs = multiLs.getNumGeometries();
		for(int i = 0; i < numLs; ++i) {
			JTSLineString2Vector((LineString)multiLs.getGeometryN(i), oriVec, lenVec, reverse);
		}
		return true;
	}
	
	/**
	 * @param geom JTS geometry to be processed
	 * @param oriVec Vector to store orientation values of edges
	 * @param lenVec Vector to store length values of edges
	 * @param useHole if holes of polygons are also considered
	 * @param reverse if the vertex order should be reversed (reversed order is always used for holes)
	 * @return 
	 */
	public static boolean JTSGeom2Vector(Geometry geom, Vector<Double> oriVec, Vector<Double> lenVec, boolean useHole, boolean reverse){
		if(geom!=null){
			String gt = geom.getGeometryType();
			if(gt.compareToIgnoreCase("POLYGON")==0){
				return JTSPolygon2Vector((Polygon)geom, oriVec, lenVec, useHole, reverse);
			}else if(gt.compareToIgnoreCase("LINESTRING")==0 || gt.compareToIgnoreCase("LINEARRING")==0){
				return JTSLineString2Vector((LineString)geom, oriVec, lenVec, reverse);
			}else if(gt.compareToIgnoreCase("MULTIPOLYGON") == 0){
				return JTSMultiPolygon2Vector((MultiPolygon)geom, oriVec, lenVec, useHole, reverse);
			}else if(gt.compareToIgnoreCase("MULTILINESTRING") == 0){
				return JTSMultiLineString2Vector((MultiLineString)geom, oriVec, lenVec, reverse);
			}else{
				return false;
			}
		}
		return false;
	}
	//
	/**
	 * @param gc
	 * @param oriVec
	 * @param lenVec
	 * @param useHole
	 * @param reverse
	 * @return
	 */
	private static boolean JTSGeomCollection2Vector(GeometryCollection gc, Vector<Double> oriVec, Vector<Double> lenVec, boolean useHole, boolean reverse){
		int numPart = gc.getNumGeometries();
		if(numPart > 0) {
			for(int i = 0; i < numPart; ++i) {
				Geometry geom = gc.getGeometryN(i);
				JTSGeom2Vector(geom, oriVec, lenVec, useHole, reverse);
			}
			return true;
		}
		return false;
	}
	//
	/**
	 * @param coords
	 * @param oriVec
	 * @param lenVec
	 * @param reverse
	 * @return
	 */
	private static boolean JTSRing2Vector(Coordinate[] coords, Vector<Double> oriVec, Vector<Double> lenVec, boolean reverse){
		Coordinate sp = coords[0], ep;
		double baseOffset = 0.0; // orientation of the longest edge, used in angle reversion
		if(reverse){
			double maxL = 0.0;
			double maxA = 0.0;
			for(int i = 1; i < coords.length;++i){
				ep = coords[i];
				double len = sp.distance(ep); 
				if(sp.equals2D(ep) || len < maxL){
					sp = ep;
					continue;
				}
				maxL = len;
				double xdiff = ep.x - sp.x;
				double ydiff = ep.y - sp.y;
				maxA = Math.atan2(ydiff, xdiff);
				sp = ep;
			}
			if(maxL > 0.0){
				baseOffset = AngleUtility.angFromTo(maxA, ANGREVTOL);
			}
		}
		
		for(int i = 1; i < coords.length;++i){
			ep = coords[i];
			double len = 0.0; 
			if(sp.equals2D(ep)){
				continue;
			}else {
				len = sp.distance(ep);
			}
			double xdiff = ep.x - sp.x;
			double ydiff = ep.y - sp.y;
			double ori = Math.atan2(ydiff, xdiff);
			double oriOff = AngleUtility.angSum(ori, baseOffset);
			if(oriOff < 0.0 && reverse){
				//ori = JTSUtility.angReverse(ori, baseOffset);
				ori = AngleUtility.angReverse(ori);
			}
			oriVec.add(ori);
			lenVec.add(len);
			sp = ep;
		}
		return true;
	}
	
	/** clustering angle values stored in a orientation vector and a length/weight vectors 
	 * @param orisVec orientation values
	 * @param lensVec weight values in the same order as orientation vector
	 * @param varTol maximum within-cluster variance ratio over total angular variance. 1.0 if not used for control. if < 1.0, will be checked even if the maximum angular width is not reached
	 * @param maxAW maximum angular width (in DEGREE) for generated clusters. >=360 if not used for control
	 * @param lenThreshold minimum ratio of projected weight sum over the sum of all projected weight for all clusters
	 * @param modulo modulo operation to be applied. 0 for none, 2 for PI and 4 for 0.5PI
	 * @param pertAngDeg perturbance angle  
	 * @return
	 */
	private static Collection<AngCluster> clusteringVectors(Vector<Double> orisVec, Vector<Double> lensVec, double varTol, double maxAW, double lenThreshold, int modulo, double pertAngDeg){
		if(pertAngDeg!=0.0){
			pertAngDeg = Math.toRadians(pertAngDeg);
		}

		int numEdges = orisVec.size();
		if(numEdges > 0){
			double[] oris = new double[numEdges];
			double[] weights = new double[numEdges];
			for(int i = 0; i < numEdges; ++i){
				oris[i] = orisVec.get(i);
				weights[i] = lensVec.get(i);
			}
			
			Collection<AngCluster> clus = AngClusterGenerator.clusteringAngleArray(oris, weights, varTol, maxAW, modulo, pertAngDeg);
			
			double[] warray = new double[clus.size()];
			double totW = 0.0; 
			for(AngCluster clu:clus){
				totW += clu.getLenSum();
			}
			double threshold = totW*lenThreshold;
			Collection<AngCluster> rlt = new Vector<AngCluster>();
			for(AngCluster clu:clus){
				if(clu.getLenSum() >= threshold){
					rlt.add(clu);
				}
			}
			return rlt;
		}
		return null;

	}
		
	/**test using OSMM Topo layer data (building polygons with attributes)
	 * @param fn_in input shapefile name (full path)
	 * @param fn_out output shapefile name for clusters represented by radial from centroid of the polygon 
	 * @param varTol maximum intra-cluster variance ratio
	 * @param modulo modulo (0 for no modulo, 2 for PI and 4 for PI/2)
	 * @param tiltAngDeg perturbance angle (0.0 for no perturbance)
	 */
	private static void testOSMM(String fn_in, String fn_out, double varTol, double maxAW, int modulo, double tiltAngDeg){
		List<Feature> feats = JUMPUtility.loadShapeFile(fn_in, false, "fid");
		FeatureSchema sch = new FeatureSchema();
		
		sch.addAttribute("GEOMETRY",  AttributeType.GEOMETRY);
		sch.addAttribute("TOID", AttributeType.STRING);
		sch.addAttribute("NUMSIDES", AttributeType.INTEGER);
		sch.addAttribute("ROUNDNESS", AttributeType.DOUBLE);
		sch.addAttribute("SQUARENESS", AttributeType.DOUBLE);

		FeatureSchema schClu = new FeatureSchema();
		schClu.addAttribute("GEOMETRY",  AttributeType.GEOMETRY);
		schClu.addAttribute("TOID", AttributeType.STRING);
		schClu.addAttribute("EDGECOUNT", AttributeType.INTEGER);
		schClu.addAttribute("DIRECTION", AttributeType.DOUBLE);
		schClu.addAttribute("VARIANCE", AttributeType.DOUBLE);
		schClu.addAttribute("WEIGHT", AttributeType.DOUBLE);
		schClu.addAttribute("LENGTH_PROJ", AttributeType.DOUBLE);

		FeatureCollection fc = new FeatureDataset(sch);
		FeatureCollection fcClu = new FeatureDataset(schClu);
		int i = 0;
		for(Feature feat:feats){
			String toid = (String) feat.getAttribute("fid");
			Geometry geom = ((Polygon)feat.getGeometry()).getExteriorRing();
			Coordinate centroid = geom.getCentroid().getCoordinate();
			double area = Math.abs(JTSUtility.areaDouble(geom.getCoordinates())*0.5);
			double perimeter = geom.getLength();
			Geometry mmbr = JTSUtility.compMinMBRGeom(geom, false);
			double areaMmbr = mmbr.getArea();
			double roundness = GeomUtility.roundness(area, perimeter);
			double squareness = area/Math.abs(JTSUtility.areaDouble(mmbr.getCoordinates())*0.5);;
			Collection<AngCluster> clus2 = clusteringJTSGeometry(geom, false, varTol, maxAW, 0.001, modulo, tiltAngDeg);
			String clutxt = "";
			for(AngCluster clu:clus2){
				Feature featClu = new BasicFeature(schClu);
				Coordinate[] coord = new Coordinate[2];
				featClu.setAttribute("TOID", toid);
				featClu.setAttribute("EDGECOUNT", clu.getMembers().size());
				double ang = clu.getAng();
				double projLen = clu.getLenSum();
				double weight = clu.getWeight();
				double variance = clu.getVariance();
				double xe = centroid.x + projLen*Math.cos(ang);
				double ye = centroid.y + projLen*Math.sin(ang);
				coord[0] = new Coordinate(centroid);
				coord[1] = new Coordinate(xe, ye);
				Geometry vec = geom.getFactory().createLineString(coord);
				featClu.setGeometry(vec);
				featClu.setAttribute("DIRECTION", ang*180.0/Math.PI); // in degree
				featClu.setAttribute("VARIANCE", variance);
				featClu.setAttribute("WEIGHT", weight);
				featClu.setAttribute("LENGTH_PROJ", projLen);
				fcClu.add(featClu);
			}

			//Geometry hull = clusterCH(clus2);

			//double areaHull = Math.abs(JTSUtility.areaDouble(hull.getCoordinates()));
			Feature feature = new BasicFeature(sch);
			//featClu.setAttribute("GEOMETRY", hull);
			feature.setAttribute("GEOMETRY", feat.getGeometry());
			feature.setAttribute("NUMSIDES", clus2.size());
			feature.setAttribute("ROUNDNESS", roundness);
			feature.setAttribute("SQUARENESS", squareness);
			feature.setAttribute("TOID", toid);

			fc.add(feature);
			i++;
			if(i%100 == 0){
				System.out.println(i + " of " + feats.size());
			}
		}
		System.out.println("exporting...");
		// export clusters
		JUMPUtility.exportShapeFile(fcClu, fn_out);
		// export polygons with additional information (number of sides, roundness, squarness) if required
		//String fn_out_main = fn_out.substring(0, fn_out.length()-4);
		//JUMPUtility.exportShapeFile(fc, fn_out_main+"_plg.shp");
		
		System.out.println("Done...");
		
	}
	/**test using shape file converted from wkt (no addition attribute)
	 * @param fn_in
	 * @param fn_out
	 * @param varTol
	 * @param modulo
	 * @param tiltAngDeg
	 */
	private static void testWKT(String fn_in, String fn_out, double varTol, double maxAW, int modulo, double tiltAngDeg){
		//List<Feature> feats = JUMPUtility.loadShapeFile("d:/temp/SX99-TA_descriptiv_contains_Building.shp", false, "TOID");
		//List<Feature> feats = JUMPUtility.loadShapeFile("d:/temp/clu_bld01.shp", false, "FID");
		List<Feature> feats = JUMPUtility.loadShapeFile(fn_in, false, "fid");
		FeatureSchema sch = new FeatureSchema();
		sch.addAttribute("GEOMETRY",  AttributeType.GEOMETRY);
		sch.addAttribute("FID", AttributeType.STRING);
		sch.addAttribute("NUMSIDES", AttributeType.INTEGER);
		sch.addAttribute("ROUNDNESS", AttributeType.DOUBLE);
		sch.addAttribute("SQUARENESS", AttributeType.DOUBLE);

		FeatureSchema schClu = new FeatureSchema();
		schClu.addAttribute("GEOMETRY",  AttributeType.GEOMETRY);
		schClu.addAttribute("FID", AttributeType.STRING);
		schClu.addAttribute("EDGECOUNT", AttributeType.INTEGER);
		schClu.addAttribute("DIRECTION", AttributeType.DOUBLE);
		schClu.addAttribute("BOUNDARY_R", AttributeType.DOUBLE);
		schClu.addAttribute("BOUNDARY_L", AttributeType.DOUBLE);
		schClu.addAttribute("ANGLEWIDTH", AttributeType.DOUBLE);
		schClu.addAttribute("VARIANCE", AttributeType.DOUBLE);
		schClu.addAttribute("WEIGHT", AttributeType.DOUBLE);
		schClu.addAttribute("LENGTH_PROJ", AttributeType.DOUBLE);

		FeatureCollection fc = new FeatureDataset(sch);//
		FeatureCollection fcClu = new FeatureDataset(schClu);
		int i = 0;
		for(Feature feat:feats){
			String fid = Integer.toString(i);
			Geometry geom = ((Polygon)feat.getGeometry()).getExteriorRing();
			Coordinate centroid = geom.getCentroid().getCoordinate();
			double area = Math.abs(JTSUtility.areaDouble(geom.getCoordinates())*0.5);
			double perimeter = geom.getLength();
			Geometry mmbr = JTSUtility.compMinMBRGeom(geom, false);
			double areaMmbr = mmbr.getArea();
			double roundness = GeomUtility.roundness(area, perimeter);
			double squareness = area/Math.abs(JTSUtility.areaDouble(mmbr.getCoordinates())*0.5);;
			Collection<AngCluster> clus2 = clusteringJTSGeometry(geom, false, varTol, maxAW, 0.001, modulo, tiltAngDeg);
			AngClusterGenerator.report(clus2);
			String clutxt = "";
			for(AngCluster clu:clus2){
				Feature featClu = new BasicFeature(schClu);
				Coordinate[] coord = new Coordinate[2];
				featClu.setAttribute("FID", fid);
				featClu.setAttribute("EDGECOUNT", clu.getMembers().size());
				double ang = clu.getAng();
				double angR = clu.getAngF().getAng();
				double angL = clu.getAngT().getAng();
				double angW = clu.getAngWidth();
				double projLen = clu.getLenSum();
				double weight = clu.getWeight();
				double variance = clu.getVariance();
				double xe = centroid.x + projLen*Math.cos(ang);
				double ye = centroid.y + projLen*Math.sin(ang);
				coord[0] = new Coordinate(centroid);
				coord[1] = new Coordinate(xe, ye);
				Geometry vec = geom.getFactory().createLineString(coord);
				featClu.setGeometry(vec);
				featClu.setAttribute("DIRECTION", ang*180.0/Math.PI); // in degree
				featClu.setAttribute("BOUNDARY_R", angR*180/Math.PI);
				featClu.setAttribute("BOUNDARY_L", angL*180/Math.PI);
				featClu.setAttribute("ANGLEWIDTH", angW*180/Math.PI);
				featClu.setAttribute("VARIANCE", variance);
				featClu.setAttribute("WEIGHT", weight);
				featClu.setAttribute("LENGTH_PROJ", projLen);
				fcClu.add(featClu);
			}
			//Geometry hull = clusterCH(clus2);
			//double areaHull = Math.abs(JTSUtility.areaDouble(hull.getCoordinates()));
			Feature feature = new BasicFeature(sch);
			feature.setAttribute("GEOMETRY", feat.getGeometry());
			feature.setAttribute("NUMSIDES", clus2.size());
			feature.setAttribute("ROUNDNESS", roundness);
			feature.setAttribute("SQUARENESS", squareness);
			feature.setAttribute("FID", fid);

			fc.add(feature);
			i++;
			if(i%100 == 0){
				System.out.println(i + " of " + feats.size());
			}
		}
		System.out.println("exporting...");
		JUMPUtility.exportShapeFile(fcClu, fn_out);
		// export polygons with additional information (number of sides, roundness, squarness) if required
		//String fn_out_main = fn_out.substring(0, fn_out.length()-4);
		//JUMPUtility.exportShapeFile(fc, fn_out_main+"_plg.shp");
		System.out.println("Done...");

	}
	//
	public static void testWKTElbow(String fn_in, String fn_out, int modulo, double tiltAngDeg){
		//List<Feature> feats = JUMPUtility.loadShapeFile("d:/temp/SX99-TA_descriptiv_contains_Building.shp", false, "TOID");
		//List<Feature> feats = JUMPUtility.loadShapeFile("d:/temp/clu_bld01.shp", false, "FID");
		List<Feature> feats = JUMPUtility.loadShapeFile(fn_in, false, "fid");
		FeatureSchema sch = new FeatureSchema();
		sch.addAttribute("GEOMETRY",  AttributeType.GEOMETRY);
		sch.addAttribute("TOID", AttributeType.STRING);
		sch.addAttribute("NUMSIDES", AttributeType.INTEGER);
		sch.addAttribute("ROUNDNESS", AttributeType.DOUBLE);
		sch.addAttribute("SQUARENESS", AttributeType.DOUBLE);

		FeatureSchema schClu = new FeatureSchema();
		schClu.addAttribute("GEOMETRY",  AttributeType.GEOMETRY);
		schClu.addAttribute("TOID", AttributeType.STRING);
		schClu.addAttribute("EDGECOUNT", AttributeType.INTEGER);
		schClu.addAttribute("DIRECTION", AttributeType.DOUBLE);
		schClu.addAttribute("BOUNDARY_R", AttributeType.DOUBLE);
		schClu.addAttribute("BOUNDARY_L", AttributeType.DOUBLE);
		schClu.addAttribute("ANGLEWIDTH", AttributeType.DOUBLE);
		schClu.addAttribute("VARIANCE", AttributeType.DOUBLE);
		schClu.addAttribute("WEIGHT", AttributeType.DOUBLE);
		schClu.addAttribute("LENGTH_PROJ", AttributeType.DOUBLE);

		FeatureCollection fc = new FeatureDataset(sch);//
		FeatureCollection fcClu = new FeatureDataset(schClu);
		int i = 0;
		for(Feature feat:feats){
			String toid = Integer.toString(i);
			Geometry geom = ((Polygon)feat.getGeometry()).getExteriorRing();
			Coordinate centroid = geom.getCentroid().getCoordinate();
			double area = Math.abs(JTSUtility.areaDouble(geom.getCoordinates())*0.5);
			double perimeter = geom.getLength();
			Geometry mmbr = JTSUtility.compMinMBRGeom(geom, false);
			double areaMmbr = mmbr.getArea();
			double roundness = GeomUtility.roundness(area, perimeter);
			double squareness = area/Math.abs(JTSUtility.areaDouble(mmbr.getCoordinates())*0.5);
			//
			Collection<AngCluster> clus2 = clusteringJTSGeometryElbow(geom, false, 0.001, modulo, tiltAngDeg);
			//
			AngClusterGenerator.report(clus2);
			String clutxt = "";
			for(AngCluster clu:clus2){
				Feature featClu = new BasicFeature(schClu);
				Coordinate[] coord = new Coordinate[2];
				featClu.setAttribute("TOID", toid);
				featClu.setAttribute("EDGECOUNT", clu.getMembers().size());
				double ang = clu.getAng();
				double angR = clu.getAngF().getAng();
				double angL = clu.getAngT().getAng();
				double angW = clu.getAngWidth();
				double projLen = clu.getLenSum();
				double weight = clu.getWeight();
				double variance = clu.getVariance();
				double xe = centroid.x + projLen*Math.cos(ang);
				double ye = centroid.y + projLen*Math.sin(ang);
				coord[0] = new Coordinate(centroid);
				coord[1] = new Coordinate(xe, ye);
				Geometry vec = geom.getFactory().createLineString(coord);
				featClu.setGeometry(vec);
				featClu.setAttribute("DIRECTION", ang*180.0/Math.PI); // in degree
				featClu.setAttribute("BOUNDARY_R", angR*180/Math.PI);
				featClu.setAttribute("BOUNDARY_L", angL*180/Math.PI);
				featClu.setAttribute("ANGLEWIDTH", angW*180/Math.PI);
				featClu.setAttribute("VARIANCE", variance);
				featClu.setAttribute("WEIGHT", weight);
				featClu.setAttribute("LENGTH_PROJ", projLen);
				fcClu.add(featClu);
			}
			//Geometry hull = clusterCH(clus2);
			//double areaHull = Math.abs(JTSUtility.areaDouble(hull.getCoordinates()));
			Feature feature = new BasicFeature(sch);
			feature.setAttribute("GEOMETRY", feat.getGeometry());
			feature.setAttribute("NUMSIDES", clus2.size());
			feature.setAttribute("ROUNDNESS", roundness);
			feature.setAttribute("SQUARENESS", squareness);
			feature.setAttribute("TOID", toid);

			fc.add(feature);
			i++;
			if(i%100 == 0){
				System.out.println(i + " of " + feats.size());
			}
		}
		System.out.println("exporting...");
		JUMPUtility.exportShapeFile(fcClu, fn_out);
		// export polygons with additional information (number of sides, roundness, squarness) if required
		//String fn_out_main = fn_out.substring(0, fn_out.length()-4);
		//JUMPUtility.exportShapeFile(fc, fn_out_main+"_plg.shp");
		System.out.println("Done...");

	}	
	// CH of angle vectors
	//
	/**
	 * @param clus
	 * @return
	 */
	private static Geometry clusterCH(Collection<AngCluster> clus){
		Coordinate[] coords = new Coordinate[clus.size()];
		int i = 0;
		for(AngCluster clu:clus){
			double ang = clu.getAng();
			double wei = clu.getLenSum();
			coords[i++] = new Coordinate(Math.cos(ang)*wei, Math.sin(ang)*wei);
		}
		GeometryFactory gf = new GeometryFactory();
		ConvexHull ch = new ConvexHull(coords, gf);
		return ch.getConvexHull();
	}
	//
	
	/** Clustering edge segments of a JTS Polygon
	 * @param geom a JTS Polygon
	 * @param reverse if the order of coordinate array is to be reversed
	 * @param varTol maximum ratio of intra-cluster variance over the total variance (1.0 if not used for control; if < 1.0, will be checked even if the cluster width is within the width threshold  
	 * @param maxAW maximum angular width (in DEGREE) of a generated cluster (not used if >=360
	 * @param lenThreshold minimum ratio of projected length sum over total length for a cluster to be included in output 
	 * @param modulo base of modulo operation (0 for none, 2 for PI and 4 for PI/2)
	 * @param tiltAngDeg perturbance angle (in DEGREE)
	 * @return
	 */
	private static Collection<AngCluster> clusteringJTSGeometry(Geometry geom, boolean reverse, double varTol, double maxAW, double lenThreshold, int modulo, double tiltAngDeg){
		if(tiltAngDeg!=0.0){
			tiltAngDeg = Math.toRadians(tiltAngDeg);
		}
		// edge orientations
		Vector<Double> orisVec = new Vector<Double>();
		// edge lengths
		Vector<Double> lensVec = new Vector<Double>();
		JTSGeom2Vector(geom, orisVec, lensVec, false, false);
		return clusteringVectors(orisVec, lensVec, varTol, maxAW, lenThreshold, modulo, tiltAngDeg);
	}
	
	/** testing elbow and silhouette, no cluster generated uet
	 * @param geom
	 * @param reverse
	 * @param lenThreshold
	 * @param modulo
	 * @param tiltAngDeg
	 * @return
	 */
	private static Collection<AngCluster> clusteringJTSGeometryElbow(Geometry geom, boolean reverse, double lenThreshold, int modulo, double pertAngDeg){
		if(pertAngDeg!=0.0){
			pertAngDeg = Math.toRadians(pertAngDeg);
		}
		// edge orientations
		Vector<Double> orisVec = new Vector<Double>();
		// edge lengths
		Vector<Double> lensVec = new Vector<Double>();
		JTSGeom2Vector(geom, orisVec, lensVec, false, false);
		
		int numEdges = orisVec.size();
		if(numEdges > 0){
			double[] oris = new double[numEdges];
			double[] weights = new double[numEdges];
			for(int i = 0; i < numEdges; ++i){
				oris[i] = orisVec.get(i);
				weights[i] = lensVec.get(i);
			}
			//
			Collection<AngCluster> clus = AngClusterGenerator.clusteringAngleArrayElbow(oris, weights, modulo, pertAngDeg);
			//
			double[] warray = new double[clus.size()];
			double totW = 0.0; 
			for(AngCluster clu:clus){
				totW += clu.getLenSum();
			}
			double threshold = totW*lenThreshold;
			Collection<AngCluster> rlt = new Vector<AngCluster>();
			for(AngCluster clu:clus){
				if(clu.getLenSum() >= threshold){
					rlt.add(clu);
				}
			}
			return rlt;
		}
		return null;

	}
	//
	public static void main(String[] args){
		//List<Feature> feats = JUMPUtility.loadShapeFile("D:/temp/ac02.shp", null, false, null);
		/*
		List<Feature> feats = JUMPUtility.loadShapeFile("D:/temp/SX99-TA_descriptiv_contains_Building.shp", "fid", false, null);		
		//clusteringJTSFeatures(feats, "D:/temp/ac02_m2t00.shp", null, 1.0, 15, 2, 0.0);
		clusteringJUMPFeatures(feats, "D:/temp/sx99rlt_10_15_0_15.shp", "fid", 1.0, 15, 0, 15.0);
		clusteringJUMPFeatures(feats, "D:/temp/sx99rlt_10_15_2_15.shp", "fid", 1.0, 15, 2, 15.0);
		clusteringJUMPFeatures(feats, "D:/temp/sx99rlt_10_15_4_15.shp", "fid", 1.0, 15, 4, 15.0);
		*/
		testWKTElbow("D:/temp/ac02.shp", "D:/temp/ac02_m2t00_elbow.shp", 0, 0.0);
		//testWKTElbow("D:/temp/ac03.shp", "D:/temp/ac03_m0t00_elbow.shp", 2, 15.0);
		//testWKTElbow("D:/temp/ac01s.shp", "D:/temp/ac01s_m2t00_elbow.shp", 2, 15.0);
		//testWKTElbow("D:/temp/ac04.shp", "D:/temp/ac04_m0t00_elbow.shp", 0, 0.0);
		/*
		testWKT("D:/temp/ac01.shp", "D:/temp/ac01_m0t00ec2.shp", 0.01, 0, 0.0);
		testWKT("D:/temp/ac01.shp", "D:/temp/ac01_m2t00ec2.shp", 0.01, 2, 0.0);
		testWKT("D:/temp/ac01.shp", "D:/temp/ac01_m2t15ec2.shp", 0.01, 2, 15.0);
		testWKT("D:/temp/ac01.shp", "D:/temp/ac01_m4t00ec2.shp", 0.01, 4, 0.0);
		testWKT("D:/temp/ac01.shp", "D:/temp/ac01_m4t15ec2.shp", 0.01, 4, 15.0);
		*/
		//testWKT("D:/temp/ac02.shp", "D:/temp/ac02_m0t00ec3.shp", 0.01, 0, 0.0);

		
	}

}
