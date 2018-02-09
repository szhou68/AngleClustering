/* Jump based IO etc.
 * 
 */
package uk.osgb.utilities;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiPoint;
import com.vividsolutions.jts.index.strtree.STRtree;
import com.vividsolutions.jump.feature.Feature;
import com.vividsolutions.jump.feature.FeatureCollection;
import com.vividsolutions.jump.feature.FeatureDataset;
import com.vividsolutions.jump.feature.FeatureSchema;
import com.vividsolutions.jump.io.DriverProperties;
import com.vividsolutions.jump.io.IllegalParametersException;
import com.vividsolutions.jump.io.ShapefileReader;
import com.vividsolutions.jump.io.ShapefileWriter;
import com.vividsolutions.jump.io.WKTReader;
import com.vividsolutions.jump.io.WKTWriter;
import com.vividsolutions.jump.io.geojson.GeoJSONReader;
import com.vividsolutions.jump.io.geojson.GeoJSONWriter;

public class JUMPUtility {
	/***************************
	 * writers
	 **************************/
	/**
	 * 
	 * @param fc
	 * @param sfn
	 */
	public static void exportGeoJSON(FeatureCollection fc, String sfn) {
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		try {
			GeoJSONWriter gjWriter = new GeoJSONWriter();
			gjWriter.write(fc, dp);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	//
	public static void exportShapeFile(FeatureCollection fc, String sfn) {
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		try {
			ShapefileWriter shpWriter = new ShapefileWriter();
			shpWriter.write(fc, dp);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	//
	public static void exportWKTFile(FeatureCollection fc, String sfn) {
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		try {
			WKTWriter wktWriter = new WKTWriter();
			wktWriter.write(fc, dp);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	//
	/******************************
	 * 
	 * Readers
	 * 
	 ******************************/
	public static FeatureCollection loadWKT(String sfn) {
		FeatureCollection featColl = null;
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		try {
			WKTReader wktRead = new WKTReader();
			featColl = wktRead.read(dp);
		} catch (Exception e) {
			System.out.println("\tReading from WKT file error " + e);
		}
		return featColl;
	}

	public static FeatureCollection loadGeoJSON(String sfn) {
		FeatureCollection featColl = null;
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		try {
			GeoJSONReader gjRead = new GeoJSONReader();
			featColl = gjRead.read(dp);
		} catch (Exception e) {
			System.out.println("\tReading from GeoJSON file error " + e);
		}
		return featColl;
	}

	//
	public static List<Feature> loadGeoJSON(String sfn, boolean rmvDuplicate, String dupKey) {
		FeatureCollection featColl = loadGeoJSON(sfn);
		TreeMap<String, Feature> featIdx = null;
		if (featColl != null) {
			if (rmvDuplicate) {
				if (dupKey != null) {
					FeatureSchema sch = featColl.getFeatureSchema();
					if (sch.hasAttribute(dupKey)) {
						featIdx = new TreeMap<String, Feature>();
						List<Feature> features = featColl.getFeatures();
						// List<Feature> featToRmv = new LinkedList<Feature>();
						Iterator<Feature> iter = features.iterator();
						while (iter.hasNext()) {
							Feature feat = iter.next();
							Geometry geom = feat.getGeometry();
							String toid = (String) feat.getAttribute(dupKey);
							Feature existFeat = null;
							if ((existFeat = (Feature) featIdx.get(toid)) == null) {
								featIdx.put(toid, feat);
							} else {// multiple copies
								Geometry existGeom = existFeat.getGeometry();
								if (existGeom.compareTo(geom) != 0) {// geometry
																		// not
																		// identical,
																		// shouldn't
																		// happen
									Geometry newgeom = existGeom.union(geom);
									existFeat.setGeometry(newgeom);
								}
							}
						}
						FeatureCollection featCollNew = new FeatureDataset(sch);
						featCollNew.addAll(featIdx.values());
						featIdx.clear();
						return featCollNew.getFeatures();
					}
				}
			}
			return featColl.getFeatures();
		} else {
			return null;
		}
	}
	public static FeatureCollection loadShapeFile(String sfn) {
		FeatureCollection featColl = null;
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		TreeMap<String, Feature> featIdx = null;
		try {
			ShapefileReader shpRead = new ShapefileReader();
			featColl = shpRead.read(dp);
		} catch (Exception e) {
			System.out.println("\tReading from shapefile to file error " + e);
		}
		return featColl;
	}
	//
	public static List<Feature> loadShapeFile(String sfn, boolean rmvDuplicate, String dupKey) {
		FeatureCollection featColl = null;
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		TreeMap<String, Feature> featIdx = null;
		try {
			ShapefileReader shpRead = new ShapefileReader();
			featColl = shpRead.read(dp);
		} catch (Exception e) {
			System.out.println("\tReading from shapefile to file error " + e);
		}
		if (featColl != null) {
			if (rmvDuplicate) {
				if (dupKey != null) {
					FeatureSchema sch = featColl.getFeatureSchema();
					if (sch.hasAttribute(dupKey)) {
						featIdx = new TreeMap<String, Feature>();
						List<Feature> features = featColl.getFeatures();
						// List<Feature> featToRmv = new LinkedList<Feature>();
						Iterator<Feature> iter = features.iterator();
						while (iter.hasNext()) {
							Feature feat = iter.next();
							Geometry geom = feat.getGeometry();
							String toid = (String) feat.getAttribute(dupKey);
							Feature existFeat = null;
							if ((existFeat = (Feature) featIdx.get(toid)) == null) {
								featIdx.put(toid, feat);
							} else {// multiple copies
								Geometry existGeom = existFeat.getGeometry();
								if (existGeom.compareTo(geom) != 0) {// geometry
																		// not
																		// identical,
																		// shouldn't
																		// happen
									Geometry newgeom = existGeom.union(geom);
									existFeat.setGeometry(newgeom);
								}
							}
						}
						FeatureCollection featCollNew = new FeatureDataset(sch);
						featCollNew.addAll(featIdx.values());
						featIdx.clear();
						return featCollNew.getFeatures();
					}
				}
			}
			return featColl.getFeatures();
		} else {
			return null;
		}
	}

	//
	public static List<Feature> loadShapeFile(String sfn, double minx, double miny, double maxx, double maxy,
			boolean rmvDuplicate, String dupKey) {
		FeatureCollection featColl = null;
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		TreeMap<String, Feature> featIdx = null;
		try {
			ShapefileReader shpRead = new ShapefileReader();
			featColl = shpRead.read(dp);
		} catch (Exception e) {
			System.out.println("\tReading from shapefile to file error " + e);
		}
		Geometry extent = null;
		if (featColl != null) {
			FeatureSchema sch = featColl.getFeatureSchema();
			FeatureCollection featCollNew = new FeatureDataset(sch);
			GeometryFactory gf = null;
			List<Feature> features = featColl.getFeatures();
			if (rmvDuplicate && dupKey != null && sch.hasAttribute(dupKey)) {
				featIdx = new TreeMap<String, Feature>();
				// List<Feature> featToRmv = new LinkedList<Feature>();
				Iterator<Feature> iter = features.iterator();
				while (iter.hasNext()) {
					Feature feat = iter.next();
					Geometry geom = feat.getGeometry();
					if (gf == null) {
						gf = geom.getFactory();
						extent = JTSUtility.rect2Geometry(minx, miny, maxx, maxy, gf);
					}
					if (extent.intersects(geom)) {
						String toid = (String) feat.getAttribute(dupKey);
						Feature existFeat = null;
						if ((existFeat = (Feature) featIdx.get(toid)) == null) {
							featIdx.put(toid, feat);
						} else {// multiple copies
							Geometry existGeom = existFeat.getGeometry();
							if (existGeom.compareTo(geom) != 0) {// geometry not
																	// identical,
																	// shouldn't
																	// happen
								Geometry newgeom = existGeom.union(geom);
								existFeat.setGeometry(newgeom);
							}
						}
					}
				}
				featCollNew.addAll(featIdx.values());
				featIdx.clear();
			} else {
				for (Feature feat : features) {
					Geometry geom = feat.getGeometry();
					if (gf == null) {
						gf = geom.getFactory();
						extent = JTSUtility.rect2Geometry(minx, miny, maxx, maxy, gf);
					}
					if (extent.intersects(geom)) {
						featCollNew.add(feat);
					}
				}
			}
			return featCollNew.getFeatures();
		} else {
			return null;
		}
	}

	public static List<Feature> loadShapeFile(String sfn, String outFn, String dupKey, double minx, double miny,
			double maxx, double maxy) {
		FeatureCollection featColl = null;
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		TreeMap<String, Feature> featIdx = null;
		try {
			ShapefileReader shpRead = new ShapefileReader();
			featColl = shpRead.read(dp);
		} catch (Exception e) {
			System.out.println("\tReading from shapefile to file error " + e);
		}
		Geometry extent = null;
		if (featColl != null) {
			FeatureSchema sch = featColl.getFeatureSchema();
			FeatureCollection featCollNew = new FeatureDataset(sch);
			GeometryFactory gf = null;
			List<Feature> features = featColl.getFeatures();
			if (dupKey != null && sch.hasAttribute(dupKey)) {
				featIdx = new TreeMap<String, Feature>();
				// List<Feature> featToRmv = new LinkedList<Feature>();
				Iterator<Feature> iter = features.iterator();
				while (iter.hasNext()) {
					Feature feat = iter.next();
					Geometry geom = feat.getGeometry();
					if (gf == null) {
						gf = geom.getFactory();
						extent = JTSUtility.rect2Geometry(minx, miny, maxx, maxy, gf);
					}
					if (extent.intersects(geom)) {
						String toid = (String) feat.getAttribute(dupKey);
						Feature existFeat = null;
						if ((existFeat = (Feature) featIdx.get(toid)) == null) {
							featIdx.put(toid, feat);
						} else {// multiple copies
							Geometry existGeom = existFeat.getGeometry();
							if (existGeom.compareTo(geom) != 0) {// geometry not
																	// identical,
																	// shouldn't
																	// happen
								Geometry newgeom = existGeom.union(geom);
								existFeat.setGeometry(newgeom);
							}
						}
					}
				}
				featCollNew.addAll(featIdx.values());
				featIdx.clear();
			} else {
				for (Feature feat : features) {
					Geometry geom = feat.getGeometry();
					if (gf == null) {
						gf = geom.getFactory();
						extent = JTSUtility.rect2Geometry(minx, miny, maxx, maxy, gf);
					}
					if (extent.intersects(geom)) {
						featCollNew.add(feat);
					}
				}

			}
			dp.set("DefaultValue", outFn);
			ShapefileWriter shpWriter = new ShapefileWriter();
			try {
				shpWriter.write(featCollNew, dp);
				System.out.println(featCollNew.size() + " features exported...");
			} catch (IllegalParametersException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return featCollNew.getFeatures();
		} else {
			return null;
		}
	}

	// load shape file, remove duplicated features and save to outFn
	public static List<Feature> loadShapeFile(String sfn, String fid_nm, boolean rmvDup, String outFn) {
		FeatureCollection featColl = null;
		DriverProperties dp = new DriverProperties();
		dp.set("DefaultValue", sfn);
		TreeMap<String, Feature> featIdx = null;
		try {
			ShapefileReader shpRead = new ShapefileReader();
			featColl = shpRead.read(dp);
		} catch (Exception e) {
			System.out.println("\tReading from shapefile to file error " + e);
		}
		if (featColl != null) {
			System.out.println(featColl.size() + " features loaded...");
			if (rmvDup) {
				featIdx = new TreeMap<String, Feature>();
				List<Feature> features = featColl.getFeatures();
				// List<Feature> featToRmv = new LinkedList<Feature>();
				Iterator<Feature> iter = features.iterator();
				while (iter.hasNext()) {
					Feature feat = iter.next();
					Geometry geom = feat.getGeometry();
					String fid = (String) feat.getAttribute(fid_nm);
					Feature existFeat = null;
					if ((existFeat = (Feature) featIdx.get(fid)) == null) {
						featIdx.put(fid, feat);
					} else {// multiple copies
						Geometry existGeom = existFeat.getGeometry();
						if (existGeom.compareTo(geom) != 0) {// geometry not
																// identical,
																// shouldn't
																// happen
							Geometry newgeom = existGeom.union(geom);
							existFeat.setGeometry(newgeom);
						}
					}
				}
				FeatureCollection featCollNew = new FeatureDataset(featColl.getFeatureSchema());
				featCollNew.addAll(featIdx.values());
				featIdx.clear();
				if (outFn != null) {
					// save to new shapefile
					dp.set("DefaultValue", outFn);
					ShapefileWriter shpWriter = new ShapefileWriter();
					try {
						shpWriter.write(featCollNew, dp);
						System.out.println(featCollNew.size() + " features exported...");
					} catch (IllegalParametersException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				return featCollNew.getFeatures();
			} else {
				return featColl.getFeatures();
			}
		} else {
			System.out.println(sfn + ": no feature loaded...exit...");
			return null;
		}
	}

	/*******************
	 * 
	 * Geometric operations on FeatureCollections
	 * 
	 ******************/
	// convex hull of a feature collection
	public static Geometry compConvexHullFC(FeatureCollection fc) {
		List<Feature> feats = fc.getFeatures();
		LinkedList<Coordinate> coords = new LinkedList<Coordinate>();
		Iterator<Feature> iter = feats.iterator();
		// graph
		while (iter.hasNext()) {
			Feature feat = iter.next();
			Geometry geom = feat.getGeometry();
			Geometry hull = geom.convexHull();
			for (Coordinate coord : hull.getCoordinates()) {
				coords.add(coord);
			}
		}
		GeometryFactory factory = new GeometryFactory();
		MultiPoint mp = factory.createMultiPoint((Coordinate[]) coords.toArray());
		feats.clear();
		return mp.convexHull();
	}
	//
	public static void main(String[] args) {
		FeatureCollection fc = loadShapeFile("d:/temp/as02.shp");
		exportGeoJSON(fc, "d:/temp/as02.json");
		System.out.println("done");
	}
}
