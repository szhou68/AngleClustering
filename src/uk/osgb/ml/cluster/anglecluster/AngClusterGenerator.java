/**
 *  utility class to provide methods for angle value clustering
 * 
 */
package uk.osgb.ml.cluster.anglecluster;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.Vector;
import java.util.Map.Entry;

import java.util.Set;
import java.util.TreeMap;

import uk.osgb.datastructures.MultiTreeMap;
import uk.osgb.utilities.AngleUtility;
import uk.osgb.utilities.GeomUtility;

public class AngClusterGenerator<T> {
	//
	// total_weight-cluster map to hold the results of clustering
	//
	MultiTreeMap<Double, AngCluster<T>> clusters = new MultiTreeMap<Double, AngCluster<T>>();
	// sum of variances
	double varianceSum = 0.0;
	// sum of weights
	double weightSum = 0.0;
	// sum of weights projected to cluster means
	double projWeightSum = 0.0;
	//
	MergeEvaluator evaluator = null;
	TerminateEvaluator terminator = null;
	//
	/**
	 * 
	 */
	public AngClusterGenerator(){
		evaluator = new MergeEvaluatorMinVar();
	}
	/**
	 * @param eval
	 */
	public AngClusterGenerator(MergeEvaluator eval){
		if(eval!=null){
			evaluator = eval;
		}else{
			eval = new MergeEvaluatorMinVar();
		}
	}
	/**
	 * @param ref
	 * @return
	 */
	public Collection<AngCluster<T>> getCluster(AngObject<T> ref){
		if(clusters!=null){
			return clusters.get(ref.getAng());
		}
		return null;
	}
	/**
	 * @return all clusters in a collection
	 */
	public Collection<AngCluster<T>> getAllClusters(){
		Collection<Entry<Double, AngCluster<T>>> entries = clusters.entrySet();
		Vector<AngCluster<T>> rlt = new Vector<AngCluster<T>>(entries.size());
		for(Entry<Double, AngCluster<T>> entry:entries){
			rlt.add(entry.getValue());
		}
		return rlt;
	}
	//
	/**
	 * @param ref
	 * @return 
	 */
	public AngCluster<T> getFloorCluster(AngObject<T> ref){
		if(clusters!=null){
			Double floorkey = clusters.floorKey(ref.getAng());
			Collection<AngCluster<T>> clus;
			if(floorkey!=null){
				clus = clusters.get(floorkey);
			}else{
				clus = clusters.lastValue();
			}
			AngCluster<T> rlt = null;
			for(AngCluster<T> clu:clus){
				if(rlt == null){
					rlt = clu;
				}else if(clu.compareTo(rlt) > 0){ // clu is lower
					rlt = clu;
				}
			}
			return rlt;
		}
		return null;
	}
	
	/**
	 * @param ref
	 * @return
	 */
	public AngCluster<T> getCeilingCluster(AngObject<T> ref){
		if(clusters!=null){
			Double ceilingkey = clusters.ceilingKey(ref.getAng());
			Collection<AngCluster<T>> clus;
			if(ceilingkey!=null){
				clus = clusters.get(ceilingkey);
			}else{
				clus = clusters.firstValue();
			}
			AngCluster<T> rlt = null;
			for(AngCluster<T> clu:clus){
				if(rlt == null){
					rlt = clu;
				}else if(clu.compareTo(rlt) < 0){ // clu is higher
					rlt = clu;
				}
			}
			return rlt;
		}
		return null;
	}
	public TerminateEvaluator getTerminator() {
		return terminator;
	}
	public void setTerminator(TerminateEvaluator terminator) {
		this.terminator = terminator;
	}
	//
	/**
	 * @return
	 */
	public double compClusterWeightSum(){
		double wSum = 0.0;
		Collection<AngCluster<T>> clus = new Vector<AngCluster<T>>();
		clusters.values(clus);
		for(AngCluster<T> clu:clus){
			wSum+=clu.compWeightSum();
		}
		return weightSum = wSum;
	}
	/**
	 * @return
	 */
	public double compClusterProjectedSum(){
		double wSum = 0.0;
		Collection<AngCluster<T>> clus = new Vector<AngCluster<T>>();
		clusters.values(clus);
		for(AngCluster<T> clu:clus){
			wSum+=clu.compProjectedWeightSum();
		}
		
		return projWeightSum = wSum;
	}
	
	/**
	 * @param cluSet
	 * @return
	 */
	public static double compWithinCluserVariation(Collection<AngCluster> cluSet){
		double vSum = 0.0;
		for(AngCluster clu:cluSet){
			vSum+=clu.compVariation();
		}
		return vSum;
	}
	/**
	 * @return
	 */
	public double compWithinClusterVariationSum(){
		double vSum = 0.0;
		Collection<AngCluster<T>> clus = new Vector<AngCluster<T>>();
		clusters.values(clus);
		for(AngCluster<T> clu:clus){
			vSum+=clu.compVariation();
		}
		return vSum;
	}
	/**
	 * @param cluSet
	 * @return
	 */
	public static double getWithinCluserVariation(Collection<AngCluster> cluSet){
		double vSum = 0.0;
		for(AngCluster clu:cluSet){
			vSum+=clu.getVariance();
		}
		return vSum;
	}
	//
	/**
	 * @return
	 */
	public double getWithinClusterVariationSum(){
		double vSum = 0.0;
		Collection<AngCluster<T>> clus = new Vector<AngCluster<T>>();
		clusters.values(clus);
		for(AngCluster<T> clu:clus){
			vSum+=clu.getVariance();
		}
		return vSum;
	}
	/**
	 * 
	 */
	public void compClusterAngMean(){
		Collection<AngCluster<T>> clus = new Vector<AngCluster<T>>();
		clusters.values(clus);
		for(AngCluster<T> clu:clus){
			clu.compClusterMean();
		}
	}
	/**
	 * 
	 */
	public void clearClusterMembers(){
		Collection<AngCluster<T>> clus = new Vector<AngCluster<T>>();
		clusters.values(clus);
		for(AngCluster<T> clu:clus){
			clu.clearAngRefs();;
		}
	}
	/** compute the from and to AngRef objects for each cluster
	 * 
	 */
	public void compFTRefs(){
		Collection<AngCluster<T>> clus = new Vector<AngCluster<T>>();
		clusters.values(clus);
		for(AngCluster<T> clu:clus){
			clu.compFromToAngles();
		}
	}
	//
	/**
	 * @param angles
	 * @param weights
	 * @param varTol
	 * @param modulo
	 * @param tiltAng
	 * @return
	 */
	public static Collection<AngCluster> clusteringAngleArray(double[] angles, double[] weights, double varTol, double maxAW, int modulo, double tiltAng){
		if(angles==null || angles.length == 0){
			return null;
		}

		int numPts = angles.length;
		double[][] pts = new double[numPts][2];
		double wSum = 1.0;
		//
		//double wSamp = 1.0/numPts; // normalised weight for each sample if individual weight not specified
		for(int i = 0; i < numPts; ++i){
			pts[i][0] = angles[i];
			pts[i][1] = weights==null?1.0:weights[i];
		}
		Vector<AngRef> angRefs = new Vector<AngRef>(numPts);
		for(int i = 0; i < numPts; ++i){
			AngRef ref = new AngRef(pts[i][0], pts[i][1]); 
			angRefs.add(ref);
		}
		Collection<AngCluster> clus = clusteringAngleRefCollectionModulo(angRefs, varTol, maxAW, tiltAng, modulo);
		return clus;
	}
	/** 
	 * @param angles
	 * @param weights
	 * @param varTol
	 * @param modulo
	 * @param tiltAng
	 * @param normWeight
	 * @return
	 */
	public static Collection<AngCluster> clusteringAngleArrayElbow(double[] angles, double[] weights, int modulo, double tiltAng){
		if(angles==null || angles.length == 0){
			return null;
		}
		int numPts = angles.length;
		double[][] pts = new double[numPts][2];
		double wSum = 1.0;
		//
		//double wSamp = 1.0/numPts; // normalised weight for each sample if individual weight not specified
		for(int i = 0; i < numPts; ++i){
			pts[i][0] = angles[i];
			pts[i][1] = weights==null?1.0:weights[i];
		}
		Vector<AngRef> angRefs = new Vector<AngRef>(numPts);
		for(int i = 0; i < numPts; ++i){
			AngRef ref = new AngRef(pts[i][0], pts[i][1]); 
			angRefs.add(ref);
		}
		Collection<AngCluster> clus = clusteringAngleRefCollectionModuloElbow(angRefs, tiltAng, modulo);
		return clus;
	}
	/** clustering angle references, controlled by remaining variance 
	 * @param angRefs input collection of angle reference
	 * @param varTol ratio of remaining variance over total variance, in [0, 1)
	 * @return a collection of angle clusters
	 */
	public static Collection<AngCluster> clusteringAngleRefCollection(Collection<AngRef> angRefs, double varTol, double maxAWDeg){
		AngClusterGenerator gen;
		int numPts = angRefs.size();
		boolean zeroVar = false;
		//Coordinate[] ls = new Coordinate[numPts];
		int clsId = -1;
		double totalVar = 0.0;
		double varThreshold = 0.0;
		gen = new AngClusterGenerator();
		totalVar = compRefColAngleVariation(angRefs, null, null);
		TerminateEvaluator te;
		if(varTol < 1.0 && maxAWDeg >= 360) {
			te = new TerminateEvaluatorVar(totalVar, varTol);
		}else {
			te = new TerminateEvaluatorMaxAW(totalVar, varTol, maxAWDeg);
		}
		gen.setTerminator(te);
		gen.clusteringAnglesHierarchical(angRefs, 0);
		return gen.getAllClusters();
	}
	/** clustering angle references, controlled by remaining variance. All references will be turned an angle prior to clustering
	 * @param angRefs
	 * @param varTol
	 * @param turnAng 
	 * @return a collection of angle clusters
	 */
	public static Collection<AngCluster> clusteringAngleRefCollection(Collection<AngRef> angRefs, double varTol, double maxAW, double turnAng, int modulo){
		if(turnAng!=0.0){
			for(AngRef ref:angRefs){
				ref.turnAng(turnAng);
			}
		}
		ArrayList<AngRef> refMod = new ArrayList<AngRef>(angRefs.size());
		for(AngRef ref:angRefs){
			refMod.add(ref.modulo(modulo));
		}
		Collection<AngCluster> clsters = clusteringAngleRefCollection(refMod, varTol, maxAW);
		if(turnAng!=0.0){// restore to original angle values
			for(AngRef ref:refMod){
				ref.turnAng(-turnAng);
			}
		}
		//
		for(AngCluster clster:clsters){
			clster.compClusterMean();
			//clster.compProjectedSum(); // this shouldn't change
		}
		return clsters;
	}
	
	/** clustering angle references, controlled by remaining variance. Modulo operation is applied to map the angle values to
	 *  [0, PI) )if modulo == 2) or [0, PI/2) (if modulo == 4). Therefore, the resulting clusters do not containing the original
	 *  angle references any more. Also two rounds of clustering are carried out. During the second one all angle references are tilted
	 *  by a given value and the result with fewer clusters will be returned. The tilting operation is for handling angles near 0 and PI 
	 *  (in case of modulo 2) or near 0, +/- PI/2 and PI
	 * @param angRefs angle references to be clustered
	 * @param varTol variance tolerance to control the clustering process. The value represents the maximum percentage of remaining variance 
	 * @param tiltAng
	 * @param modulo 1 for original, 2 for modulo PI and 4 for modulo PI/2
	 * @return
	 */
	public static Collection<AngCluster> clusteringAngleRefCollectionModulo(Collection<AngRef> angRefs, double varTol, double maxAW, double tiltAng, int modulo){
		if(modulo>1){
			// do modulo
			ArrayList<AngRef> refMod = new ArrayList<AngRef>(angRefs.size());
			for(AngRef ref:angRefs){
				refMod.add(ref.modulo(modulo));
			}
			// find the longest vector and use it as base?
			// do clustering
			Collection<AngCluster> clsters1 = clusteringAngleRefCollection(refMod, varTol, maxAW); 
			// tilt then modulo and clustering
			Collection<AngCluster> clsters2 = clusteringAngleRefCollection(angRefs, varTol, maxAW, tiltAng, modulo);
			if(clsters1.size() <= clsters2.size()){
				return clsters1;
			}else{
				return clsters2;
			}
		}else{
			return clusteringAngleRefCollection(angRefs, varTol, maxAW);
		}
	}
	/**cluster angle reference collection using hierarchical method, merge and termination are controlled by evaulators supplied in generator
	 * resulting clusters are stored in this.clusters
	 * @param angRefs
	 * @param numClusters maximum number of clusters to be generated; if <=0, it will be decided by varTol 
	 * @param varTol ratio of within-cluster variation over total variation
	 * @return
	 */
	private boolean clusteringAnglesHierarchical(Collection<AngRef<T>> angRefs, int numClusters){
		if(angRefs.size() > 0){
			// generate a cluster for each individual angle value
			TreeMap<Double, AngCluster<T>> clus = new TreeMap<Double, AngCluster<T>>();
			for(AngRef<T> ref:angRefs){
				double ang = ref.getAng();
				AngCluster<T> clu = clus.get(ang);
				if(clu!=null){// cluster exists
					clu.addAngleRef(ref);
				}else{ // new cluster
					clu= new AngCluster<T>();
					clu.addAngleRef(ref); // mean angle initialised here 
					clus.put(ref.getAng(), clu);
				}
			}
			// entry set
			Set<Map.Entry<Double, AngCluster<T>>> cluEntrySet = clus.entrySet();
			int cnt = 0;
			int numClu = cluEntrySet.size();
			
			for(Map.Entry<Double, AngCluster<T>> cluEntry:cluEntrySet){
				AngCluster clu = cluEntry.getValue();
				clu.compFromToAngles();
				clu.compClusterMean();
				clu.compWeightSum();
				clu.compProjectedWeightSum();
			}
			double totalVar = compRefColAngleVariation(angRefs);
			varianceSum = totalVar;
			if(numClu == 1 || numClu <= numClusters){// only 1 cluster or no more than maximum number of clusters to be generated
				for(Map.Entry<Double, AngCluster<T>> cluEntry :cluEntrySet){
					clusters.put(cluEntry.getKey(), cluEntry.getValue());
				}
				return true;
			}
			//
			TreeSet<AngCluster<T>> cluSet = new TreeSet<AngCluster<T>>(); // index of current clusters
			// more than one cluster, link them up to form a cyclic chain
			AngCluster firstClu = null;
			AngCluster prevClu = null;
			for(Map.Entry<Double, AngCluster<T>> cluEntry:cluEntrySet){// ordered by weighted mean of cluster
				cnt++;
				AngCluster clu = cluEntry.getValue();
				cluSet.add(clu);
				if(cnt == 1){// first cluster
					firstClu = clu;
					prevClu = clu;
				}else{// following clusters
					clu.setPrev(prevClu);
					prevClu.setNext(clu);
					prevClu = clu;
					if(cnt==numClu){ // last cluster
						clu.setNext(firstClu);
						firstClu.setPrev(clu);
					}
				}
			}
			// merge test result, in ascending order 
			TreeSet<MergeEvalRlt> rltSet = new TreeSet<MergeEvalRlt>();
			TreeMap<AngCluster, MergeEvalRlt> rltIdx = new TreeMap<AngCluster, MergeEvalRlt>();
			for(AngCluster clu:cluSet){
				MergeEvalRlt rlt = evaluator.evaluate(clu, clu.getNext());
				rltSet.add(rlt);
				rltIdx.put(rlt.getCluster(), rlt);
			}
			//System.out.println("total Variance:" + totalVar);
			//
			while(cluSet.size() > 0){
				//double inVarSum = this.compWithinCluserVariation(cluSet);
				//System.out.println(Double.toString(inVarSum) + "\t"+ cluSet.size());
				MergeEvalRlt minRlt = null; 
				boolean done = false;
				// checkf for goal:
				if(numClusters > 0 && cluSet.size() <=numClusters){// by number of clusters, 
					done = true;
				}else if(cluSet.size() == 1){
					done = true;
				}else{
					while(!rltSet.isEmpty()){
						minRlt = rltSet.pollFirst();
						if(terminator.terminate(cluSet, minRlt)){
							if(rltSet.isEmpty()){
								done = true;
							}else{
								rltIdx.remove(minRlt.getCluster(), minRlt);
								continue;
							}
						}else{
							break;
						}
					}
					if(minRlt == null){
						done = true;
					}
				}
				if(done){
					for(AngCluster clu:cluSet){
						clusters.put(clu.getWeight(), clu);
					}
					return true;
				}else{
					// merge
					AngCluster clu = minRlt.getCluster();
					AngCluster cluNext = clu.getNext();

					if(clu!=cluNext){
						cluSet.remove(clu);
						cluSet.remove(cluNext);
						MergeEvalRlt rltNext = rltIdx.get(cluNext);
						rltSet.remove(rltNext);

						rltIdx.remove(clu);
						rltIdx.remove(cluNext);

						clu.merge(cluNext);

						cluSet.add(clu);
						cluNext = clu.getNext();
						// reset affected merge evaluation after merge
						if(clu!=cluNext){
							MergeEvalRlt rltNew = evaluator.evaluate(clu, cluNext);
							rltSet.add(rltNew);
							rltIdx.put(clu, rltNew);
						}
						AngCluster cluPrev = clu.getPrev();
						if(cluPrev!= clu){
							MergeEvalRlt rltPrev = rltIdx.get(cluPrev);
							rltSet.remove(rltPrev);
							rltIdx.remove(cluPrev);
							rltPrev = evaluator.evaluate(cluPrev, clu);
							rltSet.add(rltPrev);
							rltIdx.put(cluPrev, rltPrev);
						}
					}

				}
			}
			return true;
		}else{
			return false;
		}
	}
/************************* 
 * 
 * Testing Elbow/Silhouette 
 * 	
 *************************/
	//
	public static Collection<AngCluster> clusteringAngleRefCollectionElbow(Collection<AngRef> angRefs){
		AngClusterGenerator gen;
		int numPts = angRefs.size();
		boolean zeroVar = false;
		//Coordinate[] ls = new Coordinate[numPts];
		int clsId = -1;
		double totalVar = 0.0;
		double varThreshold = 0.0;
		gen = new AngClusterGenerator();
		totalVar = compRefColAngleVariation(angRefs, null, null);
		
		TerminateEvaluator te = new TerminateEvaluatorVar(totalVar, 1.0); // dummy terminator
		
		gen.setTerminator(te);
		gen.clusteringAnglesElbow(angRefs);
		
		return gen.getAllClusters();
	}
	
	public static Collection<AngCluster> clusteringAngleRefCollectionElbow(Collection<AngRef> angRefs, double turnAng, int modulo){
		if(turnAng!=0.0){
			for(AngRef ref:angRefs){
				ref.turnAng(turnAng);
			}
		}
		ArrayList<AngRef> refMod = new ArrayList<AngRef>(angRefs.size());
		for(AngRef ref:angRefs){
			refMod.add(ref.modulo(modulo));
		}
		Collection<AngCluster> clsters = clusteringAngleRefCollectionElbow(refMod);
		if(turnAng!=0.0){// restore to original angle values
			for(AngRef ref:refMod){
				ref.turnAng(-turnAng);
			}
		}
		//
		for(AngCluster clster:clsters){
			clster.compClusterMean();
			//clster.compProjectedSum(); // this shouldn't change
		}
		return clsters;
	}

	
	/**
	 * @param angRefs
	 * @param tiltAng
	 * @param modulo
	 * @return
	 */
	public static Collection<AngCluster> clusteringAngleRefCollectionModuloElbow(Collection<AngRef> angRefs, double tiltAng, int modulo){
		if(modulo>1){
			// do modulo
			ArrayList<AngRef> refMod = new ArrayList<AngRef>(angRefs.size());
			for(AngRef ref:angRefs){
				refMod.add(ref.modulo(modulo));
			}
			// do clustering
			Collection<AngCluster> clsters1 = clusteringAngleRefCollectionElbow(refMod); 
			// tilt then modulo and clustering
			Collection<AngCluster> clsters2 = clusteringAngleRefCollectionElbow(angRefs, tiltAng, modulo);
			if(clsters1.size() <= clsters2.size()){
				return clsters1;
			}else{
				return clsters2;
			}
		}else{// no modulo
			return clusteringAngleRefCollectionElbow(angRefs);
		}
	}

	/** for testing elbow method, no cluster is outputted yet
	 * @param angRefs
	 * @return
	 */
	private boolean clusteringAnglesElbow(Collection<AngRef> angRefs){
		if(angRefs.size() > 0){
			// generate a cluster for each individual angle value
			TreeMap<Double, AngCluster<T>> clus = new TreeMap<Double, AngCluster<T>>();
			for(AngRef<T> ref:angRefs){
				double ang = ref.getAng();
				AngCluster<T> clu = clus.get(ang);
				if(clu!=null){// cluster exists
					clu.addAngleRef(ref);
				}else{ // new cluster
					clu= new AngCluster<T>();
					clu.addAngleRef(ref); // mean angle initialised here 
					clus.put(ref.getAng(), clu);
				}
			}
			// entry set
			Set<Map.Entry<Double, AngCluster<T>>> cluEntrySet = clus.entrySet();
			int cnt = 0;
			int numClu = cluEntrySet.size();
			
			for(Map.Entry<Double, AngCluster<T>> cluEntry:cluEntrySet){
				AngCluster clu = cluEntry.getValue();
				clu.compFromToAngles();
				clu.compClusterMean();
				clu.compWeightSum();
				clu.compProjectedWeightSum();
			}
			double totalVar = compRefColAngleVariation(angRefs);
			varianceSum = totalVar;
			TreeSet<AngCluster> cluSet = new TreeSet<AngCluster>(); // index of current clusters
			// more than one cluster, link them up to form a cyclic chain
			AngCluster firstClu = null;
			AngCluster prevClu = null;
			for(Map.Entry<Double, AngCluster<T>> cluEntry:cluEntrySet){
				cnt++;
				AngCluster clu = cluEntry.getValue();
				cluSet.add(clu);
				if(cnt == 1){// first cluster
					firstClu = clu;
					prevClu = clu;
				}else{// following clusters
					clu.setPrev(prevClu);
					prevClu.setNext(clu);
					prevClu = clu;
					if(cnt==numClu){ // last cluster
						clu.setNext(firstClu);
						firstClu.setPrev(clu);
					}
				}
			}
			// merge test result, in ascending order 
			TreeSet<MergeEvalRlt> rltSet = new TreeSet<MergeEvalRlt>();
			TreeMap<AngCluster, MergeEvalRlt> rltIdx = new TreeMap<AngCluster, MergeEvalRlt>();
			for(AngCluster clu:cluSet){
				MergeEvalRlt rlt = evaluator.evaluate(clu, clu.getNext());
				rltSet.add(rlt);
				rltIdx.put(rlt.getCluster(), rlt);
			}
			System.out.println("total Variance:" + totalVar);
			//
			Double[] varArray = new Double[cluSet.size()];
			Double[] silArray = new Double[cluSet.size()];
			
			for(int i = 0; i < varArray.length; ++i){
				varArray[i] = getWithinCluserVariation(cluSet);
				//silArray[i] = AngClusterGenerator.computeAverageSilhouette(cluSet);
				silArray[i] = AngClusterGenerator.computeAverageSilhouetteWeighted(cluSet);
				//report(cluSet);
				MergeEvalRlt minRlt = rltSet.pollFirst();
				if(minRlt!=null){
					// merge
					AngCluster clu = minRlt.getCluster();
					AngCluster cluNext = clu.getNext();

					if(clu!=cluNext){
						cluSet.remove(clu);
						cluSet.remove(cluNext);
						MergeEvalRlt rltNext = rltIdx.get(cluNext);
						rltSet.remove(rltNext);

						rltIdx.remove(clu);
						rltIdx.remove(cluNext);

						clu.merge(cluNext);

						cluSet.add(clu);
						cluNext = clu.getNext();
						// reset affected merge evaluation after merge
						if(clu!=cluNext){
							MergeEvalRlt rltNew = evaluator.evaluate(clu, cluNext);
							rltSet.add(rltNew);
							rltIdx.put(clu, rltNew);
						}
						AngCluster cluPrev = clu.getPrev();
						if(cluPrev!= clu){
							MergeEvalRlt rltPrev = rltIdx.get(cluPrev);
							rltSet.remove(rltPrev);
							rltIdx.remove(cluPrev);
							rltPrev = evaluator.evaluate(cluPrev, clu);
							rltSet.add(rltPrev);
							rltIdx.put(cluPrev, rltPrev);
						}
					}
				}
			}
			// find the optimial number of clusters
			System.out.println("Total within-cluster variance and Average Silhouette: ");
			for(int i = 0; i < varArray.length; ++i){
				System.out.println((varArray.length - i) + ", " + varArray[i] + ", " + silArray[i]);
			}
			//System.out.println("symmetric difference quotient: ");
			double sf = varArray.length-1;
			System.out.println("difference: ");
			for(int i = 1; i < varArray.length - 1; ++i){
				double x0 = (i-1);
				double y0 = varArray[i-1] * sf;
				double x1 = i;
				double y1 = varArray[i]*sf;
				double x2 = i+1;
				double y2 = varArray[i+1]*sf;
				double area2 = GeomUtility.triAreaDouble2(x0, y0, x1, y1, x2, y2);
				double baseLen = Math.sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0));
				double offset = area2 / baseLen;
				//System.out.println((varArray.length - i) + ", " + offset);
				//double der = (varArray[i+1] - varArray[i-1])*0.5;
				//System.out.println((varArray.length - i) + ", " + der);
			}
			//
			return true;
		}else{
			return false;
		}
	}
/*
 *  angle mean
 * 
 */
	public static double compRefColAngleMeanArithmetic(Collection<? extends AngObject> refs, AngObject angF, AngObject angT){
		if(refs.isEmpty()){
			return Double.MAX_VALUE;// not initilaised
		}
		if(refs.size() == 1){
			return refs.iterator().next().getAng();
		}
		double awSum = 0.0;
		double wSum = 0.0;
		double angBase = angF.getAng();
		for(AngObject ref:refs){
			double ang = ref.getAng();
			double angDiff = AngleUtility.angFromTo(angBase, ref.getAng());
			double weight = ref.getWeight();
			awSum+=weight*angDiff;
			wSum+=weight;
		}
		return AngleUtility.angSum(awSum / wSum, angBase);
	}
	//
	public static double compRefColAngleMean(Collection<? extends AngObject> refs, AngObject angF, AngObject angT){
		if(refs.isEmpty()){
			return Double.MAX_VALUE;// not initilaised
		}
		if(refs.size() == 1){
			return refs.iterator().next().getAng();
		}
		if(angF!=null && angT!=null){
			if(angF.getAng() == angT.getAng()){
				return angF.getAng();
			}
		}
		double wSum = 0.0; // total weights
		for(AngObject ref:refs){
			wSum+=ref.getWeight();
		}
		double x = 0.0;
		double y = 0.0;
		for(AngObject ref:refs){
			double w = ref.getWeight();
			double ang = ref.getAng();
			x+=w*Math.cos(ang);
			y+=w*Math.sin(ang);
		}
		if(x==0.0 && y == 0.0){
			if(angF!=null && angT!=null){
				return compRefColAngleMeanArithmetic(refs, angF, angT);
			}else{
				return Double.NaN; // not defined
			}
		}else{
			return Math.atan2(y, x); // don't really need total weight
		}
	}
	public static double compRefColAngleMean(Collection<? extends AngObject> refs){
		return  compRefColAngleMean(refs, null, null);
	}
	//
	// circular variance from "Directional Statistics
	//
	public static double compRefColAngleVariation(Collection<? extends AngObject> refs, AngObject angF, AngObject angT){
		if(refs.isEmpty()){
			return 0.0;// not initilaised
		}
		if(refs.size() == 1){
			return 0.0; // no variance
		}
		double wSum = 0.0; // total weights
		for(AngObject ref:refs){
			wSum+=ref.getWeight();
		}
		double x = 0.0;
		double y = 0.0;
		for(AngObject ref:refs){
			double w = ref.getWeight();
			double ang = ref.getAng();
			x+=w*Math.cos(ang);
			y+=w*Math.sin(ang);
		}
		// 
		double r = Math.sqrt(x*x + y*y) / wSum;
		if(r>1.0){
			r = 1.0;
		}else if(r<0.0){
			r = 0.0;
		}
		return 1.0 - r;
	}
	public static double compRefColAngleVariation(Collection<? extends AngObject> refs){
		return compRefColAngleVariation(refs, null, null);
	}
/*
 * 
 *  Silhouette
 * 	
 */
	// this is a modified definition, taking the weights into consideration (implicitly)
	// for each point in a cluster, compute distance to the weighted cluster mean as ai, 
	// and distance to the weighted mean of the nearest cluster as bi, then compute si
	// the average si of all points is returned
	//
	private double computeSilhouetteAverage(){
		double sumSi = 0.0; // sum of silhouette 
		double numSamp = 0; // number of samples
		Collection<AngCluster<T>> clus = new Vector<AngCluster<T>>();
		clusters.values(clus);
		for(AngCluster<T> clu:clus){// for each cluster
			double sumSiClu = 0.0;
			Collection<AngRef<T>> refs = clu.getMembers();
			numSamp+=refs.size();
			for(AngRef ref:refs){
				double ai = ref.angDistance(clu);
				double bi = 10.0; // > PI
				for(AngCluster<T> clu2:clus){
					if(clu2!=clu){
						double angDiff = ref.angDistance(clu2);
						if(angDiff < bi){
							bi = angDiff;
						}
					}
				}
				double si = (bi - ai) / Math.max(bi, ai);
				sumSi+=si;
				sumSiClu+=si;
			}
//			System.out.println(refs.size()+", "+sumSiClu/refs.size());
		}
		double sumSiAve = sumSi/numSamp;
//		System.out.println("average silhouette:" + sumSiAve);
		return sumSiAve;
	}
	
	private static double computeAverageSilhouette(Collection<AngCluster> clus){
		double sumSi = 0.0; // sum of silhouette 
		double numSamp = 0; // number of samples
		if(clus.size() == 1){
			return 0.0; //?
		}
		for(AngCluster clu:clus){// for each cluster
			//System.out.println(clu.getAng());
			double sumSiClu = 0.0;
			Collection<AngRef> refs = clu.getMembers();
			numSamp+=refs.size();
			for(AngRef ref:refs){
				double ai = ref.angDistance(clu);
				double bi = Double.MAX_VALUE; // > PI
				double si;
				if(refs.size() == 1){
					si = 1.0; // arbitrary def (the original paper uses 0.0)
				}else if(ai==0.0){
					si = 1.0;
				}else{
					for(AngCluster clu2:clus){
						if(clu2 == clu.prev || clu2 == clu.next){
							double angDiff = ref.angDistance(clu2);
							boolean inCCW = AngleUtility.containsAng(ref.getAng(), clu.getAng(), clu2.getAng());
							if(inCCW){
								angDiff = AngleUtility.angFromTo(ref.getAng(), clu2.getAng());
							}else{
								angDiff = AngleUtility.angFromTo(clu2.getAng(), ref.getAng());
							}
							if(angDiff < bi){
								bi = angDiff;
							}
						}
					}
					si = (bi - ai) / Math.max(bi, ai);
				}
				sumSi+=si;
				sumSiClu+=si;
			}
			//System.out.println(refs.size()+", "+sumSiClu/refs.size());
		}
		double sumSiAve = sumSi/numSamp;
		//System.out.println("average silhouette:" + sumSiAve);
		return sumSiAve;
	}
	//
	private static double computeAverageSilhouetteWeighted(Collection<AngCluster> clus){
		double sumSi = 0.0; // sum of silhouette 
		double numSamp = 0; // number of samples
		if(clus.size() == 1){
			return 0.0; //?
		}
		for(AngCluster clu:clus){// for each cluster
			//System.out.println(clu.getAng());
			double sumSiClu = 0.0;
			Collection<AngRef> refs = clu.getMembers();
			double aveWeight = clu.getWeight() / refs.size();
			numSamp+=refs.size();
			for(AngRef ref:refs){
				double angDiff01 = ref.angDistance(clu);
				//double ai = angDiff01*aveWeight;
				// double ai = angDiff01*angDiff01/clu.getWeight()
				double bi = Double.MAX_VALUE; // > PI
				double si;
				if(refs.size() == 1){
					si = 1.0; // arbitrary def (the original paper uses 0.0)
				}else if(angDiff01 == 0){
					si = 1.0;
				}else{
					// should exclude the weight of current ref?
					//double ai = angDiff01/Math.sqrt(clu.getWeight());
					double ai = angDiff01*angDiff01/ clu.getWeight();
					//
					for(AngCluster clu2:clus){
						if(clu2!=clu){
							double aveWeight2 = clu2.getWeight() / clu2.getMembers().size();
							double angDiff = ref.angDistance(clu2);
							boolean inCCW = AngleUtility.containsAng(ref.getAng(), clu.getAng(), clu2.getAng());
							if(inCCW){
								angDiff = AngleUtility.angFromTo(ref.getAng(), clu2.getAng());
							}else{
								angDiff = AngleUtility.angFromTo(clu2.getAng(), ref.getAng());
							}
							//double biNew = angDiff*aveWeight2;
							//double biNew = angDiff / Math.sqrt(clu2.getWeight());
							double biNew = angDiff * angDiff / clu2.getWeight();
							if(biNew < bi){
								bi = biNew;
							}
						}
					}
					si = (bi - ai) / Math.max(bi, ai);
				}
				sumSi+=si;
				sumSiClu+=si;
			}
			//System.out.println(refs.size()+", "+sumSiClu/refs.size());
		}
		double sumSiAve = sumSi/numSamp;
		//System.out.println("average silhouette:" + sumSiAve);
		return sumSiAve;
	}
/*
 * 
 * debug info
 * 	
 */
	//
	public void report(int i){
		Collection<AngCluster<T>> clus = new Vector<AngCluster<T>>();
		clusters.values(clus);
		System.out.println("Iteration "+i);
		double wVar = 0.0;
		for(AngCluster<T> clu:clus){
			wVar+=clu.getVariance();
			System.out.println(clu.getId() + " numRef: " + clu.getMembers().size() + "- ang: " + Math.toDegrees(clu.getAng()) + " weight: " + clu.getWeight() + " projected length: " + clu.getLenSum() + " variance: " + clu.getVariance());
		}
		computeSilhouetteAverage();
		System.out.println("Totl variation: "+varianceSum + " within-cluster variation ratio: " + wVar/varianceSum);
	}
	
	public static void report(Collection<AngCluster> clus){
		double wVar = 0.0;
		System.out.println("*** \nNumber of Clusters: " + clus.size());
		for(AngCluster clu:clus){
			wVar+=clu.getVariance();
			System.out.println(clu.getId() + " numRef: " + clu.getMembers().size() + " - mean ang: " + Math.toDegrees(clu.getAng()) + " weight: " + clu.getWeight() + " projected length: " + clu.getLenSum() + " variance: " + clu.getVariance());
		}
		System.out.println("total within-cluster variance: " + wVar);
		double aveSih = computeAverageSilhouette(clus);
		System.out.println("Average Silhouette: " + aveSih);
		
	}
	//
	
	//
	public static void main(String[] args){
//		List<Feature> feats = JUMPUtility.loadShapeFile("d:/temp/ac01s.shp", false, "fid");
		
		double a90 = Math.PI*0.5;
		double a120 = Math.PI*2/3;
		double a180 = Math.PI;
		double a240 = -a120;
		double a270 = -a90;
		
		AngRef[] refs = new AngRef[12];
		
		refs[0] = new AngRef(0.11, 4.0);
		refs[4] = new AngRef(0.0, 3.0);
		refs[8] = new AngRef(-0.01, 2.0);
		refs[1] = new AngRef(a90, 3.0);
		refs[5] = new AngRef(a90, 2.0);
		refs[9] = new AngRef(a90, 1.0);
		refs[2] = new AngRef(a180, 2.0);
		refs[6] = new AngRef(a180, 1.0);
		refs[10] = new AngRef(a180, 3.0);
		refs[3] = new AngRef(a270, 3.0);
		refs[7] = new AngRef(a270, 2.0);
		refs[11] = new AngRef(a270, 1.0);
		
		refs[0] = new AngRef(0.0, 4.0);
		refs[4] = new AngRef(a120, 3.0);
		refs[8] = new AngRef(0.0, 2.0);
		refs[1] = new AngRef(a120, 3.0);
		refs[5] = new AngRef(a240, 3.0);
		refs[9] = new AngRef(a90, 1.0);
		refs[2] = new AngRef(a240, 3.0);
		refs[6] = new AngRef(a180, 1.0);
		refs[10] = new AngRef(a180, 3.0);
		refs[3] = new AngRef(0.0, 4.0);
		refs[7] = new AngRef(a270, 2.0);
		refs[11] = new AngRef(a270, 1.0);
		
		List<AngRef> angVec = new ArrayList<AngRef>(12);
		for(int i = 0; i < 12; ++i){
			angVec.add(refs[i]);
		}
		AngClusterGenerator gen = new AngClusterGenerator();
		gen.setTerminator(new TerminateEvaluatorVar(angVec, 1.0));
		gen.clusteringAnglesElbow(angVec);
		System.out.println("Done...");
/*
		for(int i = 1; i <= 12; ++i){
			AngClusterGenerator gen = new AngClusterGenerator();
			System.out.println("****** Cluster Number: "+i + "******");
			gen.setTerminator(new TerminateEvaluatorVar(angVec, 0.0));
			gen.clusteringAnglesHierarchical(angVec, i);
			gen.report(i);
			System.out.println("Done...");
		}
		//
		for(int i = 0; i < 12; ++i){
			AngClusterGenerator gen = new AngClusterGenerator();
			
			System.out.println("****** within-cluster variation ratio: "+ i*0.01 + "******");
			//gen.angleClusterKMean(angVec, 4);
			gen.setTerminator(new TerminateEvaluatorVar(angVec, i*0.01));
			gen.clusteringAnglesHierarchical(angVec, 0);
			gen.report(i);
			System.out.println("Done...");
		}
*/		
	}
}
