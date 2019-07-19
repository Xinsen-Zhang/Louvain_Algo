package com.graphsee.cypher.alao.distribute.louvain;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import javax.cache.Cache.Entry;

import org.cltech.crayfish.Crayfish;
import org.cltech.crayfish.CrayfishCache;
import org.cltech.crayfish.CrayfishDataStreamer;
import org.cltech.crayfish.cache.CacheMode;
import org.cltech.crayfish.cache.CacheWriteSynchronizationMode;
import org.cltech.crayfish.configuration.CacheConfiguration;

import com.graphsee.v5.core.CrayfishManager;

public class VOSClusteringTechniqueDisV2 {
	protected NetworkDis network;
	protected ClusteringDis clustering;
	protected double resolution;
	private static Crayfish crayfish = CrayfishManager.Instance.getCrayfish();
	private static int SPECIALNODE = 1000;

	public VOSClusteringTechniqueDisV2(NetworkDis network, double resolution) {
		this.network = network;
		clustering = new ClusteringDis(network.nNodes);
		clustering.initSingletonClustersV2();
		this.resolution = resolution;
	}

	public boolean runLouvainAlgorithm() {
		boolean update = true;
		VOSClusteringTechniqueDisV2 vOSClusteringTechnique;
		if (network.nNodes == 1)
			return false;
		
		while(update) {
			update = runLocalMovingAlgorithmV2();
			vOSClusteringTechnique = new VOSClusteringTechniqueDisV2(network.createReducedNetworkV2(clustering), resolution);
			clustering.mergeClustersV2(vOSClusteringTechnique.clustering);
		}
		
		System.out.println("==== is over ====");
		return update;
	}

	// 模块度计算
	public boolean runLocalMovingAlgorithmV2() {
		System.out.println("----> 开始移动  ----->");
		if (network.nNodes == 1)
			return false;

		Collection<String> cacheString = crayfish.cacheNames();
		// 设置簇的信息
		String clisterWeight = "clusterWeight"; // 存储，每一点的信息
		String clusterMachine = "clusterMachine";

		if (cacheString.contains(clisterWeight)) {
			crayfish.cache(clisterWeight).clear();
		} else {
			CacheConfiguration<Integer, Double> cacheConfig = new CacheConfiguration<>();
			cacheConfig.setName(clisterWeight);
			cacheConfig.setCacheMode(CacheMode.REPLICATED);
			cacheConfig.setBackups(0);
			cacheConfig.setWriteSynchronizationMode(CacheWriteSynchronizationMode.FULL_SYNC);
			crayfish.getOrCreateCache(cacheConfig);
		}

		if (cacheString.contains(clusterMachine)) {
			crayfish.cache(clusterMachine).clear();
		} else {
			crayfish.createCache(clusterMachine);
		}

		//同步
		boolean isResult = true, updateAll = false;
		//直到都不更新为止
		int iterator = 0;
		
		while(isResult) {
			// Parallel local clustering
			System.out.println("");
			Collection<Boolean> values = crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
				
				double maxQualityFunction, qualityFunction, reduceWeight, totWeight, weight;
				double[] edgeWeightPerCluster;
				boolean[] nodeIsLocal;
				boolean update = false, oldBestIsLocal = false;
				int bestCluster, id, k, l, nNeighboringClusters, t = 0, nodes = network.nNodes;
				//上一轮点的权重
				String nodeWeight = network.nodeWeight;
				String cahceNeighbor = network.neighbor;
				String cacheClusterName = clustering.strClustersName;
				
				
				CrayfishCache<Integer, Double> nodeWeightCache = crayfish.cache(nodeWeight);
				// 原来的数据    弄好最后更新该cache
				CrayfishCache<Integer, Integer> clusterCache = crayfish.cache(cacheClusterName);
				CrayfishDataStreamer<Integer, Integer> clusterDataStreamer = crayfish.dataStreamer(cacheClusterName);
				clusterDataStreamer.allowOverwrite(true);
				
				int key, clusterId;
				Double oldValue;
				
				// 申请一轮空间 要移动簇的到邻居的空间
				int[] neighboringCluster = new int[nodes - 1];
				
				//点的权重值  这个不改
				CrayfishCache<Integer, List<LouvainData>> louvainValue = crayfish.cache(cahceNeighbor);
				
				//每个簇的权重存储在本地上   可以通过cluster拿取
				Map<Integer, Double> clusterWeightAll = new HashMap<>();
				//每个点的簇号
				Map<Integer, Integer> nodeCluster = new HashMap<>();
				
				edgeWeightPerCluster = new double[nodes];
				nodeIsLocal = new boolean[nodes];
				
				//本机的点
				Set<Integer> localIds = new HashSet<>();
				Iterator<Entry<Integer, Double>> iteratorNodeWeightIterator = nodeWeightCache.localEntries().iterator();
				while (iteratorNodeWeightIterator.hasNext()) {
					Entry<Integer, Double> item = iteratorNodeWeightIterator.next();
					localIds.add(item.getKey());
				}
				
				//更新上轮的点info
				//每轮重新初始化TOT cluster值
				iteratorNodeWeightIterator = nodeWeightCache.iterator();
				while (iteratorNodeWeightIterator.hasNext()) {
					Entry<Integer, Double> item = iteratorNodeWeightIterator.next();
					key = item.getKey();
					clusterId = clusterCache.localPeek(key);
					
					oldValue = clusterWeightAll.get(clusterId);
					if(oldValue == null) {
						clusterWeightAll.put(clusterId, item.getValue());
					}else {
						clusterWeightAll.put(clusterId, oldValue + item.getValue());
					}
					nodeCluster.put(key, clusterId);
				}
				
				// 迭代本地点进行计算法
				Iterator<Entry<Integer, List<LouvainData>>> localNeighbor = louvainValue.localEntries().iterator();
				while (localNeighbor.hasNext()) {
					Entry<Integer, List<LouvainData>> value = localNeighbor.next();
					id = value.getKey();

					// 这个值不会变化
					weight = nodeWeightCache.localPeek(id);
					clusterId = nodeCluster.get(id);

					List<LouvainData> lists = value.getValue();
					// 耗时在去其他簇在邻居

					nNeighboringClusters = 0;
					
					//计算 w u->c 的权重值
					for (LouvainData data : lists) {
						int elementId = data.getElementId();
						l = nodeCluster.get(elementId);
						if (edgeWeightPerCluster[l] == 0) {
							neighboringCluster[nNeighboringClusters] = l;
							
							if(localIds.contains(elementId)) {
								nodeIsLocal[nNeighboringClusters] = true;
							}else {
								nodeIsLocal[nNeighboringClusters] = false;
							}
							nNeighboringClusters++;
						}
						edgeWeightPerCluster[l] += data.getValue();
					}

					reduceWeight = clusterWeightAll.get(clusterId) - weight;

					bestCluster = Integer.MAX_VALUE;
					maxQualityFunction = 0;
					for (k = 0; k < nNeighboringClusters; k++) {
						l = neighboringCluster[k];
						// 邻居和自己的权重相同
						if (l == clusterId) {
							totWeight = reduceWeight;
						} else {
							totWeight = clusterWeightAll.get(l);
						}

						qualityFunction = edgeWeightPerCluster[l] - weight * totWeight * resolution;
						
						//**** 优选 小标签原则 优化本地的 优选个数大于1的
						if (qualityFunction >= maxQualityFunction) {
							if(qualityFunction == maxQualityFunction) {
								//本地的优新更新  优化本地的
								if(nodeIsLocal[l]) {
									
									//bestCluster 原值不是本地直接更新
									if(!oldBestIsLocal) {
										bestCluster = l;
									}else {
										if(l < bestCluster) {
											bestCluster = l;
										}
									}
									oldBestIsLocal = true;
								}else {
									//非本地的更新 小标签原则
									if(!oldBestIsLocal && (l < bestCluster)) {
										bestCluster = l;
										oldBestIsLocal = false;
									}
								}
								
							}else {
								bestCluster = l;
								
								if(nodeIsLocal[l]) {
									oldBestIsLocal = true;
								}else {
									oldBestIsLocal = false;
								}
							}
							maxQualityFunction = qualityFunction;
						}
						edgeWeightPerCluster[l] = 0;
					}
					
					//通过簇号判断是否更新
					if (bestCluster == Integer.MAX_VALUE)
						bestCluster = clusterId;

					//更新
					if (bestCluster != clusterId) {
						clusterWeightAll.put(clusterId, reduceWeight);
						clusterWeightAll.put(bestCluster, clusterWeightAll.get(bestCluster) + weight);
						nodeCluster.put(id, bestCluster);
						// 设置为1
						update = true;
					}

					if (t++ % 3000 == 0) {
						System.out.print(".");
					}
				}
				//Synchronize Vhigh states   这一步暂时不采用
				
				//Synchronize ghost vertex community states  因为这个put慢的原因吗
				iteratorNodeWeightIterator = nodeWeightCache.localEntries().iterator();
				int temp = 0;
				while (iteratorNodeWeightIterator.hasNext()) {
					Entry<Integer, Double> item = iteratorNodeWeightIterator.next();
					key = item.getKey();
					clusterDataStreamer.addData(key, nodeCluster.get(key));
					//更新到本地修改到cache中
					temp++;
				}
				System.out.println("update is " + temp);
				//Calculate partial modularity 可以先不用				
				return update;
			});
			
			isResult = false;
			//同步了
			for (Boolean id : values) {
				if (id != false) {
					isResult = true;
					updateAll = true;
					break;
				}
			}
			System.out.println("k is " + (iterator++));
		}

		return updateAll;
	}

	
	// 模块度计算
	public boolean runLocalMovingAlgorithmV2NoUpdate() {
		System.out.println("----> 开始移动  ----->");
		if (network.nNodes == 1)
			return false;

		Collection<String> cacheString = crayfish.cacheNames();
		// 设置簇的信息
		String clisterWeight = "clusterWeight"; // 存储，每一点的信息
		String clusterMachine = "clusterMachine";

		if (cacheString.contains(clisterWeight)) {
			crayfish.cache(clisterWeight).clear();
		} else {
			CacheConfiguration<Integer, Double> cacheConfig = new CacheConfiguration<>();
			cacheConfig.setName(clisterWeight);
			cacheConfig.setCacheMode(CacheMode.REPLICATED);
			cacheConfig.setBackups(0);
			cacheConfig.setWriteSynchronizationMode(CacheWriteSynchronizationMode.FULL_SYNC);
			crayfish.getOrCreateCache(cacheConfig);
		}

		if (cacheString.contains(clusterMachine)) {
			crayfish.cache(clusterMachine).clear();
		} else {
			crayfish.createCache(clusterMachine);
		}

		//同步
		boolean isResult = true, updateAll = false;
		//直到都不更新为止
		while(isResult) {
			
			// Parallel local clustering with delegates
			Collection<Boolean> values = crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
				
				double maxQualityFunction, qualityFunction, reduceWeight, totWeight, weight;
				double[] edgeWeightPerCluster;
				boolean update = false;
				int bestCluster, id, k, l, nNeighboringClusters, t = 0, nodes = network.nNodes;
				//上一轮点的权重
				String nodeWeight = network.nodeWeight;
				String cahceNodeWeightName = network.nodeWeight;
				String strNeighbor = network.neighbor;
				String cacheClusterName = clustering.strClustersName;
				
				
				CrayfishCache<Integer, Double> nodeWeightCache = crayfish.cache(cahceNodeWeightName);
				// 原来的数据    弄好最后更新该cache
				CrayfishCache<Integer, Integer> clusterCache = crayfish.cache(cacheClusterName);
				
				int key, clusterId;
				Double oldValue;
				
				// 申请一轮空间 要移动簇的到邻居的空间
				int[] neighboringCluster = new int[nodes - 1];
				
				//点的权重值  这个不改
				CrayfishCache<Integer, Double> nodeWeightTemp = crayfish.cache(nodeWeight);
				CrayfishCache<Integer, List<LouvainData>> louvainValue = crayfish.cache(strNeighbor);
				
				//每个簇的权重存储在本地上   可以通过cluster拿取
				Map<Integer, Double> clusterWeightAll = new HashMap<>();
				//每个点的簇号
				Map<Integer, Integer> nodeCluster = new HashMap<>();
				edgeWeightPerCluster = new double[nodes];
				//更新上轮的点info
				//每轮重新初始化TOT cluster值
				Iterator<Entry<Integer, Double>> iteratorNodeWeightIterator = nodeWeightCache.localEntries().iterator();
				while (iteratorNodeWeightIterator.hasNext()) {
					Entry<Integer, Double> item = iteratorNodeWeightIterator.next();
					key = item.getKey();
					clusterId = clusterCache.localPeek(key);
					
					oldValue = clusterWeightAll.get(clusterId);
					if(oldValue == null) {
						clusterWeightAll.put(clusterId, item.getValue());
					}else {
						clusterWeightAll.put(clusterId, oldValue + item.getValue());
					}
					nodeCluster.put(key, clusterId);
				}
				
				

				//计算 w u->c 的权重值
				Map<Integer, Map<Integer,Double>> inC = new HashMap<>();
				Iterator<Entry<Integer, List<LouvainData>>> localNeighbor = louvainValue.localEntries().iterator();
				while (localNeighbor.hasNext()) {
					Entry<Integer, List<LouvainData>> value = localNeighbor.next();
					List<LouvainData> lists = value.getValue();
					// 耗时在去其他簇在邻居
					nNeighboringClusters = 0;
					//计算 w u->c 的权重值
					for (LouvainData data : lists) {
						int elementId = data.getElementId();
						l = nodeCluster.get(elementId);
						if (edgeWeightPerCluster[l] == 0) {
							neighboringCluster[nNeighboringClusters] = l;
							nNeighboringClusters++;
						}
						edgeWeightPerCluster[l] += data.getValue();
					}
					
					Map<Integer,Double> in = new HashMap<>(nNeighboringClusters, 1);
					for (k = 0; k < nNeighboringClusters; k++) {
						in.put(neighboringCluster[k], edgeWeightPerCluster[k]);
					}
					
					inC.put(value.getKey(), in);
				}
				
				// 迭代本地点进行计算法
				localNeighbor = louvainValue.localEntries().iterator();
				while (localNeighbor.hasNext()) {
					Entry<Integer, List<LouvainData>> value = localNeighbor.next();
					id = value.getKey();

					// 这个值不会变化
					weight = nodeWeightTemp.localPeek(id);
					clusterId = nodeCluster.get(id);

					//List<LouvainData> lists = value.getValue();
					
					reduceWeight = clusterWeightAll.get(clusterId) - weight;
					
					bestCluster = Integer.MAX_VALUE;
					maxQualityFunction = 0;
					
					Map<Integer, Double> in = inC.get(id);
					//计算法Q的增益，  计算数据根据原始的in
					for(int clusterIdin : in.keySet()) {
						if(clusterIdin == clusterId) {
							totWeight = reduceWeight; 
						}else {
							totWeight = clusterWeightAll.get(clusterIdin);
						}
						
						qualityFunction = edgeWeightPerCluster[clusterIdin] - weight * totWeight * resolution;
						
						//**** 优选 小标签原则 优化本地的 优选个数大于1的
						if (qualityFunction >= maxQualityFunction) {
							if(qualityFunction == maxQualityFunction) {
								//本地的优新更新  优化本地的
								if(isLoalCluster(clusterIdin)) {
									//bestCluster 不是本地直接更新
									if(!isLoalCluster(bestCluster)) {
										bestCluster = clusterIdin;
									}
									
									if(clusterIdin < bestCluster) {
										bestCluster = clusterIdin;
									}
								}else {
									//非本地的更新 小标签原则
									if(clusterIdin < bestCluster) {
										bestCluster = clusterIdin;
									}
								}
							}else {
								bestCluster = clusterIdin;
							}
							maxQualityFunction = qualityFunction;
						}
					}

					
					//通过簇号判断是否更新
					if (bestCluster == Integer.MAX_VALUE)
						bestCluster = clusterId;

					//更新
					if (bestCluster != clusterId) {
						//clusterWeightAll.put(clusterId, reduceWeight);
						//clusterWeightAll.put(bestCluster, clusterWeightAll.get(bestCluster) + weight);
						nodeCluster.put(id, bestCluster);
						// 设置为1
						update = true;
					}

					if (t++ % 3000 == 0) {
						System.out.print(".");
					}
				}
				//Synchronize Vhigh states   这一步暂时不采用
				
				System.out.println("所有点遍历完成");
				
				//Synchronize ghost vertex community states
				iteratorNodeWeightIterator = nodeWeightCache.localEntries().iterator();
				while (iteratorNodeWeightIterator.hasNext()) {
					Entry<Integer, Double> item = iteratorNodeWeightIterator.next();
					key = item.getKey();
					
					if(nodeCluster.get(key) == null) {
						System.out.println("===  error === " + key);
					}
					clusterCache.put(key, nodeCluster.get(key));
					//更新到本地修改到cache中
				}
				//Calculate partial modularity 可以先不用				
				return update;
			});
			
			isResult = false;
			//同步了
			for (Boolean id : values) {
				if (id != false) {
					isResult = true;
					updateAll = true;
					break;
				}
			}
		}

		return updateAll;
	}

	
	public boolean isLoalCluster(int l) {
		boolean res = true;
		return res;
	}

	// 模块度计算
	public double calcQualityFunctionV2() {
		double qualityFunction = 0;

		String strNeighbor = network.neighbor;

		// 每台机子计算其值
		Collection<Double> qualityFunctionSeg = crayfish.compute(CrayfishManager.Instance.getDBClusterGroup())
				.broadcast(() -> {
					CrayfishCache<Integer, Integer> clusterTemp = crayfish.cache(clustering.strClustersName);

					CrayfishCache<Integer, List<LouvainData>> louvainValue = crayfish.cache(strNeighbor);
					Iterator<Entry<Integer, List<LouvainData>>> localNeighbor = louvainValue.localEntries().iterator();
					int id, key;

					double qualityFunctionSum = 0;
					while (localNeighbor.hasNext()) {
						Entry<Integer, List<LouvainData>> entry = localNeighbor.next();
						List<LouvainData> lists = entry.getValue();
						key = entry.getKey();
						id = clusterTemp.localPeek(key);
						for (LouvainData data : lists) {
							if (clusterTemp.localPeek(data.getElementId()) == id)
								qualityFunctionSum += data.getValue();
						}
					}

					return qualityFunctionSum;
				});

		for (double value : qualityFunctionSeg) {
			qualityFunction += value;
		}

		qualityFunction += network.totalEdgeWeightSelfLinks;

		String clisterWeight = "clusterWeight";
		crayfish.cache(clisterWeight).clear();
		// 可以优化 从上一轮获取
		crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			// 存储，每一点的信息
			CrayfishCache<Integer, Double> clusterWeightCache = crayfish.cache(clisterWeight);

			CrayfishCache<Integer, Double> clusterWeightCacheTemp = crayfish.cache(network.nodeWeight);
			CrayfishCache<Integer, Integer> clusteTemp = crayfish.cache(clustering.strClustersName);
			Iterator<Entry<Integer, Double>> clusterIterator = clusterWeightCacheTemp.localEntries().iterator();
			int key, clusterId;
			while (clusterIterator.hasNext()) {
				Entry<Integer, Double> entry = clusterIterator.next();
				key = entry.getKey();

				clusterId = clusteTemp.localPeek(key);

				clusterWeightCache.invoke(clusterId, (entrys, args) -> {
					double temp = (double) args[0];
					if (entrys.getValue() == null) {
						entrys.setValue(temp);
					} else {
						entrys.setValue(entrys.getValue() + temp);
					}
					return null;
				}, entry.getValue());
			}

		});

		Collection<Double> qualityFunctionSums = crayfish.compute(CrayfishManager.Instance.getDBClusterGroup())
				.broadcast(() -> {
					CrayfishCache<Integer, Double> tempCache = crayfish.cache(clisterWeight);
					Iterator<Entry<Integer, Double>> iterator = tempCache.localEntries().iterator();
					double sum = 0, clusterWeightValue = 0;
					while (iterator.hasNext()) {
						Entry<Integer, Double> entry = iterator.next();
						clusterWeightValue = entry.getValue();
						sum += clusterWeightValue * clusterWeightValue * resolution;
					}

					return sum;
				});

		for (Double value : qualityFunctionSums) {
			qualityFunction -= value;
		}

		qualityFunction /= 2 * network.getTotalEdgeWeightV2() + network.totalEdgeWeightSelfLinks;

		System.out.println("current qualityFunction is " + qualityFunction);
		return qualityFunction;
	}

	// 模块度计算
	public boolean runLocalMovingAlgorithmEdge(Random random) {
		System.out.println("----> 开始移动  ----->");
		int i;
		if (network.nNodes == 1)
			return false;

		Collection<String> cacheString = crayfish.cacheNames();
		// 设置簇的信息
		String clisterWeight = "clusterWeight"; // 存储，每一点的信息
		String nStableNodesValue = "nStableNodes";
		String clusterMachine = "clusterMachine";

		if (cacheString.contains(clisterWeight)) {
			crayfish.cache(clisterWeight).clear();
		} else {
			CacheConfiguration<Integer, Double> cacheConfig = new CacheConfiguration<>();
			cacheConfig.setName(clisterWeight);
			cacheConfig.setCacheMode(CacheMode.REPLICATED);
			cacheConfig.setBackups(0);
			cacheConfig.setWriteSynchronizationMode(CacheWriteSynchronizationMode.FULL_SYNC);
			crayfish.getOrCreateCache(cacheConfig);
		}

		if (cacheString.contains(nStableNodesValue)) {
			crayfish.cache(nStableNodesValue).clear();
		} else {
			crayfish.createCache(nStableNodesValue);
		}

		if (cacheString.contains(clusterMachine)) {
			crayfish.cache(clusterMachine).clear();
		} else {
			crayfish.createCache(clusterMachine);
		}

		Collection<Boolean> size = crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			String cahceNodeWeightName = network.nodeWeight;
			String cacheClusterName = clustering.strClustersName;

			// 原来的数据
			CrayfishCache<Integer, Double> nodeWeightCache = crayfish.cache(cahceNodeWeightName);
			CrayfishCache<Integer, Integer> clusterCache = crayfish.cache(cacheClusterName);

			// 设置簇的信息
			CrayfishCache<Integer, Double> nodeWeightCacheNew = crayfish.cache(clisterWeight);

			int key, clusterId;
			Iterator<Entry<Integer, Double>> iteratorNodeWeightIterator = nodeWeightCache.localEntries().iterator();
			while (iteratorNodeWeightIterator.hasNext()) {
				Entry<Integer, Double> item = iteratorNodeWeightIterator.next();
				key = item.getKey();

				clusterId = clusterCache.localPeek(key);
				nodeWeightCacheNew.invoke(clusterId, (entrys, args) -> {
					double temp = (double) args[0];
					if (entrys.getValue() == null) {
						entrys.setValue(temp);
					} else {
						entrys.setValue(entrys.getValue() + temp);
					}
					return null;
				}, item.getValue());

			}

			// 设置 nStableNodes值
			CrayfishCache<Integer, Integer> nStableNodesValueCache = crayfish.cache(nStableNodesValue);
			CrayfishCache<Integer, Integer> clusterMachineCache = crayfish.cache(clusterMachine);
			int affinity = CrayfishManager.Instance.getAffinity();
			nStableNodesValueCache.put(affinity, 0);
			clusterMachineCache.put(affinity, 1);
			return true;
		});

		int machineSize = size.size();

		CrayfishCache<Integer, Integer> clusterMachineCache1 = crayfish.cache(clusterMachine);
		clusterMachineCache1.put(SPECIALNODE, machineSize);

		// 开始迭代计算
		Collection<Boolean> values = crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			double maxQualityFunction, qualityFunction, weightCacheTemp, calcWeightCacheTemp;
			double[] edgeWeightPerCluster;
			boolean update = false;
			int bestCluster, j, k, l, nNeighboringClusters, t = 0, affinity = CrayfishManager.Instance.getAffinity(),
					nodes = network.nNodes, machineSizeTemp = machineSize;

			int stableNodes = 1;
			String nodeWeight = network.nodeWeight;

			edgeWeightPerCluster = new double[nodes];
			int[] neighboringCluster = new int[nodes - 1];

			// clustering.cluster 4个cache
			CrayfishCache<Integer, Double> nodeWeightTemp = crayfish.cache(nodeWeight);
			CrayfishCache<Integer, Integer> clusterTemp = crayfish.cache(clustering.strClustersName);
			CrayfishCache<Integer, Double> clusterWeightCacheTemp = crayfish.cache(clisterWeight);
			CrayfishCache<Integer, Integer> nStableNodesValueCacheTemp = crayfish.cache(nStableNodesValue);
			CrayfishCache<Integer, Integer> clusterMachineCache = crayfish.cache(clusterMachine);

			String strNeighbor = network.neighbor;
			CrayfishCache<Integer, List<LouvainData>> louvainValue = crayfish.cache(strNeighbor);
			int localSize = louvainValue.localSize();

			boolean isContinue = true;

			System.out.printf("-------- 计算开始迭代 %d ---------\n", network.nNodes);
			while (isContinue) {
				// 同步操作
				int nodeSize = 0;

				do {
					nodeSize = 0;
					Iterator<Entry<Integer, Integer>> segIterator = clusterMachineCache.iterator();
					while (segIterator.hasNext()) {
						Entry<Integer, Integer> item = segIterator.next();
						if (item.getKey() != SPECIALNODE) {
							nodeSize += item.getValue();
						}
					}
					Thread.sleep(1000);
				} while (nodeSize != machineSizeTemp);

				// 三台机都通过
				nodeSize = clusterMachineCache.get(SPECIALNODE);

				if (nodeSize != machineSizeTemp) {
					clusterMachineCache.invoke(SPECIALNODE, (entery, args) -> {
						entery.setValue(machineSizeTemp);
						return null;
					});
				}

				while (clusterMachineCache.get(SPECIALNODE) != machineSizeTemp) {
					Thread.sleep(1000);
				}
				// 拿到每台机子的stableValue值
				Iterator<Entry<Integer, Integer>> stableIteror = nStableNodesValueCacheTemp.iterator();

				nodeSize = 0;
				while (stableIteror.hasNext()) {
					Entry<Integer, Integer> item = stableIteror.next();
					nodeSize += item.getValue();
				}

				if (nodeSize >= nodes) {
					isContinue = false;
				}

				clusterMachineCache.put(affinity, 0);
				// 更新本机的
				do {
					nodeSize = 0;
					Iterator<Entry<Integer, Integer>> segIterator = clusterMachineCache.iterator();
					while (segIterator.hasNext()) {
						Entry<Integer, Integer> item = segIterator.next();
						if (item.getKey() != SPECIALNODE) {
							nodeSize += item.getValue();
						}
					}
					Thread.sleep(1000);
				} while (nodeSize != 0);

				nodeSize = clusterMachineCache.get(SPECIALNODE);
				if (nodeSize != 0) {
					clusterMachineCache.invoke(SPECIALNODE, (entery, args) -> {
						entery.setValue(0);
						return null;
					});
				}

				Iterator<Entry<Integer, List<LouvainData>>> localNeighbor = louvainValue.localEntries().iterator();
				while (localNeighbor.hasNext() && isContinue) {
					Entry<Integer, List<LouvainData>> value = localNeighbor.next();
					j = value.getKey();

					// 这个值不会变化
					double weight = nodeWeightTemp.localPeek(j);
					int clusterId = clusterTemp.localPeek(j);

					List<LouvainData> lists = value.getValue();
					// 耗时在去其他簇在邻居

					nNeighboringClusters = 0;
					for (LouvainData data : lists) {
						int elementId = data.getElementId();
						l = clusterTemp.localPeek(elementId);
						if (edgeWeightPerCluster[l] == 0) {
							neighboringCluster[nNeighboringClusters] = l;
							nNeighboringClusters++;
						}
						edgeWeightPerCluster[l] += data.getValue();
					}

					weightCacheTemp = clusterWeightCacheTemp.localPeek(clusterId) - weight;

					bestCluster = -1;
					maxQualityFunction = 0;
					for (k = 0; k < nNeighboringClusters; k++) {
						l = neighboringCluster[k];
						// 邻居和自己的权重相同
						if (l == clusterId) {
							calcWeightCacheTemp = weightCacheTemp;
						} else {
							calcWeightCacheTemp = clusterWeightCacheTemp.localPeek(l);
						}

						qualityFunction = edgeWeightPerCluster[l] - weight * calcWeightCacheTemp * resolution;
						if (qualityFunction > maxQualityFunction) {
							bestCluster = l;
							maxQualityFunction = qualityFunction;
						}
						edgeWeightPerCluster[l] = 0;
					}

					if (bestCluster == -1)
						bestCluster = clusterId;

					if (bestCluster != clusterId) {

						clusterWeightCacheTemp.invoke(clusterId, (entery, args) -> {
							entery.setValue(entery.getValue() - weight);
							return null;
						});

						clusterWeightCacheTemp.invoke(bestCluster, (entery, args) -> {
							entery.setValue(entery.getValue() + weight);
							return null;
						});
						// clusterWeightCacheTemp.put(clusterId, clusterWeightCacheTemp.get(clusterId) -
						// weight);

						// clusterWeightCacheTemp.put(bestCluster,
						// clusterWeightCacheTemp.get(bestCluster) + weight);

						clusterTemp.put(j, bestCluster);
						// 设置为1
						stableNodes = 1;
						update = true;

					} else {
						if ((++stableNodes) >= localSize) {
							// 结束整个程序
							System.out.println("\ncurrent t is " + t + " update is " + update + " local size is "
									+ localSize + " \n ");
							break;
						}
					}

					if (t++ % 3000 == 0) {
						System.out.print(".");
					}

					// t的数量 本机迭代结束
					if (t % localSize == 0) {
						break;
					}
				}

				nStableNodesValueCacheTemp.invoke(affinity, (entery, args) -> {
					entery.setValue((int) args[0]);
					return null;
				}, stableNodes);

				// 修改Machine为通过
				clusterMachineCache.invoke(affinity, (entery, args) -> {
					entery.setValue(1);
					return null;
				});
				// 遍历完一轮 需要
			}

			return update;
		});

		boolean isResult = false;
		for (Boolean id : values) {
			if (id != false) {
				isResult = true;
				break;
			}
		}

		int[] newCluster = new int[network.nNodes];

		CrayfishDataStreamer<Integer, Integer> clusterStreamer = crayfish.dataStreamer(clustering.strClustersName);
		CrayfishCache<Integer, Integer> clusterTemp = crayfish.cache(clustering.strClustersName);

		int[] newClusterTemp = new int[network.nNodes];
		int[] perNode = new int[network.nNodes];
		// 有cluster的
		for (i = 0; i < network.nNodes; i++) {
			newClusterTemp[i] = clusterTemp.localPeek(i);
			perNode[newClusterTemp[i]] += 1;
		}
		// 清空数据
		clusterTemp.clear();

		clustering.nClusters = 0;

		for (i = 0; i < network.nNodes; i++)
			if (perNode[i] > 0) {
				newCluster[i] = clustering.nClusters;
				clustering.nClusters++;
			}

		perNode = null;
		int clusterId;
		for (i = 0; i < network.nNodes; i++) {
			clusterId = newCluster[newClusterTemp[i]];
			clusterStreamer.addData(i, clusterId);
		}

		clusterStreamer.flush();
		clusterStreamer.close();

		return isResult;
	}

	public ClusteringDis getClustering() {
		return clustering;
	}

}
