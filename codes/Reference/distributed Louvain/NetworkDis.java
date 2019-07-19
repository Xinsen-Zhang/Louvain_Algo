package com.graphsee.cypher.alao.distribute.louvain;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.cache.Cache.Entry;
import org.cltech.crayfish.Crayfish;
import org.cltech.crayfish.CrayfishCache;
import org.cltech.crayfish.CrayfishDataStreamer;

import com.graphsee.v5.core.CrayfishManager;

public class NetworkDis {

	protected int nNodes;
	// 备份的邻居
	private static Crayfish crayfish = CrayfishManager.Instance.getCrayfish();

	protected double totalEdgeWeightSelfLinks;
	protected String onlyEdge;
	protected String neighbor;
	// protected double[] nodeWeight;
	// 保存nodeweight
	protected String nodeWeight;
	protected String STRNODEWEIGHT = "_nodeWeight";

	private static CrayfishCache<Integer, List<LouvainData>> nNeighbors2;
	private static CrayfishCache<Integer, Integer> clusterTemp;

	public NetworkDis(int nodeSize, String nEIGHBOR2) {
		this.nNodes = nodeSize;
		System.out.println("NetworkDis node is " + nNodes);
		neighbor = nEIGHBOR2;
		totalEdgeWeightSelfLinks = 0;
		this.nodeWeight = getTotalEdgeWeightPerNodeV2();
	}

	private NetworkDis() {

	}

	/*
	 * public NetworkDis createReducedNetwork(ClusteringDis clustering) { NetworkDis
	 * reducedNetwork;
	 * 
	 * reducedNetwork = new NetworkDis();
	 * 
	 * reducedNetwork.nNodes = clustering.nClusters;
	 * 
	 * reducedNetwork.totalEdgeWeightSelfLinks = totalEdgeWeightSelfLinks;
	 * 
	 * String[] name = neighbor.split("_"); int current = Integer.parseInt(name[1]);
	 * 
	 * reducedNetwork.neighbor = name[0] + "_" + (current + 1);
	 * reducedNetwork.nodeWeight = reducedNetwork.neighbor + STRNODEWEIGHT;
	 * 
	 * 
	 * Collection<String> cacheString = crayfish.cacheNames(); //创建cache //
	 * nNeighborsTemp 临时的 if(cacheString.contains(reducedNetwork.nodeWeight)) {
	 * crayfish.cache(reducedNetwork.nodeWeight).clear(); }else {
	 * crayfish.createCache(reducedNetwork.nodeWeight); }
	 * 
	 * if(cacheString.contains(reducedNetwork.neighbor)) {
	 * crayfish.cache(reducedNetwork.neighbor).clear(); }else {
	 * crayfish.createCache(reducedNetwork.neighbor); }
	 * 
	 * //设置邻居 Collection<Double> totalEdgeWeights =
	 * crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(()
	 * -> { int clusterNum = clustering.nClusters, i, n, k = 0; String
	 * netWorkNeighbor = reducedNetwork.neighbor; String netWorknodeWeight =
	 * reducedNetwork.nodeWeight;
	 * 
	 * //创建新的cache CrayfishDataStreamer<Integer, List<LouvainData>> neighborStreamer
	 * = crayfish.dataStreamer(netWorkNeighbor); CrayfishDataStreamer<Integer,
	 * Double> nodeWeightStreamer = crayfish.dataStreamer(netWorknodeWeight);
	 * 
	 * int[] reducedNetworkNeighbor2 = new int[clusterNum];
	 * 
	 * double[] reducedNetworkEdgeWeight2 = new double[clusterNum];
	 * 
	 * 
	 * CrayfishCache<Integer, List<LouvainData>> nNeighbors2 =
	 * crayfish.cache(neighbor); CrayfishCache<Integer, Double> nodeWeight2 =
	 * crayfish.cache(nodeWeight); CrayfishCache<Integer, Integer> clusterTemp =
	 * crayfish.cache(clustering.strClustersName);
	 * 
	 * int length = clusterNum, elementId; double nodeWeightValue = 0, totalEdgeSum
	 * = 0, value;
	 * 
	 * for (i = 0; i < length; i++) { Set<Integer> neighbor = new HashSet<>();
	 * nodeWeightValue = 0; //每台机子只计算自己的 //Iterator<Entry<Integer,
	 * List<LouvainData>>> neighborIterator = nNeighbors2.localEntries().iterator();
	 * Iterator<Entry<Integer, List<LouvainData>>> neighborIterator =
	 * nNeighbors2.iterator(); while(neighborIterator.hasNext()) { Entry<Integer,
	 * List<LouvainData>> itemData = neighborIterator.next(); Integer key =
	 * itemData.getKey(); //if (clusterTemp.localPeek(key) == i) { if
	 * (clusterTemp.get(key) == i) { List<LouvainData> lists = itemData.getValue();
	 * //nodeWeightValue += nodeWeight2.localPeek(key); nodeWeightValue +=
	 * nodeWeight2.get(key); k = 0; for (LouvainData item : lists) { elementId =
	 * item.getElementId(); value = item.getValue(); //n =
	 * clusterTemp.localPeek(elementId); n = clusterTemp.get(elementId); if (n != i)
	 * { if (reducedNetworkEdgeWeight2[n] == 0) { reducedNetworkNeighbor2[k] = n;
	 * k++; neighbor.add(n); } reducedNetworkEdgeWeight2[n] += value; } else
	 * totalEdgeSum += value; } } }
	 * 
	 * List<LouvainData> lists = new ArrayList<>(); for (Integer id : neighbor) {
	 * lists.add(new LouvainData(id, reducedNetworkEdgeWeight2[id]));
	 * reducedNetworkEdgeWeight2[id] = 0; }
	 * 
	 * neighborStreamer.addData(i, lists); nodeWeightStreamer.addData(i,
	 * nodeWeightValue); }
	 * 
	 * neighborStreamer.flush(); neighborStreamer.close();
	 * 
	 * nodeWeightStreamer.flush(); nodeWeightStreamer.close();
	 * 
	 * return totalEdgeSum; });
	 * 
	 * 
	 * for(double item : totalEdgeWeights) { reducedNetwork.totalEdgeWeightSelfLinks
	 * = item; }
	 * 
	 * System.out.println("reducedNetwork.totalEdgeWeightSelfLinks is " +
	 * reducedNetwork.totalEdgeWeightSelfLinks);
	 * 
	 * return reducedNetwork; }
	 */

	public NetworkDis createReducedNetworkV2(ClusteringDis clustering) {
		System.out.println("--- create network ---");
		NetworkDis reducedNetwork;

		reducedNetwork = new NetworkDis();

		reducedNetwork.nNodes = clustering.nClusters;

		reducedNetwork.totalEdgeWeightSelfLinks = totalEdgeWeightSelfLinks;

		String[] name = neighbor.split("_");
		int current = Integer.parseInt(name[1]);

		reducedNetwork.neighbor = name[0] + "_" + (current + 1);
		reducedNetwork.nodeWeight = reducedNetwork.neighbor + STRNODEWEIGHT;

		Collection<String> cacheString = crayfish.cacheNames();
		// 创建cache
		// nNeighborsTemp 临时的
		if (cacheString.contains(reducedNetwork.nodeWeight)) {
			crayfish.cache(reducedNetwork.nodeWeight).clear();
		} else {
			crayfish.createCache(reducedNetwork.nodeWeight);
		}

		if (cacheString.contains(reducedNetwork.neighbor)) {
			crayfish.cache(reducedNetwork.neighbor).clear();
		} else {
			crayfish.createCache(reducedNetwork.neighbor);
		}

		//为了拼合的弄的
		CrayfishCache<Integer, List<String>> neighborStringCacheTotal;
		String neighborString = "neighborString";
		if (cacheString.contains(neighborString)) {
			neighborStringCacheTotal = crayfish.cache(neighborString);
			neighborStringCacheTotal.clear();
		} else {
			neighborStringCacheTotal = crayfish.createCache(neighborString);
		}

		// 设置邻居
		Collection<Double> totalEdgeWeights = crayfish.compute(CrayfishManager.Instance.getDBClusterGroup())
				.broadcast(() -> {
					int n, length = clustering.nClusters;

					CrayfishCache<Integer, Double> nodeWeightCache = crayfish.cache(reducedNetwork.nodeWeight);
					CrayfishCache<Integer, List<String>> neighborStringCache = crayfish.cache(neighborString);

					nNeighbors2 = crayfish.cache(neighbor);
					clusterTemp = crayfish.cache(clustering.strClustersName);

					CrayfishCache<Integer, Double> nodeWeight2 = crayfish.cache(nodeWeight);
					double totalEdgeSum = 0;

					// 每台机子只计算自己的
					Iterator<Entry<Integer, List<LouvainData>>> neighborIterator = nNeighbors2.localEntries().iterator();
					Map<Integer, Map<Integer, Double>> neighborValue = new HashMap<>();

					Map<Integer, Double> nodeValue = new HashMap<>(length, 1);

					int clusterId, elementId;
					double value;
					while (neighborIterator.hasNext()) {
						Entry<Integer, List<LouvainData>> itemData = neighborIterator.next();
						Integer key = itemData.getKey();
						clusterId = clusterTemp.localPeek(key);

						Map<Integer, Double> neighbors = neighborValue.get(clusterId);
						if (neighbors == null) {
							neighbors = new HashMap<>(length - 1, 1);
						}

						List<LouvainData> lists = itemData.getValue();

						if (nodeValue.containsKey(clusterId)) {
							nodeValue.put(clusterId, nodeValue.get(clusterId) + nodeWeight2.localPeek(key));
						} else {
							nodeValue.put(clusterId, nodeWeight2.localPeek(key));
						}

						for (LouvainData item : lists) {
							elementId = item.getElementId();
							value = item.getValue();
							n = clusterTemp.localPeek(elementId);
							if (n != clusterId) {
								if (neighbors.containsKey(n)) {
									neighbors.put(n, neighbors.get(n) + value);
								} else {
									neighbors.put(n, value);
								}
							} else
								totalEdgeSum += value;
						}
						neighborValue.put(clusterId, neighbors);
					}

					// 设置nodeweight
					for (int nid : nodeValue.keySet()) {
						double nodeweight = nodeValue.get(nid);
						nodeWeightCache.invoke(nid, (entery, args) -> {
							if (entery.getValue() == null) {
								entery.setValue(nodeweight);
							} else {
								entery.setValue(entery.getValue() + nodeweight);
							}
							return null;
						});
					}

					// 设置邻居

					for (int id : neighborValue.keySet()) {
						Map<Integer, Double> neighbors = neighborValue.get(id);
						StringBuilder sb = new StringBuilder();
						
						for (int key : neighbors.keySet()) {
							sb.append(key + ":" + neighbors.get(key) + ",");
						}

						String strRes = sb.toString();
						neighborStringCache.invoke(id, (entery, args) -> {
							List<String> res = entery.getValue();
							if (res == null) {
								res = new ArrayList<>();
								res.add(strRes);
								entery.setValue(res);
							} else {
								res.add(strRes);
								entery.setValue(res);
							}
							return null;
						});
					}
					
					return totalEdgeSum;
				});
		
		
		//拼合邻居
		crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			String netWorkNeighbor = reducedNetwork.neighbor;
			// 创建新的cache
			CrayfishDataStreamer<Integer, List<LouvainData>> neighborStreamer = crayfish.dataStreamer(netWorkNeighbor);

			Iterator<Entry<Integer, List<String>>> neighborStringIterator = neighborStringCacheTotal.localEntries()
					.iterator();
			while (neighborStringIterator.hasNext()) {
				Entry<Integer, List<String>> entry = neighborStringIterator.next();
				int i = entry.getKey(), key;
				double nodevalue;
				List<String> value = entry.getValue();
				Map<Integer, Double> maps = new HashMap<>();
				for (String str : value) {
					String[] item = str.split(",");
					
					if(item[0].equals(""))
						continue;
					
					for (String data : item) {
						String[] dataValue = data.split(":");
						key = Integer.parseInt(dataValue[0]);
						nodevalue = Double.parseDouble(dataValue[1]);
						Double temp = maps.get(key);
						if (temp == null) {
							maps.put(key, nodevalue);
						} else {
							maps.put(key, temp + nodevalue);
						}
					}
				}

				List<LouvainData> lists = new ArrayList<>();
				for (int elementId : maps.keySet()) {
					lists.add(new LouvainData(elementId, maps.get(elementId)));
				}
				neighborStreamer.addData(i, lists);
			}
			
			neighborStreamer.flush();
			neighborStreamer.close();
		});

		neighborStringCacheTotal.clear();
		
		for (double item : totalEdgeWeights) {
			reducedNetwork.totalEdgeWeightSelfLinks = item;
		}
		System.out.println("reducedNetwork.totalEdgeWeightSelfLinks " + reducedNetwork.totalEdgeWeightSelfLinks);

		return reducedNetwork;
	}

	public double[] getTotalEdgeWeightPerNode() {
		double[] totalEdgeWeightPerNode;

		totalEdgeWeightPerNode = new double[nNodes];
		String cache_nodeWeight = neighbor + "_nodeWeight";

		CrayfishCache<Integer, Double> nodeWeight;
		Collection<String> cacheString = crayfish.cacheNames();

		if (cacheString.contains(cache_nodeWeight)) {
			nodeWeight = crayfish.cache(cache_nodeWeight);
		} else {
			nodeWeight = crayfish.createCache(cache_nodeWeight);
		}

		crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			CrayfishCache<Integer, List<LouvainData>> edgeWeight = crayfish.cache(neighbor);
			Iterator<Entry<Integer, List<LouvainData>>> iterator = edgeWeight.localEntries().iterator();
			// 写数据
			CrayfishDataStreamer<Integer, Double> nodeWeightStreamer = crayfish.dataStreamer(cache_nodeWeight);

			while (iterator.hasNext()) {
				Entry<Integer, List<LouvainData>> entry = iterator.next();
				Integer key = entry.getKey();
				List<LouvainData> value = entry.getValue();
				double sum = 0;

				for (LouvainData item : value) {
					sum += item.getValue();
				}
				nodeWeightStreamer.addData(key, sum);
			}

			nodeWeightStreamer.flush();
			nodeWeightStreamer.close();
		});

		long startTime = System.currentTimeMillis();
		Iterator<Entry<Integer, Double>> iterator = nodeWeight.iterator();
		while (iterator.hasNext()) {
			Entry<Integer, Double> entry = iterator.next();
			totalEdgeWeightPerNode[entry.getKey()] = entry.getValue();
		}
		long endTime = System.currentTimeMillis();
		System.out.println("迭代一轮点集合所花时间 " + (endTime - startTime) + " ms");

		nodeWeight.clear();
		nodeWeight.close();
		nodeWeight.destroy();

		return totalEdgeWeightPerNode;
	}

	/**
	 * 可以移到初始化的部分
	 * 
	 * @return
	 */
	public String getTotalEdgeWeightPerNodeV2() {

		String cache_nodeWeight = neighbor + STRNODEWEIGHT;
		Collection<String> cacheString = crayfish.cacheNames();

		if (cacheString.contains(cache_nodeWeight)) {
			crayfish.cache(cache_nodeWeight).clear();
		} else {
			crayfish.createCache(cache_nodeWeight);
		}

		crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			CrayfishCache<Integer, List<LouvainData>> edgeWeight = crayfish.cache(neighbor);
			Iterator<Entry<Integer, List<LouvainData>>> iterator = edgeWeight.localEntries().iterator();
			// 写数据
			CrayfishDataStreamer<Integer, Double> nodeWeightStreamer = crayfish.dataStreamer(cache_nodeWeight);

			while (iterator.hasNext()) {
				Entry<Integer, List<LouvainData>> entry = iterator.next();
				Integer key = entry.getKey();
				List<LouvainData> value = entry.getValue();
				double sum = 0;

				for (LouvainData item : value) {
					sum += item.getValue();
				}
				nodeWeightStreamer.addData(key, sum);
			}
			nodeWeightStreamer.flush();
			nodeWeightStreamer.close();
		});

		// 不摧毁
		// nodeWeight.clear();
		// nodeWeight.close();
		// nodeWeight.destroy();

		return cache_nodeWeight;
	}

	public double getTotalEdgeWeightV2() {
		double res = 0;
		Collection<Double> sums = crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			CrayfishCache<Integer, Double> nodeWeightCache = crayfish.cache(this.nodeWeight);
			Iterator<Entry<Integer, Double>> iterator = nodeWeightCache.localEntries().iterator();
			// 写数据
			double sum = 0;

			while (iterator.hasNext()) {
				Entry<Integer, Double> entry = iterator.next();
				sum += entry.getValue();
			}
			return sum;
		});

		for (double item : sums) {
			res += item;
		}
		
		System.out.println("sum 计算得到 " + res);
		
		return res / 2;
	}

	/*
	 * 换成分布式的 public double getTotalEdgeWeight() { double sum = 0;
	 * System.out.println("---------------------------\n"); for (int i =
	 * nodeWeight.length - 1; i >= 0; i--) { sum += nodeWeight[i]; } return sum / 2;
	 * }
	 */
}