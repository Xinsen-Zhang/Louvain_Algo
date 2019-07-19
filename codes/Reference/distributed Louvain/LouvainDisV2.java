package com.graphsee.cypher.alao.distribute.louvain;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import javax.cache.Cache.Entry;

import org.cltech.crayfish.Crayfish;
import org.cltech.crayfish.CrayfishCache;
import org.cltech.crayfish.CrayfishDataStreamer;
import org.cltech.crayfish.binary.BinaryObject;
import org.cltech.crayfish.cluster.ClusterNode;

import com.graphsee.cypher.algo.helper.Helper;
import com.graphsee.v5.algo.utils.GraphHelper;
import com.graphsee.v5.algo.utils.Pools;
import com.graphsee.v5.core.CrayfishManager;
import com.graphsee.v5.core.GraphManager;
import com.graphsee.v5.core.element.GsEdgeDirection;
import com.graphsee.v5.core.entry.NeighborValue;
import com.graphsee.v5.core.exception.GsGraphException;
import com.graphsee.v5.core.graph.GsCompleteGraph;
import com.graphsee.v5.query.NeighborIterator;

public class LouvainDisV2 {
	private static GsCompleteGraph graph;
	private final String graphName;
	private static int ITERATIONS;
	private static String VERTEX;
	private static String NEIGHBOR = "nNeighbors";
	private static int NODESIZE;
	private static CrayfishDataStreamer<Integer, Integer> nodeSize_streamer;
	private static Crayfish crayfish = CrayfishManager.Instance.getCrayfish();
	private static AtomicInteger EDGENUM = new AtomicInteger(0);

	public LouvainDisV2(GsCompleteGraph graph, int nIterations, String vertex, String relation) {
		// TODO Auto-generated constructor stub
		this.graphName = graph.getName();
		ITERATIONS = nIterations;
		VERTEX = vertex;
	}

	public Map<String, String> compute(boolean isPrim, boolean isWrite) throws GsGraphException {
		boolean update;
		int j;
		double resolution2, modularity, maxModularity;

		NetworkDis network;
		VOSClusteringTechniqueDisV2 voSClusteringTechnique;
		ClusteringDis clusteringDis = null;

		// 创建id 映射
		network = init2();

		System.out.println("---初始化完成---");
		resolution2 = 1.0 / (2 * network.getTotalEdgeWeightV2() + network.totalEdgeWeightSelfLinks);
		
		System.out.println("----- 计算法EdgeWeight ------");
		voSClusteringTechnique = new VOSClusteringTechniqueDisV2(network, resolution2);

		j = 0;
		update = true;
		maxModularity = Double.NEGATIVE_INFINITY;
		int clustering = Integer.MAX_VALUE, tNodes = 0;

		do {
			
			long startTime = System.currentTimeMillis();
			update = voSClusteringTechnique.runLouvainAlgorithm();
			
			modularity = voSClusteringTechnique.calcQualityFunctionV2();
			
			long endTime = System.currentTimeMillis();
			System.out.println("测试得到 所话时间 " + (endTime - startTime));

			Collection<String> cacheString = crayfish.cacheNames();
			for (String item : cacheString) {
				if (item.contains(NEIGHBOR + "_") && (!item.contains(NEIGHBOR + "_0"))) {
					crayfish.cache(item).close();
					crayfish.cache(item).destroy();
				}
			}
			
			tNodes = voSClusteringTechnique.clustering.nClusters;
			if (clustering > tNodes) {
				clustering = tNodes;
			} else {
				update = false;
			}
			
			j++;
			System.out.println("stop clustering is " + clustering);
		} while ((j < ITERATIONS) && update);

		if (modularity > maxModularity) {
			clusteringDis = voSClusteringTechnique.getClustering();
			maxModularity = modularity;
			System.out.println(" 最终大小,   模块度值:  " + maxModularity);
			maxModularity = modularity;
		}

		if (isPrim) {
			return resKey(clusteringDis);
		}

		if (isWrite) {
			return resKeyMulitMachine(clusteringDis);
		}

		return res(clusteringDis);
	}

	private Map<String, String> res(ClusteringDis clustering) {

		Map<String, String> maps = new HashMap<>();
		int i, nNodes;
		nNodes = clustering.getNNodes();
		clustering.orderClustersByNNodes();
		CrayfishCache<Integer, Long> ids = crayfish.cache("ids");
		for (i = 1; i < nNodes; i++) {
			maps.put(ids.get(i) + "", Integer.toString(clustering.getCluster(i)));
		}

		ids.clear();
		ids.destroy();

		return maps;
	}

	private Map<String, String> resKeyMulitMachine(ClusteringDis clustering) throws GsGraphException {
		crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			Map<String, String> res = new HashMap<>();
			clustering.orderClustersByNNodes();
			CrayfishCache<Integer, Long> ids = crayfish.cache("ids");
			CrayfishCache<Integer, Integer> clusterCache = crayfish.cache(clustering.strClustersName);

			Iterator<Entry<Integer, Long>> idsIterator = ids.localEntries().iterator();
			while (idsIterator.hasNext()) {
				Entry<Integer, Long> entry = idsIterator.next();
				entry.getValue();
				try {
					res.put(graph.getPrimaryKey(entry.getValue()),Integer.toString(clusterCache.localPeek(entry.getKey())));
				} catch (GsGraphException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			System.out.println("该机 点数 :" + res.size());
			
			if (res.keySet().size() != 0)
				GraphHelper.writeFile(res,  "Louvain");
		});

		crayfish.cache("ids").destroy();

		return new HashMap<>();
	}

	private Map<String, String> resKey(ClusteringDis clustering) throws GsGraphException {

		Map<String, String> maps = new HashMap<>();
		int i, nNodes;
		nNodes = clustering.getNNodes();
		clustering.orderClustersByNNodes();
		CrayfishCache<Integer, Long> ids = crayfish.cache("ids");
		for (i = 1; i < nNodes; i++) {
			maps.put(graph.getPrimaryKey(ids.get(i)), Integer.toString(clustering.getCluster(i)));
		}

		ids.clear();
		ids.destroy();

		return maps;
	}

	public NetworkDis init2() {

		NetworkDis network = null;
		long startTime = System.currentTimeMillis();
		createIdsMulit();
		long endTime = System.currentTimeMillis();
		System.out.println("花费时间 : " + (endTime -startTime));

		System.out.println("点索引创建成功！！！");
		// 测试****

		Collection<String> cacheString = crayfish.cacheNames();
		CrayfishCache<Integer, List<LouvainData>> neighbor;
		
		if (cacheString.contains(NEIGHBOR + "_0")) {
			neighbor = crayfish.cache(NEIGHBOR + "_0");
			neighbor.clear();
		} else {
			neighbor = crayfish.createCache(NEIGHBOR + "_0");
		}

		crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {

			CrayfishCache<Long, Integer> idToLong = crayfish.cache("idToLong");
			CrayfishDataStreamer<Integer, List<LouvainData>> nNeighborsStreamer = crayfish.dataStreamer(NEIGHBOR + "_0");

			ExecutorService pool = Pools.DEFAULT;
			int threadCount = Runtime.getRuntime().availableProcessors() * 20;
			List<Callable<Boolean>> work = new ArrayList<>();

			Iterator<Entry<Long, Integer>> iter_vertex = idToLong.localEntries().iterator();
			for (int count = 0; count < threadCount; count++) {
				work.add(new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception {
						// TODO Auto-generated method stub
						while (true) {
							Entry<Long, Integer> node;

							synchronized (iter_vertex) {
								if (iter_vertex.hasNext()) {
									node = iter_vertex.next();
								} else {
									break;
								}
							}

							if (EDGENUM.incrementAndGet() % 4000 == 0)
								System.out.print(".");

							long nodeID = node.getKey();
							int indexF = node.getValue();
							NeighborIterator iteor = graph.getNeighborByIndex(nodeID, GsEdgeDirection.Any, null, null, false);
							
							int toId;
							Map<Integer, Integer> segment = new HashMap<>();
							while (iteor.hasNext()) {
								NeighborValue value = iteor.next();
								toId = idToLong.get(value.getVertexId());
								if (segment.containsKey(toId)) {
									segment.put(toId, segment.get(toId) + 1);
								} else {
									segment.put(toId, 1);
								}
							}

							List<LouvainData> lists = new ArrayList<>();
							for (Integer key : segment.keySet()) {
								lists.add(new LouvainData(key, segment.get(key)));
							}
							// dataStream 直接写
							nNeighborsStreamer.addData(indexF, lists);
						}
						return true;
					}
				});
			}

			try {
				System.out.println("执行");
				List<Future<Boolean>> results = pool.invokeAll(work);
				for (Future<Boolean> result : results) {
					try {
						if (result.get() == false)
							System.out.println("执行失败");
					} catch (ExecutionException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println("执行完成");

			nNeighborsStreamer.flush();
			nNeighborsStreamer.close();
		});

		CrayfishCache<Long, Integer> idToLong = crayfish.cache("idToLong");
		idToLong.close();
		idToLong.destroy();
		System.out.println("销毁 idToLong cache success");

		neighbor.put(0, new ArrayList<>());
		// 这里改为多机的进行

		network = new NetworkDis(NODESIZE + 1, NEIGHBOR + "_0");

		return network;
	}
	
	public void createIdsMulit() {
		System.out.println("清除上了cache");
		Crayfish crayfish = CrayfishManager.Instance.getCrayfish();

		CrayfishCache<Long, Integer> idToLongCache;
		CrayfishCache<Integer, Long> ids;
		CrayfishCache<Integer, Integer> nodeSize;

		
		Collection<String> cacheString = crayfish.cacheNames();
		if (cacheString.contains("idToLong")) {
			idToLongCache = crayfish.cache("idToLong");
			idToLongCache.clear();
		} else {
			idToLongCache = Helper.createIntegerCache("idToLong");
		}

		if (cacheString.contains("ids")) {
			ids = crayfish.cache("ids");
			ids.clear();
		} else {
			ids = crayfish.createCache("ids");
		}

		if (cacheString.contains("nodeSize")) {
			nodeSize = crayfish.cache("nodeSize");
			nodeSize.clear();
		} else {
			nodeSize = crayfish.createCache("nodeSize");
		}

		System.out.println("cache 清除成功！！！");
		crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			graph = (GsCompleteGraph) GraphManager.Instance.getGraph(graphName);
			nodeSize_streamer = crayfish.dataStreamer("nodeSize");
			if (graphName == null) {
				System.out.println("=== graphName === 序列化失败！！！ ");
			}

			int affinity = CrayfishManager.Instance.getAffinity();
			int sum = 0;
			if (VERTEX == null) {
				Set<Integer> setv = graph.meta().getVertexType();
				for (Integer typeIndex : setv) {
					sum += graph.stores().propertyStore().binaryCache(typeIndex).localSize();
				}
			} else {
				int typeIndex = graph.meta().getIndexByVertexTypeName(VERTEX);
				sum += graph.stores().propertyStore().binaryCache(typeIndex).localSize();
			}
			nodeSize_streamer.addData(affinity, sum);
			nodeSize_streamer.flush();
			nodeSize_streamer.close();
			System.out.println("affinity 为 " + affinity + " nodesize is " + sum);
		});

		crayfish.compute(CrayfishManager.Instance.getDBClusterGroup()).broadcast(() -> {
			CrayfishDataStreamer<Integer, Long> ids_streamer = crayfish.dataStreamer("ids");
			CrayfishDataStreamer<Long, Integer> idToLong_streamer = crayfish.dataStreamer("idToLong");
			CrayfishCache<Integer, Integer> nodeSizeT = crayfish.cache("nodeSize");

			int affinity = CrayfishManager.Instance.getAffinity();
			Iterator<ClusterNode> clusters = CrayfishManager.Instance.getDBClusterGroup().nodes().iterator();

			int startFlag = 0;
			int allNode = 0;

			while (clusters.hasNext()) {
				ClusterNode node = clusters.next();
				int nodeAffient = Integer.parseInt(node.attribute(CrayfishManager.AFFINITY_ATTRIBUTE_KEY));
				int localNode = nodeSizeT.get(nodeAffient);
				if (nodeAffient < affinity) {
					startFlag += localNode;
				}
				allNode += localNode;
			}
			System.out.println("startFlag is " + startFlag);

			NODESIZE = allNode;

			long start = System.currentTimeMillis();
			//这边赫曼
			if (VERTEX == null) {
				Set<Integer> setv = graph.meta().getVertexType();
				for (Integer typeIndex : setv) {
					Iterator<Entry<Long, BinaryObject>> iter_vertexs = graph.stores().propertyStore().binaryCache(typeIndex).localEntries().iterator();
					while (iter_vertexs.hasNext()) {
						Entry<Long, BinaryObject> item = iter_vertexs.next();
						Long key = item.getKey();
						idToLong_streamer.addData(key, ++startFlag);
						ids_streamer.addData(startFlag, key);
					}
				}
				
			} else {
				int typeIndex = graph.meta().getIndexByVertexTypeName(VERTEX);
				Iterator<Entry<Long, BinaryObject>> iter_vertexs = graph.stores().propertyStore().binaryCache(typeIndex)
						.localEntries().iterator();
				while (iter_vertexs.hasNext()) {
					Entry<Long, BinaryObject> item = iter_vertexs.next();
					Long key = item.getKey();
					idToLong_streamer.addData(key, ++startFlag);
					ids_streamer.addData(startFlag, key);
				}
			}
			
			long end = System.currentTimeMillis();
			System.out.println("spend TIme is " + (end - start));

			ids_streamer.flush();
			ids_streamer.close();
			idToLong_streamer.flush();
			idToLong_streamer.close();
			
		});
	}
}
