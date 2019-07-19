package Main;

import java.io.Serializable;
import java.util.Random;
 
public class Network implements Cloneable, Serializable {
	private static final long serialVersionUID = 1;
 
	private int nNodes; // 网络的节点数
	private int[] firstNeighborIndex; // 去重
	private int[] neighbor; // 有几条edge长度就有多长
	private double[] edgeWeight; // 长度与neighbor数组相同
 
	private double[] nodeWeight;  // 节点的权重, 长度与nNodes相同
	private int nClusters;  // nClusters表示簇的数量
	private int[] cluster;  // 用来保存簇的信息, cluster[i] 表示the cluster id of node_i
 
	private double[] clusterWeight;  // cluster的权重
	private int[] nNodesPerCluster;  // 每个cluster保存的节点的数量
	private int[][] nodePerCluster;  // 保存了每个cluster中存在的节点的id
	private boolean clusteringStatsAvailable;  // 表示是否需要计算cluster的信息
 
//	用数组保存的neighbor是什么样的
//	cluster的数据长什么样
//	第i个node对应的cluster标号
	
	public Network(int nNodes, int[] firstNeighborIndex, int[] neighbor, double[] edgeWeight, double[] nodeWeight) {
//		构造器
		this(nNodes, firstNeighborIndex, neighbor, edgeWeight, nodeWeight, null);
	}
 
	public Network(int nNodes, int[] firstNeighborIndex, int[] neighbor, double[] edgeWeight, double[] nodeWeight,
			int[] cluster) {
		int i, nEdges;
		this.nNodes = nNodes;
		this.firstNeighborIndex = firstNeighborIndex;
		this.neighbor = neighbor;
		if (edgeWeight == null) {
			nEdges = neighbor.length;
			this.edgeWeight = new double[nEdges];
			for (i = 0; i < nEdges; i++)
				this.edgeWeight[i] = 1;
		} else
			this.edgeWeight = edgeWeight;
 
		if (nodeWeight == null) {
			this.nodeWeight = new double[nNodes];
			for (i = 0; i < nNodes; i++)
				this.nodeWeight[i] = 1;
		} else
			this.nodeWeight = nodeWeight;
	}
 
	public int getNNodes() {
		return nNodes;
	}
 
	public int getNEdges() {
		return neighbor.length;
	}
 
	public double getTotalEdgeWeight() // 计算总边权
	{
		double totalEdgeWeight;
		int i;
		totalEdgeWeight = 0;
		for (i = 0; i < neighbor.length; i++)
			totalEdgeWeight += edgeWeight[i];
		return totalEdgeWeight;
	}
 
	public double[] getEdgeWeights() {
		return edgeWeight;
	}
 
	public double[] getNodeWeights() {
		return nodeWeight;
	}
 
	public int getNClusters() {
		return nClusters;
	}
 
	public int[] getClusters() {
		return cluster;
	}
 
	public void initSingletonClusters() {
		int i;
		nClusters = nNodes;
		cluster = new int[nNodes];
		for (i = 0; i < nNodes; i++)
			cluster[i] = i;
 
		deleteClusteringStats();
	}
 
	public void mergeClusters(int[] newCluster) {
//		todo: mergeClusters的原理是什么(暂时用不上)
		int i, j, k;
		if (cluster == null)
			return;
		i = 0;
		for (j = 0; j < nNodes; j++) {
			k = newCluster[cluster[j]]; // 节点j的簇在新的簇里面的id
			if (k > i)
				i = k;
			cluster[j] = k;
		}
		nClusters = i + 1;
		deleteClusteringStats();
	}
 
	public Network getReducedNetwork() {
		/*
		* 生成一个压缩的网络
		* nNodes:
		* 	数组的长度是 nClusters;
		* todo: 为什么fistNrighborIndex 的长度是 nCluster+1
		* nodeWeight:
		*	数据的长度是 nCluster
		* todo: 什么是reduce1, reduce2
		*/
		double[] reducedNetworkEdgeWeight1, reducedNetworkEdgeWeight2;
		int i, j, k, l, m, reducedNetworkNEdges1, reducedNetworkNEdges2;
		int[] reducedNetworkNeighbor1, reducedNetworkNeighbor2;
		Network reducedNetwork;
		if (cluster == null)
			return null;
		if (!clusteringStatsAvailable)
			calcClusteringStats();
 
		// 实例化一个压缩的网络
		reducedNetwork = new Network();
		reducedNetwork.nNodes = nClusters;
		reducedNetwork.firstNeighborIndex = new int[nClusters + 1];
		reducedNetwork.nodeWeight = new double[nClusters];
		
		reducedNetworkNeighbor1 = new int[neighbor.length];
		reducedNetworkEdgeWeight1 = new double[edgeWeight.length];
		
		/*
		* 节点i表示压缩之后网络的节点, 对应于压缩前的簇i
		* reducedNetworkNeighbor2: 在压缩的网络中, 节点i相连的neighbor数组. 不考虑自连接的情况
		* reducedNetworkEdgeWeight2[k]: 在压缩的网络中, 节点i与节点k所连接的边的权重
		*/
		reducedNetworkNeighbor2 = new int[nClusters - 1]; 
		reducedNetworkEdgeWeight2 = new double[nClusters]; 
		reducedNetworkNEdges1 = 0;
		
		for (i = 0; i < nClusters; i++) { // 开始遍历cluster, 也就是压缩后网络的节点
			reducedNetworkNEdges2 = 0; // reducedNetworkNEdges2: 在压缩的网络中, 与节点i相连接的边的计数
			for (j = 0; j < nodePerCluster[i].length; j++) {
				// 遍历 cluster_i 上面的节点
				k = nodePerCluster[i][j]; // k是簇i中第j个节点的id
				for (l = firstNeighborIndex[k]; l < firstNeighborIndex[k + 1]; l++) {
					// neighbor[l] 是与id为k的节点相邻的节点
					m = cluster[neighbor[l]]; // m是k的在l位置的邻居节点所属的簇id
					if (m != i) { // m 和 i对应不同的簇的时候 => 压缩的网络的节点m和节点i相连接
						if (reducedNetworkEdgeWeight2[m] == 0) {
							// 首次发现压缩之后的网络上, 当前节点i与第m个节点相连的时候
							reducedNetworkNeighbor2[reducedNetworkNEdges2] = m; // 与当前节点相连接的第reducedNetworkNEdges2个节点是m
							reducedNetworkNEdges2++; // 计数+1
						}
						// 把原来第l条边的权重赋给压缩后的网络的当前节点
						reducedNetworkEdgeWeight2[m] += edgeWeight[l]; 
					}
				}
				reducedNetwork.nodeWeight[i] += nodeWeight[k]; // 将簇i中的节点的权重累加到簇i上
			}

			// 遍历在压缩前的网络中, 与簇i相连接的其他的簇
			for (j = 0; j < reducedNetworkNEdges2; j++) { // 遍历在压缩的网络中, 与节点i相连的节点
				reducedNetworkNeighbor1[reducedNetworkNEdges1 + j] = reducedNetworkNeighbor2[j];
				reducedNetworkEdgeWeight1[reducedNetworkNEdges1 + j] = reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]];
						// todo: 这个地方为什么要置未0?
						// 反复使用变量, 需要置空
				reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]] = 0;
			}
			reducedNetworkNEdges1 += reducedNetworkNEdges2; // 在压缩之前, 与簇i相连接的簇的数量(不包括本身). 并且在数值上做累计
			reducedNetwork.firstNeighborIndex[i + 1] = reducedNetworkNEdges1; // todo: 这个地方的firstNeighborIndex 不是很明白是怎么回事
		}
		reducedNetwork.neighbor = new int[reducedNetworkNEdges1];
		reducedNetwork.edgeWeight = new double[reducedNetworkNEdges1];
		System.arraycopy(reducedNetworkNeighbor1, 0, reducedNetwork.neighbor, 0, reducedNetworkNEdges1);
		System.arraycopy(reducedNetworkEdgeWeight1, 0, reducedNetwork.edgeWeight, 0, reducedNetworkNEdges1);
		return reducedNetwork;
	}
	
	public double calcQualityFunction(double resolution) {
		double qualityFunction, totalEdgeWeight;
		int i, j, k;
		if (cluster == null)
			return Double.NaN;
		if (!clusteringStatsAvailable)
			calcClusteringStats();
		qualityFunction = 0;
		totalEdgeWeight = 0;
		for (i = 0; i < nNodes; i++) {
			j = cluster[i]; //  节点i的cluster
			for (k = firstNeighborIndex[i]; k < firstNeighborIndex[i + 1]; k++) {
				// 节点i相连的节点在neighbor中的索引值
				if (cluster[neighbor[k]] == j) // 属于同一个簇
					qualityFunction += edgeWeight[k]; // 簇内的边的权重的累计
				totalEdgeWeight += edgeWeight[k]; // 总权重
			}
		}
		for (i = 0; i < nClusters; i++)
			qualityFunction -= clusterWeight[i] * clusterWeight[i] * resolution;
		qualityFunction /= totalEdgeWeight;
		return qualityFunction;
	}
	
	public boolean runLocalMovingAlgorithm(double resolution) {
		return runLocalMovingAlgorithm(resolution, new Random());
	}
	public boolean runLocalMovingAlgorithm(double resolution, Random random) {
		boolean update;
		double maxQualityFunction, qualityFunction;
		double[] clusterWeight, edgeWeightPerCluster;
		int bestCluster, i, j, k, l, nNeighboringClusters, nStableNodes, nUnusedClusters;
		int[] neighboringCluster, newCluster, nNodesPerCluster, nodeOrder, unusedCluster;
		if ((cluster == null) || (nNodes == 1))
			return false;
		update = false;
		clusterWeight = new double[nNodes];
		nNodesPerCluster = new int[nNodes];
		for (i = 0; i < nNodes; i++) {
			/*
			  *       使用clusterWeight和nNodesPerCluster来保存移动之后的簇, 用这个来更新网络的信息
			 *  clusterWeight: 节点weight的累计. 保存节点i所在的cluster的weight
			 *  nNodesPerCluster: 节点数量的统计
			 */
			clusterWeight[cluster[i]] += nodeWeight[i];
			nNodesPerCluster[cluster[i]]++;
		}
		// todo: unUsedCluster是做什么用的
		nUnusedClusters = 0;
		unusedCluster = new int[nNodes];
		for (i = 0; i < nNodes; i++) 
			if (nNodesPerCluster[i] == 0) { // 如果cluster没有节点, 逐次存储内容
				unusedCluster[nUnusedClusters] = i;
				nUnusedClusters++;
			}
		nodeOrder = new int[nNodes];
		for (i = 0; i < nNodes; i++)
			nodeOrder[i] = i;
		for (i = 0; i < nNodes; i++) { // 随机交换n次nodeOrder的顺序
			j = random.nextInt(nNodes);
			k = nodeOrder[i];
			nodeOrder[i] = nodeOrder[j];
			nodeOrder[j] = k;
		}
		
		// todo: edgeWeightPerCluster和neighboringCluster的作用 
		edgeWeightPerCluster = new double[nNodes];
		neighboringCluster = new int[nNodes - 1];
		nStableNodes = 0; // 稳定点的计数
		i = 0;
		do {
			j = nodeOrder[i];
			// 随机选择节点j进行增益的计算
			nNeighboringClusters = 0; // 相邻的簇的计数
			for (k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++) { // 与随机节点j相连接的neighbor的index
				l = cluster[neighbor[k]]; // 节点j相连接的簇的id
				if (edgeWeightPerCluster[l] == 0) { // 首次扫描这个相邻的簇
					neighboringCluster[nNeighboringClusters] = l;
					nNeighboringClusters++;
				}
				edgeWeightPerCluster[l] += edgeWeight[k]; // 累加计算这个簇的权重, \ Sigma_{i, in}
			}
			// 将节点j从原来的簇中移除
			clusterWeight[cluster[j]] -= nodeWeight[j];
			nNodesPerCluster[cluster[j]]--;
			if (nNodesPerCluster[cluster[j]] == 0) { // 如果移除这个节点之后, 簇为空
				unusedCluster[nUnusedClusters] = cluster[j];
				nUnusedClusters++;
			}
			bestCluster = -1;
			maxQualityFunction = 0;
			for (k = 0; k < nNeighboringClusters; k++) {
				l = neighboringCluster[k]; // l是节点j相邻的cluster
				qualityFunction = edgeWeightPerCluster[l] - nodeWeight[j] * clusterWeight[l] * resolution;
				if ((qualityFunction > maxQualityFunction)
						|| ((qualityFunction == maxQualityFunction) && (l < bestCluster))) { // todo: 后面那条规则的意思?
					bestCluster = l;
					maxQualityFunction = qualityFunction;
				}
				edgeWeightPerCluster[l] = 0;
			}
			if (maxQualityFunction == 0) { // 等于0意味着不需要进行簇的转移
				bestCluster = unusedCluster[nUnusedClusters - 1];
				nUnusedClusters--;
			}
			// 前面先减去是为了后面能够统一操作
			clusterWeight[bestCluster] += nodeWeight[j];
			nNodesPerCluster[bestCluster]++;
			if (bestCluster == cluster[j])
				nStableNodes++;
			else {
				cluster[j] = bestCluster;
				nStableNodes = 1;
				update = true;
			}
			i = (i < nNodes - 1) ? (i + 1) : 0; // 不停地循环nodeOrder中的顺序
		} while (nStableNodes < nNodes); // 优化步骤是直到所有的点都稳定下来才结束
		newCluster = new int[nNodes];
		nClusters = 0; // 更新clusters
		for (i = 0; i < nNodes; i++)
			if (nNodesPerCluster[i] > 0) {
				// 簇内有节点
				newCluster[i] = nClusters;
				nClusters++;
			}
		for (i = 0; i < nNodes; i++)
			cluster[i] = newCluster[cluster[i]];
		deleteClusteringStats();
		return update;
	}
	
	public boolean runLouvainAlgorithm(double resolution) {
		return runLouvainAlgorithm(resolution, new Random());
	}
	
	public boolean runLouvainAlgorithm(double resolution, Random random) {
		boolean update, update2;
		Network reducedNetwork;
		if ((cluster == null) || (nNodes == 1))
			return false;
		update = runLocalMovingAlgorithm(resolution, random);
		if (nClusters < nNodes) {
			reducedNetwork = getReducedNetwork();
			reducedNetwork.initSingletonClusters();
			update2 = reducedNetwork.runLouvainAlgorithm(resolution, random);
			if (update2) {
				update = true;
				mergeClusters(reducedNetwork.getClusters());
			}
		}
		deleteClusteringStats();
		return update;
	}
	
	private Network() {
	}
	
	private void calcClusteringStats() {
		/*
		 * 
		 */
		int i, j;
		clusterWeight = new double[nClusters];
		nNodesPerCluster = new int[nClusters];
		nodePerCluster = new int[nClusters][];
//		第一步: 对所有的点进行扫描
//		第一步: 计算每个cluster的weight
//		第一步: 每个cluster含有的节点的计数+1
		for (i = 0; i < nNodes; i++) {
			clusterWeight[cluster[i]] += nodeWeight[i];
			nNodesPerCluster[cluster[i]]++;
		}
//		第二步: 对所有的cluster进行扫描
//		第二步: 对每个cluster初始化一个对应存储的节点数量的数组
//		第二步: 每个cluster对应的节点数重置为0
		for (i = 0; i < nClusters; i++) {
			nodePerCluster[i] = new int[nNodesPerCluster[i]];
			nNodesPerCluster[i] = 0;
		}
//		第三步: 扫描所有的节点, 用来更新每个cluster中含有的节点的index
		for (i = 0; i < nNodes; i++) {
			j = cluster[i]; // node_i 对应的 cluster index为j
			nodePerCluster[j][nNodesPerCluster[j]] = i;
			nNodesPerCluster[j]++;
		}
		clusteringStatsAvailable = true;
	}
	
	private void deleteClusteringStats() {
		/*
		 * 
		 */
		clusterWeight = null;
		nNodesPerCluster = null;
		nodePerCluster = null;
		clusteringStatsAvailable = false;
	}
}