package singleMachine;

import java.io.Serializable;
import java.util.Random;

/**
 * 
* @Title: Network.java
* @Package singleMachine
* @Description: TODO(用一句话描述该文件做什么)
* @author zhang
* @date 2019年7月21日 下午7:09:43
* @version V1.0
 */

public class Network implements Cloneable, Serializable {
	/**
	 * 用来表示网络的数据结构.
	 * e.g. 数据类型如下:
	 * 0: {3,4}
	 * 1: {3}
	 * 3: {0,1}
	 * 4: {0}
	 * nNodes				4
	 * firstNeighborIndex	[0,2,3,5]
	 * neighbor				[3,4,3,0,1,0]
	 * edgeWeight			[1,1,1,1,1,1]
	 */
	private static final long serialVersionUID = 1;
	
	private int nNodes; 
	private int [] firstNeighborIndex; 
	private int [] neighbor;
	private double [] edgeWeight;
	private double [] nodeWeight;
	
	private int nClusters;
	private int [] cluster; // cluster[i] 表示节点i的簇id
	
	private double [] clusterWeight;
	private int [] nNodesPerCluster;
	private int [][] nodePerCluster;
	private boolean cluteringStatsAvailable;
	
	private Network () {}
	
	public Network (int nNodes, int [] firstNeighborIndex, int [] neighbor,
			double [] edgeWeight, double [] nodeWeight) {
		this(nNodes, firstNeighborIndex, neighbor, edgeWeight, nodeWeight, null);
	}
	
	public Network (int nNodes, int [] firstNeighborIndex, int [] neighbor,
			double [] edgeWeight, double [] nodeWeight, int [] cluster) {
		int i, nEdge;
		this.nNodes = nNodes;
		this.firstNeighborIndex = firstNeighborIndex;
		this.neighbor = neighbor;
		if (edgeWeight == null) {
			nEdge = neighbor.length;
			this.edgeWeight = new double[nEdge];
			for (i = 0; i < nEdge; i++) {
				this.edgeWeight[i] = 1;
			}
		} else {
				this.edgeWeight = edgeWeight;
			}
		
		if (nodeWeight == null) {
			this.nodeWeight = new double[nNodes];
			for (i = 0; i < nNodes; i++) {
				this.nodeWeight[i] = 1;
			}
		} else {
			this.nodeWeight = nodeWeight;
		}
		
	}

	public void initSingletonClusters () {
		int i;
		this.nClusters = this.nNodes;
		this.cluster = new int[this.nClusters];
		for (i = 0; i < this.nClusters; i++) {
			cluster[i] = i;
		}
		
		deleteClusteringStats();
	}
	
	public Network getReducedNetwork () {
		Network reducedNetwork;
		double [] reducedNetworkEdgeWeight1, reducedNetworkEdgeWeight2;
		int [] reducedNetworkNeighbor1, reducedNetworkNeighbor2;
		int reducedNetworkNEdges1, reducedNetworkNEdges2;
		int i, j, k, m, index;
		
		if (this.cluster == null) {
			return null;
		}
		if (!this.cluteringStatsAvailable) {
			this.calcClusteringStats();
		}
		
		reducedNetwork = new Network();
		reducedNetwork.nNodes = this.nClusters;
		reducedNetwork.firstNeighborIndex = new int[this.nClusters + 1];
		reducedNetwork.nodeWeight = new double [this.nClusters];
		
		reducedNetworkEdgeWeight2 = new double [this.nClusters];
		reducedNetworkNeighbor2 = new int [this.nClusters - 1];
		
		reducedNetworkNeighbor1 = new int [this.neighbor.length];
		reducedNetworkEdgeWeight1 = new double [this.edgeWeight.length];
		
		
		reducedNetworkNEdges1 = 0;
		
		for (i = 0; i < this.nClusters; i++) {
			reducedNetworkNEdges2 = 0;
			for (j = 0; j < this.nNodesPerCluster[i]; j++) {
				k = this.nodePerCluster[i][j]; // 簇i中有节点k
				for(index = this.firstNeighborIndex[k]; index < this.firstNeighborIndex[k+1]; index ++) {
					m = cluster[this.neighbor[index]];
					if (m != i) {
						if (reducedNetworkEdgeWeight2[m] == 0) {
							reducedNetworkNeighbor2[reducedNetworkNEdges2] = m;
							reducedNetworkNEdges2 += 1;
						}
						reducedNetworkEdgeWeight2[m] += this.edgeWeight[index];
					}
				}
				reducedNetwork.nodeWeight[i] += this.nodeWeight[k];
			}
			
			for (j = 0; j < reducedNetworkNEdges2; j++) {
				reducedNetworkNeighbor1[reducedNetworkNEdges1 + j] = reducedNetworkNeighbor2[j];
				reducedNetworkEdgeWeight1[reducedNetworkNEdges1 + j] = reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]];
				reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]] = 0;
			}
			reducedNetworkNEdges1 += reducedNetworkNEdges2;
			reducedNetwork.firstNeighborIndex[i+1] = reducedNetworkNEdges1;
		}
		reducedNetwork.neighbor = new int [reducedNetworkNEdges1];
		reducedNetwork.edgeWeight = new double [reducedNetworkNEdges1];
		
		System.arraycopy(reducedNetworkNeighbor1, 0, reducedNetwork.neighbor, 0, reducedNetworkNEdges1);
		System.arraycopy(reducedNetworkEdgeWeight1, 0, reducedNetwork.edgeWeight, 0, reducedNetworkNEdges1);
		return reducedNetwork;
	}

	public boolean runLocalMovingAlgorithm (double resolution) {
		/**
		 * @author zhang
		 * @param resulution : 1/m
		 */
		return this.runLocalMovingAlgorithm(resolution, new Random());
	}
	
	public boolean runLocalMovingAlgorithm (double resolution, Random random) {
		/**
		 * deltaQ = k_{i,in} - 
		 */
		int i, j, k, l, bestCluster, nUnusedClusters, nStableNodes, nNeighboringClusters, nClusters;
		boolean update;
		double maxDeltaQuality, deltaQuality;
		
		double [] clusterWeight, edgeWeightPerCluster;
		int [] nodeOrder, nNodesPerCluster, unusedCluster, neighboringCluster, newCluster;
		
		if (this.cluster == null || this.nNodes == 1) {
			return false;
		}
		clusterWeight = new double [this.nNodes];
		nNodesPerCluster = new int [this.nNodes];
		
		for (i = 0; i < this.nNodes; i++) {
			clusterWeight[this.cluster[i]] += this.nodeWeight[i];
			nNodesPerCluster[this.cluster[i]] ++;
		}
		
		update = false;
		// 顺序列表的产生
		nodeOrder = new int [this.nNodes];
		for (i = 0; i < this.nNodes; i++) {
			nodeOrder[i] = i;
		}
		for (i = 0; i < this.nNodes; i++) {
			j = random.nextInt(nNodes);
			k = nodeOrder[i];
			nodeOrder[i] = nodeOrder[j];
			nodeOrder[j] = k;
		}
		// 初始化统计不含节点的cluster的信息
		nUnusedClusters = 0;
		unusedCluster = new int [this.nNodes]; 
		for (i = 0; i < this.nNodes; i++) {
			if (nNodesPerCluster[i] == 0) {
				unusedCluster[nUnusedClusters] = i;
				nUnusedClusters += 1;
			}
		}
		
		edgeWeightPerCluster = new double [this.nNodes];
		neighboringCluster = new int [this.nNodes - 1];
		nStableNodes = 0;
		i = 0;
		
		do {
			j = nodeOrder[i];
			nNeighboringClusters = 0;
			for (k = this.firstNeighborIndex[j]; k < this.firstNeighborIndex[j + 1]; k++) {
				l = this.cluster[this.neighbor[k]];
				if (edgeWeightPerCluster[l] == 0) {
					neighboringCluster[nNeighboringClusters] = l;
					nNeighboringClusters += 1;
				}
				edgeWeightPerCluster[l] += edgeWeight[k];
			}
			clusterWeight[this.cluster[j]] -= this.nodeWeight[j];
			nNodesPerCluster[cluster[j]]--;
			if (nNodesPerCluster[cluster[j]] == 0) {
				unusedCluster[nUnusedClusters] = cluster[j];
				nUnusedClusters += 1;
			}
			bestCluster = -1;
			maxDeltaQuality = 0;
			deltaQuality = -1;
			for (k = 0; k < nNeighboringClusters; k++) {
				l = neighboringCluster[k];
				deltaQuality = edgeWeightPerCluster[l] - this.nodeWeight[j] * clusterWeight[l] * resolution;
				if ((deltaQuality > maxDeltaQuality) || (deltaQuality == maxDeltaQuality && l < bestCluster)) {
					bestCluster = l;
					maxDeltaQuality = deltaQuality;
				}
				edgeWeightPerCluster[l] = 0;
				if (maxDeltaQuality == 0) { // 当前节点不需要移动进入其它的簇
					nUnusedClusters--;
					bestCluster = unusedCluster[nUnusedClusters];
				}
				clusterWeight[bestCluster] += this.nodeWeight[j];
				nNodesPerCluster[bestCluster]++;
				if (bestCluster == cluster[j]) {
					nStableNodes++;
				} else {
					cluster[j] = bestCluster;
					nStableNodes = 1;
					update = true;
				}
				i = (i < this.nNodes - 1) ? i + 1 : 0;
			}
		} while (nStableNodes < this.nNodes);
		
		newCluster = new int [this.nNodes];
		nClusters = 0;
		for (i = 0; i < this.nNodes; i++) {
			if (nNodesPerCluster[i] > 0) {
				newCluster[i] = nClusters;
				nClusters++;
			}
		}
		
		for (i = 0; i < this.nNodes; i++) {
			this.cluster[i] = newCluster[cluster[i]];
		}
		this.deleteClusteringStats();
		return update;
	}
	
	public void calcClusteringStats () {
		int i, k;
		this.clusterWeight = new double [this.nClusters];
		this.nNodesPerCluster = new int [this.nClusters];
		this.nodePerCluster = new int [this.nClusters][];
		
		for (i = 0; i < this.nNodes; i++) {
			this.clusterWeight[cluster[i]] += this.nodeWeight[i];
			this.nNodesPerCluster[cluster[i]] += 1;
		}
		
		for (i = 0; i < this.nClusters; i++) {
			this.nNodesPerCluster[i] = 0;
		}
		
		for (i = 0; i < this.nNodes; i++) {
			k = this.cluster[i];
			this.nodePerCluster[k][this.nNodesPerCluster[k]] = i;
			this.nNodesPerCluster[k] += 1;
		}
		
		this.cluteringStatsAvailable = true;
	};
	
	public void deleteClusteringStats () {};

 	public int getNnodes () {
		return this.nNodes;
	}
	
	public int getNEdges () {
		return this.neighbor.length;
	}

	public double getTotalEdgeWeight () {
		double totalWeight;
		int i;
		totalWeight = 0;
		for (i = 0; i < this.neighbor.length; i++) {
			totalWeight += this.edgeWeight[i];
		}
		return totalWeight;
	}

	public double getTotalNodeWeight () {
		double totalWeight;
		int i;
		totalWeight = 0;
		for (i = 0; i < this.nNodes; i++) {
			totalWeight += this.nodeWeight[i];
		}
		return totalWeight;
	}
	
	public double [] getNodeWeight () {
		return this.nodeWeight;
	}
	
	public double [] getEdgeWeight () {
		return this.edgeWeight;
	}

	public int getNCluster () {return this.nClusters;}
	
	public int [] getClusters () {return this.cluster;}
}
