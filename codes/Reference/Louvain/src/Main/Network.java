package Main;

import java.io.Serializable;
import java.util.Random;
 
public class Network implements Cloneable, Serializable {
	private static final long serialVersionUID = 1;
 
	private int nNodes; // ����Ľڵ���
	private int[] firstNeighborIndex; // ȥ�������, �б����ӵĽڵ�(idС���Ǹ��ڵ�)
	
	private int[] neighbor; // �м���edge���Ⱦ��ж೤
	private double[] edgeWeight; // ������neighbor������ͬ
 
	private double[] nodeWeight;  // �ڵ��Ȩ��, ������nNodes��ͬ
	private int nClusters;  // nClusters��ʾ�ص�����
	private int[] cluster;  // ��������ص���Ϣ, cluster[i] ��ʾthe cluster id of node_i
 
	private double[] clusterWeight;  // cluster��Ȩ��
	private int[] nNodesPerCluster;  // ÿ��cluster����Ľڵ������
	private int[][] nodePerCluster;  // ������ÿ��cluster�д��ڵĽڵ��id
	private boolean clusteringStatsAvailable;  // ��ʾ�Ƿ���Ҫ����cluster����Ϣ
 
//	�����鱣���neighbor��ʲô����
//	cluster�����ݳ�ʲô��
//	��i��node��Ӧ��cluster���
	
	public Network(int nNodes, int[] firstNeighborIndex, int[] neighbor, double[] edgeWeight, double[] nodeWeight) {
//		������
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
 
	public double getTotalEdgeWeight() // �����ܱ�Ȩ
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
//		todo: mergeClusters��ԭ����ʲô
		int i, j, k;
		if (cluster == null)
			return;
		i = 0;
		for (j = 0; j < nNodes; j++) {
			k = newCluster[cluster[j]];
			if (k > i)
				i = k;
			cluster[j] = k;
		}
		nClusters = i + 1;
		deleteClusteringStats();
	}
 
	public Network getReducedNetwork() {
		/*
		* ����һ��ѹ��������
		* nNodes:
		* 	����ĳ����� nClusters;
		* todo: ΪʲôfistNrighborIndex �ĳ����� nCluster+1
		* nodeWeight:
		*	���ݵĳ����� nCluster
		* todo: ʲô��reduce1, reduce2
		*/
		double[] reducedNetworkEdgeWeight1, reducedNetworkEdgeWeight2;
		int i, j, k, l, m, reducedNetworkNEdges1, reducedNetworkNEdges2;
		int[] reducedNetworkNeighbor1, reducedNetworkNeighbor2;
		Network reducedNetwork;
		if (cluster == null)
			return null;
		if (!clusteringStatsAvailable)
			calcClusteringStats();
 
		// ʵ����һ��ѹ��������
		reducedNetwork = new Network();
		reducedNetwork.nNodes = nClusters;
		reducedNetwork.firstNeighborIndex = new int[nClusters + 1];
		reducedNetwork.nodeWeight = new double[nClusters];
		
		reducedNetworkNeighbor1 = new int[neighbor.length];
		reducedNetworkEdgeWeight1 = new double[edgeWeight.length];
		
		// ���������ע����, �ڵ�i��ʾ��ǰ�ڵ�
		reducedNetworkNeighbor2 = new int[nClusters - 1]; // reducedNetworkNeighbor2: ��ѹ����������, �ڵ�i������neighbor����
		reducedNetworkEdgeWeight2 = new double[nClusters]; // reducedNetworkEdgeWeight2[k]: ��ѹ����������, �ڵ�i��ڵ�k�����ӵıߵ�Ȩ��
		reducedNetworkNEdges1 = 0;
		
		for (i = 0; i < nClusters; i++) {
			reducedNetworkNEdges2 = 0; // reducedNetworkNEdges2: ��ѹ����������, ��ڵ�i�����ӵıߵļ���
			// cluster_i
			for (j = 0; j < nodePerCluster[i].length; j++) {
				// ���� cluster_i ����Ľڵ�
				k = nodePerCluster[i][j]; // k�Ǵ�i�е�j���ڵ��id
				for (l = firstNeighborIndex[k]; l < firstNeighborIndex[k + 1]; l++) {
					// neighbor[l] ����idΪk�Ľڵ����ڵĽڵ�
					m = cluster[neighbor[l]]; // m��k����lλ�õ��ھӽڵ������Ĵ�id
					if (m != i) { // m �� i��Ӧ��ͬ�Ĵص�ʱ�� => ѹ��������Ľڵ�m�ͽڵ�i������
						if (reducedNetworkEdgeWeight2[m] == 0) {
							// �״η���ѹ��֮���������, ��ǰ�ڵ�i���m���ڵ�������ʱ��
							reducedNetworkNeighbor2[reducedNetworkNEdges2] = m; // �뵱ǰ�ڵ������ӵĵ�reducedNetworkNEdges2���ڵ���m
							reducedNetworkNEdges2++; // ����+1
						}
						// ��ԭ����l���ߵ�Ȩ�ظ���ѹ���������ĵ�ǰ�ڵ�
						reducedNetworkEdgeWeight2[m] += edgeWeight[l]; 
					}
				}
				reducedNetwork.nodeWeight[i] += nodeWeight[k]; // ����i�еĽڵ��Ȩ���ۼӵ���i��
			}
			for (j = 0; j < reducedNetworkNEdges2; j++) { // ������ѹ����������, ��ڵ�i�����Ľڵ�
				reducedNetworkNeighbor1[reducedNetworkNEdges1 + j] = reducedNetworkNeighbor2[j];
				reducedNetworkEdgeWeight1[reducedNetworkNEdges1
						+ j] = reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]];
				reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]] = 0;
			}
			reducedNetworkNEdges1 += reducedNetworkNEdges2; // ѹ����������, ��ڵ�i�����ӵĽڵ������
			reducedNetwork.firstNeighborIndex[i + 1] = reducedNetworkNEdges1;
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
			j = cluster[i];
			for (k = firstNeighborIndex[i]; k < firstNeighborIndex[i + 1]; k++) {
				if (cluster[neighbor[k]] == j)
					qualityFunction += edgeWeight[k];
				totalEdgeWeight += edgeWeight[k];
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
			clusterWeight[cluster[i]] += nodeWeight[i];
			nNodesPerCluster[cluster[i]]++;
		}
		nUnusedClusters = 0;
		unusedCluster = new int[nNodes];
		for (i = 0; i < nNodes; i++)
			if (nNodesPerCluster[i] == 0) {
				unusedCluster[nUnusedClusters] = i;
				nUnusedClusters++;
			}
		nodeOrder = new int[nNodes];
		for (i = 0; i < nNodes; i++)
			nodeOrder[i] = i;
		for (i = 0; i < nNodes; i++) {
			j = random.nextInt(nNodes);
			k = nodeOrder[i];
			nodeOrder[i] = nodeOrder[j];
			nodeOrder[j] = k;
		}
		edgeWeightPerCluster = new double[nNodes];
		neighboringCluster = new int[nNodes - 1];
		nStableNodes = 0;
		i = 0;
		do {
			j = nodeOrder[i];
			nNeighboringClusters = 0;
			for (k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++) {
				l = cluster[neighbor[k]];
				if (edgeWeightPerCluster[l] == 0) {
					neighboringCluster[nNeighboringClusters] = l;
					nNeighboringClusters++;
				}
				edgeWeightPerCluster[l] += edgeWeight[k];
			}
			clusterWeight[cluster[j]] -= nodeWeight[j];
			nNodesPerCluster[cluster[j]]--;
			if (nNodesPerCluster[cluster[j]] == 0) {
				unusedCluster[nUnusedClusters] = cluster[j];
				nUnusedClusters++;
			}
			bestCluster = -1;
			maxQualityFunction = 0;
			for (k = 0; k < nNeighboringClusters; k++) {
				l = neighboringCluster[k];
				qualityFunction = edgeWeightPerCluster[l] - nodeWeight[j] * clusterWeight[l] * resolution;
				if ((qualityFunction > maxQualityFunction)
						|| ((qualityFunction == maxQualityFunction) && (l < bestCluster))) {
					bestCluster = l;
					maxQualityFunction = qualityFunction;
				}
				edgeWeightPerCluster[l] = 0;
			}
			if (maxQualityFunction == 0) {
				bestCluster = unusedCluster[nUnusedClusters - 1];
				nUnusedClusters--;
			}
			clusterWeight[bestCluster] += nodeWeight[j];
			nNodesPerCluster[bestCluster]++;
			if (bestCluster == cluster[j])
				nStableNodes++;
			else {
				cluster[j] = bestCluster;
				nStableNodes = 1;
				update = true;
			}
			i = (i < nNodes - 1) ? (i + 1) : 0;
		} while (nStableNodes < nNodes); // �Ż�������ֱ�����еĵ㶼�ȶ������Ž���
		newCluster = new int[nNodes];
		nClusters = 0;
		for (i = 0; i < nNodes; i++)
			if (nNodesPerCluster[i] > 0) {
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
//		��һ��: �����еĵ����ɨ��
//		��һ��: ����ÿ��cluster��weight
//		��һ��: ÿ��cluster���еĽڵ�ļ���+1
		for (i = 0; i < nNodes; i++) {
			clusterWeight[cluster[i]] += nodeWeight[i];
			nNodesPerCluster[cluster[i]]++;
		}
//		�ڶ���: �����е�cluster����ɨ��
//		�ڶ���: ��ÿ��cluster��ʼ��һ����Ӧ�洢�Ľڵ�����������
//		�ڶ���: ÿ��cluster��Ӧ�Ľڵ�������Ϊ0
		for (i = 0; i < nClusters; i++) {
			nodePerCluster[i] = new int[nNodesPerCluster[i]];
			nNodesPerCluster[i] = 0;
		}
//		������: ɨ�����еĽڵ�, ��������ÿ��cluster�к��еĽڵ��index
		for (i = 0; i < nNodes; i++) {
			j = cluster[i]; // node_i ��Ӧ�� cluster indexΪj
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