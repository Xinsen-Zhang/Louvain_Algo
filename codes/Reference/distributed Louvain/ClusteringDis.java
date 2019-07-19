package com.graphsee.cypher.alao.distribute.louvain;

import java.util.Arrays;

import org.cltech.crayfish.Crayfish;
import org.cltech.crayfish.CrayfishCache;
import org.cltech.crayfish.CrayfishDataStreamer;
import org.cltech.crayfish.cache.CacheMode;
import org.cltech.crayfish.cache.CacheWriteSynchronizationMode;
import org.cltech.crayfish.configuration.CacheConfiguration;

import com.graphsee.v5.core.CrayfishManager;

public class ClusteringDis {
	protected int nNodes;
	protected int nClusters;
	protected String strClustersName;
	protected String STRCLUSTERS = "Cluster";
	private static Crayfish crayfish = CrayfishManager.Instance.getCrayfish();

//	public ClusteringDis(int nNodes) {
//		this.nNodes = nNodes;
//		cluster = new int[nNodes];
//		nClusters = 1;
//	}
	
	//修改Cluster
	public ClusteringDis(int nNodes) {
		this.nNodes = nNodes;
		strClustersName = STRCLUSTERS + nNodes;
		nClusters = 1;
	}
	
	public void initSingletonClustersV2() {
		int i;
		CacheConfiguration<Integer, Integer> cacheConfig = new CacheConfiguration<>();
		cacheConfig.setName(strClustersName);
		cacheConfig.setCacheMode(CacheMode.REPLICATED);
		cacheConfig.setBackups(0);
		cacheConfig.setWriteSynchronizationMode(CacheWriteSynchronizationMode.FULL_SYNC);
		crayfish.getOrCreateCache(cacheConfig);
		CrayfishDataStreamer<Integer, Integer> clusterStreamer = crayfish.dataStreamer(strClustersName);
		for (i = 0; i < nNodes; i++)
			clusterStreamer.addData(i, i);
		
		clusterStreamer.flush();
		clusterStreamer.close();
		
		nClusters = nNodes;
	}

	public int getNNodes() {
		return nNodes;
	}

	public void orderClustersByNNodes() {
		class ClusterNNodes implements Comparable<ClusterNNodes> {
			public int cluster;
			public int nNodes;

			public ClusterNNodes(int cluster, int nNodes) {
				this.cluster = cluster;
				this.nNodes = nNodes;
			}

			@Override
			public int compareTo(ClusterNNodes clusterNNodes) {
				return (clusterNNodes.nNodes > nNodes) ? 1 : ((clusterNNodes.nNodes < nNodes) ? -1 : 0);
			}
		}

		ClusterNNodes[] clusterNNodes;
		int i;
		int[] newCluster, nNodesPerCluster;

		nNodesPerCluster = getNNodesPerCluster();
		clusterNNodes = new ClusterNNodes[nClusters];
		for (i = 0; i < nClusters; i++)
			clusterNNodes[i] = new ClusterNNodes(i, nNodesPerCluster[i]);

		Arrays.sort(clusterNNodes);

		newCluster = new int[nClusters];
		i = 0;
		do {
			newCluster[clusterNNodes[i].cluster] = i;
			i++;
		} while ((i < nClusters) && (clusterNNodes[i].nNodes > 0));
		nClusters = i;
		
		CrayfishCache<Integer, Integer> clusterLocal = crayfish.cache(strClustersName);
		for (i = 0; i < nNodes; i++) {
			clusterLocal.put(i, newCluster[clusterLocal.localPeek(i)]);
		}
	}

	public int getCluster(int node) {
		return (int) crayfish.cache(strClustersName).localPeek(node);
	}

	public int[] getNNodesPerCluster() {
		int i;
		int[] nNodesPerCluster;
		
		CrayfishCache<Integer, Integer> clusterLocal = crayfish.cache(strClustersName);
		
		nNodesPerCluster = new int[nClusters];
		for (i = 0; i < nNodes; i++)
			nNodesPerCluster[clusterLocal.localPeek(i)]++;
		return nNodesPerCluster;
	}

//	public void mergeClusters(ClusteringDis clustering) {
//		int i;
//		for (i = 0; i < nNodes; i++)
//			cluster[i] = clustering.cluster[cluster[i]];
//		nClusters = clustering.nClusters;
//	}
	
	public void mergeClustersV2(ClusteringDis clustering) {
		int i;
		CrayfishDataStreamer<Integer, Integer> dataStreamer = crayfish.dataStreamer(strClustersName);
		dataStreamer.allowOverwrite(true);
		
		CrayfishCache<Integer, Integer> valueClusterNew = crayfish.cache(strClustersName);
		CrayfishCache<Integer, Integer> valueCluster = crayfish.cache(clustering.strClustersName);
		
		System.out.println("拼合簇 --->");
		for (i = 0; i < nNodes; i++) {
			dataStreamer.addData(i, valueCluster.localPeek(valueClusterNew.localPeek(i)));
		}
		
		dataStreamer.flush();
		dataStreamer.close();
		nClusters = clustering.nClusters;
	}
	
}
