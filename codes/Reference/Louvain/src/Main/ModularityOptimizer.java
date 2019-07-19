package Main;

import java.io.BufferedReader;  
import java.io.BufferedWriter;  
import java.io.FileReader;  
import java.io.FileWriter;  
import java.io.IOException;  
import java.util.Arrays;  
import java.util.Random;  
  
public class ModularityOptimizer  
{  
    public static void main(String[] args) throws IOException  
    {  
        boolean  update;  
         
        double modularity, maxModularity, resolution, resolution2;  
        int algorithm, i, j, modularityFunction, nClusters, nIterations, nRandomStarts;  
        int[] cluster;  
        double beginTime, endTime;  
        Network network;  
        Random random;  
        String inputFileName, outputFileName;  
  
          
        inputFileName = "/home/eason/Louvain/facebook_combined.txt";  
        outputFileName = "/home/eason/Louvain/answer.txt";  
        modularityFunction = 1;  
        resolution = 1.0;  
        algorithm = 1;  
        nRandomStarts = 10;  
        nIterations = 3;  
                
        System.out.println("Modularity Optimizer version 1.2.0 by Ludo Waltman and Nees Jan van Eck");  
        System.out.println();  
                              
        System.out.println("Reading input file...");  
        System.out.println();  
          
        network = readInputFile(inputFileName, modularityFunction);  
        
        System.out.format("Number of nodes: %d%n", network.getNNodes());  
        System.out.format("Number of edges: %d%n", network.getNEdges() / 2);  
        System.out.println();  
        System.out.println("Running " + ((algorithm == 1) ? "Louvain algorithm" : ((algorithm == 2) ? "Louvain algorithm with multilevel refinement" : "smart local moving algorithm")) + "...");  
        System.out.println();  
          
        resolution2 = ((modularityFunction == 1) ? (resolution / network.getTotalEdgeWeight()) : resolution);  
  
        beginTime = System.currentTimeMillis();  
        cluster = null;  
        nClusters = -1;  
        maxModularity = Double.NEGATIVE_INFINITY;  
        random = new Random(100);  
        for (i = 0; i < nRandomStarts; i++)  
        {  
            if (nRandomStarts > 1)  
                System.out.format("Random start: %d%n", i + 1);  
  
            network.initSingletonClusters();       //网络初始化，每个节点一个簇  
  
            j = 0;  
            update = true;  
            do  
            {  
                if (nIterations > 1)  
                    System.out.format("Iteration: %d%n", j + 1);  
  
                if (algorithm == 1)  
                    update = network.runLouvainAlgorithm(resolution2, random);  
           
                j++;  
  
                modularity = network.calcQualityFunction(resolution2);  
  
                if (nIterations > 1)  
                    System.out.format("Modularity: %.4f%n", modularity);  
            }  
            while ((j < nIterations) && update);  
            if (modularity > maxModularity)     {  
             
                cluster = network.getClusters();  
                nClusters = network.getNClusters();  
                maxModularity = modularity;  
            }  
  
            if (nRandomStarts > 1)  
            {  
                if (nIterations == 1)  
                    System.out.format("Modularity: %.8f%n", modularity);  
                System.out.println();  
            }  
        }  
        endTime = System.currentTimeMillis();  
  
          
        if (nRandomStarts == 1)  
        {  
            if (nIterations > 1)  
                System.out.println();  
            System.out.format("Modularity: %.8f%n", maxModularity);  
        }  
        else  
            System.out.format("Maximum modularity in %d random starts: %f%n", nRandomStarts, maxModularity);  
        System.out.format("Number of communities: %d%n", nClusters);  
        System.out.format("Elapsed time: %.4f seconds%n", (endTime - beginTime) / 1000.0);  
        System.out.println();  
        System.out.println("Writing output file...");  
        System.out.println();  
         
  
        writeOutputFile(outputFileName, cluster);  
    }  
  
    private static Network readInputFile(String fileName, int modularityFunction) throws IOException  
    {  
        BufferedReader bufferedReader;  
        double[] edgeWeight1, edgeWeight2, nodeWeight;  
        int i, j, nEdges, nLines, nNodes;  
        int[] firstNeighborIndex, neighbor, nNeighbors, node1, node2;  
        Network network;  
        String[] splittedLine;  
  
        bufferedReader = new BufferedReader(new FileReader(fileName));  
  
        nLines = 0;  
        while (bufferedReader.readLine() != null)  
            nLines++;  
  
        bufferedReader.close();  
  
        bufferedReader = new BufferedReader(new FileReader(fileName));  
  
        node1 = new int[nLines];  
        node2 = new int[nLines];  
        edgeWeight1 = new double[nLines];  
        i = -1;  
        bufferedReader.readLine() ;
        for (j = 1; j < nLines; j++)  
        {  
            splittedLine = bufferedReader.readLine().split(" ");  
            node1[j] = Integer.parseInt(splittedLine[0]);  
            if (node1[j] > i)  
                i = node1[j];  
            node2[j] = Integer.parseInt(splittedLine[1]);  
            if (node2[j] > i)  
                i = node2[j];  
            edgeWeight1[j] = (splittedLine.length > 2) ? Double.parseDouble(splittedLine[2]) : 1;  
        }  
        nNodes = i + 1;  
  
        bufferedReader.close();  
  
        nNeighbors = new int[nNodes];  
        for (i = 0; i < nLines; i++)  
            if (node1[i] < node2[i])  
            {  
                nNeighbors[node1[i]]++;  
                nNeighbors[node2[i]]++;  
            }  
  
        firstNeighborIndex = new int[nNodes + 1];  
        nEdges = 0;  
        for (i = 0; i < nNodes; i++)  
        {  
            firstNeighborIndex[i] = nEdges;  
            nEdges += nNeighbors[i];  
        }  
        firstNeighborIndex[nNodes] = nEdges;  
  
        neighbor = new int[nEdges];  
        edgeWeight2 = new double[nEdges];  
        Arrays.fill(nNeighbors, 0);  
        for (i = 0; i < nLines; i++)  
            if (node1[i] < node2[i])  
            {  
                j = firstNeighborIndex[node1[i]] + nNeighbors[node1[i]];  
                neighbor[j] = node2[i];  
                edgeWeight2[j] = edgeWeight1[i];  
                nNeighbors[node1[i]]++;  
                j = firstNeighborIndex[node2[i]] + nNeighbors[node2[i]];  
                neighbor[j] = node1[i];  
                edgeWeight2[j] = edgeWeight1[i];  
                nNeighbors[node2[i]]++;  
            }  
  
     
        {  
            nodeWeight = new double[nNodes];  
            for (i = 0; i < nEdges; i++)  
                nodeWeight[neighbor[i]] += edgeWeight2[i];  
            network = new Network(nNodes, firstNeighborIndex, neighbor, edgeWeight2, nodeWeight);  
        }  
        
  
        return network;  
    }  
  
    private static void writeOutputFile(String fileName, int[] cluster) throws IOException  
    {  
        BufferedWriter bufferedWriter;  
        int i;  
  
        bufferedWriter = new BufferedWriter(new FileWriter(fileName));  
  
        for (i = 0; i < cluster.length; i++)  
        {  
            bufferedWriter.write(Integer.toString(cluster[i]));  
            bufferedWriter.newLine();  
        }  
  
        bufferedWriter.close();  
    }  
}  
