package com.graphsee.cypher.alao.distribute.louvain;

public class LouvainData {
	private int elementId;
	private double value;
	
	public LouvainData(int elementId, double value) {
		// TODO Auto-generated constructor stub
		this.elementId = elementId;
		this.value = value;
	}

	public int getElementId() {
		return elementId;
	}

	public void setElementId(int elementId) {
		this.elementId = elementId;
	}

	public double getValue() {
		return value;
	}

	public void setValue(double value) {
		this.value = value;
	}

}
