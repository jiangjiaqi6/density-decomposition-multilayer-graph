/*
* File: Graph.h
* Description: Defines the Graph class, which serves as a foundational component
*              for representing a single layer within a multilayer graph structure.
* This file is part of the multilayer graph project.
* Author: Jiaqi Jiang
* Created Date: 2025-01-14
* Last Modified: []
*/

#ifndef GRAPH_H
#define GRAPH_H

#include "../utilities/Defines.h"
#include "../utilities/Utility.h"
#include "../utilities/Time.h"
#include "../utilities/Tools.h"
#include "../utilities/Log.h"


typedef struct edge{
	ui src{0};
	ui end{0};
	edge() = default;
    edge(ui s_, ui e_):src(s_),end(e_){};
	bool operator<(const edge &e) const {
        return src < e.src || (src == e.src && end < e.end);
    }
	// 定义等于运算符，用于去重
    bool operator==(const edge &e) const {
        return src == e.src && end == e.end;
    }
}edge;

typedef struct di_edge{
	ui src{0};
	ui end{0};
	ui to{0};
	di_edge() = default;
    di_edge(ui s_, ui e_):src(s_),end(e_){};

}di_edge;


typedef struct Set
{
	ui* nodes;
	bool* in;
	ui size, E_number, capacity;
	bool have_found;
	Set() {}
	Set(ui sz) { size = 0; capacity = sz; nodes = (ui*)malloc(sz * sizeof(ui)); in = (bool*)malloc(sz * sizeof(bool)); memset(in, 0, sz * sizeof(bool)); E_number = 0; }
	void reset(ui sz) { size = 0; capacity = sz; nodes = (ui*)malloc(sz * sizeof(ui)); in = (bool*)malloc(sz * sizeof(bool)); memset(in, 0, sz * sizeof(bool)); E_number = 0; }
	void free_memory() { free(nodes), free(in); }
	void push(ui x) { nodes[size++] = x; in[x] = true; }
	void refresh() { size = 0; memset(in, 0, capacity * sizeof(bool)); E_number = 0;}
	void copy(Set tmp){
		size = tmp.size;
		capacity = tmp.capacity;
		E_number = tmp.E_number;
		have_found = tmp.have_found;
		memcpy(nodes,tmp.nodes,sizeof(ui)*capacity);
		memcpy(in,tmp.in,sizeof(bool)*capacity);
	}
}Set;

class Graph {
private:
	ui n; //#nodes of the graph
	ui m; //#edges of the graph
	ui maxDeg; //max degree

	ui *pstart; //start positions of neighbors of vertices in the array "edges"
	ui *pend; //end positions of neighbors of vertices in the array "edges"
	ui *edgeList; //concatenation of neighbors of all vertices
	ui *eid; // id of each edge
	ui *deg; // the out degree of each vertex
	ui *abs_deg; // the degree of each vertex
	di_edge *edges;

	string name;
	std::ofstream fout;
	Clock graph_clock;

	//output file
	// std::vector<std::pair<ui,double>> out_res;  
	std::vector<std::pair<ui,ui>> out_res;

public:
	Graph();
	~Graph() ;
	void BuildFromEdgeLst(vector<edge>& edge_buf, uint num_of_vtx, uint num_of_edge, string data_name);
	ui core_decomposition(ui* core);
	void extract_subgraph_from_original();
	ui * get_pstart(){
		return pstart;
	}
	ui * get_pend(){
		return pend;
	}
	ui * get_eid(){
		return eid;
	}
	ui * get_deg(){
		return deg;
	}
	ui * get_abs_deg(){
		return abs_deg;
	}
	ui * get_edgeList(){
		return edgeList;
	}
	di_edge * get_edges(){
		return edges;
	}
	string get_data_name(){
		return name;
	}
	ui get_n(){
		return n;
	}
	ui get_m(){
		return m;
	}
	
private:
	
	void greedy_and_reduce(ui *peel_sequence, ui *core, ui mid = 0);
	void greedy_and_reduce_plus(ui *peel_sequence, ui *core, ui input_mid = 0);

};

#endif /* GRAPH_H_ */