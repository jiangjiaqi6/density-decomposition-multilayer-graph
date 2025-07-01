/*
* File: Graph.h
* Description: Implements the Graph class, which provides core functionalities.
*              This implementation file defines methods declared in Graph.h and handles
*              many operations on graph. 
* This file is part of the multilayer graph project.
* Author: Jiaqi Jiang
* Created Date: 2025-01-14
*/

#include "Graph.h"

Graph::Graph():graph_clock("graph"){
    n = m = maxDeg = 0;
    pstart = pend = edgeList = eid = deg = nullptr;
    edges = nullptr;
}

Graph::~Graph() {
    delete[] pstart;
    delete[] pend;
    delete[] edgeList;
    delete[] eid;
    delete[] deg;
	delete[] abs_deg;
    delete[] edges;
}

// Using CSR to store graph.
// void Graph::BuildFromEdgeLst(edge *edge_buf, uint num_of_vtx, uint num_of_edge) {
//     log_info(graph_clock.Count("num_of_vtx: %u, num_of_edge: %u",num_of_vtx,num_of_edge));

//     uint i, j, pj;

//     n = num_of_vtx;
//     pstart = new ui[n]();
//     pend = new ui[n]();
//     deg = new ui[n]();
//     edgeList = new ui[num_of_edge]();
//     eid = new ui[num_of_edge]();
//     edges = new di_edge[num_of_edge/2];

//     std::sort(edge_buf, edge_buf + num_of_edge);

//     pstart[0] = 0;
// 	ui idx = 0;
// 	for(ui i = 0; i < n; i++) {
// 		pend[i] = pstart[i];
//         if(idx < num_of_edge && edge_buf[idx].src == i){
//             edgeList[pend[i]++] = edge_buf[idx++].end;
//             while(idx < num_of_edge && edge_buf[idx].src == i) {
//                 if(edge_buf[idx].end == edge_buf[idx-1].end) idx++;
//                 else
//                     edgeList[pend[i] ++] = edge_buf[idx++].end;
//             }
//         }
// 		if(i < n-1)
// 			pstart[i+1] = pend[i];
// 		maxDeg = (pend[i] - pstart[i]) > maxDeg ? (pend[i] - pstart[i]) : maxDeg;
// 	}

//     for(ui i = 0; i < n; i++)
//         m += (pend[i] - pstart[i]);
//     log_info(graph_clock.Count("m: %u",m));

// }


void Graph::BuildFromEdgeLst(vector<edge> &edge_buf, uint num_of_vtx, uint num_of_edge, string data_name) {
    name = data_name;
    uint i, j, pj;
    n = num_of_vtx;
    pstart = new ui[n]();
    pend = new ui[n]();
    deg = new ui[n]();
	abs_deg = new ui[n]();
    edgeList = new ui[num_of_edge]();
    eid = new ui[num_of_edge]();
    edges = new di_edge[num_of_edge/2];

    sort(edge_buf.begin(), edge_buf.end());
    edge_buf.erase(unique(edge_buf.begin(), edge_buf.end()), edge_buf.end());

    num_of_edge = edge_buf.size();

    pstart[0] = 0;
	ui idx = 0;
	for(ui i = 0; i < n; i++) {
		pend[i] = pstart[i];
        while(idx < num_of_edge && edge_buf[idx].src == i) {
            edgeList[pend[i] ++] = edge_buf[idx++].end;
        }
		if(i < n-1)
			pstart[i+1] = pend[i];
		maxDeg = (pend[i] - pstart[i]) > maxDeg ? (pend[i] - pstart[i]) : maxDeg;
        deg[i] = pend[i] - pstart[i];
	}

    log_info(graph_clock.Count("num_of_vtx: %u, num_of_edge: %u",num_of_vtx,num_of_edge/2));


    // set eid

    idx = 0;
	for(ui i = 0; i < n; i++){
		for(ui j = pstart[i]; j < pend[i]; j++){
			ui v = edgeList[j];
			if(v < i) continue;
			eid[j] = idx;
			eid[pstart[v]++] = idx;
			idx++;
		}
	}
	for(ui i = 0; i < n; i++){
		pstart[i] = pend[i] - deg[i];
		deg[i] = 0;
	}

    for(ui i = 0; i < n; i++){
        m += (pend[i] - pstart[i]);
    }
}

//
ui Graph::core_decomposition(ui *core) {
	ui *degree = new ui[n];
	for(ui i = 0;i < n;i ++) degree[i] = pend[i]-pstart[i];

	ui *rid = new ui[n]();
	ui *id = new ui[n]();
    // ui *core = new ui[n]();
	memset(id, 0, sizeof(ui)*n);
	memset(deg, 0, sizeof(ui)*n);
	for(ui i = 0;i < n;i ++) ++ id[degree[i]];
	for(ui i = 1;i < n;i ++) id[i] += id[i-1];

	for(ui i = 0;i < n;i ++) rid[i] = -- id[degree[i]];
	for(ui i = 0;i < n;i ++) id[rid[i]] = i;

	ui *degree_start = new ui[n+1];
	for(ui i = 0, j = 0;i <= n;i ++) {
		while(j < n&&degree[id[j]] < i) ++ j;
		degree_start[i] = j;
	}

	ui max_core = 0, min_core = 0;
	for(ui i = 0;i < n;i ++) {
		ui u = id[i];
		// assert(degree_start[degree[u]] == i);
		if(degree[u] > max_core) {
			max_core = degree[u];
		}
		core[u] = max_core;
		++ degree_start[degree[u]];
		if(degree[u] == 0){
			// printf("u: %u\n",u);
			continue;
		} 

		degree_start[degree[u]-1] = degree_start[degree[u]];
		for(ui j = pstart[u];j < pend[u];j ++) if(rid[edgeList[j]] > i) {
			ui v = edgeList[j];
			ui pos1 = degree_start[degree[v]], pos2 = rid[v];
			std::swap(id[pos1], id[pos2]);
			rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
			assert(degree[v]>=degree[u]);
			++ degree_start[degree[v]];
			-- degree[v];
			deg[u]++;
			edges[eid[j]].src = u;
			edges[eid[j]].end = v;
			edges[eid[j]].to = v;
		}
	}

	delete[] degree;
	delete[] degree_start;
	delete[] rid;

    delete[] id;
	
    // delete[] core;
	// log_info(graphClock_.Count("core decomposition, max_core: %u, min_core: %u",max_core,min_core));

	return max_core;
}