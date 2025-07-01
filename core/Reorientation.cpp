/*
* File: Reorientation.cpp
* Description: 
    Implements the methods and algorithms defined in Reorientation.h to execute 
    network flow computations on reoriented networks. The implementation focuses 
    on identifying dense subgraphs with integer density k, leveraging advanced flow
    techniques in a single layer graph.
* This file is part of the multilayer graph project.
* Author: Jiaqi Jiang
* Created Date: 2025-01-17
* Last Modified: []
*/

#include "Reorientation.h"

ReorienNetwork::ReorienNetwork(){
    pa = max_pa = app_pa = test_value = 0;
	iter_time = 0;
	network_flow_time = 0;
	update_graph_time = 0;
	reconfig_info_time = 0;
	peel_core_time = 0;

    pstart = nullptr;
    pend = nullptr;
    edgeList = nullptr;
    deg = nullptr;
	abs_deg = nullptr;
    eid = nullptr;
    edges = nullptr;
    
}
ReorienNetwork::~ReorienNetwork(){
    delete[] p;
    delete[] dist;
    delete[] cur;
	delete[] high_outdegree;
    delete[] cur_vis;
	delete[] cur_vis_copy;
    delete[] r;
	delete[] core;
	// delete[] R;
}

void ReorienNetwork::use_graph(bool *vis_info){
	memcpy(cur_vis,vis_info,sizeof(bool)*n);
}

void ReorienNetwork::check_density(ui integer_threshold){
	ui tmp_m = 0;
	ui tmp_v = 0;
	for(ui i = 0; i < n; i++){
		ui u = i;
		if(cur_vis[u]) continue;
		tmp_v++;
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(cur_vis[v] || u > v) continue;
			tmp_m++;
		}
	}
	if(tmp_v == 0)
		printf("error: vetex zero\n");
	else
		printf("integer_threshold: %u, density: %lf\n",integer_threshold,1.0*tmp_m/tmp_v);
}

// if the minimum deg of subgraph less k, it must be computed;
// else it is not be computed
bool ReorienNetwork::compute_lower_bound(ui k, vector<ui> veck)
{
	adjust_deg();
	ui min_out_deg = n;
	ui c_n = 0, c_m = 0;
	unordered_set<ui> deg_less_k;
	unordered_set<ui> to_remove;
	for(ui i = 0; i < n; i++){
		if(cur_vis[i]) continue;
		min_out_deg = std::min(min_out_deg,deg[i]);
		c_n++;
		c_m += deg[i];
		if(deg[i] < k) deg_less_k.insert(i);
	}
	ui inter_den = std::ceil(c_m*1.0/c_n);
	// if(veck[0] == 0 && veck[1] == 0 && veck[2] == 1 && veck[3] == 2 && veck[4] == 1){
	// 	printf("k: %u, deg_less_k size: %u, inter_den: %u, min_out_deg: %u\n",k,deg_less_k.size(),inter_den,min_out_deg);
	// }
	if(min_out_deg >= k) return true;

	// if(veck[0] == 1 && veck[1] == 0 && veck[2] == 0 && veck[3] == 0 && veck[4] == 2){
	// 	printf("k: %u, deg_less_k size: %u, inter_den: %u\n",k,deg_less_k.size(),inter_den);
	// 	for(auto u:deg_less_k)
	// 		printf("u: %u, deg: %u\n",u,deg[u]);
	// 	string ss = "test_out.csv";
	// 	auto out = ofstream(ss);
	// 	for(ui u = 0; u < n; u++){
	// 		if(cur_vis[u]) continue;
	// 		for(ui j = pstart[u]; j < pend[u]; j++){
	// 			ui v = edgeList[j];
	// 			if(cur_vis[v] || v < u) continue;
	// 			di_edge die = edges[eid[j]];
	// 			ui from = die.to == die.src ? die.end : die.src;
	// 			printf("from: %u, to: %u\n",from,die.to);
	// 			out<<from<<","<<die.to<<endl;
	// 		}
	// 	}
	// 	out.close();
	// }

	for(auto u : deg_less_k){
		queue<ui> nbr;
		unordered_map<ui,ui> path;
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(cur_vis[v]) continue;
			di_edge die = edges[eid[j]];
			if(die.to == u){
				nbr.push(v);
				path[v] = eid[j];
			}
		}
		if(nbr.size() == 0) return false;
		bool flag = true;

		while(!nbr.empty()){
			ui uu = nbr.front();
			nbr.pop();
			if(deg[uu] > k){
				
				while(uu != u){
					deg[uu]--;
					di_edge& die = edges[path[uu]];
					ui from = die.to == die.end ? die.src : die.end;

					ui to = die.to;
					deg[to]++;
					uu = to;
					die.to = from;
				}
				if(uu == u) {
					flag = false;
					break;
				}
			}
			else if(deg[uu] == k){
				flag = false;
				break;
			}
			for(ui j = pstart[uu]; j < pend[uu]; j++){
				ui v = edgeList[j];
				if(cur_vis[v] || path.find(v) != path.end()) continue;
				di_edge die = edges[eid[j]];
				if(die.to == uu){
					nbr.push(v);
					path[v] = eid[j];
				}
			}
		}
		if(!flag){
			to_remove.insert(u);  
		}
	}

	for(auto u : to_remove) {
		deg_less_k.erase(u);
	}
	
	
		
	if(deg_less_k.empty()) return true;
	else{
		return false;
	}
		
}

ui ReorienNetwork::core_decomposition_incremental(){
	ui *degree = new ui[n]();
	memset(core,0,sizeof(ui)*n);
	for(ui u = 0;u < n;u ++){
		if(cur_vis[u]){
			continue;
		}
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(cur_vis[v] || v < u) continue;
			degree[u]++;
			degree[v]++;
		}
	}

	ui *rid = new ui[n]();
	ui *id = new ui[n]();
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
			assert(cur_vis[v] == false);
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
	
	// log_info(graphClock_.Count("core decomposition, max_core: %u, min_core: %u",max_core,min_core));

	return max_core;
}

void ReorienNetwork::initial(MultilayerGraph & mg, ui dim_){
    dim = dim_;
    n = mg.GetN();
    p = new ui[n]();
    dist = new ui[n]();
    cur = new ui[n]();
	high_outdegree = new ui[n]();
    cur_vis = new bool[n]();
	cur_vis_copy = new bool[n]();
    r = new ui[n]();
	core = new ui[n]();

    auto &g = mg.GetGraph(dim);
    ui maxcore = g.core_decomposition(core);
	max_pa = maxcore;
	string ss = mg.GetDataPath();
	if (ss.back() == '/') {
        ss = ss.substr(0, ss.length() - 1);
    }
	size_t pos = ss.find_last_of('/');
    ss = ss.substr(pos + 1);
    cout<< "name: " << ss << ", node: " << g.get_n() << ", edge: " << g.get_m()/2 <<", layer: " << dim << ", maxcore: " << maxcore << endl;
	set_graph_name(ss);
    m = g.get_m();
    pstart = g.get_pstart();
    pend = g.get_pend();
    edgeList = g.get_edgeList();
    deg = g.get_deg();
	abs_deg = g.get_abs_deg();
    eid = g.get_eid();
    edges = g.get_edges();
}


void ReorienNetwork::dispaly_density_every_layer(string ss){
	size_t last_slash_pos = ss.find_last_of('/');
    size_t last_dot_pos = ss.find_last_of('.');
    std::string path = ss.substr(0, last_slash_pos + 1); // 包含最后的 '/'
    std::string file_name = ss.substr(last_slash_pos + 1, last_dot_pos - last_slash_pos - 1); // 不包含扩展名

    ss = path + file_name + "_density_result.txt";
    // cout<<ss<<endl;

    auto out = ofstream(ss);
    ui *id = new ui[n]();
	for(ui i = 0; i < n; i++)id[i] = i;
	std::sort(id,id+n,[&](ui t1, ui t2){return r[t1] < r[t2];});
    for(ui u = 0; u < n; u++){
        ui i = id[u];
        // if(r[i] == 0) continue;
        out << "r[" << i <<"] = " << r[i] <<endl;
    }
    delete[] id;
    out.close();
}

void ReorienNetwork::compute_each_layer_sub_edge(ReorienNetwork *rn, vector<ui>& veck, ui layer, bool* vis){
	ui n = rn[0].get_n();
    ui* pstart_;
    ui* pend_;
    ui* edgeList_;
	for(ui i = 0; i < layer; i++){
		pstart_ = rn[i].get_pstart();
		pend_ = rn[i].get_pend();
		edgeList_ = rn[i].get_edgeList();
		ui c = 0;
		for(ui j = 0; j < n; j++){
			if(vis[j]) continue;
			for(ui x = pstart_[j]; x < pend_[j]; x++){
				ui v = edgeList_[x];
				if(vis[v]) continue;
				if(v < j) continue;
				c++;
			}
		}
		veck.push_back(c);
	}
}

void ReorienNetwork::peel_k_core_multilayer(ReorienNetwork *rn, vector<ui> veck, ui layer){
	ui n = rn[0].get_n();
    ui* pstart_;
    ui* pend_;
    ui* edgeList_;
    ui* abs_deg_;
	ui* deg_;
    bool* cur_vis_;
    queue<ui> q;
    bool is_change = true;
	vector<ui> del_store;
    bool flag = false;
    // if(veck[0] == 2 && veck[1] == 0 && veck[2] == 0)
    //     flag = true;
    ui debug_c = 0;
    while(is_change){
        is_change = false;
        for(ui i = 0; i < layer; i++){
            ui key = veck[i];
            if(key == 0) continue;
            pstart_ = rn[i].get_pstart();
            pend_ = rn[i].get_pend();
            edgeList_ = rn[i].get_edgeList();
            abs_deg_ = rn[i].get_abs_deg();
            cur_vis_ = rn[i].get_cur_vis();
            
            for(ui j = 0; j < n; j++){
                if(cur_vis_[j]) continue;
                if(abs_deg_[j] < key){
                    q.push(j);
                    cur_vis_[j] = true;
                }
            }
            if(q.empty())continue;
            is_change = true;
            while(!q.empty()){
                ui j = q.front();
                q.pop();
				del_store.push_back(j);
                for(ui x = pstart_[j]; x < pend_[j]; x++){
                    ui v = edgeList_[x];
                    if(cur_vis_[v]) continue;
                    assert(abs_deg_[v]>0);
                    abs_deg_[v]--;
                    if(abs_deg_[v] < key){
                        q.push(v);
                        cur_vis_[v] = true; 
                    }
                }
                for(ui c = 0; c < layer; c++){
                    // if(c != i && veck[c])
					if(c != i)
                        rn[c].peel_k_core_singlelayer(j);
                }   
            }   
        }
    }

	for(ui i = 0; i < n; i++)
		if(!cur_vis[i])
			debug_c++;
	if(debug_c == 0)return;

	ui *eid_;
	di_edge *edges_;
	for(ui i = 0; i < del_store.size(); i++){
		ui u = del_store[i];
		for(ui x = 0; x < layer; x++){
			pstart_ = rn[x].get_pstart();
            pend_ = rn[x].get_pend();
            edgeList_ = rn[x].get_edgeList();
            abs_deg_ = rn[x].get_abs_deg();
			deg_ = rn[x].get_subgraph_d();
            cur_vis_ = rn[x].get_cur_vis();
			eid_ = rn[x].get_eid();
			edges_ = rn[x].get_edges();
			for(ui j = pstart_[u]; j < pend_[u]; j++){
				ui v = edgeList_[j];
				if(cur_vis_[v])continue;
				di_edge& de = edges_[eid_[j]];
				if(de.to == u){
					de.to = v;
					assert(deg_[v]>0);
					deg_[v]--;
					deg_[u]++;
				}
			}
		}
		
	}
}


void ReorienNetwork::update_kdense_subgraph_high_dim_bottom_up(ReorienNetwork *rn, ui index, vector<ui> veck, ui Z_rank,ui &num_of_call, Set& Z)
{
	ui layer = index;
	bool *remove = new bool[n]();
	memcpy(remove,cur_vis,sizeof(bool)*n);
	vector<ui> tmp_veck;
	tmp_veck  = veck;
	tmp_veck.push_back(Z_rank);
	for(ui i = 0; i < layer; i++){
		rn[i].use_graph(remove);
        rn[i].adjust_deg(); 
    }

	timeCount t;
	t.StartTime();
	ui before_c = 0, after_c = 0;
	memset(abs_deg,0,sizeof(ui)*n);
	for(ui i = 0; i < n; i++){
		if(cur_vis[i]) continue;
		before_c++;
		for(ui j = pstart[i]; j < pend[i]; j++){
			ui v = edgeList[j];
			if(cur_vis[v] || v < i) continue;
			abs_deg[i]++;
			abs_deg[v]++;
		}
	}
	assert(tmp_veck.size() == index+1);
	peel_k_core_multilayer(rn,tmp_veck,index+1);
	
	for(ui i = 0; i < n; i++)
		if(!cur_vis[i])
			after_c++;
	if(after_c == 0){
		Z.refresh();
		return;
	}
	t.EndTime();
	peel_core_time += t.QueryTime();

	// printf("before_c: %u, after_c: %u\n",before_c,after_c);
	
	vector<ui> varVertex;
	bool* cur_vis_;
    bool exist = true, tmp_flag = false;
	vector<ui> isContain, change;

	// if(index == 2 && veck[0] == 9 && veck[1] >= 1 && Z_rank == 5)
	// 	tmp_flag = true;
    ui iter = 0, i, c_v = 0, old_z_size = Z.size, isContainOldSize = 0;
	if(tmp_flag ){
		printf("before Z_rank: %u, graph size: %u\n",Z_rank,Z.size);
		check_density(Z_rank);
	}
        

    while(true){
        bool all_unchanged = true;
		iter++;
		isContain.clear();
        for(i = 0; i < layer; i++){      
			num_of_call++;   
			c_v = 0;   
			cur_vis_ = rn[i].get_cur_vis(); 
			rn[i].find_k_dense_sub_plus(isContain,change,veck[i],i);
			Z.refresh(); 
			for(ui j = 0; j < isContain.size(); j++)
				Z.push(isContain[j]);
			
			if(tmp_flag){
				printf("iter: %u, i : %u, graph size: %u, isContainOldSize: %u\n",iter,i,Z.size,isContainOldSize);
				rn[i].check_density(Z_rank);
			}
			
			timeCount t;
			t.StartTime();
			for(ui jj = 0; jj < change.size(); jj++){
				ui j = change[jj];
				remove[j] = true;
				// cur_vis_[j] = true;
				assert(cur_vis_[j] == true);
				high_outdegree[c_v++] = j;
			}

			t.EndTime();
			reconfig_info_time += t.QueryTime();

			for(ui j = 0; j <= index; j++){
				if(j == i) continue;
				if(j < index)
					rn[j].update_graph_incremental(high_outdegree,c_v,Z);
				else{
					update_graph_incremental(high_outdegree,c_v,Z);
				}
			}
			if(Z.size == 0){
                exist = false;
                break;
            }	
            
            if(iter == 1 && i == 0){
                isContainOldSize = Z.size;
            }
            else{
                if(Z.size != isContainOldSize){
                    all_unchanged = false;
                    isContainOldSize = Z.size;
                }     
            }
        }
		if(Z.size == 0){
			exist = false;
			break;
		}
		if(tmp_flag)
			printf("before Z.size: %u, isContainOldSize: %u\n",Z.size,isContainOldSize);
		if(old_z_size != isContainOldSize)
		{
			timeCount t;
			t.StartTime();
			memset(r,0,sizeof(ui)*n);
			Z.refresh();
			if(!ReTestBoUp(Z_rank)){
				for(ui j = 0; j < n; j++)
					if(!cur_vis[j])
						Z.push(j);
			}
			else{
				if(!Z.size){
					exist = false;
					break;
				}
			}

			t.EndTime();
			network_flow_time += t.QueryTime();

			t.StartTime();

			c_v = 0;
			for(ui k = 0; k < n; k++){
				if(remove[k]) continue;
				if(Z.in[k]) continue;
				remove[k] = true;
				// cur_vis[k] = true;
				assert(cur_vis[k] == true);
				high_outdegree[c_v++] = k;
			}
			if(tmp_flag){
				printf("after Z.size: %u, isContainOldSize: %u\n",Z.size,isContainOldSize);
				check_density(Z_rank);
			}
				

			t.EndTime();
			reconfig_info_time += t.QueryTime();

			for(ui j = 0; j < index; j++)
				rn[j].update_graph_incremental(high_outdegree,c_v,Z);

			if(Z.size != isContainOldSize){
				all_unchanged = false;
				isContainOldSize= Z.size;
				old_z_size = Z.size;
			} 
			if(!Z.size){
				exist = false;
				break;
			}
		}
        if(all_unchanged)break;
        if(!exist) break;
    }

	if(tmp_flag)
        printf("index: %u, varVertex size: %u, final graph size: %u\n",index,varVertex.size(),Z.size);


	delete[] remove;
}






void ReorienNetwork::update_kdense_subgraph_high_dim_incremental(ReorienNetwork *rn, ui index, ui &num_of_iter, vector<ui> veck, ui X_rank, ui Y_rank, ui Z_rank, Set& Z)
{
	ui layer = index;
	Set& X = R[X_rank], & Y = R[Y_rank];
	bool *remove = new bool[n]();
	for(ui i = 0; i < n; i++) remove[i] = true;

	for (ui p = 0; p < Z.size; p++){
		remove[Z.nodes[p]] = false;
	}
	memcpy(cur_vis,remove,sizeof(bool)*n);


	
	vector<ui> tmp_veck;
	tmp_veck  = veck;
	tmp_veck.push_back(Z_rank);
	for(ui i = 0; i < layer; i++){
		rn[i].use_graph(remove);
        rn[i].adjust_deg(); 
    }

	timeCount t;
	t.StartTime();
	ui before_c = 0, after_c = 0;
	memset(abs_deg,0,sizeof(ui)*n);
	for(ui i = 0; i < n; i++){
		if(cur_vis[i]) continue;
		before_c++;
		for(ui j = pstart[i]; j < pend[i]; j++){
			ui v = edgeList[j];
			if(cur_vis[v] || v < i) continue;
			abs_deg[i]++;
			abs_deg[v]++;
		}
	}
	assert(tmp_veck.size() == index+1);
	peel_k_core_multilayer(rn,tmp_veck,index+1);
	
	for(ui i = 0; i < n; i++)
		if(!cur_vis[i])
			after_c++;
	if(after_c == 0){
		Z.refresh();
		return;
	}
	t.EndTime();
	peel_core_time += t.QueryTime();
	
	vector<ui> varVertex;
	bool* cur_vis_;


    bool exist = true, tmp_flag = false;
	// if(index == 2 && veck[0] == 1 && veck[1] == 1)
	// 	tmp_flag = true;
    ui iter = 0, i, c_v = 0, old_z_size = Z.size, isContainOldSize = 0;
	if(tmp_flag)
        printf("X_rank: %u, Y_rank: %u, Z_rank: %u, graph size: %u\n",X_rank,Y_rank,Z_rank,Z.size);

    while(true){
        bool all_unchanged = true;
		iter++;
        for(i = 0; i < layer; i++){         
			c_v = 0;   
			cur_vis_ = rn[i].get_cur_vis();  
			num_of_iter++;
			rn[i].find_k_dense_sub_incremental(veck[i],X,X_rank,Y,Y_rank,Z);
			if(tmp_flag)
                printf("iter: %u, i : %u, graph size: %u, isContainOldSize: %u\n",iter,i,Z.size,isContainOldSize);
			
			timeCount t;
			t.StartTime();

			for(ui j = 0; j < n; j++){
				if(remove[j]) continue;
				if(Z.in[j]) continue;
				varVertex.push_back(j);
				remove[j] = true;
				cur_vis_[j] = true;
				high_outdegree[c_v++] = j;
			}

			t.EndTime();
			reconfig_info_time += t.QueryTime();

			for(ui j = 0; j <= index; j++){
				if(j == i) continue;
				if(j < index)
					rn[j].update_graph_incremental(high_outdegree,c_v,Z);
				else{
					update_graph_incremental(high_outdegree,c_v,Z);
				}
			}
			if(Z.size == 0){
                exist = false;
                break;
            }	
            
            if(iter == 1 && i == 0){
                isContainOldSize = Z.size;
            }
            else{
                if(Z.size != isContainOldSize){
                    all_unchanged = false;
                    isContainOldSize = Z.size;
                }     
            }
        }
		if(Z.size == 0){
			exist = false;
			break;
		}
		if(tmp_flag)
			printf("before Z.size: %u, isContainOldSize: %u\n",Z.size,isContainOldSize);
		if(old_z_size != isContainOldSize)
		{
			timeCount t;
			t.StartTime();

			test_value--;
			ui high_num = 0;
			num_of_iter++;
			while (DinicBFSDivAndCon(X, Y, high_num))
			{
				for (ui p = 0; p < Y.size; p++)
					cur[Y.nodes[p]] = pstart[Y.nodes[p]];
				for (ui i = 0; i < high_num; i++)
				{
					p[high_outdegree[i]] = INF + 1;
					DinicDFS(high_outdegree[i]);
				}
			}
			test_value++;
			
			Z.refresh();
			update_density_div_and_con_incremental(X, Y, Z);

			t.EndTime();
			network_flow_time += t.QueryTime();

			t.StartTime();

			c_v = 0;
			for(ui k = 0; k < n; k++){
				if(remove[k]) continue;
				if(Z.in[k]) continue;
				remove[k] = true;
				cur_vis[k] = true;
				high_outdegree[c_v++] = k;
				varVertex.push_back(k);
			}
			if(tmp_flag)
				printf("after Z.size: %u, isContainOldSize: %u\n",Z.size,isContainOldSize);

			t.EndTime();
			reconfig_info_time += t.QueryTime();

			for(ui j = 0; j < index; j++)
				rn[j].update_graph_incremental(high_outdegree,c_v,Z);

			if(Z.size != isContainOldSize){
				all_unchanged = false;
				isContainOldSize= Z.size;
				old_z_size = Z.size;
			} 
			if(!Z.size){
				exist = false;
				break;
			}
		}
	

        if(all_unchanged)break;
        if(!exist) break;
    }

	// if(after_c != before_c)
	// 	printf("before_c: %u, after_c: %u, final graph size: %u, iter: %u\n",before_c,after_c,Z.size,iter);

	if(tmp_flag)
        printf("index: %u, varVertex size: %u, final graph size: %u\n",index,varVertex.size(),Z.size);

	timeCount t1;
	t1.StartTime();
	
	for(ui i = 0; i < varVertex.size(); i++){
		ui u = varVertex[i];
		deg[u] = 0;
		r[u] = Y_rank;
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(cur_vis_copy[v]) continue;
			di_edge& ed = edges[eid[j]];
			
			if(Z.in[v]){
				deg[u]++;
				ed.to = v;
			}
			else{
				if(ed.to == v)
					deg[u]++;
			}
		}
	}

	t1.EndTime();
	reconfig_info_time += t1.QueryTime();

	for(ui j = 0; j < index; j++)
		rn[j].adjust_deg_incremental(varVertex,cur_vis_copy,Z);

	delete[] remove;
}




void ReorienNetwork::update_kdense_subgraph_2dim_incremental(ReorienNetwork *rn, ui index, vector<ui> veck, ui X_rank, ui Y_rank, ui Z_rank, Set& Z)
{
	ui layer = index;
	Set& X = R[X_rank], & Y = R[Y_rank];
	bool *remove = new bool[n]();
	for(ui i = 0; i < n; i++) remove[i] = true;

	for (ui p = 0; p < Z.size; p++){
		remove[Z.nodes[p]] = false;
	}
	memcpy(cur_vis,remove,sizeof(bool)*n);

	for(ui i = 0; i < layer; i++){
		rn[i].use_graph(remove);
        rn[i].adjust_deg(); 
    }

	vector<ui> isContain, isContainOld, change;
	vector<ui> varVertex;

    bool* cur_vis_ = rn[0].get_cur_vis();
    ui* pstart_;
    ui* pend_;
    ui* edgeList_;
    ui* abs_deg_;
    bool exist = true, tmp_flag = false;
    ui iter = 1, i, c_v = 0, isContainSize,isContainOldSize;


    while(true){
        bool all_unchanged = true;
		c_v = 0;

        for(i = 0; i < layer; i++){              
			rn[i].find_k_dense_sub_incremental(veck[i],X,X_rank,Y,Y_rank,Z);
			if(Z.size == 0){
                exist = false;
                break;
            }
			for(ui j = 0; j < n; j++){
				if(remove[j]) continue;
				if(Z.in[j]) continue;
				remove[j] = true;
				cur_vis[j] = true;
				cur_vis_[j] = true;
				high_outdegree[c_v++] = j;
			}

			if(tmp_flag)
                printf("iter: %u, i : %u, graph size: %u\n",iter,i,Z.size);
			
            
            if(iter == 1 && i == 0){
                isContainOldSize = isContainSize;
                all_unchanged = false;
            }
            else{
                if(Z.size != isContainOldSize){
                    all_unchanged = false;
                    isContainOldSize = Z.size;
                }     
            }
        }
		if(tmp_flag)
			printf("before Z.size: %u, Z.num: %u\n",Z.size,Z.E_number);
		bool is_update = false;
		if(c_v) is_update = true;
		for (ui p = 0; p < c_v; p++){
			ui u = high_outdegree[p];
			varVertex.push_back(u);
			// u is in Y
			assert(Y.in[u] == true);
			if(Y.in[u] != true){
				printf("X_rank: %u, Y_rank: %u, Z_rank: %u, u: %u, r[u]: %u, deg[u]: %u\n",X_rank,Y_rank,Z_rank,u,r[u],deg[u]);
			}
			for(ui j = pstart[u]; j < pend[u]; j++){
				ui v = edgeList[j];
				if(Z.in[v]){
					di_edge de = edges[eid[j]];
					ui from = de.to == de.src ? de.end : de.src;
					if(from == v){
						assert(deg[v]>0);
						deg[v]--;
					}
				}
			}
		}
		if(is_update){

			test_value--;
			ui high_num = 0;
			while (DinicBFSDivAndCon(X, Y, high_num))
			{
				for (ui p = 0; p < Y.size; p++)
					cur[Y.nodes[p]] = pstart[Y.nodes[p]];
				for (ui i = 0; i < high_num; i++)
				{
					p[high_outdegree[i]] = INF + 1;
					DinicDFS(high_outdegree[i]);
				}
			}
			test_value++;
			
			Z.refresh();
			update_density_div_and_con_incremental(X, Y, Z);
			c_v = 0;
			for(ui k = 0; k < n; k++){
				if(remove[k]) continue;
				if(Z.in[k]) continue;
				remove[k] = true;
				cur_vis[k] = true;
				cur_vis_[k] = true;
				high_outdegree[c_v++] = k;
				varVertex.push_back(k);
			}
			if(tmp_flag)
				printf("after Z.size: %u, Z.num: %u\n",Z.size,Z.E_number);

			if(!Z.size){
				exist = false;
				break;
			}

			// rn[0].update_graph(remove,veck[0]);
			rn[0].update_graph_incremental(high_outdegree,c_v,Z);

			if(Z.size != isContainOldSize){
				all_unchanged = false;
				isContainOldSize= Z.size;
			} 
		}
	

        if(all_unchanged)break;
        if(!exist) break;
        iter++;
    }

	if(tmp_flag)
        printf("index: %u, varVertex size: %u, Z.size: %u\n",index,varVertex.size(),Z.size);


	for(ui i = 0; i < varVertex.size(); i++){
		ui u = varVertex[i];
		deg[u] = 0;
		r[u] = Y_rank;
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(cur_vis_copy[v]) continue;
			di_edge& ed = edges[eid[j]];
			
			if(Z.in[v]){
				deg[u]++;
				ed.to = v;
			}
			else{
				if(ed.to == v)
					deg[u]++;
			}
		}
	}
	rn[0].adjust_deg_incremental(varVertex,cur_vis_copy,Z);

	delete[] remove;
}

void ReorienNetwork::adjust_deg_incremental(vector<ui> varVertex, bool* vis, Set& Z){
	timeCount t;
	t.StartTime();
	for(ui i = 0; i < Z.size; i++){
		ui u = Z.nodes[i];
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(vis[v]) continue;
			if(Z.in[v]) continue;
			di_edge& ed = edges[eid[j]];
			ed.to = u;
		}
	}
	// for(ui i = 0; i < varVertex.size(); i++){
	// 	ui u = varVertex[i];
	// 	// deg[u] = 0;
	// 	for(ui j = pstart[u]; j < pend[u]; j++){
	// 		ui v = edgeList[j];
	// 		// if(cur_vis_copy[v]) continue;
	// 		di_edge& ed = edges[eid[j]];
			
	// 		if(Z.in[v]){
	// 			// deg[u]++;
	// 			ed.to = v;
	// 		}
	// 		// else{
	// 		// 	if(ed.to == v)
	// 		// 		deg[u]++;
	// 		// }
	// 	}
	// }
	t.EndTime();
	update_graph_time += t.QueryTime();
}



void ReorienNetwork::update_kdense_subgraph_incremental(ReorienNetwork *rn, ui index, vector<ui> veck, ui X_rank, ui Y_rank, ui Z_rank, Set& Z)
{
	ui layer = index;
	Set& X = R[X_rank], & Y = R[Y_rank];
	bool *remove = new bool[n]();
	for(ui i = 0; i < n; i++) remove[i] = true;

	for (ui p = 0; p < Z.size; p++){
		remove[Z.nodes[p]] = false;
	}

	for(ui i = 0; i < layer; i++){
		rn[i].use_graph(remove);
        rn[i].adjust_deg(); 
    }

	vector<ui> isContain, isContainOld, change;
	vector<ui> varVertex;

    bool* cur_vis_;
    ui* pstart_;
    ui* pend_;
    ui* edgeList_;
    ui* abs_deg_;
    bool exist = true, tmp_flag = false;
    ui iter = 1, i, c_v = 0, isContainSize,isContainOldSize;


    while(true){
        bool all_unchanged = true;
        isContain.clear();
        for(i = 0; i < layer; i++){              
            // rn[i].find_k_dense_sub(isContain,veck[i],i);   // Optimize according to the order of element values in the veck vector
            rn[i].find_k_dense_sub_plus(isContain,change,veck[i],i); 

            // if(tmp_flag)
            //     printf("iter: %u, i : %u, before graph size: %u, change size: %u\n",iter,i,isContain.size(),change.size());
            // c_v = 0;
            // if(change.size() && isContain.size()){
            //     queue<ui> Q;
            //     // ui* abs_deg = rn[i].get_abs_deg();
            //     // bool* cur_vis = rn[i].get_cur_vis();
            //     for(ui j = 0; j < change.size(); j++){
            //         ui u = change[j];
            //         Q.push(u);
            //         // remove[u] = true;
            //     }
            //     while(!Q.empty()){
            //         ui u = Q.front();
            //         Q.pop();
            //         if(remove[u]) continue;
            //         c_v++;
            //         remove[u] = true;
            //         for(ui j = 0; j < layer; j++){
            //             abs_deg_ = rn[j].get_abs_deg();
            //             cur_vis_ = rn[j].get_cur_vis();
            //             pstart_ = rn[j].get_pstart();
            //             pend_ = rn[j].get_pend();
            //             edgeList_ = rn[j].get_edgeList();
            //             for(ui k = pstart_[u]; k < pend_[u]; k++){
            //                 ui v = edgeList_[k];       
            //                 if(remove[v])continue;
            //                 assert(abs_deg_[v]<n);
            //                 abs_deg_[v]--;
            //                 if(abs_deg_[v] < veck[j]){
            //                     Q.push(v);
            //                     // remove[v] = true;
            //                 }
            //             }
            //         }
            //     }
            //     isContain.clear();
            //     bool no_subgraph = true;
            //     for(ui j = 0; j < n; j++)
            //         if(!remove[j])
            //         {
            //             no_subgraph = false;
            //             break;
            //         }
            //     if(no_subgraph){
            //         exist = false;
            //         break;
            //     }
            //     memcpy(rn[i].get_cur_vis(),remove,n*sizeof(bool));
            //     for(ui ii = 0; ii < n; ii++)
            //         if(!remove[ii])
            //             isContain.push_back(ii);
            // }
            
			for(ui j = 0; j < change.size(); j++)
				remove[change[j]] = true;

			if(tmp_flag)
                printf("iter: %u, i : %u, graph size: %u, c_v: %u\n",iter,i,isContain.size(),c_v);

            if(!isContain.size()){
                exist = false;
                break;
            }
			
			if(i+1 != layer)
			{
				rn[(i+1)%layer].update_graph(rn[i].get_cur_vis(),veck[(i+1)%layer]);
			}	

            
            if(iter == 1 && i == 0){
                isContainOldSize = isContainSize;
                all_unchanged = false;
            }
            else{
                if(isContain.size() != isContainOldSize){
                    all_unchanged = false;
                    isContainOldSize = isContain.size() ;
                }     
            }
        }
		if(tmp_flag)
			printf("before Z.size: %u, Z.num: %u\n",Z.size,Z.E_number);
		bool is_update = false;
		for (ui p = 0; p < Z.size; p++){
			ui u = Z.nodes[p];
			if(remove[u])
			{
				is_update = true;
				varVertex.push_back(u);
				Z.in[u] = false;
				cur_vis[u] = true;
				// u is in Y
				assert(Y.in[u] == true);
				if(Y.in[u] != true){
					printf("X_rank: %u, Y_rank: %u, Z_rank: %u, u: %u, r[u]: %u, deg[u]: %u\n",X_rank,Y_rank,Z_rank,u,r[u],deg[u]);
				}
				for(ui j = pstart[u]; j < pend[u]; j++){
					ui v = edgeList[j];
					if(Z.in[v]){
						di_edge de = edges[eid[j]];
						ui from = de.to == de.src ? de.end : de.src;
						if(from == v){
							assert(deg[v]>0);
							deg[v]--;
						}
					}
				}
			}
		}
		memcpy(cur_vis,remove,sizeof(bool)*n);	
		if(is_update){

			test_value--;
			ui high_num = 0;
			while (DinicBFSDivAndCon(X, Y, high_num))
			{
				for (ui p = 0; p < Y.size; p++)
					cur[Y.nodes[p]] = pstart[Y.nodes[p]];
				for (ui i = 0; i < high_num; i++)
				{
					p[high_outdegree[i]] = INF + 1;
					DinicDFS(high_outdegree[i]);
				}
			}
			test_value++;
			
			Z.refresh();
			
			update_density_div_and_con_incremental(X, Y, Z);
			for(ui k = 0; k < n; k++){
				if(remove[k]) continue;
				if(Z.in[k]) continue;
				varVertex.push_back(k);
			}
			if(tmp_flag)
				printf("after Z.size: %u, Z.num: %u\n",Z.size,Z.E_number);

			if(!Z.size){
				exist = false;
				break;
			}
			for(ui k = 0; k < n; k++) remove[k] = true;
			for(ui k = 0; k < Z.size; k++){
				remove[Z.nodes[k]] = false;
			}
			rn[0].update_graph(remove,veck[0]);

			if(Z.size != isContainOldSize){
				all_unchanged = false;
				isContainOldSize = Z.size;
				// isContainOld.clear();
				// for(ui k = 0; k < Z.size; k++){
				// 	isContainOld.push_back(Z.nodes[k]);
				// }
			} 
		}
	

        if(all_unchanged)break;
        if(!exist) break;
        iter++;
    }

	if(tmp_flag)
        printf("index: %u, varVertex size: %u, Z.size: %u\n",index,varVertex.size(),Z.size);


	for(ui i = 0; i < varVertex.size(); i++){
		ui u = varVertex[i];
		deg[u] = 0;
		r[u] = Y_rank;
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(cur_vis_copy[v]) continue;
			di_edge& ed = edges[eid[j]];
			
			if(Z.in[v]){
				deg[u]++;
				ed.to = v;
			}
			else{
				if(ed.to == v)
					deg[u]++;
			}
		}
	}


	delete[] remove;
}

//there is optimization based on Z, not on Y
void ReorienNetwork::update_density_div_and_con_incremental(Set& X, Set& Y, Set& Z){

	queue<ui> Q;
	bool* vis = new bool[n]();
	for (ui p = 0; p < Y.size; p++)
	{
		ui i = Y.nodes[p];
		if (X.in[i] || cur_vis[i]) continue;
		vis[i] = false;
		if (deg[i] >= test_value)
		{
			r[i] = max(r[i], test_value);
			vis[i] = true;
			Q.push(i);
		}
	}

	while(!Q.empty()){
		ui x = Q.front();
		Q.pop();
		Z.push(x);
		for (ui i = pstart[x]; i < pend[x]; i++)
		{
			ui v = edgeList[i];
			
			if (vis[v] || X.in[v] || cur_vis[v])  continue;
			ui id = eid[i];
			di_edge& ne = edges[id];
			if (ne.to == x ) continue;
			ui to = ne.to;
			vis[to] = true;
			Q.push(to);
			r[to] = max(r[to],test_value);
		}
	}
	for (ui p = 0; p < X.size; p++)
		Z.push(X.nodes[p]);
	
	delete[] vis;
}

void ReorienNetwork::update_density_div_and_con(ui X_rank, ui Y_rank, ui Z_rank)
{
	Set& X = R[X_rank], & Y = R[Y_rank], & Z = R[Z_rank];
	queue<ui> Q;
	bool* vis = (bool*)malloc(n * sizeof(bool));
	for (ui p = 0; p < Y.size; p++)
	{
		ui i = Y.nodes[p];
		if (X.in[i] || cur_vis[i]) continue;
		vis[i] = false;
		if (deg[i] >= test_value)
		{
			r[i] = max(r[i], test_value);
			vis[i] = true;
			Q.push(i);
		}
	}

	while(!Q.empty()){
		ui x = Q.front();
		Q.pop();
		Z.push(x);
		for (ui i = pstart[x]; i < pend[x]; i++)
		{
			ui v = edgeList[i];
			ui id = eid[i];
			di_edge& ne = edges[id];
			if (ne.to == x ) continue;
			if (vis[v] || X.in[v] || cur_vis[v])  continue;
			ui to = ne.to;
			vis[to] = true;
			Q.push(to);
			r[to] = max(r[to],test_value);
		}
	}
	for (ui p = 0; p < X.size; p++)
		Z.push(X.nodes[p]);
	Z.E_number = 0;
	for (ui p = 0; p < Z.size; p++){
		Z.E_number += deg[Z.nodes[p]];
	}
	

	delete[] vis;
}

void ReorienNetwork::update_density_of_each_vertex_basic(){
	queue<ui> Q;
	int pre_num = 0;
	for (ui i = 0; i < n; i++)
	{
		cur_vis_copy[i] = true;
		if (!cur_vis[i] && deg[i] >= test_value){
			Q.push(i);
			r[i] = test_value;
		}
	}

	// bool *cur_vis_copy = new bool[n]();
	// for(ui i = 0; i < n; i++) cur_vis_copy[i] = true;
	
	while(!Q.empty()){
		ui x = Q.front();
		Q.pop();
		cur_vis_copy[x] = false;
		for (ui i = pstart[x]; i < pend[x]; i++)
		{
			ui v = edgeList[i];
			if (cur_vis[v]) continue;
			ui id = eid[i];
			di_edge& ne = edges[id];
			if (ne.to == x ) continue;
			
			ui to = ne.to;
			if(r[to] == test_value) continue;
			Q.push(to);
			r[to] = test_value;
		}
	}
	// printf("test_val: %u, c_v: %u, c_m: %u, density: %lf, connect_last_layer: %u\n",test_value,c_v,c_m,c_m*1.0/c_v,connect_ver);
	memcpy(cur_vis,cur_vis_copy,sizeof(bool)*n);	
}


void ReorienNetwork::bottom_up_density_decom()
{
	memset(r,0,sizeof(ui)*n);
    ui now_rank = 1;
    while(!ReTestBoUp(now_rank)){
        now_rank++;
    }
	for(ui i = 0; i < n; i++){
		cur_vis[i] = false;
	}
	// printf("now_rank: %u\n",now_rank);
}

void ReorienNetwork::top_down_density_decom()
{

}

bool ReorienNetwork::compute_R_density(bool *vis_info, ui k){
	memset(deg,0,sizeof(ui)*n);
	for(ui i = 0; i < n; i++){
		if(vis_info[i]) continue;
		for(ui j = pstart[i]; j < pend[i]; j++){
			ui v = edgeList[j];
			if(vis_info[v] || v < i) continue;
			di_edge die = edges[eid[j]];
			ui from = die.to == die.end ? die.src : die.end;
			assert(from == v || from == i);
			deg[from]++;
		}
	}
	ui c_n = 0, c_m = 0;
	for(ui i = 0; i < n; i++){
		if(vis_info[i]) continue;
		c_n++;
		c_m += deg[i];
	}
	if(c_n == 0) return true;
	ui int_den = ceil(c_m*1.0/c_n);
	if(int_den < k)
		return true;
	return false; 

}


void ReorienNetwork::extract_R_sub(ui k, bool* vis){

	for(ui i = 0; i < n; i++){
		vis[i] = true;
		if(!cur_vis[i] && r[i] >= k){
			vis[i] = false;
		}
	}
}

void ReorienNetwork::filter_subgraphs(vector<vector<bool>>& filter_sub, vector<vector<ui>>& filter_vec, vector<ui>& initial_vec){
	// vector<bool> root(n,false);
	// filter_sub.push_back(root);
	ui max_d = get_max_d();
	std::vector<std::vector<bool>> tmp_filter_sub((max_d+1), std::vector<bool>(n, true));
	for(ui i = 0; i < n; i++){
		if(cur_vis[i]) continue;
		ui den = r[i];
		for(ui j = 0; j <= den; j++){
			tmp_filter_sub[j][i] = false;
		}
	}
	for(ui i = 0; i < max_d+1; i++){
		filter_sub.push_back(tmp_filter_sub[i]);
		vector<ui> tmp_filter_vec = initial_vec;
		tmp_filter_vec.push_back(i);
		filter_vec.push_back(tmp_filter_vec);
	}
		
	// std::swap(filter_sub, tmp_filter_sub);
}


void ReorienNetwork::divide_and_conquer_density_decom_plus(){
	ui max_d = get_max_d();
	memset(r,0,sizeof(ui)*n);
	// printf("Degeneracy: %u\n",max_d);
	R = (Set*)malloc((max_d + 2) * sizeof(Set));
	for (ui i = 0; i <= max_d + 1; i++)
	{
		R[i].reset(n);
		R[i].have_found = false;
	}
	//R[0] is V
	ui tmp_m = 0;
	for (ui i = 0; i < n; i++){
		if(!cur_vis[i]){
			R[0].push(i);
			tmp_m += deg[i];
		}
			
	}
		
	R[0].E_number = tmp_m;
	R[0].have_found = true;
	//R[max_d+1] is empty set
	R[max_d + 1].have_found = true;
	//begin recursion
	recursive_function(max_d+1, 0);

	delete[] R;
}

void ReorienNetwork::divide_and_conquer_density_decom_incremental_space_efficient(ReorienNetwork *rn, ui index_sub, ui index, ui layer, ui& ub, vector<ui> veck,
	 ui &num_densub, ui &computed_densub, ui &num_of_iter,
	vector<vector<ui>>& veck_set, vector<ui>& veck_set_to_size, vector<vector<bool>>& filter_sub, vector<vector<ui>>& filter_veck, vector<vector<ui>>& edge_set)
{
	if(index_sub == 0)
	{
		timeCount t;
		t.StartTime();
		max_pa = core_decomposition_incremental(); 
		t.EndTime();
		peel_core_time+= t.QueryTime();
	}
	else{
		max_pa = ub;
	}
	ui max_d = max_pa;
	// printf("Degeneracy: %u, index_sub: %u, layer: %u, ub: %u\n",max_d,index_sub,index,ub);
	if(max_d == 0) {
		
		vector<ui> tmp_filter_veck = veck;
		tmp_filter_veck.push_back(0);
		
		if(index == layer-1)
		{
			veck_set.push_back(tmp_filter_veck);
			ui c = 0;
			vector<ui> tmp_vector;
			for(ui i = 0; i < n; i++){
				if(!cur_vis[i]){
					c++;
					// tmp_vector.push_back(i);
				}
			}
			veck_set_to_size.push_back(c);
			num_densub++;		
			// compute_each_layer_sub_edge(rn,tmp_vector,layer,cur_vis);
			// edge_set.push_back(tmp_vector);
			return;
		}



		vector<bool> tmp_filter_sub(n,false);
		for(ui i = 0; i < n; i++){
			if(!cur_vis[i]) continue;
			tmp_filter_sub[i] = cur_vis[i];
		}
		filter_sub.push_back(tmp_filter_sub);
		filter_veck.push_back(tmp_filter_veck);

		ub = max_pa;
		return ;
	}
	timeCount t;
	t.StartTime();

	max_pa = 0;
	memset(r,0,sizeof(ui)*n);
	R = (Set*)malloc((max_d + 2) * sizeof(Set));
	for (ui i = 0; i <= max_d + 1; i++)
	{
		R[i].reset(n);
		R[i].have_found = false;
	}
	//R[0] is V
	ui tmp_m = 0;
	for (ui i = 0; i < n; i++){
		if(!cur_vis[i]){
			R[0].push(i);
			tmp_m += deg[i];
		}		
	}
		
	R[0].E_number = tmp_m;
	R[0].have_found = true;
	//R[max_d+1] is empty set
	R[max_d + 1].have_found = true;
	//begin recursion

	t.EndTime();
	reconfig_info_time += t.QueryTime();

	recursive_function_incremental(rn, index_sub, index, veck, computed_densub, iter_time, max_d+1, 0);

	t.StartTime();
	ub = max_pa;

	for(ui i = 0; i <= max_d; i++){
		if(R[i].size){
			vector<ui> tmp_filter_veck = veck;
			tmp_filter_veck.push_back(i);
			
			if(index == layer-1){
				num_densub++;
				veck_set.push_back(tmp_filter_veck);
				veck_set_to_size.push_back(R[i].size);
				// vector<ui> tmp_vector;
				// for(ui x = 0; x < R[i].size; x++)
				// 	tmp_vector.push_back(R[i].nodes[x]);
				// edge_set.push_back(tmp_vector);
				continue;
			}
			

			vector<bool> tmp_filter_sub(n,false);
			for(ui j = 0; j < n; j++){
				if(R[i].in[j])
					tmp_filter_sub[j] = false;
				else
					tmp_filter_sub[j] = true;
			}
			filter_sub.push_back(tmp_filter_sub);
			filter_veck.push_back(tmp_filter_veck);
		}
	}
	for (ui i = 0; i <= max_d + 1; i++)
	{
		R[i].free_memory();
	}
	free(R);
	t.EndTime();
	reconfig_info_time += t.QueryTime();
}



void ReorienNetwork::bottom_up_density_decom_incremental(ReorienNetwork *rn, ui index_sub, ui index, vector<ui> veck, ui &num_densub,ui &num_computed_densub,
	ui &num_of_call,
	vector<vector<ui>>& veck_set, vector<ui>& veck_set_to_size, vector<vector<bool>>& filter_sub, vector<vector<ui>>& filter_veck)
{
	if(index >= 2 && veck[index-1] == 0){
		max_pa = core_decomposition_incremental();  
	}
	ui max_d = max_pa;
	// printf("Degeneracy: %u\n",max_d);
	
	{
		num_densub++;
		vector<ui> tmp_filter_veck = veck;
		tmp_filter_veck.push_back(0);
		filter_veck.push_back(tmp_filter_veck);
		veck_set.push_back(tmp_filter_veck);
		ui c = 0;
		for(ui i = 0; i < n; i++)
			if(!cur_vis[i])
				c++;
		veck_set_to_size.push_back(c);
		vector<bool> tmp_filter_sub(n,false);
		for(ui i = 0; i < n; i++)
			tmp_filter_sub[i] = cur_vis[i];
		filter_sub.push_back(tmp_filter_sub);
		
		if(max_d == 0) 
			return ;
	}
	max_pa = 0;
	memset(r,0,sizeof(ui)*n);
	Set Z(n);
	ui now_rank = 1;
    while(!ReTestBoUp(now_rank)){
		num_computed_densub++;
		for(ui j = 0; j < n; j++)
			if(!cur_vis[j])
				Z.push(j);
		num_of_call++;
		update_kdense_subgraph_high_dim_bottom_up(rn,index,veck,now_rank,num_of_call,Z);
		if(Z.size == 0)
			break;
		else{
			num_densub++;
			vector<ui> tmp_filter_veck = veck;
			tmp_filter_veck.push_back(now_rank);
			filter_veck.push_back(tmp_filter_veck);
			veck_set.push_back(tmp_filter_veck);
			veck_set_to_size.push_back(Z.size);
			vector<bool> tmp_filter_sub(n,false);
			for(ui i = 0; i < n; i++)
				tmp_filter_sub[i] = cur_vis[i];
			filter_sub.push_back(tmp_filter_sub);
		}
		// printf("Z.rank: %u, Z.size: %u\n",now_rank,Z.size);
		
        now_rank++;
		Z.refresh();
    }
	max_pa = now_rank;
	// printf("now_rank: %u\n",now_rank);
}

void ReorienNetwork::divide_and_conquer_density_decom_incremental(ReorienNetwork *rn, ui index_sub, ui index, vector<ui> veck, ui &num_densub,
	vector<vector<ui>>& veck_set, vector<ui>& veck_set_to_size, vector<vector<bool>>& filter_sub, vector<vector<ui>>& filter_veck)
{
	if(index >= 2 && veck[index-1] == 0){
		max_pa = core_decomposition_incremental();  
	}
	ui max_d = max_pa;
	// printf("Degeneracy: %u\n",max_d);
	if(max_d == 0) {
		num_densub++;
		vector<ui> tmp_filter_veck = veck;
		tmp_filter_veck.push_back(0);
		filter_veck.push_back(tmp_filter_veck);
		veck_set.push_back(tmp_filter_veck);
		ui c = 0;
		for(ui i = 0; i < n; i++)
			if(!cur_vis[i])
				c++;
		veck_set_to_size.push_back(c);
		vector<bool> tmp_filter_sub(n,false);
		for(ui i = 0; i < n; i++)
			tmp_filter_sub[i] = cur_vis[i];
		filter_sub.push_back(tmp_filter_sub);
		
		return ;
	}
	max_pa = 0;

	memset(r,0,sizeof(ui)*n);
	
	R = (Set*)malloc((max_d + 2) * sizeof(Set));
	for (ui i = 0; i <= max_d + 1; i++)
	{
		R[i].reset(n);
		R[i].have_found = false;
	}
	//R[0] is V
	ui tmp_m = 0;
	for (ui i = 0; i < n; i++){
		if(!cur_vis[i]){
			R[0].push(i);
			tmp_m += deg[i];
		}		
	}
		
	R[0].E_number = tmp_m;
	R[0].have_found = true;
	//R[max_d+1] is empty set
	R[max_d + 1].have_found = true;
	//begin recursion

	recursive_function_incremental(rn, index_sub, index, veck, num_densub,tmp_m, max_d+1, 0);

	for(ui i = 0; i <= max_d; i++){
		if(R[i].size){
			num_densub++;
			vector<ui> tmp_filter_veck = veck;
			// for(ui j = 0; j < veck.size(); j++)
			// 	printf("index :%u, i: %u, veck[%u]: %u\n",index,i,j,veck[j]);
			tmp_filter_veck.push_back(i);
			filter_veck.push_back(tmp_filter_veck);
			veck_set.push_back(tmp_filter_veck);
			veck_set_to_size.push_back(R[i].size);

			vector<bool> tmp_filter_sub(n,false);
			for(ui j = 0; j < n; j++){
				if(R[i].in[j])
					tmp_filter_sub[j] = false;
				else
					tmp_filter_sub[j] = true;
			}
				
			filter_sub.push_back(tmp_filter_sub);
		}
	}
	for (ui i = 0; i <= max_d + 1; i++)
	{
		R[i].free_memory();
	}
	free(R);
}

void ReorienNetwork::divide_and_conquer_density_decom()
{
	timeCount t;
	t.StartTime();
	ui max_d = get_max_d();
	memset(r,0,sizeof(ui)*n);
	// printf("Degeneracy: %u\n",max_d);
	R = (Set*)malloc((max_d + 2) * sizeof(Set));
	for (ui i = 0; i <= max_d + 1; i++)
	{
		R[i].reset(n);
		R[i].have_found = false;
	}
	//R[0] is V
	for (ui i = 0; i < n; i++)
		R[0].push(i);
	R[0].E_number = m;
	R[0].have_found = true;
	//R[max_d+1] is empty set
	R[max_d + 1].have_found = true;
	//begin recursion
	recursive_function(max_d+1, 0);
	delete[] R;

	t.EndTime();
	network_flow_time += t.QueryTime();
}

bool ReorienNetwork::DinicBFSDivAndCon(Set& X, Set& Y, ui &high_num)
{
	high_num = 0;
	ui t_dist = INF;
	// X is corresponding to the most densest subgraph than Y
	queue<ui> Q;
	for (ui p = 0; p < Y.size; p++)
	{
		ui i = Y.nodes[p];
		dist[i] = INF;
		if (X.in[i]||cur_vis[i]) continue;
		if (deg[i] > test_value)
		{
			dist[i] = 0;
			Q.push(i);
			high_outdegree[high_num++] = i;
		}
	}

	while (!Q.empty())
	{
		ui x = Q.front();
		Q.pop();
		for (ui i = pstart[x]; i < pend[x]; i++)
		{
			ui v = edgeList[i];
			if (cur_vis[v]) continue;
			di_edge& ne = edges[eid[i]];
			if (ne.to == x) continue;
			ui to = ne.to; 
			if (dist[to] != INF || X.in[to]) continue;
			dist[to] = dist[x] + 1;
			if (deg[to] < test_value) {
				t_dist = dist[x] + 1;
			}
			Q.push(to);
		}
	}
	return t_dist != INF;
}




void ReorienNetwork::ReTestDivAndCon_incremental(ReorienNetwork *rn, ui index, vector<ui> veck, ui& computed_num_densub, ui& num_of_iter, ui Z_rank)
{
	bool flag = false;
	// if(index == 2 && veck[0] == 2 && veck[1] == 0)
	// 	flag = true;
	//pruning: if R[tem] have been found, then return directly
	if (R[Z_rank].have_found)
		return;
	//find R[Z_rank] within the range R[Y_rank] - R[X_rank], now find the integer X_rank and Y_rank

	computed_num_densub++;

	timeCount t;
    t.StartTime();

	ui X_rank, Y_rank;
	X_rank = Z_rank + 1;
	while (!R[X_rank].have_found)
		X_rank++;
	Y_rank = Z_rank - 1;
	while (!R[Y_rank].have_found)
		Y_rank--;
	Set& X = R[X_rank], & Y = R[Y_rank], & Z = R[Z_rank];

	test_value = Z_rank - 1;
	ui high_num = 0;

	Set tmp_Z(n);
	// memcpy(cur_vis_copy,cur_vis,sizeof(bool)*n);

	// printf("X_rank: %u, Y_rank: %u\n",X_rank,Y_rank);
	num_of_iter++;
	while (DinicBFSDivAndCon( X, Y, high_num))
	{
		for (ui p = 0; p < Y.size; p++)
			cur[Y.nodes[p]] = pstart[Y.nodes[p]];
		for (ui i = 0; i < high_num; i++)
		{
			p[high_outdegree[i]] = INF + 1;
			DinicDFS(high_outdegree[i]);
		}
	}
	test_value++;
	// update_density_div_and_con(X_rank, Y_rank, Z_rank);
	update_density_div_and_con_incremental(X, Y, tmp_Z);

	t.EndTime();
    network_flow_time += t.QueryTime();

	if(flag)
		printf("before update veck[0]: %u, Z_rank: %u, Z.size: %u\n",veck[0],Z_rank,tmp_Z.size);
	// update_kdense_subgraph_incremental(rn, index, veck, X_rank, Y_rank, Z_rank, tmp_Z);
	// update_kdense_subgraph_2dim_incremental(rn, index, veck, X_rank, Y_rank, Z_rank, tmp_Z);
	update_kdense_subgraph_high_dim_incremental(rn, index, num_of_iter, veck, X_rank, Y_rank, Z_rank, tmp_Z);
	if(flag)
		printf("after update veck[0]: %u, Z_rank: %u, Z.size: %u\n",veck[0],Z_rank,tmp_Z.size);

	t.StartTime();	

	if(tmp_Z.size != 0){
		max_pa = max(max_pa,Z_rank);
		
	}
	R[Z_rank].have_found = true;

	Z.E_number = 0;
	for(ui i = 0; i < tmp_Z.size; i++){
		Z.push(tmp_Z.nodes[i]);
		Z.E_number += deg[tmp_Z.nodes[i]];
	}

	memcpy(cur_vis,cur_vis_copy,sizeof(bool)*n);
	tmp_Z.free_memory();

	t.EndTime();
	reconfig_info_time += t.QueryTime();
}



void ReorienNetwork::ReTestDivAndCon(ui Z_rank)
{
	//pruning: if R[tem] have been found, then return directly
	if (R[Z_rank].have_found)
		return;
	//find R[Z_rank] within the range R[Y_rank] - R[X_rank], now find the integer X_rank and Y_rank
	ui X_rank, Y_rank;
	X_rank = Z_rank + 1;
	while (!R[X_rank].have_found)
		X_rank++;
	Y_rank = Z_rank - 1;
	while (!R[Y_rank].have_found)
		Y_rank--;
	Set& X = R[X_rank], & Y = R[Y_rank];
	test_value = Z_rank - 1;
	ui high_num = 0;
	// ui *high_outdegree = new ui[n]();
    // std::cout << "X_rank: " << X_rank << ", Z_rank: "<< Z_rank << ", Y_rank: " << Y_rank << ", Thread ID: " << this_id << std::endl;
	// printf("X_rank: %u, Y_rank: %u\n",X_rank,Y_rank);
	while (DinicBFSDivAndCon(X, Y, high_num))
	{
		for (ui p = 0; p < Y.size; p++)
			cur[Y.nodes[p]] = pstart[Y.nodes[p]];
		for (ui i = 0; i < high_num; i++)
		{
			p[high_outdegree[i]] = INF + 1;
			DinicDFS(high_outdegree[i]);
		}
	}
	test_value++;
	update_density_div_and_con(X_rank, Y_rank, Z_rank);
	R[Z_rank].have_found = true;


	// delete[] high_outdegree;
}


void ReorienNetwork::recursive_function(ui X_rank, ui Y_rank){
	Set& X = R[X_rank], & Y = R[Y_rank];
	//get the number of edges
	ui XY_E = R[Y_rank].E_number - R[X_rank].E_number;
	//define variables of the binary search
	ui l = Y_rank, r = X_rank, l_layer = Y_rank, r_layer = X_rank;
	//first find the R[i], such that i is the maximum number satisfying |XY_E - E(R[i])| < |XY_E|/2
	while (r > l)
	{
		ui mid = (r + l + 1) / 2;
		ReTestDivAndCon(mid);
		if (R[Y_rank].E_number - R[mid].E_number < XY_E / 2)
			l = mid, l_layer = mid;
		else
			r = mid - 1, r_layer = mid;
	}
	ui k = l;

	//begin deeper recursion
	// cout << "X_rank, k, Y_rank, time: " << X_rank << ", " << k << ", " << Y_rank << "\n";
	if (k - Y_rank >= 2 && R[Y_rank].size != R[k].size)
		recursive_function(k, Y_rank);
	k++;
	if (X_rank - k >= 2 && R[X_rank].size != R[k].size)
		recursive_function(X_rank, k);
	return;
}

void ReorienNetwork::recursive_function_incremental(ReorienNetwork *rn, ui index_sub, ui index, vector<ui> veck, ui &computed_num_densub, ui & num_of_iter, ui X_rank, ui Y_rank)
{
	bool flag = false;
	// if(index == 1)
		// flag = true;
	
	Set& X = R[X_rank], & Y = R[Y_rank];
	//get the number of edges
	ui XY_E = R[Y_rank].E_number - R[X_rank].E_number;
	//define variables of the binary search
	ui l = Y_rank, r = X_rank;

	//using binary to find the R[i], such that i is the maximum number satisfying |XY_E - E(R[i])| < |XY_E|/2
	// while (r > l)
	// {
	// 	ui mid = (r + l + 1) / 2;
	// 	ReTestDivAndCon_incremental(rn,index,veck,num_densub,mid);
	// 	if(flag){
	// 		printf("r: %u, l: %u, mid: %u, R[mid].size: %u\n",r,l,mid,R[mid].size);
	// 	}
	// 	if (R[Y_rank].E_number - R[mid].E_number < XY_E / 2)
	// 		l = mid;
	// 	else
	// 		r = mid - 1;
	// }
	// ui k = l;

	//no using binary method
	memcpy(cur_vis_copy,cur_vis,sizeof(bool)*n);
	ui mid = (r + l + 1) / 2;
	if(index_sub == 0){
		timeCount t;
		t.StartTime();
		for(ui u = 0; u < rn[0].get_n(); u++){
			if(!cur_vis[u] && core[u] < mid)
			{
				cur_vis[u] = true;
				for(ui j = pstart[u]; j < pend[u]; j++){
					ui v = edgeList[j];
					if(cur_vis[v] || core[v] < mid) continue;
					di_edge& ed = edges[eid[j]];
					if(ed.to == u){
						deg[v]--;
						ed.to = v;
						deg[u]++;
					}
				}
			}
		}
		t.EndTime();
		peel_core_time += t.QueryTime();
	}
	ReTestDivAndCon_incremental(rn,index,veck,computed_num_densub,num_of_iter,mid);

	memcpy(cur_vis,cur_vis_copy,sizeof(bool)*n);
	ui k = mid;


	//begin deeper recursion

	if (k - Y_rank >= 2 && R[Y_rank].size != R[k].size)
		recursive_function_incremental(rn,index_sub,index,veck,computed_num_densub,num_of_iter, k, Y_rank);
	else if(k - Y_rank >= 2 && R[Y_rank].size == R[k].size){
		timeCount t;
		t.StartTime();

		for(ui i = Y_rank+1; i < k; i++){
			if(R[i].have_found) continue;
			R[i].copy(R[Y_rank]);
			// for(ui j = 0; j < R[Y_rank].size; j++){
			// 	R[i].push(R[Y_rank].nodes[j]);
			// }
			// R[i].have_found = true;
		}

		t.EndTime();
		reconfig_info_time += t.QueryTime();
	}

	// k++;
	if (X_rank - k >= 2 && R[X_rank].size != R[k].size)
		recursive_function_incremental(rn,index_sub,index,veck,computed_num_densub,num_of_iter, X_rank, k);
	else if(X_rank - k >= 2 && R[X_rank].size == R[k].size){

		timeCount t;
		t.StartTime();

		for(ui i = k+1; i < X_rank; i++){
			if(R[i].have_found) continue;
			R[i].copy(R[k]);
			// for(ui j = 0; j < R[k].size; j++){
			// 	R[i].push(R[k].nodes[j]);
			// }
			// R[i].have_found = true;
		}

		t.EndTime();
		reconfig_info_time += t.QueryTime();
	}
	return;
}

void ReorienNetwork::adjust_deg() {
	timeCount t;
    t.StartTime();

	memset(deg,0,sizeof(ui)*n);
	memset(abs_deg,0,sizeof(ui)*n);
	for(ui i = 0; i < n; i++){
		if(cur_vis[i]) continue;
		for(ui j = pstart[i]; j < pend[i]; j++){
			ui v = edgeList[j];
			if(cur_vis[v] || v < i) continue;
			di_edge die = edges[eid[j]];
			ui from = die.to == die.end ? die.src : die.end;
			assert(from == v || from == i);
			deg[from]++;
			abs_deg[die.src]++;
			abs_deg[die.end]++;
		}
	}
	t.EndTime();
    update_graph_time += t.QueryTime();

}


void ReorienNetwork::peel_k_core_singlelayer(ui j){
	cur_vis[j] = true;
	for(ui x = pstart[j]; x < pend[j]; x++){
		ui v = edgeList[x];
		if(cur_vis[v]) continue;
		assert(abs_deg[v]>0);
		abs_deg[v]--;
	}
}


void ReorienNetwork::peel_prune_vertex(ui k){
	memset(deg,0,sizeof(ui)*n);
	for(ui i = 0; i < n; i++){
		if(cur_vis[i]) continue;
		for(ui j = pstart[i]; j < pend[i]; j++){
			ui v = edgeList[j];
			if(cur_vis[v] || v < i) continue;
			deg[v]++;
			deg[i]++;
		}
	}
	queue<ui> to_remove;
	for(ui i = 0; i < n; i++){
		if(cur_vis[i]) continue;
		if(deg[i] < k) {
			to_remove.push(i);
			cur_vis[i] = true;
		}
	}
	// if(to_remove.size())
	// 	printf("to_remove.size: %u\n",to_remove.size());
	while(!to_remove.empty()){
		ui u = to_remove.front();
		to_remove.pop();
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(cur_vis[v]) continue;
			assert(deg[v]>0);
			deg[v]--;
			if(deg[v] < k) {
				to_remove.push(v);
				cur_vis[v] = true;
			}
		}
	}
}


void ReorienNetwork::find_k_dense_sub_incremental(ui tem, Set& X, ui X_rank, Set& Y, ui Y_rank, Set& Z)
{
	timeCount t;
    t.StartTime();

	if(tem == 0){
		Z.refresh(); 
		for(ui i = 0; i < n; i++){
			if(cur_vis[i])continue;
			Z.push(i);
		}
		return;
	}
	    
    assert(tem>0);
    test_value = tem-1;

	ui high_num = 0;
	while (DinicBFSDivAndCon(X, Y, high_num))
	{
		for (ui p = 0; p < Y.size; p++)
			cur[Y.nodes[p]] = pstart[Y.nodes[p]];
		for (ui i = 0; i < high_num; i++)
		{
			p[high_outdegree[i]] = INF + 1;
			DinicDFS(high_outdegree[i]);
		}
	}
	test_value++;
	
	Z.refresh();
	update_density_div_and_con_incremental(X, Y, Z);

	t.EndTime();
    network_flow_time += t.QueryTime();
}



void ReorienNetwork::find_k_dense_sub_plus(vector<ui>& isContainVer, vector<ui>& change, ui tem, ui layer_){
	timeCount t;
    t.StartTime();
    isContainVer.clear();
	change.clear();
    if(tem == 0){
        for(ui i = 0; i < n; i++){
            if(!cur_vis[i])
                isContainVer.push_back(i);
        }

		#ifdef DEBUG
		string ss = "/home/jjq/research/NetworkFlow/multilayer_graph/density_decomposition/dataset/"+graph_name+"/" + to_string(layer_) +  "_dense_graph.txt";
		auto out = ofstream(ss);
		for(ui u = 0; u < n; u++)
		{
			if(cur_vis[u]) continue;
			for(ui j = pstart[u]; j < pend[u]; j++){
				ui v = edgeList[j];
				if(v < u) continue;
				if(cur_vis[v]) continue;
				out << u << "," << v << endl;
			}
		}
		out.close();
		#endif
        return;
    }
	    
    assert(tem>0);
    test_value = tem-1;
	unordered_set<ui> nodes_to_decrease;
	memset(dist,0,sizeof(ui)*n);
	memset(cur,0,sizeof(ui)*n);
	for (unsigned i = 0; i < n; i++)
	{
		if (!cur_vis[i] && deg[i] > test_value) //vertex with outdegree greater than test_value
		{
			nodes_to_decrease.insert(i);
		}
	}
	while (DinicBFS(nodes_to_decrease))
	{
		memcpy(cur, pstart, n * sizeof(ui));
		vector<ui> to_erase;
		for (auto i : nodes_to_decrease)
		{
			p[i] = INF + 1;
			DinicDFS(i);
			if (deg[i] <= test_value){
                to_erase.push_back(i);
            }
		}
		for (auto i : to_erase)
			nodes_to_decrease.erase(i);
	}
    test_value++;

    memset(r,0,sizeof(ui)*n);

    queue<ui> Q;
	for (ui i = 0; i < n; i++)
	{
		cur_vis_copy[i] = true;
		if (!cur_vis[i] && deg[i] >= test_value){
			Q.push(i);
			r[i] = test_value;
            isContainVer.push_back(i);
		}
	}

	
	while(!Q.empty()){
		ui x = Q.front();
		Q.pop();
		cur_vis_copy[x] = false;
		for (ui i = pstart[x]; i < pend[x]; i++)
		{
			ui v = edgeList[i];
			ui id = eid[i];
			di_edge& ne = edges[id];
			if (ne.to == x ) continue;
			if (cur_vis[v]) continue;
			ui to = ne.to;
			if(r[to] == test_value) continue;
			Q.push(to);
			r[to] = test_value;
            isContainVer.push_back(to);
		}
	}

	for(ui i = 0; i < n; i++){
		if(cur_vis[i]) continue;
		if(cur_vis_copy[i]) change.push_back(i);
	}

    memcpy(cur_vis,cur_vis_copy,sizeof(bool)*n);	

	t.EndTime();
    network_flow_time += t.QueryTime();
}


void ReorienNetwork::find_k_dense_sub(vector<ui>& isContainVer, ui tem, ui layer_){
	timeCount t;
    t.StartTime();
    isContainVer.clear();
    if(tem == 0){
        for(ui i = 0; i < n; i++){
            if(!cur_vis[i])
                isContainVer.push_back(i);
        }
		#ifdef DEBUG
		string ss = "/home/jjq/research/NetworkFlow/multilayer_graph/density_decomposition/dataset/"+graph_name+"/" + to_string(layer_) +  "_dense_graph.txt";
		auto out = ofstream(ss);
		for(ui u = 0; u < n; u++)
		{
			if(cur_vis[u]) continue;
			for(ui j = pstart[u]; j < pend[u]; j++){
				ui v = edgeList[j];
				if(v < u) continue;
				if(cur_vis[v]) continue;
				out << u << "," << v << endl;
			}
		}
		out.close();
		#endif
        return;
    }
	    
    assert(tem>0);
    test_value = tem-1;
	unordered_set<ui> nodes_to_decrease;
	memset(dist,0,sizeof(ui)*n);
	memset(cur,0,sizeof(ui)*n);
	for (unsigned i = 0; i < n; i++)
	{
		if (!cur_vis[i] && deg[i] > test_value) //vertex with outdegree greater than test_value
		{
			nodes_to_decrease.insert(i);
		}
	}
	while (DinicBFS(nodes_to_decrease))
	{
		memcpy(cur, pstart, n * sizeof(ui));
		vector<ui> to_erase;
		for (auto i : nodes_to_decrease)
		{
			p[i] = INF + 1;
			DinicDFS(i);
			if (deg[i] <= test_value){
                to_erase.push_back(i);
            }
		}
		for (auto i : to_erase)
			nodes_to_decrease.erase(i);
	}
    test_value++;

    memset(r,0,sizeof(ui)*n);

    queue<ui> Q;
	for (ui i = 0; i < n; i++)
	{
		cur_vis_copy[i] = true;
		if (!cur_vis[i] && deg[i] >= test_value){
			Q.push(i);
			r[i] = test_value;
            isContainVer.push_back(i);
		}
	}

	
	while(!Q.empty()){
		ui x = Q.front();
		Q.pop();
		cur_vis_copy[x] = false;
		for (ui i = pstart[x]; i < pend[x]; i++)
		{
			ui v = edgeList[i];
			ui id = eid[i];
			di_edge& ne = edges[id];
			if (ne.to == x ) continue;
			if (cur_vis[v]) 
			{
				continue;
			}
			ui to = ne.to;
			if(r[to] == test_value) continue;
			Q.push(to);
			r[to] = test_value;
            isContainVer.push_back(to);
		}
	}

    memcpy(cur_vis,cur_vis_copy,sizeof(bool)*n);	
	// delete[] cur_vis_copy;

	t.EndTime();
    network_flow_time += t.QueryTime();

	#ifdef DEBUG
    string ss = "/home/jjq/research/NetworkFlow/multilayer_graph/density_decomposition/dataset/"+graph_name+"/" + to_string(layer_) +  "_dense_graph.txt";
    auto out = ofstream(ss);
    for(ui u = 0; u < n; u++)
    {
        if(cur_vis[u]) continue;
        for(ui j = pstart[u]; j < pend[u]; j++){
            ui v = edgeList[j];
            if(v < u) continue;
            if(cur_vis[v]) continue;
            out << u << "," << v << endl;
        }
    }
    out.close();
	#endif
}


void ReorienNetwork::update_graph_incremental(ui *val, ui c, Set& Z){
	timeCount t;
	t.StartTime();
	for(ui i = 0; i < c; i++){
		ui u = val[i];
		cur_vis[u] = true;
		for(ui j = pstart[u]; j < pend[u]; j++){
			ui v = edgeList[j];
			if(Z.in[v]){
				di_edge de = edges[eid[j]];
				ui from = de.to == de.src ? de.end : de.src;
				if(from == v){
					assert(deg[v]>0);
					deg[v]--;
				}
			}
		}
	}

	t.EndTime();
    update_graph_time += t.QueryTime();
}

void ReorienNetwork::update_graph(bool *vis, ui k)
{
	timeCount t;
    t.StartTime();
    memcpy(cur_vis,vis,sizeof(bool)*n);
	// if(k > 1){
	// 	peel_prune_vertex(k);
	// }
	memset(deg,0,sizeof(ui)*n);

	// memset(abs_deg,0,sizeof(ui)*n);
    
	for (ui i = 0; i < n; i++) {
		if(cur_vis[i]) continue;
		ui ele = i;
		for(ui j = pstart[ele]; j < pend[ele]; j++){
			ui e_id = eid[j];
			ui v = edgeList[j];
			if(v < ele) continue;
			di_edge tmp = edges[e_id];
			if(!cur_vis[v]){
				ui from = tmp.to == tmp.end ? tmp.src : tmp.end;
				assert(from == v || from == ele);
				deg[from]++;
				// abs_deg[tmp.src]++;
				// abs_deg[tmp.end]++;
			}
		}
    }

	t.EndTime();
    update_graph_time += t.QueryTime();
}



bool ReorienNetwork::ReTestBoUp(ui tem)
{
    test_value = tem-1;
	unordered_set<ui> nodes_to_decrease;
	memset(dist,0,sizeof(ui)*n);
	memset(cur,0,sizeof(ui)*n);
	for (unsigned i = 0; i < n; i++)
	{
		if (!cur_vis[i] && deg[i] > test_value) //vertex with outdegree greater than test_value
		{
			nodes_to_decrease.insert(i);
		}
	}
	if(nodes_to_decrease.size() == 0)
		return true;
	while (DinicBFS(nodes_to_decrease))
	{
		memcpy(cur, pstart, n * sizeof(ui));
		vector<ui> to_erase;
		for (auto i : nodes_to_decrease)
		{
			p[i] = INF + 1;
			DinicDFS(i);
			if (deg[i] <= test_value){
                to_erase.push_back(i);
            }
		}
		for (auto i : to_erase)
			nodes_to_decrease.erase(i);
	}
    test_value++;
	update_density_of_each_vertex_basic();

	return get_max_d() < test_value;
}


bool ReorienNetwork::DinicBFS(unordered_set<ui>& nodes_to_decrease)
{
	if (nodes_to_decrease.empty()) return false;	
	bool can_reach_low_indegree = false;
	queue<ui> Q;
	for (unsigned i = 0; i < n; i++)
		dist[i] = INF;
	for (auto i : nodes_to_decrease)
	{
		Q.push(i);
		dist[i] = 0;
	}
	while (!Q.empty())
	{
		ui x = Q.front();
		Q.pop();
		for (ui i = pstart[x]; i < pend[x]; i++)
		{
			ui v = edgeList[i];
			if (cur_vis[v]) continue;
			di_edge& ne = edges[eid[i]];
			if (ne.to == x) continue;
			ui to = ne.to; 
			if (dist[to] != INF) continue;
			dist[to] = dist[x] + 1;
			if (deg[to] < test_value) {
				can_reach_low_indegree = true;
				continue;
			}
			Q.push(to);
		}
	}

	return can_reach_low_indegree;
}

bool ReorienNetwork::DinicDFS(ui x)
{
	if (deg[x] < test_value)
	{
		assert(p[x]<2*m);
		ui from = edges[p[x]].to == edges[p[x]].src ? edges[p[x]].end : edges[p[x]].src;
		deg[edges[p[x]].to]++;
		assert(deg[from]>0);
		deg[from]--;
		edges[p[x]].to = from;
		return true;
	}
	for (ui& i = cur[x]; i < pend[x]; i++)
	{
		di_edge& ne = edges[eid[i]];
		if (ne.to == x) continue;
		ui to = ne.to;
		if (dist[to] != dist[x] + 1 || cur_vis[to]) continue;
		p[to] = eid[i];
		if (DinicDFS(to))
		{
			if (p[x] >= INF){
				if(deg[x] == test_value) {
					// printf("x: %u, test_val: %u\n",x,test_value);
					return true;
				}
				continue;
			} 
			ui from = edges[p[x]].to == edges[p[x]].src ? edges[p[x]].end : edges[p[x]].src;
			deg[edges[p[x]].to]++;
			deg[from]--;
			edges[p[x]].to = from;
			return true;
		}
	}
	return false;
}


ui ReorienNetwork::get_max_d()
{
	ui maxd = 0;
	for (ui i = 0; i < n; i++)
		if(!cur_vis[i])
			maxd = maxd > deg[i] ? maxd : deg[i];
	return maxd;
}

ui ReorienNetwork::get_max_r()
{
	ui maxd = 0;
	for (ui i = 0; i < n; i++)
		if(!cur_vis[i])
			maxd = maxd > r[i] ? maxd : r[i];
	return maxd;
}

void ReorienNetwork::takeSnapshot(ReorienNetworkSnapshot& snapshot) {
	snapshot.cur_vis_ = new bool[n]();
	snapshot.deg_ = new ui[n]();
	snapshot.edges_ = new di_edge[m];
	memcpy(snapshot.cur_vis_,cur_vis,sizeof(bool)*n);
	memcpy(snapshot.deg_,deg,sizeof(ui)*n);
	memcpy(snapshot.edges_,edges,sizeof(di_edge)*m);
}

void ReorienNetwork::restoreSnapshot(ReorienNetworkSnapshot& snapshot) {
	memcpy(cur_vis,snapshot.cur_vis_,sizeof(bool)*n);
	memcpy(deg,snapshot.deg_,sizeof(ui)*n);
	memcpy(edges,snapshot.edges_,sizeof(di_edge)*m);
	delete[] snapshot.deg_;
	delete[] snapshot.cur_vis_;
	delete[] snapshot.edges_;
}

double ReorienNetwork::com_global_cluster_coefficient(){
    double total_cc = 0.0;
    int valid_nodes = 0;

    for(int u = 0; u < n; ++u){
        if(cur_vis[u]) continue;

        // 构建邻居集合
        vector<int> neighbors;
        for(ui j = pstart[u]; j < pend[u]; ++j){
            int v = edgeList[j];
            if(cur_vis[v]) continue;
            neighbors.push_back(v);
        }

        int deg = neighbors.size();
        if(deg < 2) continue;

        // 统计三角形数量（邻居两两之间是否连接）
        int triangle_count = 0;
        unordered_set<int> neighbor_set(neighbors.begin(), neighbors.end());
        for(int i = 0; i < deg; ++i){
            int a = neighbors[i];
            for(ui j = pstart[a]; j < pend[a]; ++j){
                int b = edgeList[j];
                if(cur_vis[b]) continue;
                if(b > a && neighbor_set.count(b)){
                    triangle_count++;
                }
            }
        }

        // 局部聚类系数
        double cc_u = (2.0 * triangle_count) / (deg * (deg - 1));
        total_cc += cc_u;
        valid_nodes++;
    }

    if(valid_nodes == 0) return 0.0;
    return total_cc / valid_nodes;
}


double ReorienNetwork::compute_conductance(const vector<int>& subgraph_nodes) {

    unordered_set<int> sub_set(subgraph_nodes.begin(), subgraph_nodes.end());
    int cut_edges = 0;     // 跨边
    int vol_S = 0;         // 子图体积

    for(int u : subgraph_nodes){
        for(ui j = pstart[u]; j < pend[u]; ++j){
            int v = edgeList[j];
            if(cur_vis[v]) {
				cut_edges++;
				continue;
			}
        }
    }

    for(int u = 0; u < n; ++u){
        if(cur_vis[u]) continue;
        vol_S += pend[u] - pstart[u];
    }

	// printf("cut_edges:%d, vol_S: %d\n",cut_edges,vol_S);

    if(vol_S == 0) return 0.0;
    return static_cast<double>(cut_edges) / vol_S;
}
