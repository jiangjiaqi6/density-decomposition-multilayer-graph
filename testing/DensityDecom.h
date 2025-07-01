

#include "../core/Reorientation.h"
#include "../utilities/MemoryUtils.h"
#include <omp.h>

static void restore(ReorienNetwork *rn, ui layer, ui old_e, CoreIndex& kci, ui** degs){
    ui* pstart, *pend, *edgeList;
    for(ui i = 0; i < layer; i++){
        pstart = rn[i].get_pstart();
        pend = rn[i].get_pend();
        edgeList = rn[i].get_edgeList();
        for (uint j = old_e; j < kci.e; j++) {
            ui v = kci.vert[j];
            for (uint l = pstart[v]; l < pend[v]; l++) {
                ui u = edgeList[l];
                degs[i][u]++;
            }
        }
    }
    kci.e = old_e;
    kci.s = old_e;
}



static void compute_k_core_set_multilayer_with_index(ReorienNetwork *rn, vector<ui> veck, ui layer, CoreIndex& kci, ui** degs){
    ui old_e;
    auto &s = kci.s;
    auto &e = kci.e;
    ui* pstart, *pend, *edgeList;
    bool* cur_vis;
    while(s < e){
        old_e = e;
        for(ui i = 0; i < layer; i++){
            pstart = rn[i].get_pstart();
            pend = rn[i].get_pend();
            edgeList = rn[i].get_edgeList();
            cur_vis = rn[i].get_cur_vis();
            for(ui j = s; j < old_e; j++){
                ui v = kci.vert[j];
                for(ui x = pstart[v]; x < pend[v]; x++){
                    ui u = edgeList[x];
                    assert(degs[i][u] > 0);
                    degs[i][u]--;
                    if(kci.pos[u] >= e){
                        if (degs[i][u] < veck[i]) {
                            kci.Remove(u);
                        }
                    }
                }
            }
        }
        s = old_e;
    }
}

static void compute_k_core_set_multilayer(ReorienNetwork *rn, vector<ui> veck, ui layer, vector<ui>& core_set)
{
    ui n = rn[0].get_n();
    ui* pstart;
    ui* pend;
    ui* edgeList;
    ui* abs_deg;
    bool* cur_vis;
    queue<ui> q;
    bool is_change = true;

    bool flag = false;
    if(veck[0] == 2 && veck[1] == 0 && veck[2] == 0)
        flag = true;
    ui debug_C = 0;
    while(is_change){
        is_change = false;
        for(ui i = 0; i < layer; i++){
            ui key = veck[i];
            if(key == 0) continue;
            pstart = rn[i].get_pstart();
            pend = rn[i].get_pend();
            edgeList = rn[i].get_edgeList();
            abs_deg = rn[i].get_abs_deg();
            cur_vis = rn[i].get_cur_vis();
            
            for(ui j = 0; j < n; j++){
                if(cur_vis[j]) continue;
                if(abs_deg[j] < key){
                    q.push(j);
                    cur_vis[j] = true;
                }
            }
            if(q.empty())continue;
            is_change = true;
            while(!q.empty()){
                ui j = q.front();
                q.pop();
                for(ui x = pstart[j]; x < pend[j]; x++){
                    ui v = edgeList[x];
                    if(cur_vis[v]) continue;
                    assert(abs_deg[v]>0);
                    abs_deg[v]--;
                    if(abs_deg[v] < key){
                        q.push(v);
                        cur_vis[v] = true; 
                    }
                }
                for(ui c = 0; c < layer; c++){
                    if(c != i && veck[c])
                        rn[c].peel_k_core_singlelayer(j);
                }
                
            }   
    
        }
        debug_C++;
    }
    cur_vis = rn[0].get_cur_vis();
    for(ui i = 0; i < n; i++)
        if(!cur_vis[i])
            core_set.push_back(i);
    
}


static vector<ui> KDenseSubNaive(ReorienNetwork *rn, vector<ui> veck, ui layer, ui& iter_time){

    for(ui i = 0; i < layer; i++)
        rn[i].adjust_deg(); 
    
    vector<ui> high_core;
    vector<ui> low_core;
    ui n = rn[0].get_n();
    bool* vis = new bool[n]();
    memcpy(vis,rn[0].get_cur_vis(),sizeof(bool)*n);
    ui c = 0;
    for(ui i = 0; i < n; i++)
        if(!vis[i])
            c++;

    timeCount t;
    t.StartTime();

    vector<ui> thres_veck(layer,0);

    for(ui i = 0; i < layer; i++){
        thres_veck[i] = veck[i];
    }
    
    compute_k_core_set_multilayer(rn,thres_veck,layer,low_core);

    t.EndTime();
    double compute_core_time = t.QueryTime();

    // for(ui i = 0; i < layer; i++){
    //     if(veck[i] == 0)
    //         thres_veck[i] = 0;
    //     else
    //         thres_veck[i] = 2*veck[i] - 1;
    // }
    // compute_k_core_set_multilayer(rn,thres_veck,layer,high_core);

    for(ui i = 0; i < layer; i++){
        // rn[i].use_graph(vis);
        rn[i].adjust_deg(); 
    }
        
    t.StartTime();

    vector<ui> isContain, isContainOld;
    bool exist = true;
    ui iter = 1, i;
    while(true){
        bool all_unchanged = true;
        isContain.clear();
        for(i = 0; i < layer; i++){  
            if(veck[i] > 0) iter_time++;       
            rn[i].find_k_dense_sub(isContain,veck[i],i);   // Optimize according to the order of element values in the veck vector
            if(!isContain.size()){
                exist = false;
                break;
            }
            rn[(i+1)%layer].update_graph(rn[i].get_cur_vis(),veck[(i+1)%layer]);
            if(iter == 1 && i == 0){
                isContainOld = isContain;
                all_unchanged = false;
            }
            else{
                if(isContain.size() != isContainOld.size()){
                    all_unchanged = false;
                    isContainOld = isContain;
                }     
            }
        }
        if(all_unchanged)break;
        if(!exist) break;
        iter++;
    }


    t.EndTime();
    double compute_kdense_time = t.QueryTime();

    delete[] vis;

    if(!exist){
        // printf("There is no dense subgraph with veck\n");
        return isContain;
    }
    else{


        // printf("vector: ");
        // for(ui i = 0; i < layer; i++)
        //     printf("%u ",veck[i]);
        // printf("\n");
        // printf("before graph size: %u, low_core.size: %u, high_core.size: %u, isContainOld.size: %u\n",c,low_core.size(),high_core.size(),isContainOld.size());

        string ss = "core_info.txt";
		auto out = ofstream(ss, std::ios::app);
        out<< "vector: ";
        for(ui i = 0; i < layer; i++)
            out<<veck[i]<<" ";
        out<<", original size: " << c << ", low_core: " << low_core.size()<< ", high_core: "<<high_core.size()<<", now size: "<<isContainOld.size()<< ", compute core time: " << compute_core_time <<"s, compute kdense time: "<< compute_kdense_time <<"s" <<endl;
		out.close();
        return isContainOld;
    }
}


static vector<ui> KDenseSub(ReorienNetwork *rn, vector<ui> veck, ui index, ui layer, ui& iter){

    // for(ui i = 0; i < layer; i++)
        // if(veck[index] > 1)
        //     rn[index].peel_prune_vertex(veck[index]);
    rn[index].adjust_deg(); 

    vector<ui> isContain, isContainOld;
    bool exist = true;
    // ui iter = 1;
    ui c = 0, i;
    while(true){
        bool all_unchanged = true;
        isContain.clear();
        if(iter == 1) i = index;
        else i = 0;
        for(; i < layer; i++){         
            rn[i].find_k_dense_sub(isContain,veck[i],i); 
            if(!isContain.size()){
                exist = false;
                break;
            }
            
            rn[(i+1)%layer].update_graph(rn[i].get_cur_vis(),veck[(i+1)%layer]);
            
            if((iter == 1 && i == index)){
                isContainOld = isContain;
                c=1;
            }
            else{
                if(isContain.size() != isContainOld.size()){
                    all_unchanged = false;
                    isContainOld = isContain;
                    c = 1;
                }
                else{
                    if(c == layer) break;
                    c++; 
                }    
            }
        }
        if(all_unchanged && c == layer)break;
        if(!exist) break;
        iter++;
    }    

    if(!exist){
        // printf("There is no dense subgraph with veck\n");
        return isContain;
    }
    else{
        return isContainOld;
    }
}

// vector is transformed to string
std::string VecToString(const std::vector<ui>& vec) {
    std::string result;
    for (ui val : vec) {
        result += std::to_string(val) + ","; // 使用 ',' 分隔元素
    }
    return result;
}

static void dfs_peel(ReorienNetwork *rn, vector<vector<ui>> &veck_set, vector<vector<ui>> &vertex_set, vector<ui> &veck_set_to_size,
vector<ui> veck, vector<bool> level_vis, CoreIndex& kci, uint** degs,
vector<ui>& iter_time, ui lev_start, ui &num_of_core, ui &num_of_computed_core , ui layer, ui n)
{
    ui old_e = kci.e;
    bool *global_vis = new bool[n]();
    for(ui i = 0; i < n; i++) global_vis[i] = level_vis[i];

    for(ui k = lev_start; k < layer; k++){

        // vector<ReorienNetworkSnapshot> snapshots(layer);
        // for (ui i = 0; i < layer; i++) {
        //     rn[i].takeSnapshot(snapshots[i]);
        // }

        for (uint j = old_e; j < kci.n; j++) {
            ui v = kci.vert[j];
            if (degs[k][v] <= veck[k]) {
                kci.Remove(v);
            }
        }
        ui len = veck_set.size();
        assert(len>0);
        veck[k]++;
        if(kci.e != old_e)
            compute_k_core_set_multilayer_with_index(rn,veck,layer,kci,degs);
        
        for(ui i = 0; i < n; i++){
            global_vis[i] = true;
            if(kci.pos[i] >= kci.e)
                global_vis[i] = false;
        }
        for(ui i = 0; i < layer; i++){
            rn[i].use_graph(global_vis);
        }
        //prune technique

        // bool lb_noless_than_vk = rn[k].compute_lower_bound(veck[k],veck);
        // if(lb_noless_than_vk){
        //     // printf("lb_noless_than_vk is true, k: %u\n",k);
        //     // for(ui i = 0; i < layer; i++)
        //     //     printf("%u ",veck[i]);
        //     // printf("\n");
        //     for(ui i = 0; i < n; i++) level_vis[i] = global_vis[i];
        //     dfs_peel(rn,veck_set,vertex_set,veck_set_to_size,veck,level_vis,k,num_of_core,num_of_computed_core,layer,n);
        // }
        // else
        {
            ui iters = 1;
            // unordered_set<ui> kds = KDenseSub(rn,veck,k,layer,iters);
            vector<ui> kds = KDenseSub(rn,veck,k,layer,iters);
            num_of_computed_core++;
            if(kds.size() != 0){
                for(ui i = 0; i < n; i++) level_vis[i] = true;
                num_of_core++;
                vector<ui> vertex_tmp;
                vertex_tmp.reserve(kds.size());
                for(auto x : kds) {
                    vertex_tmp.push_back(x);
                    level_vis[x] = false;
                }
                ui tmp_e = kci.e;
                for(ui x = 0; x < n; x++){
                    if(global_vis[x])continue;
                    if(level_vis[x]){
                        kci.Remove(x);
                    }
                }
                if(kci.e != tmp_e)
                    compute_k_core_set_multilayer_with_index(rn,veck,layer,kci,degs);
                // sort(vertex_tmp.begin(),vertex_tmp.end());
                vertex_set.push_back(vertex_tmp);
                veck_set.push_back(veck);
                veck_set_to_size.push_back(kds.size());
                iter_time.push_back(iters);
                dfs_peel(rn,veck_set,vertex_set,veck_set_to_size,veck,level_vis,kci,degs,iter_time,k,num_of_core,num_of_computed_core,layer,n);
            }

        }
        veck[k]--;
        restore(rn,layer,old_e,kci,degs);
        // for (ui i = 0; i < layer; i++) {
        //     rn[i].restoreSnapshot(snapshots[i]);
        // }
    }

    delete[] global_vis;
}



static ui compute_sub_satisfy_vec(ui index, vector<ui> &filter_vec, ReorienNetwork* rn)
{
    bool *vis = new bool[rn[index].get_n()]();
    ui max_d = rn[index].get_max_r();
    while(true){
        if(max_d == 0){
            rn[index].extract_R_sub(max_d,vis);
            break;
        }
        else{
            rn[index].extract_R_sub(max_d,vis);
            bool flag = true;
            for(ui i = 0; i < filter_vec.size(); i++){
                ui test_d = filter_vec[i];
                if(rn[i].compute_R_density(vis,test_d))//There is no subgraph with a density of test_d on rn[i]
                {
                    flag = false;
                    break;
                }
            }
            if(flag){
                break;
            }
        }
        max_d--;
    }
    delete[] vis;
    return max_d;
}


static void process_bottom_up_high_dimen(vector<vector<bool>> &filter_sub, vector<vector<ui>> &filter_vec, ui index, ReorienNetwork* rn, ui& num_dense_sub,
    ui& num_of_computed_dense_sub, ui &num_of_call, vector<vector<ui>>& veck_set, vector<ui>& veck_set_to_size){
    ui n = rn[index].get_n();
    bool *arr_vis = new bool[n]();
    vector<vector<bool>> filter_sub_new;
    vector<vector<ui>> filter_vec_new;
    log_info(MultilayerGraph::GetMulGraphClock().Count("index: %u, filter_sub.size: %u",index,filter_sub.size()));

    veck_set.clear();
    veck_set_to_size.clear();

    for(ui i = 0; i < filter_sub.size(); i++){
        vector<bool> tmp_filter_sub = filter_sub[i];
        vector<ui> tmp_filter_vec = filter_vec[i];
        assert(index == tmp_filter_vec.size());
        ui c = 0, num_dense = 0;
        for(ui j = 0; j < n; j++) {
            arr_vis[j] = tmp_filter_sub[j];
            if(!tmp_filter_sub[j])c++;
        }

        rn[index].use_graph(arr_vis);
        rn[index].adjust_deg();
        rn[index].bottom_up_density_decom_incremental(rn,i,index,tmp_filter_vec,num_dense_sub,num_of_computed_dense_sub,num_of_call,veck_set,veck_set_to_size,filter_sub_new,filter_vec_new);

        num_dense_sub += num_dense;
    }
    delete[] arr_vis;
    swap(filter_sub_new,filter_sub);
    swap(filter_vec_new,filter_vec);
    
}


static void process_high_dimen_incremental(vector<vector<bool>> &filter_sub, vector<vector<ui>> &filter_vec, ui index, ReorienNetwork* rn, ui& num_dense_sub,
    vector<vector<ui>>& veck_set, vector<ui>& veck_set_to_size){
    ui n = rn[index].get_n();
    bool *arr_vis = new bool[n]();
    vector<vector<bool>> filter_sub_new;
    vector<vector<ui>> filter_vec_new;
    log_info(MultilayerGraph::GetMulGraphClock().Count("index: %u, filter_sub.size: %u",index,filter_sub.size()));

    veck_set.clear();
    veck_set_to_size.clear();

    for(ui i = 0; i < filter_sub.size(); i++){
        vector<bool> tmp_filter_sub = filter_sub[i];
        vector<ui> tmp_filter_vec = filter_vec[i];
        assert(index == tmp_filter_vec.size());
        ui c = 0, num_dense = 0;
        for(ui j = 0; j < n; j++) {
            arr_vis[j] = tmp_filter_sub[j];
            if(!tmp_filter_sub[j])c++;
        }
        // printf("veck \n");
        // for(ui j = 0; j < tmp_filter_vec.size(); j++)
        //     printf("%u ",tmp_filter_vec[j]);
        // printf("\nsize: %u\n",c);

        rn[index].use_graph(arr_vis);
        rn[index].adjust_deg();
        rn[index].divide_and_conquer_density_decom_incremental(rn,i,index,tmp_filter_vec,num_dense,veck_set,veck_set_to_size,filter_sub_new,filter_vec_new);
        num_dense_sub += num_dense;
    }
    delete[] arr_vis;
    swap(filter_sub_new,filter_sub);
    swap(filter_vec_new,filter_vec);
    
}


//Hirarchical progressive

static void HPDD(MultilayerGraph &mg){
    double mem;
    log_info(MultilayerGraph::GetMulGraphClock().Count("Hierarchical progressive method starts"));
    ui layer = mg.GetLayerNumber();
    ui num_of_dense_sub = 0, num_of_computed_dense_sub = 0, num_of_call = 0; 

    // Record the dense subgraphs corresponding to different vectors
    vector<vector<ui>> veck_set;
    vector<vector<ui>> vertex_set;
    vector<ui> veck_set_to_size;
    
    ReorienNetwork *rn = new ReorienNetwork[layer];
    for(ui i = 0; i < layer; i++){
        rn[i].initial(mg,i);
    }

    vector<vector<bool>> filter_subs;
    vector<vector<ui>> filter_vec;
    vector<ui> initial_vec;

    vector<ui> sub_num_each_layer;
    rn[0].bottom_up_density_decom();
    assert(rn[0].get_max_d() == rn[0].get_max_r());
    rn[0].filter_subgraphs(filter_subs,filter_vec,initial_vec);

    num_of_computed_dense_sub += rn[0].get_max_d();

    sub_num_each_layer.push_back(rn[0].get_max_d());
    
    for(ui i = 1; i < layer; i++)
    {
        num_of_dense_sub = 0;
        process_bottom_up_high_dimen(filter_subs,filter_vec,i,rn,num_of_dense_sub,num_of_computed_dense_sub,num_of_call,veck_set,veck_set_to_size);
    }
    

    mem = GetPeakRSSInMB();

    // string ss = mg.GetDataPath() + "dense_subgraph_res_hpdd.txt";
    // auto out = ofstream(ss);
    // for(ui i = 0; i < veck_set.size(); i++){
    //     out << "[";
    //     for(ui j = 0; j < layer; j++){
    //         if(j == layer-1)
    //             out << veck_set[i][j];
    //         else
    //             out << veck_set[i][j] << ",";
    //     }     
    //     out << "]: " << veck_set_to_size[i] << endl;
    // }
    // out.close();

    double total_flow_time = 0, total_update_graph_time = 0;
    for(ui i = 0; i < layer; i++){
        total_flow_time += rn[i].get_network_flow_time();
        total_update_graph_time += rn[i].get_update_graph_time();
        // printf("%u layer, number of sub: %u\n",i,sub_num_each_layer[i]);
    }

    log_info(MultilayerGraph::GetMulGraphClock().Count("The number of dense subgraph: %u, the number of computed dense subgraph: %u",num_of_dense_sub,num_of_computed_dense_sub));
    log_info(MultilayerGraph::GetMulGraphClock().Count("Memory: %lf MB, the number of times SLBF is called %u",mem,num_of_call));
    log_info(MultilayerGraph::GetMulGraphClock().Count("Network flow time: %lf s, update graph time: %lf s",total_flow_time,total_update_graph_time));
    log_info(MultilayerGraph::GetMulGraphClock().Count("HPDD method ends!"));
    
    delete[] rn;

}


static void Recursive(ReorienNetwork *rn, bool* vis, ui index_sub, ui index, ui layer, ui& ub,  vector<ui> veck, ui &num_dense, ui &computed_densub, ui& num_of_iter, vector<vector<ui>>& veck_set, vector<ui>& veck_set_to_size, vector<vector<ui>>& edge_set){
    ui n = rn[0].get_n();
    vector<vector<bool>> filter_sub;
    vector<vector<ui>> filter_vec;
    vector<ui> initial_vec = veck;
    rn[index].use_graph(vis);
    rn[index].adjust_deg();
    rn[index].divide_and_conquer_density_decom_incremental_space_efficient(rn,index_sub,index,layer,ub,initial_vec,num_dense,computed_densub,num_of_iter,veck_set,veck_set_to_size,filter_sub,filter_vec,edge_set);
    ui tmp_ub = ub;
    // printf("veck size: %u\n",veck.size());
    // for(ui i = 0; i < veck.size(); i++)
    //     printf("%u ",veck[i]);
    // printf("\n");
    // printf("layer: %u, max_pa: %u, subgraphs size: %u\n",index,ub,filter_sub.size());
    if(index+1 < layer){
        for(ui i = 0; i < filter_sub.size(); i++){
            vector<bool> tmp_filter_sub = filter_sub[i];
            for(ui j = 0; j < n; j++)
                vis[j] = tmp_filter_sub[j];
            vector<ui> tmp_filter_vec = filter_vec[i]; 
            Recursive(rn,vis,i,index+1,layer,ub,tmp_filter_vec,num_dense,computed_densub,num_of_iter,veck_set,veck_set_to_size,edge_set);
        }
    }
    ub = tmp_ub;
}


//Space efficient method based divide and conquer
static void SpaceEfficientDivAndCon(MultilayerGraph &mg){
    double mem;
    log_info(MultilayerGraph::GetMulGraphClock().Count("Space efficient method starts"));
    ui layer = mg.GetLayerNumber();
    ui num_of_dense_sub = 0;
    ui num_of_computed_dense_sub = 0; 
    ui iter = 0;

    // Record the dense subgraphs corresponding to different vectors
    vector<vector<ui>> veck_set;
    vector<vector<ui>> vertex_set;
    vector<ui> veck_set_to_size;
    vector<ui> ini_veck;

    vector<vector<ui>> edge_set;

    ReorienNetwork *rn = new ReorienNetwork[layer];
    for(ui i = 0; i < layer; i++){
        rn[i].initial(mg,i);
    }
    bool *vis = new bool[rn[0].get_n()]();
    ui ub = rn[0].get_max_d();

    Recursive(rn,vis,0,0,layer,ub,ini_veck,num_of_dense_sub,num_of_computed_dense_sub,iter,veck_set,veck_set_to_size,vertex_set);


    delete[] vis;
    
    mem = GetPeakRSSInMB();


    string ss = mg.GetDataPath() + "dense_subgraph_res_space_efficient_include_vertices.txt";
    auto out = ofstream(ss);
    for(ui i = 0; i < veck_set.size(); i++){
        out << "[";
        for(ui j = 0; j < layer; j++){
            if(j == layer-1)
                out << veck_set[i][j];
            else
                out << veck_set[i][j] << ",";
        }     
        out << "]: " << veck_set_to_size[i] << endl;
        // out << "], vertex:[";
        // for(ui j = 0; j < vertex_set[i].size(); j++){
        //     if(j == vertex_set[i].size()-1)
        //         out << vertex_set[i][j];
        //     else
        //         out << vertex_set[i][j] << ",";
        // }  
        // out << "]: " << veck_set_to_size[i] << endl;
    }
    out.close();

    double total_flow_time = 0, total_update_graph_time = 0, total_reconfig_info_time = 0, total_peel_core_time = 0;
    for(ui i = 0; i < layer; i++){
        total_flow_time += rn[i].get_network_flow_time();
        total_update_graph_time += rn[i].get_update_graph_time();
        total_reconfig_info_time += rn[i].get_reconfig_info_time();
        total_peel_core_time += rn[i].get_peel_core_time();
    }

    log_info(MultilayerGraph::GetMulGraphClock().Count("The number of dense subgraph: %u, the number of computed dense subgraphs: %u",num_of_dense_sub,num_of_computed_dense_sub));
    log_info(MultilayerGraph::GetMulGraphClock().Count("Memory: %lf MB, the call of slfbe is: ",mem, iter));
    log_info(MultilayerGraph::GetMulGraphClock().Count("Network flow time: %lf s, update graph time: %lf s, reconfig info time: %lf s, peel core time: %lf s",total_flow_time,total_update_graph_time,total_reconfig_info_time,total_peel_core_time));
    log_info(MultilayerGraph::GetMulGraphClock().Count("Space efficient method ends!"));
    
    delete[] rn;
}


static vector<ui> IntersecOfTwo(const vector<ui>& A, const vector<ui>& B) {
    vector<ui> result;
    size_t i = 0, j = 0;
    while (i < A.size() && j < B.size()) {
        if (A[i] < B[j]) {
            i++;
        } else if (A[i] > B[j]) {
            j++;
        } else {
            // A[i] == B[j]
            result.push_back(A[i]);
            i++;
            j++;
        }
    }
    return result;
}

static void ComputeFatherIntersec(vector<vector<ui>> &vertex_set,bool *local_vis, unordered_map<string,ui>& vector2id, vector<string>& father_set, ui n){
    for(ui i = 0; i < n; i++) local_vis[i] = true;
    vector<ui> pos;
    for(ui i = 0; i < father_set.size(); i++){
        pos.push_back(vector2id[father_set[i]]);
    }
    vector<ui> result = vertex_set[pos[0]];
    for (ui i = 1; i < pos.size(); i++) {
        result = IntersecOfTwo(result, vertex_set[pos[i]]);
        if (result.empty()) break;
    }
    for(ui i = 0; i < result.size(); i++)
        local_vis[result[i]] = false;
}





static void VSHT(MultilayerGraph &mg){
    double mem;
    log_info(MultilayerGraph::GetMulGraphClock().Count("VSHT starts"));
    ui layer = mg.GetLayerNumber();
    ui n = mg.GetN();
    ui num_of_den_sub = 0, num_of_computed_den_sub = 0, iter_time = 0;
    unordered_set<string> visited;
    unordered_map<string,ui> vector2id;
    unordered_map<string,vector<string>> father_set;

    bool* global_vis = new bool[mg.GetN()]();
    bool* local_vis = new bool[mg.GetN()]();

    vector<vector<ui>> veck_set;
    vector<vector<ui>> vertex_set;
    vector<ui> veck_set_to_size;

    ReorienNetwork *rn = new ReorienNetwork[layer];
    for(ui i = 0; i < layer; i++){
        rn[i].initial(mg,i);
    }

    vector<ui> start_veck(layer,0);
    veck_set.push_back(start_veck);
    vector<ui> vertex_tmp;
    for(ui i = 0; i < mg.GetN(); i++)
        vertex_tmp.push_back(i);
    
    // record the different dense subgraphs
    vertex_set.push_back(vertex_tmp);
    veck_set_to_size.push_back(mg.GetN());
    num_of_den_sub +=1;

    queue<vector<ui>> que;
    for(ui i = 0; i < layer; i++){
        vector<ui> tmp(start_veck);
        tmp[i]++;
        que.push(tmp);
        visited.insert(VecToString(tmp)); 
    }

    while(!que.empty()){
        vector<ui> veck = que.front();
        que.pop();

        for(ui i = 0; i < layer; i++){
            rn[i].use_graph(global_vis);
        }
        string str_veck = VecToString(veck);
        if(father_set[str_veck].size() > 0){
            ComputeFatherIntersec(vertex_set,local_vis,vector2id,father_set[str_veck],mg.GetN());
            for(ui i = 0; i < layer; i++)
                rn[i].use_graph(local_vis);
        }
        else{
            memset(local_vis,0,sizeof(bool)*n);
            bool *cur_vis;
            ui *core;
            for(ui i = 0; i < layer; i++){
                cur_vis = rn[i].get_cur_vis();
                core = rn[i].get_core();
                for(ui j = 0; j < n; j++){
                    if(local_vis[j]) continue;
                    if(cur_vis[j] || core[j] < veck[i]) local_vis[j] = true;
                }
            }
            for(ui i = 0; i < layer; i++)
                rn[i].use_graph(local_vis);
        }

        // compute the corresponding dense subgraph
        vector<ui> kds = KDenseSubNaive(rn,veck,layer,iter_time);
        num_of_computed_den_sub++;
        if(kds.size() != 0){

            num_of_den_sub++;
            vector<ui> vertex_tmp;
            vertex_tmp.reserve(kds.size());
            for(auto x : kds) vertex_tmp.push_back(x);
            sort(vertex_tmp.begin(),vertex_tmp.end());
            vertex_set.push_back(vertex_tmp);
            veck_set.push_back(veck);
            veck_set_to_size.push_back(kds.size());

            
            vector2id[str_veck] = vertex_set.size()-1;
            for(ui i = 0; i < layer; i++){
                vector<ui> tmp(veck);
                tmp[i]++;
                string str_tmp = VecToString(tmp);
                father_set[str_tmp].push_back(str_veck);
                if(visited.find(str_tmp) == visited.end())  
                {
                    que.push(tmp);
                    visited.insert(str_tmp); 
                }
                    
            }
        }

    }

    delete[] global_vis;
    delete[] local_vis;

    //output all dense subgraph
    // string ss = mg.GetDataPath() + "dense_subgraph_output_distinct.txt";
    // auto out = ofstream(ss);

    // for(ui i = 0; i < veck_set.size(); i++){
    //     vector<ui> dense_vector = veck_set[i];
    //     vector<ui> node_list = vertex_set[i];
    //     std::string vector_str = vectorToTupleString(dense_vector);
    //     std::string node_str = vectorToCommaString(node_list);

    //     out << vector_str << "\t" << veck_set_to_size[i] << "\t" << node_str << "\n";
    // }


    // for(ui i = 0; i < veck_set.size(); i++){
    //     out << "[";
    //     for(ui j = 0; j < layer; j++){
    //         if(j == layer-1)
    //             out << veck_set[i][j];
    //         else
    //             out << veck_set[i][j] << ",";
    //     }     
    //     out << "]: " << veck_set_to_size[i] << endl;
    //     // for(ui j = 0; j < vertex_set[i].size(); j++)
    //     // {
    //     //     if(j == vertex_set[i].size()-1)
    //     //         out << vertex_set[i][j] << endl;
    //     //     else
    //     //         out << vertex_set[i][j] <<" ";
    //     // }
    // }
    // out.close();

    double total_flow_time = 0, total_update_graph_time = 0;
    for(ui i = 0; i < layer; i++){
        total_flow_time += rn[i].get_network_flow_time();
        total_update_graph_time += rn[i].get_update_graph_time();
    }

    log_info(MultilayerGraph::GetMulGraphClock().Count("The number of dense sub: %u, the number of computed dense sub: %u",num_of_den_sub,num_of_computed_den_sub));
    mem = GetPeakRSSInMB();
    log_info(MultilayerGraph::GetMulGraphClock().Count("Memory: %lf MB, the number of times SLBF is called: %u",mem,iter_time));
    log_info(MultilayerGraph::GetMulGraphClock().Count("Network flow time: %lf s, update graph time: %lf s, iter times: %u",total_flow_time,total_update_graph_time,iter_time));
    log_info(MultilayerGraph::GetMulGraphClock().Count("Naive method ends!"));
    delete[] rn;
}

