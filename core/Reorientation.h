/*
* File: Reorientation.h
* Description: 
    Defines the necessary classes and methods for performing network flow operations
    on a reoriented network, aiming to identify dense subgraphs with an integer density
    of k in a single graph.
    This file declares the relevant variables and function interfaces for network flow computations.
* This file is part of the multilayer graph project.
* Author: Jiaqi Jiang
* Created Date: 2025-01-17
* Last Modified: []
*/

#ifndef REORIENTATION_H
#define REORIENTATION_H
#include "../graph/MultilayerGraph.h"

struct ReorienNetworkSnapshot {
    ui *deg_;
    bool *cur_vis_;
    di_edge *edges_;
};



class ReorienNetwork{

public:
    ReorienNetwork();
    ~ReorienNetwork();
    void initial(MultilayerGraph & mg, ui dim);
    void bottom_up_density_decom();
    void top_down_density_decom();
    void divide_and_conquer_density_decom();
    void divide_and_conquer_density_decom_plus();

    void check_density(ui int_thres);


    void bottom_up_density_decom_incremental(ReorienNetwork *rn, ui index_sub, ui index, vector<ui> veck, ui &num_densub, ui &num_computed_densub, 
        ui &num_of_call,
        vector<vector<ui>>& res_vec, vector<ui>& res_subsize, vector<vector<bool>>& filter_sub, vector<vector<ui>>& filter_veck);

    void divide_and_conquer_density_decom_incremental(ReorienNetwork *rn, ui index_sub, ui index, vector<ui> veck, ui &num_densub, 
        vector<vector<ui>>& res_vec, vector<ui>& res_subsize, vector<vector<bool>>& filter_sub, vector<vector<ui>>& filter_veck);

    void divide_and_conquer_density_decom_incremental_space_efficient(ReorienNetwork *rn, ui index_sub, ui index, ui layer, ui & ub, vector<ui> veck, ui &num_densub, ui &computed_densub,ui &num_of_iter, 
        vector<vector<ui>>& res_vec, vector<ui>& res_subsize, vector<vector<bool>>& filter_sub, vector<vector<ui>>& filter_veck, vector<vector<ui>>& edge_set);


    bool ReTestBoUp(ui tem);
    void ReTestDivAndCon(ui Z_rank);
    void ReTestDivAndCon_incremental(ReorienNetwork *rn, ui index, vector<ui> veck, ui& num_densub, ui &num_of_iter, ui Z_rank);

    void update_density_of_each_vertex_basic();
    void update_density_div_and_con(ui X_rank, ui Y_rank, ui Z_rank);
    void update_density_div_and_con_incremental(Set& X, Set& Y, Set& tmp_Z);
    void peel_k_core_multilayer(ReorienNetwork *rn, vector<ui> veck, ui layer);
    void compute_each_layer_sub_edge(ReorienNetwork *rn, vector<ui>& veck, ui layer, bool* vis);

    void update_graph(bool *vis, ui k);
    void update_graph_incremental(ui *val, ui c, Set& Z);
    void update_kdense_subgraph_incremental(ReorienNetwork *rn, ui index, vector<ui> veck, ui X_rank, ui Y_rank, ui Z_rank, Set& tmp_Z);
    void update_kdense_subgraph_2dim_incremental(ReorienNetwork *rn, ui index, vector<ui> veck, ui X_rank, ui Y_rank, ui Z_rank, Set& tmp_Z);
    void update_kdense_subgraph_high_dim_incremental(ReorienNetwork *rn, ui index, ui &num_of_iter, vector<ui> veck, ui X_rank, ui Y_rank, ui Z_rank, Set& Z);

    void update_kdense_subgraph_high_dim_bottom_up(ReorienNetwork *rn, ui index, vector<ui> veck, ui Z_rank, ui &num_of_call, Set& Z);


    void use_graph(bool *vis_info);
    void adjust_deg();
    void adjust_deg_incremental(vector<ui> varVertex, bool* vis, Set& Z);
    void peel_prune_vertex(ui k);
    void peel_k_core_singlelayer(ui v);
    bool compute_lower_bound(ui value, vector<ui> veck);
    bool compute_R_density(bool *vis_info, ui k);

	ui get_max_d();
    ui get_max_r();
    ui get_max_pa(){return max_pa;}
    ui core_decomposition_incremental();
    void set_graph_name(string ss_){graph_name = ss_;}
    string get_graph_name(){return graph_name;}

	
	bool DinicBFS(std::unordered_set<ui>& nodes_to_decrease);
	bool DinicDFS(ui x);

    bool DinicBFSDivAndCon(Set& X, Set& Y, ui &high_num);
    void recursive_function(ui X_rank, ui Y_rank);
    void recursive_function_incremental(ReorienNetwork *rn, ui index_sub, ui index, vector<ui> veck, ui &num_densub, ui &num_of_iter, ui X_rank, ui Y_rank);

    void find_k_dense_sub(vector<ui>& isContainVer, ui tem, ui layer_);
    void find_k_dense_sub_plus(vector<ui>& isContainVer, vector<ui>& change, ui tem, ui layer_);
    void find_k_dense_sub_incremental(ui tem, Set& X, ui X_rank, Set& Y, ui Y_rank, Set& Z);

    ui get_n() {return n;}

    bool* get_cur_vis(){return cur_vis;}
    ui* get_subgraph_r(){return r;}
    ui* get_subgraph_d(){return deg;}
    ui* get_core(){return core;}
    ui* get_abs_deg(){return abs_deg;}
    ui* get_pstart(){return pstart;}
    ui* get_pend(){return pend;}
    ui* get_edgeList(){return edgeList;}
    ui* get_eid(){return eid;}
    di_edge* get_edges(){return edges;}


    double get_network_flow_time(){return network_flow_time;}
    double get_update_graph_time(){return update_graph_time;}
    double get_reconfig_info_time(){return reconfig_info_time;}
    double get_peel_core_time(){return peel_core_time;}
    ui get_iter_time(){return iter_time;}

    void dispaly_density_every_layer(string ss);
    void extract_R_sub(ui k, bool* vis);

    void filter_subgraphs(vector<vector<bool>>& filter_sub, vector<vector<ui>>& filter_vec, vector<ui>& initial_vec);

    void takeSnapshot(ReorienNetworkSnapshot& snapshot);
    void restoreSnapshot(ReorienNetworkSnapshot& snapshot);

    double com_global_cluster_coefficient();
    double compute_conductance(const vector<int>& subgraph_nodes);

private:
    ui dim; // external parameter of dimension
    ui n;
    ui m;
    ui pa; //pseudoarboricity
	ui max_pa; //maximum pseudoarboricity
	ui app_pa; //approximate pseudoarboricity
    string graph_name;

    ui iter_time;  //compute the vector dense subgraph in multilayer graph

    Set* R{nullptr}; //Subgraph corresponding to R[i] with density no less than i

	//integer density of each vertex
    ui *r{nullptr};

    //core value
    ui *core{nullptr};
    
    //pseudoarboricity
	ui *p{nullptr};
	ui test_value;

	//flow 
	ui *dist{nullptr};
	ui *cur{nullptr};
    ui *high_outdegree{nullptr};
    bool *cur_vis{nullptr};
    bool *cur_vis_copy{nullptr};


    //graph
    ui *pstart; //start positions of neighbors of vertices in the array "edges"
	ui *pend; //end positions of neighbors of vertices in the array "edges"
	ui *edgeList; //concatenation of neighbors of all vertices
	ui *eid; // id of each edge
	ui *deg; // the out degree of each vertex
    ui *abs_deg;  //the degree of each vertex
	di_edge *edges;


    double peel_core_time;
    double network_flow_time;
    double update_graph_time;
    double reconfig_info_time;


	std::vector<ui> res;
	std::unordered_set<ui> isContainDense; // used in recording integer_dense_sub

};

#endif