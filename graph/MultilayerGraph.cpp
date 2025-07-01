/*
* File: MultilayerGraph.h
* Description: 
*   Defines the MultilayerGraph class, which represents a collection of multiple graph layers.
*   This calss manages individual Graph instances as its layers and provides methods to query and manipulate 
*   the entire multilayer graph structure.
* Author: Jiaqi Jiang
* Created Date: 2025-01-14
* Last Modified: []
*/
#include "MultilayerGraph.h"

// MultilayerGraph::MultilayerGraph():mul_graph_clock("multilayer_graph"){}

Clock MultilayerGraph::mul_graph_clock("multilayer_graph");

MultilayerGraph::~MultilayerGraph(){
    delete[] graph_layers;
    delete[] order;
}

void MultilayerGraph::LoadFromFile(const string &input_path){
    vector<string> graph_files;
    vector<vector<edge>> edge_buf;
    vector<uint> edge_size;
    input_data_path = input_path + "/";
    map_file = input_data_path + "vtx_map.txt";
    unordered_map<ll_ui, ui> vtx2id;
    std::basic_ofstream<char> map_file_out;
    GetGraphFile(input_data_path, graph_files);
    n_layers = graph_files.size();
    graph_layers = new Graph[n_layers];
    
    edge_buf.resize(n_layers);
    edge_size.resize(n_layers);
    map_file_out = ofstream(map_file);
    for(ui i = 0; i < n_layers; i++){
        edge_size[i] = LoadLayer(input_data_path+graph_files[i], edge_buf[i], vtx2id, map_file_out);
    }
    map_file_out.close();
    n = vtx2id.size();
    for(ui i = 0; i < n_layers; i++){
        graph_layers[i].BuildFromEdgeLst(edge_buf[i],n,edge_size[i],graph_files[i]);
    }

}

ui MultilayerGraph::LoadLayer(const string &graph_file, vector<edge> &edge_buf, unordered_map<ll_ui, ui> &vtx2id,
                                std::basic_ofstream<char> &map_file_out){
    log_info(mul_graph_clock.Count("graph_file: %s",graph_file.c_str()));
    uint uid, vid;
    uint edge_buf_size, num_of_vtx, num_of_edge;
    ll_ui u, v;
    edge_buf.clear();
    std::basic_ifstream<char> graph_in;
    graph_in = ifstream(graph_file);
    num_of_vtx = (uint) vtx2id.size();
    FILE *file = fopen(graph_file.c_str(),"r");
    char line[200];


    size_t found = graph_file.find("dataset"); // 找到 "dataset" 的位置
    // 从 "dataset" 后的部分开始提取
    std::string result = graph_file.substr(found + 7); // 7 是 "dataset" 的长度

    // 去掉路径末尾的文件扩展名 ".txt"
    if (result.find(".txt") != std::string::npos) {
        result = result.substr(0, result.find(".txt"));
    }
    else if (result.find(".csv") != std::string::npos) {
        result = result.substr(0, result.find(".csv"));
    }

    string ss = "/home/jjq/research/NetworkFlow/multilayer_graph/density_decomposition/dataset"+result+"_newid.csv";
    auto out = ofstream(ss);
    cout<<ss<<endl;

    map_file_out<<"original_vtx,now_vtx"<<endl;
    

    while(fgets(line,200,file)){
		ui i = 0;
        ui a = 0;
		while (line[i] < '0' || line[i]>'9') i++;
		while (line[i] >= '0' && line[i] <= '9') a = a * 10 + line[i] - '0', i++;
        ui u = a;
        ui b = 0;
		while (line[i] < '0' || line[i]>'9') i++;
		while (line[i] >= '0' && line[i] <= '9') b = b * 10 + line[i] - '0', i++;
        ui v = b;
        if (u != v) {
            auto iter1 = vtx2id.find(u);
            if (iter1 != vtx2id.end()) {
                uid = iter1->second;
            } else {
                uid = num_of_vtx++;
                vtx2id.emplace(u, uid);
                map_file_out << u << "," << uid << endl;
            }

            auto iter2 = vtx2id.find(v);
            if (iter2 != vtx2id.end()) {
                vid = iter2->second;
            } else {
                vid = num_of_vtx++;
                vtx2id.emplace(v, vid);
                map_file_out << v << "," << vid << endl;
            }

            out<<uid<<" "<<vid<<endl;
            {
                edge_buf.pb(edge(uid, vid));
                edge_buf.pb(edge(vid, uid));
            }
        }
    }


    out.close();
    graph_in.close();
    return edge_buf.size();
}

void MultilayerGraph::GetGraphFile(const string &graph_path, vector<string> &graph_files){
    struct stat buffer{};
    string &&conf_file = graph_path + "mlg.conf";

    if (stat((conf_file).c_str(), &buffer)) {
        // no configure file provided
        DIR *dir;
        dirent *env;
        if ((dir = opendir(graph_path.c_str()))) {
            while ((env = readdir(dir))) {
                graph_files.emplace_back(env->d_name);
            }
            closedir(dir);
        }

    } else {
        // configure file provided
        string line;
        auto fin = ifstream(conf_file);
        while (fin.peek() != EOF) {
            getline(fin, line);
            if (!line.empty()) {
                graph_files.emplace_back(line);
            }
        }
        fin.close();
    }    
}


void MultilayerGraph::PrintStatistics(){
    ui max_m = 0, sum_m = 0, m;
    for (ui i = 0; i < n_layers; i++){
        m = graph_layers[i].get_m() >> 1;
        max_m = max(max_m,m);
        sum_m += m;
    }
    // std::cout << "|L| = " << n_layers << ", |V| = " << n << ", |E| = " << sum_m << ", max_|E_i| = " << max_m << std::endl;
    log_info(MultilayerGraph::GetMulGraphClock().Count("|L| = %u, |V| = %u, |E| = %u, max_|E_i| = %u",n_layers,n,sum_m,max_m));
}