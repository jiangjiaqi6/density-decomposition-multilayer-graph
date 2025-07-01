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

#ifndef MULTILAYERGRAPH_H
#define MULTILAYERGRAPH_H

#include "Graph.h"

class MultilayerGraph{
public:
    MultilayerGraph() = default;
    ~MultilayerGraph();

    void LoadFromFile(const string &input_path);
    void PrintStatistics();

    void LoadId2VtxMap(ll_ui * id2vtx);
    void LoadVtx2IdMap(unordered_map<ll_ui, uint> &vtx2id);

    static MultilayerGraph* Load(const string &file);

    [[nodiscard]] inline ui GetLayerNumber() const {
        return n_layers;
    }

    [[nodiscard]] inline ui GetN() const {
        return n;
    }

    [[nodiscard]] inline string GetDataPath() const {
        return input_data_path;
    }

    inline Graph &GetGraph(ui i) {
        return graph_layers[i];
    }

    // for test only
    [[nodiscard]] inline ui* GetOrder() const {
        return order;
    }
    static Clock& GetMulGraphClock() {
        return mul_graph_clock;
    }

private:
    Graph *graph_layers{nullptr};
    ui *order{nullptr};

    ui n_layers{0};
    ui n{0};
    string input_data_path;

    string map_file;
    static Clock mul_graph_clock;

    static void GetGraphFile(const string &graph_path, vector<string> &graph_files);
    // static ui LoadLayer(const string &graph_file, edge *&edge_buf, unordered_map<ll_ui, ui> &vtx2id,
    //                      std::basic_ofstream<char> &map_file_out);

    static ui LoadLayer(const string &graph_file, vector<edge> &edge_buf, unordered_map<ll_ui, ui> &vtx2id,
                         std::basic_ofstream<char> &map_file_out);

};

#endif