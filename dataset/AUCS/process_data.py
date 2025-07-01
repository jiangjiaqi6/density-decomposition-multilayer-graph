import pandas as pd

def process_multilayer_graph(input_file):
    """
    读取单个 txt 文件中的多层图数据，进行节点重新编号和层划分。
    
    :param input_file: 包含多层图数据的 txt 文件路径。
    :return: 一个包含重新编号节点和层划分的字典，以及重新编号的节点映射。
    """
    # 读取数据
    data = pd.read_csv(input_file, header=None, names=["Node1", "Node2", "Layer"])
    
    # 获取唯一节点并重新编号
    unique_nodes = pd.unique(data[["Node1", "Node2"]].values.ravel("K"))
    node_mapping = {node: idx for idx, node in enumerate(unique_nodes, start=1)}
    
    # 使用映射替换原数据中的节点
    data["Node1"] = data["Node1"].map(node_mapping)
    data["Node2"] = data["Node2"].map(node_mapping)
    
    # 按层划分数据
    layers = data.groupby("Layer")
    layer_files = {}
    for layer, layer_data in layers:
        layer_file = f"{layer}_layer.txt"
        layer_data[["Node1", "Node2"]].to_csv(layer_file, sep="\t", index=False, header=False)
        layer_files[layer] = layer_file

    return layer_files, node_mapping


# 示例用法
input_file = "aucs_edgelist.txt"  # 包含原始示例数据的文件
layer_files, node_mapping = process_multilayer_graph(input_file)

print("各层文件生成完成：")
for layer, file in layer_files.items():
    print(f"Layer {layer}: {file}")

print("\n节点映射：")
print(node_mapping)
