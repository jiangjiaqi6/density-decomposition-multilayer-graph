def read_txt_file(filename):
    """读取txt文件，将内容解析成set集合，键值对存为字符串"""
    data_set = set()
    with open(filename, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if line:  # 确保非空
                data_set.add(line)
    return data_set

def compute_sets(file1, file2):
    """计算两个文件的交集和差集"""
    set1 = read_txt_file(file1)
    set2 = read_txt_file(file2)

    # 计算交集
    intersection = set1 & set2

    # 计算差集
    difference1 = set1 - set2  # file1 独有
    difference2 = set2 - set1  # file2 独有

    return intersection, difference1, difference2

def write_results(filename, data):
    """将结果写入到文件"""
    with open(filename, 'w', encoding='utf-8') as f:
        for item in sorted(data):  # 排序后写入
            f.write(item + '\n')

if __name__ == "__main__":
    # 输入两个文件路径
    file1 = "dense_subgraph_res_naive.txt"
    file2 = "dense_subgraph_res_divandcon.txt"

    # 计算交集和差集
    intersection, difference1, difference2 = compute_sets(file1, file2)

    # 输出结果
    print(f"交集元素个数: {len(intersection)}")
    print(f"file1 独有元素个数: {len(difference1)}")
    print(f"file2 独有元素个数: {len(difference2)}")

    # 可选：写入文件
    write_results("intersection.txt", intersection)
    write_results("difference1.txt", difference1)
    write_results("difference2.txt", difference2)

    print("结果已写入 intersection.txt, difference1.txt, difference2.txt")
