import csv

# 输入和输出文件路径
input_file = 'coauthor.txt'
output_file = 'output.csv'

# 存储符合条件的数据
data = []

# 读取 txt 文件
with open(input_file, 'r') as f:
    for line in f:
        # 去除空格和换行符，并按制表符或空格分割
        parts = line.strip().split()
        if len(parts) == 2:
            num1, num2 = int(parts[0]), int(parts[1])
            # 只有第一个元素小于第二个元素时，才保留
            if num1 < num2:
                data.append([num1, num2])

# 将数据写入 CSV 文件
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['src', 'des'])  # 添加表头
    writer.writerows(data)

print("已成功保存为 CSV 文件:", output_file)
