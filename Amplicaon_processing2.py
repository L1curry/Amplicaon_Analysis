#!/usr/bin/env python3
# 导入所需模块
import os                           # 处理文件路径
import sys                          # 系统相关操作
import shutil                       # 文件与目录操作
import subprocess                   # 调用外部命令
import pandas as pd                 # 数据处理
import logging                      # 日志记录
import argparse                     # 命令行参数解析
from Bio import SeqIO

# 设置日志记录，日志写入文件 "amplicon_processing.log"，记录级别为 INFO
logging.basicConfig(
    filename="amplicon_processing.log",   # 日志文件名
    level=logging.INFO,                     # 记录级别
    format="%(asctime)s - %(levelname)s - %(message)s"  # 日志格式
)

# 定义一个通用的输入验证函数，循环提示用户输入直到满足要求或输入“exit”退出程序
def get_valid_input(prompt, validator=None, error_msg="输入不符合要求，请重新输入。"):
    """
    循环提示用户输入，直到满足 validator 条件或用户输入“exit”退出程序。
    
    参数:
        prompt (str): 提示信息
        validator (function): 一个函数，输入字符串，返回 True/False 表示是否有效
        error_msg (str): 当输入无效时显示的错误提示
    返回:
        str: 有效的输入
    """
    while True:
        value = input(prompt).strip()  # 提示输入并去除两端空白字符
        if value.lower() == "exit":   # 如果用户输入 "exit" 则退出程序
            sys.exit("用户选择结束运行。")
        if validator is None or validator(value):  # 如果不需要验证或验证通过则返回输入
            return value
        else:
            print(error_msg)  # 提示错误信息

# 检查指定工具是否在系统路径中可用，如不可用则提示用户输入工具所在目录
def check_tool(tool_name):
    """
    检查指定工具是否在系统路径中可用。
    如果不可用，提示用户输入安装路径。
    
    参数:
        tool_name (str): 工具名称
    返回:
        str: 工具的完整路径
    """
    tool_path = shutil.which(tool_name)  # 在系统路径中查找工具
    if tool_path:
        logging.info(f"{tool_name} found at {tool_path}")  # 记录日志
        return tool_path
    else:
        # 循环提示用户输入工具所在目录
        while True:
            dir_input = input(f"{tool_name} not found. Please input its installation directory (absolute path): ")
            if dir_input.strip().lower() == "exit":
                sys.exit("用户选择结束运行。")
            # 将输入路径转换为绝对路径，并移除单引号
            dir_path = os.path.abspath(dir_input.strip("'"))
            tool_full_path = os.path.join(dir_path, tool_name)  # 拼接完整工具路径
            # 检查该路径是否为可执行文件
            if os.path.isfile(tool_full_path) and os.access(tool_full_path, os.X_OK):
                logging.info(f"{tool_name} found at {tool_full_path}")
                return tool_full_path
            else:
                print("Tool not found or not executable. Please try again.")

# 定义运行外部命令的函数，并将命令输出、错误记录到日志中
def run_command(cmd, description=""):
    """
    执行外部命令并检查是否成功。
    
    参数:
        cmd (list): 命令及其参数列表
        description (str): 命令的描述，用于日志记录
    """
    logging.info(f"Running: {description}")  # 记录描述
    logging.info(f"Command: {' '.join(map(str, cmd))}")  # 记录完整命令
    try:
        # 调用外部命令，并检查是否执行成功
        result = subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        logging.info(f"{description} completed successfully.")  # 记录成功信息
        logging.debug(f"Output: {result.stdout}")  # 记录输出（DEBUG级别）
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in {description}: {e.stderr}")  # 记录错误信息
        raise

# 解析命令行参数，要求输入输入目录、元数据文件、输出目录及线程数
def parse_args():
    """
    解析命令行参数。
    
    返回:
        argparse.Namespace: 解析后的参数对象
    """
    parser = argparse.ArgumentParser(
        description="扩增子测序数据处理脚本",
        epilog="""
        该脚本用于处理扩增子测序数据并生成OTU表。处理流程包括：
        1. 样本解复用（双端模式）
        2. 合并解复用后的双端测序文件
        3. 质量过滤
        4. 去重复
        5. 聚类
        6. 嵌合体检测
        7. OTU表生成
        8. 可选的SINTAX分类注释
        使用示例: python script.py -i input_dir -m metadata.tsv -o output_dir -t 4
        """
    )
    # 添加输入目录参数
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="包含所有测序文件的输入目录（绝对或相对路径）。"
    )
    # 添加元数据文件参数
    parser.add_argument(
        "-m", "--metadata_file",
        required=True,
        help="元数据文件的路径（TSV格式，包含run_id、sample_id、引物和文件信息）。"
    )
    # 添加输出目录参数
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="结果输出目录（将自动创建）。"
    )
    # 添加线程数参数
    parser.add_argument(
        "-t", "--threads",
        type=int,
        required=True,
        help="使用的线程/核心数（整数，例如4）。"
    )
    return parser.parse_args()

########################################
# 新功能 11：OTU序列重新标记与进一步聚类（可选）
########################################
def relabel_and_cluster(tools, threads, output_dir):
    print("\n========== 11. OTU序列重新标记与进一步聚类 ==========")
    # 定义7-OTU目录路径，并进入该目录
    otu_dir = os.path.join(output_dir, "7-OTU")
    orig_dir = os.getcwd()                      # 保存当前工作目录
    os.chdir(otu_dir)                           # 切换到7-OTU目录
    # 非嵌合体文件路径（6-chimera目录下的 otus_nochim.fasta）
    nonchim_file = os.path.join(os.path.dirname(otu_dir), "6-chimera", "otus_nochim.fasta")
    # 定义重命名后临时输出文件名
    temp_otus = "otus.temp.fasta"
    # 构造 vsearch 命令，对非嵌合体序列进行重命名（加上 OTU_ 前缀）
    cmd_relabel = [
        tools["vsearch"],
        "--fastx_filter", nonchim_file,
        "--sizein", "--sizeout",
        "--fasta_width", "0",
        "--relabel", "OTU_",
        "--threads", str(threads),
        "--fastaout", temp_otus
    ]
    # 运行命令并记录日志
    run_command(cmd_relabel, "OTU重命名")
    
    # 提示用户是否进行第二次 UPARSE3 聚类
    choice = ""
    while choice not in ["1", "2"]:
        choice = input("是否进行第二次UPARSE3聚类？1. Yes; 2. No: ").strip()
    
    if choice == "1":
        # 如果选择进行聚类，提示用户输入相似性阈值
        up_identity = ""
        while True:
            up_identity = input("请输入UPARSE3聚类使用的相似性阈值（0~1，如0.97）：").strip()
            try:
                val = float(up_identity)
                if 0 <= val <= 1:
                    break
            except:
                pass
            print("输入不合法，请重新输入！")
        # 定义最终聚类结果文件名
        centroids = "otus.fasta"
        # 构造 vsearch 聚类命令，使用临时重命名文件进行聚类
        cmd_cluster = [
            tools["vsearch"],
            "--cluster_size", temp_otus,
            "--id", up_identity,
            "--strand", "plus",
            "--threads", str(threads),
            "--sizein", "--sizeout",
            "--fasta_width", "0",
            "--centroids", centroids
        ]
        run_command(cmd_cluster, "第二次UPARSE3聚类")
        
        # 定义所有去重复序列文件路径（4-dereplicate目录下的 all_samples_derep.fasta）
        all_derep_file = os.path.join(os.path.dirname(otu_dir), "4-dereplicate", "all_samples_derep.fasta")
        # 构造 vsearch 映射命令，将所有去重复序列映射到聚类结果上生成 OTU 表（otutab.txt）
        cmd_map = [
            tools["vsearch"],
            "--usearch_global", all_derep_file,
            "--db", centroids,
            "--id", up_identity,
            "--strand", "plus",
            "--threads", str(threads),
            "--sizein", "--sizeout",
            "--fasta_width", "0",
            "--qmask", "none",
            "--dbmask", "none",
            "--otutabout", "otutab.txt"
        ]
        run_command(cmd_map, "OTU映射")
        # 删除临时文件
        if os.path.exists(temp_otus):
            os.remove(temp_otus)
        print("第二轮聚类完成，结果在 otus.fasta 和 otutab.txt 中。")
    elif choice == "2": # 功能2可能会报错
        # 如果不进行第二次聚类，则直接映射，构造命令，使用默认相似性阈值 0.97
        cmd_map2 = [
            tools["vsearch"],
            "--usearch_global", os.path.join(os.path.dirname(otu_dir), "4-dereplicate", "all_samples_derep.fasta"),
            "--db", temp_otus,
            "--id", "0.97",
            "--strand", "plus",
            "--threads", str(threads),
            "--sizein", "--sizeout",
            "--fasta_width", "0",
            "--qmask", "none",
            "--dbmask", "none",
            "--otutabout", "otutab.temp.txt"
        ]
        run_command(cmd_map2, "OTU映射（直接映射）")
        # 将临时文件重命名为最终文件
        os.rename("otus.temp.fasta", "otus.fasta")
        os.rename("otutab.temp.txt", "otutab.txt")
        print("直接映射完成，结果在 otus.fasta 和 otutab.txt 中。")
    
    os.chdir(orig_dir)  # 返回原工作目录

########################################
# 新功能 12：稀释曲线绘制
########################################
def plot_rarefaction_curve(tools, threads, output_dir):
    print("\n========== 12. 稀释曲线绘制 ==========")
    # 定义7-OTU目录路径并进入该目录
    otu_dir = os.path.join(output_dir, "7-OTU")
    orig_dir = os.getcwd()      # 保存当前目录
    os.chdir(otu_dir)           # 切换目录
    try:
        cmd = "awk '{print $2}' ../../barcoding_corrected.txt > ../../id.sample" # 生成样本ID文件，注意文件存放位置
        subprocess.run(cmd, shell=True, check=True)
        pass
    except:
        print("Error！barcoding_corrected.txt未找到，或其它原因")
    # 构造 usearch 命令，计算稀释曲线，将结果保存到 rare.txt
    cmd_rarefy = [
        "Rarefy_OTUtab.R", "otutab.txt",
        "../../id.sample"
    ]
    run_command(cmd_rarefy, "计算稀释曲线")
    # 假定样本ID存放在上一级目录的 id.sample 文件中
    id_sample_file = os.path.join(os.path.dirname(output_dir), "id.sample")
    sample_ids = []
    if os.path.exists(id_sample_file):
        with open(id_sample_file, "r") as f:
            # 读取每一行作为样本ID
            sample_ids = [line.strip() for line in f if line.strip()]
    else:
        print(f"未找到{id_sample_file}，无法绘制稀释曲线。")
        os.chdir(orig_dir)
        return
    # 尝试读取稀释曲线数据文件 rare.txt，文件以制表符分隔
    try:
        rare_df = pd.read_csv("rare.txt", sep="\t")
    except Exception as e:
        print("读取稀释曲线数据失败：", e)
        os.chdir(orig_dir)
        return
    # 整理数据，每个样本取出 richness 和对应OTU数，并合并为一张表
    all_data = []
    for sid in sample_ids:
        if sid not in rare_df.columns:
            print(f"警告：样本 {sid} 不在 rare.txt 中！")
            continue
        df_sample = rare_df[['richness', sid]].copy()  # 取 richness 和该样本列
        df_sample.columns = ['X', 'num_otus']          # 重命名列名
        df_sample['sample'] = sid                      # 添加样本标识
        all_data.append(df_sample)
    if not all_data:
        print("没有找到有效的样本数据，退出绘图。")
        os.chdir(orig_dir)
        return
    # 合并所有样本数据
    all_df = pd.concat(all_data, ignore_index=True)
    # 导入 matplotlib 绘图
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 6))    # 设置图形尺寸
    # 按样本分组绘制曲线，每个样本一条曲线，点标记为圆圈
    for sid, group in all_df.groupby('sample'):
        plt.plot(group['X'], group['num_otus'], marker='o', label=sid)
    plt.xlabel("X")               # 设置X轴标签
    plt.ylabel("num_otus")        # 设置Y轴标签
    plt.title("Rarefaction curve")  # 设置图形标题
    plt.legend(title="sample", loc="upper left")  # 添加图例
    plt.ylim(bottom=0)            # Y轴最小值为0
    plt.tight_layout()            # 自动调整布局
    plt.savefig("rarefaction_curve.pdf")  # 保存图形为 PDF
    plt.close()                   # 关闭绘图
    print("稀释曲线绘图完成，文件保存为 rarefaction_curve.pdf")
    os.chdir(orig_dir)            # 返回原目录

########################################
# 新功能 13：OTU低丰度数据过滤（可选）
########################################

def filter_low_abundance_otus(tools, threads, output_dir):
    """
    OTU低丰度数据过滤替代实现：
    1. 切换到 output_dir/7-OTU 目录；
    2. 交互提示是否进行低丰度过滤，若选择“否”则直接返回；
    3. 交互获取最低计数阈值（filter_min_count）和最低频率阈值（filter_min_freq）；
    4. 读取 otutab.txt 文件，针对每个样本（列）过滤低丰度数据（count < filter_min_count 或相对频率 < filter_min_freq，则置为0），结果写入 otutab.filter.txt；
    5. 从过滤后的 OTU 表中提取 OTU ID（第一列，每个 ID 后附加分号），保存到临时文件 list.t1；
    6. 读取 otus.fasta 文件中的所有序列，筛选序列 ID 中包含 list.t1 中任一OTU ID的序列，保存到 otus.filter.fasta；
    7. 将未通过过滤的 OTU ID 写入 list.filter，并删除临时文件 list.t1；
    8. 返回原工作目录。
    """
    print("\n========== 13. OTU低丰度数据过滤 ==========")
    # 进入 7-OTU 目录
    otu_dir = os.path.join(output_dir, "7-OTU")
    orig_dir = os.getcwd()  # 保存当前工作目录
    os.chdir(otu_dir)       # 切换到 7-OTU 目录

    # 交互提示是否进行低丰度过滤
    choice = ""
    while choice not in ["1", "2"]:
        choice = input("是否需要过滤低丰度OTU（counts小于阈值设为0）：1. Yes; 2. No: ").strip()
    if choice == "2":
        print("不进行低丰度过滤。")
        os.chdir(orig_dir)
        return

    # 提示输入最低计数阈值（必须为正整数）
    filter_min_count = None
    while True:
        try:
            filter_min_count = int(input("请输入最低计数阈值（如50）：").strip())
            if filter_min_count > 0:
                break
        except:
            pass
        print("输入不合法，请输入正整数！")

    # 提示输入最低频率阈值（可以是小数）
    filter_min_freq = None
    while True:
        try:
            filter_min_freq = float(input("请输入最低频率阈值（如0.001）：").strip())
            break
        except:
            print("输入不合法，请输入数字！")
    print(f"过滤参数：最低计数 = {filter_min_count}，最低频率 = {filter_min_freq}")

    # -------------------------------
    # Step1. OTU表过滤：利用 pandas 对 otutab.txt 进行低丰度数据过滤
    # -------------------------------
    # 读取 OTU 表，假设以制表符分隔，第一列为 OTU ID
    try:
        otu_df = pd.read_csv("otutab.txt", sep="\t", index_col=0)
    except Exception as e:
        print("读取 otutab.txt 失败：", e)
        os.chdir(orig_dir)
        return

    # 过滤规则：对每个样本（列），若单个OTU计数小于 filter_min_count 或占该样本总计数的比例小于 filter_min_freq，则将该计数置为0
    otu_df_filtered = otu_df.copy()
    for sample in otu_df.columns:
        total = otu_df[sample].sum()
        # 过滤：先判断是否满足最低计数和最低频率要求，否则置为0
        otu_df_filtered[sample] = otu_df[sample].apply(lambda x: x if (x >= filter_min_count and (x / total) >= filter_min_freq) else 0)

    # 将过滤后的 OTU 表写入文件 otutab.filter.txt
    otu_df_filtered.to_csv("otutab.filter.txt", sep="\t")
    print("过滤后的 OTU 表已保存到 otutab.filter.txt")

    # -------------------------------
    # Step2. 提取过滤后的 OTU ID 到临时文件 list.t1
    # -------------------------------
    otu_ids = []
    with open("otutab.filter.txt", "r") as f:
        header = f.readline()  # 跳过表头
        for line in f:
            fields = line.strip().split("\t")
            if fields:
                # 在OTU ID后加上分号，模仿原函数行为
                otu_ids.append(fields[0] + ";")
    with open("list.t1", "w") as f:
        for otu in otu_ids:
            f.write(otu + "\n")
    print("过滤后的 OTU ID 已写入临时文件 list.t1")

    # -------------------------------
    # Step3. 筛选通过过滤的 OTU序列
    # -------------------------------
    # 读取 otus.fasta 文件中的所有序列，存入字典
    all_otus = {}
    try:
        for record in SeqIO.parse("otus.fasta", "fasta"):
            all_otus[record.id] = record
    except Exception as e:
        print("读取 otus.fasta 失败：", e)
        os.chdir(orig_dir)
        return

    # 筛选通过过滤的 OTU：如果序列 ID 中包含临时文件 list.t1 中任一 OTU ID（含分号）则保留
    passed_ids = set()
    for otu in otu_ids:
        for rid in all_otus:
            if otu in rid:
                passed_ids.add(rid)

    # 将通过过滤的 OTU 序列写入 otus.filter.fasta 文件
    with open("otus.filter.fasta", "w") as f:
        for rid in passed_ids:
            SeqIO.write(all_otus[rid], f, "fasta")
    print(f"共有 {len(passed_ids)} 个OTU通过过滤，结果保存在 otus.filter.fasta")

    # -------------------------------
    # Step4. 记录未通过过滤的 OTU ID
    # -------------------------------
    failed_ids = set(all_otus.keys()) - passed_ids
    with open("list.filter", "w") as f:
        for rid in failed_ids:
            f.write(rid + "\n")
    print("未通过过滤的OTU ID保存在 list.filter")

    # 删除临时文件 list.t1
    try:
        os.remove("list.t1")
    except Exception as e:
        print("删除临时文件 list.t1 失败：", e)

    # 返回原工作目录
    os.chdir(orig_dir)

# 原有 parse_args 函数（与前面定义的 parse_args 相同）
def parse_args():
    """
    解析命令行参数。
    
    返回:
        argparse.Namespace: 解析后的参数对象
    """
    parser = argparse.ArgumentParser(
        description="扩增子测序数据处理脚本",
        epilog="""
        该脚本用于处理扩增子测序数据并生成OTU表。处理流程包括：
        1. 样本解复用（双端模式）
        2. 合并解复用后的双端测序文件
        3. 质量过滤
        4. 去重复
        5. 聚类
        6. 嵌合体检测
        7. OTU表生成
        8. 可选的SINTAX分类注释
        使用示例: python script.py -i input_dir -m metadata.tsv -o output_dir -t 4
        """
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="包含所有测序文件的输入目录（绝对或相对路径）。"
    )
    parser.add_argument(
        "-m", "--metadata_file",
        required=True,
        help="元数据文件的路径（TSV格式，包含run_id、sample_id、引物和文件信息）。"
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="结果输出目录（将自动创建）。"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        required=True,
        help="使用的线程/核心数（整数，例如4）。"
    )
    return parser.parse_args()

# 主函数入口
def main():
    # 解析命令行参数
    args = parse_args()
    input_dir = os.path.abspath(args.input_dir)         # 获取输入目录的绝对路径
    metadata_file = os.path.abspath(args.metadata_file)   # 获取元数据文件的绝对路径
    output_dir = os.path.abspath(args.output_dir)         # 获取输出目录的绝对路径
    threads = args.threads                                # 获取线程数

    # 检查输入目录和元数据文件是否存在
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"输入目录 {input_dir} 不存在。")
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"元数据文件 {metadata_file} 不存在。")
    os.makedirs(output_dir, exist_ok=True)                # 创建输出目录（如果不存在）

    # 检查所需工具，并保存各工具的完整路径到字典中
    tools = {
        "cutadapt": check_tool("cutadapt"),
        "vsearch": check_tool("vsearch"),
        "usearch": check_tool("usearch"),
        "seqkit": check_tool("seqkit"),
        "csvtk": check_tool("csvtk")
    }

    # 读取元数据文件（假设无表头，列顺序为：run_id, sample_id, forward_primer, reverse_primer, forward_file, reverse_file）
    logging.info("正在加载元数据...")
    metadata = pd.read_csv(metadata_file, sep='\t', header=None)
    metadata.columns = ['run_id', 'sample_id', 'forward_primer', 'reverse_primer', 'forward_file', 'reverse_file']

    # Step 1: 样本解复用（双端模式）——使用 cutadapt 同时处理正向和反向原始数据
    demultiplex_dir = os.path.join(output_dir, "1-demultiplex")
    os.makedirs(demultiplex_dir, exist_ok=True)           # 创建解复用输出目录
    for _, row in metadata.iterrows():
        sample_id = row['sample_id']                      # 获取样本ID
        barcode_for = row['forward_primer']               # 获取正向引物
        barcode_rev = row['reverse_primer']               # 获取反向引物
        forward_file = os.path.join(input_dir, row['forward_file'])  # 正向文件完整路径
        reverse_file = os.path.join(input_dir, row['reverse_file'])  # 反向文件完整路径
        output_forward = os.path.join(demultiplex_dir, f"{sample_id}.R1.fastq")  # 输出正向文件路径
        output_reverse = os.path.join(demultiplex_dir, f"{sample_id}.R2.fastq")  # 输出反向文件路径
        
        # 检查正反向文件是否存在，不存在则跳过该样本
        if not os.path.exists(forward_file) or not os.path.exists(reverse_file):
            logging.warning(f"跳过 {sample_id} - 正向或反向文件未找到。")
            continue
        
        # 构造 cutadapt 命令，同时处理正向和反向数据
        cmd = [
            tools["cutadapt"],
            "-g", f"^{barcode_for}",
            "-G", f"^{barcode_rev}",
            "-o", output_forward,
            "-p", output_reverse,
            "-j", str(threads),
            "--discard-untrimmed",
            "-e", "0.1",
            forward_file, reverse_file
        ]
        # 运行 cutadapt 命令进行解复用
        run_command(cmd, f"样本解复用 {sample_id}")

    # Step 2: 合并双端测序文件（对解复用后的数据进行合并）
    merge_dir = os.path.join(output_dir, "2-merge")
    os.makedirs(merge_dir, exist_ok=True)                # 创建合并输出目录
    for _, row in metadata.iterrows():
        sample_id = row['sample_id']                      # 获取样本ID
        demux_forward = os.path.join(demultiplex_dir, f"{sample_id}.R1.fastq")  # 解复用正向文件路径
        demux_reverse = os.path.join(demultiplex_dir, f"{sample_id}.R2.fastq")  # 解复用反向文件路径
        merge_out = os.path.join(merge_dir, f"{sample_id}.merged.fastq")        # 合并后输出文件路径
        
        # 如果解复用文件不存在则跳过该样本
        if not os.path.exists(demux_forward) or not os.path.exists(demux_reverse):
            logging.warning(f"跳过 {sample_id} - 解复用文件未找到。")
            continue
        
        # 构造 vsearch 合并命令，将正向和反向文件合并为一个 merged 文件
        cmd = [
            tools["vsearch"], "--fastq_mergepairs", demux_forward,
            "--reverse", demux_reverse,
            "--threads", str(threads),
            "--fastqout", merge_out,
            "--fastq_eeout"
        ]
        run_command(cmd, f"合并样本 {sample_id} 的双端测序文件")

    # Step 3: 质量过滤
    quality_dir = os.path.join(output_dir, "3-quality")
    os.makedirs(quality_dir, exist_ok=True)             # 创建质量过滤输出目录
    # 交互提示选择 PCR 产物长度类型：1 表示在一定范围内，2 表示固定为一个或多个值
    marker_type = get_valid_input(
        "PCR 产物长度类型: 1. 在一定范围内(default 200~400bp); 2. 固定为一个/多个值: ",
        validator=lambda x: x in ["1", "2"],
        error_msg="请输入 1 或 2。"
    )
    if marker_type == "1":
        # 如果选择在一定范围内，提示输入最小和最大长度
        min_length = get_valid_input(
            "请输入目标 reads 的最小长度: ",
            validator=lambda x: x.isdigit() and int(x) > 0,
            error_msg="请输入正整数。"
        )
        max_length = get_valid_input(
            "请输入目标 reads 的最大长度: ",
            validator=lambda x: x.isdigit() and int(x) >= int(min_length),
            error_msg="请输入大于或等于最小长度的正整数。"
        )
        # 对每个样本进行质量过滤
        for sample_id in metadata['sample_id']:
            input_file = os.path.join(merge_dir, f"{sample_id}.merged.fastq")
            output_file = os.path.join(quality_dir, f"{sample_id}.filtered.fasta")
            if not os.path.exists(input_file):
                logging.warning(f"跳过 {sample_id} - 输入文件 {input_file} 未找到。")
                continue
            cmd = [
                tools["vsearch"], "--fastx_filter", input_file,
                "--fastaout", output_file,
                "--fastq_maxee", "1.0",
                "--fastq_maxee_rate", "0.01",
                "--fastq_minlen", min_length,
                "--fastq_maxlen", max_length,
                "--fastq_maxns", "0",
                "--fasta_width", "0",
                "--threads", str(threads)
            ]
            run_command(cmd, f"质量过滤 {sample_id}")
    else:
        # 如果选择固定长度，定义验证函数确保输入为正整数（可多个，用空格分隔）
        def validate_lengths(s):
            tokens = s.split()
            return all(token.isdigit() and int(token) > 0 for token in tokens) and len(tokens) > 0

        lengths_input = get_valid_input(
            "请输入目标标记的一个或多个长度值（用空格分隔）: ",
            validator=validate_lengths,
            error_msg="请输入一个或多个正整数，用空格分隔。"
        )
        lengths = lengths_input.split()
        # 对每个样本先提取固定长度序列，再进行质量过滤
        for sample_id in metadata['sample_id']:
            input_file = os.path.join(merge_dir, f"{sample_id}.merged.fastq")
            temp_file = os.path.join(quality_dir, f"{sample_id}.temp.fastq")
            output_file = os.path.join(quality_dir, f"{sample_id}.filtered.fasta")
            if not os.path.exists(input_file):
                logging.warning(f"跳过 {sample_id} - 输入文件 {input_file} 未找到。")
                continue
            # 对每个固定长度进行提取，写入临时文件
            with open(temp_file, 'w') as temp_f:
                for length in lengths:
                    cmd = [
                        tools["seqkit"], "seq", "-j", str(threads),
                        "-m", length, "-M", length, input_file
                    ]
                    subprocess.run(cmd, stdout=temp_f, check=True)
            # 使用 vsearch 对提取后的临时文件进行质量过滤
            cmd = [
                tools["vsearch"], "--fastx_filter", temp_file,
                "--fastaout", output_file,
                "--fastq_maxee", "1.0",
                "--fastq_maxee_rate", "0.01",
                "--fastq_maxns", "0",
                "--fasta_width", "0",
                "--threads", str(threads)
            ]
            run_command(cmd, f"质量过滤 {sample_id}")
            os.remove(temp_file)  # 删除临时文件

    # Step 4: 去重复
    dereplicate_dir = os.path.join(output_dir, "4-dereplicate")
    os.makedirs(dereplicate_dir, exist_ok=True)           # 创建去重复输出目录
    all_derep_file = os.path.join(dereplicate_dir, "all_samples_derep.fasta")
    with open(all_derep_file, 'w') as all_f:
        # 对每个样本进行去重复处理
        for sample_id in metadata['sample_id']:
            input_file = os.path.join(quality_dir, f"{sample_id}.filtered.fasta")
            output_file = os.path.join(dereplicate_dir, f"{sample_id}.derep.fasta")
            if not os.path.exists(input_file):
                logging.warning(f"跳过 {sample_id} - 输入文件 {input_file} 未找到。")
                continue
            cmd = [
                tools["vsearch"], "--derep_fulllength", input_file,
                "--strand", "plus",
                "--output", output_file,
                "--sizeout",
                "--relabel", f"{sample_id}.",
                "--fasta_width", "0",
                "--threads", str(threads)
            ]
            run_command(cmd, f"样本内去重复 {sample_id}")
            with open(output_file, 'r') as f:
                all_f.write(f.read())

    # Step 5: 聚类
    cluster_dir = os.path.join(output_dir, "5-cluster")
    os.makedirs(cluster_dir, exist_ok=True)               # 创建聚类输出目录
    cluster_method = get_valid_input(
        "请输入聚类策略: 1. UPARSE3; 2. UNOISE3: ",
        validator=lambda x: x in ["1", "2"],
        error_msg="请输入 1 或 2。"
    )
    otu_file = os.path.join(cluster_dir, "otus.fasta")
    if cluster_method == "1":
        # 使用 UPARSE3 聚类，设定相似性阈值为 0.97
        cmd = [
            tools["vsearch"], "--cluster_smallmem", all_derep_file,
            "--id", "0.97",
            "--usersort",
            "--centroids", otu_file,
            "--strand", "plus",
            "--sizein", "--sizeout",
            "--fasta_width", "0",
            "--threads", str(threads)
        ]
        run_command(cmd, "使用 UPARSE3 进行 OTU 聚类")
    else:
        # 使用 UNOISE3 进行 ASV 聚类
        cmd = [
            tools["vsearch"], "--unoise3", all_derep_file,
            "--centroids", otu_file,
            "--usersort",
            "--sizein", "--sizeout",
            "--fasta_width", "0",
            "--threads", str(threads)
        ]
        run_command(cmd, "使用 UNOISE3 进行 ASV 聚类")

    # Step 6: 嵌合体检测
    chimera_dir = os.path.join(output_dir, "6-chimera")
    os.makedirs(chimera_dir, exist_ok=True)               # 创建嵌合体检测输出目录
    nochim_file = os.path.join(chimera_dir, "otus_nochim.fasta")
    chimera_method = get_valid_input(
        "嵌合体检测方法: 1. de novo; 2. 参考数据库: ",
        validator=lambda x: x in ["1", "2"],
        error_msg="请输入 1 或 2。"
    )
    if chimera_method == "1":
        # 使用 de novo 嵌合体检测
        cmd = [
            tools["vsearch"], "--uchime_denovo", otu_file,
            "--nonchimeras", nochim_file,
            "--sizein",
            "--sizeout",
            "--usersort",
            "--fasta_width", "0",
            "--threads", str(threads)
        ]
        run_command(cmd, "执行 de novo 嵌合体检测")
    else:
        # 使用参考数据库进行嵌合体检测，交互提示参考数据库路径
        ref_db = get_valid_input(
            "请输入参考数据库路径: ",
            validator=lambda x: os.path.exists(x),
            error_msg="文件不存在，请输入正确的路径。"
        )
        cmd = [
            tools["vsearch"], "--uchime_ref", otu_file,
            "--db", ref_db,
            "--usersort",
            "--nonchimeras", nochim_file,
            "--sizein", "--sizeout",
            "--fasta_width", "0",
            "--threads", str(threads)
        ]
        run_command(cmd, "执行参考数据库嵌合体检测")

    # Step 7: OTU 表生成
    otu_dir = os.path.join(output_dir, "7-OTU")
    os.makedirs(otu_dir, exist_ok=True)                   # 创建7-OTU目录
    otu_table = os.path.join(otu_dir, "otu_table.txt")
    cmd = [
        tools["vsearch"], "--usearch_global", all_derep_file,
        "--db", nochim_file,
        "--usersort",
        "--id", "0.97",
        "--otutabout", otu_table,
        "--strand", "plus",
        "--sizein",
        "--threads", str(threads)
    ]
    run_command(cmd, "生成 OTU 表")
    
    # 插入新功能 11、12、13：OTU序列重新标记与进一步聚类、稀释曲线绘制、低丰度OTU过滤
    relabel_and_cluster(tools, threads, output_dir)       # 调用新功能 11
    plot_rarefaction_curve(tools, threads, output_dir)      # 调用新功能 12
    filter_low_abundance_otus(tools, threads, output_dir)   # 调用新功能 13

    # Step 8: 分类（可选）
    classify = get_valid_input(
        "是否执行 SINTAX 分类注释？(yes/no): ",
        validator=lambda x: x.lower() in ["yes", "no"],
        error_msg="请输入 yes 或 no。"
    )
    if classify.lower() == "yes":
        sintax_dir = os.path.join(output_dir, "8-SINTAX")
        os.makedirs(sintax_dir, exist_ok=True)           # 创建SINTAX分类输出目录
        sintax_out = os.path.join(sintax_dir, "otus_sintax.txt")
        # 提示用户输入SINTAX参考数据库路径
        ref_db = get_valid_input(
            "请输入分类参考数据库路径: ",
            validator=lambda x: os.path.exists(x),
            error_msg="文件不存在，请输入正确的路径。"
        )
        cmd = [
            tools["vsearch"], "--sintax", nochim_file,
            "--db", ref_db,
            "--sintax_cutoff", "0.8",
            "--usersort",
            "--strand", "both",
            "--tabbedout", sintax_out,
            "--threads", str(threads)
        ]
        run_command(cmd, "执行 SINTAX 分类注释")

    # 日志记录处理完成，并打印输出目录
    logging.info("扩增子测序数据处理完成！")
    print("处理完成！输出结果位于:", output_dir)

# 程序主入口，捕捉异常并记录日志
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(f"发生错误: {e}")
        print(f"发生错误: {e}")
