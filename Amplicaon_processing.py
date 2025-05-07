#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
import pandas as pd
import logging
import argparse

# 设置日志记录
logging.basicConfig(
    filename="amplicon_processing.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

def get_valid_input(prompt, validator=None, error_msg="输入不符合要求，请重新输入。"):
    """
    循环提示用户输入，直到满足 validator 条件或用户输入“exist”退出程序。
    
    参数:
        prompt (str): 提示信息
        validator (function): 一个函数，输入字符串，返回 True/False 表示是否有效
        error_msg (str): 当输入无效时显示的错误提示
    返回:
        str: 有效的输入
    """
    while True:
        value = input(prompt).strip()
        if value.lower() == "exist":
            sys.exit("用户选择结束运行。")
        if validator is None or validator(value):
            return value
        else:
            print(error_msg)

def check_tool(tool_name):
    """
    检查指定工具是否在系统路径中可用。
    如果不可用，提示用户输入安装路径。
    
    参数:
        tool_name (str): 工具名称
    返回:
        str: 工具的完整路径
    """
    tool_path = shutil.which(tool_name)
    if tool_path:
        logging.info(f"{tool_name} found at {tool_path}")
        return tool_path
    else:
        while True:
            dir_input = input(f"{tool_name} not found. Please input its installation directory (absolute path): ")
            if dir_input.strip().lower() == "exist":
                sys.exit("用户选择结束运行。")
            dir_path = os.path.abspath(dir_input.strip("'"))
            tool_full_path = os.path.join(dir_path, tool_name)
            if os.path.isfile(tool_full_path) and os.access(tool_full_path, os.X_OK):
                logging.info(f"{tool_name} found at {tool_full_path}")
                return tool_full_path
            else:
                print("Tool not found or not executable. Please try again.")

def run_command(cmd, description=""):
    """
    执行外部命令并检查是否成功。
    
    参数:
        cmd (list): 命令及其参数列表
        description (str): 命令的描述，用于日志记录
    """
    logging.info(f"Running: {description}")
    # 将 cmd 列表中的所有元素转换为字符串再拼接
    logging.info(f"Command: {' '.join(map(str, cmd))}")
    try:
        result = subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        logging.info(f"{description} completed successfully.")
        logging.debug(f"Output: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in {description}: {e.stderr}")
        raise

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

def main():
    # 解析命令行参数
    args = parse_args()
    input_dir = os.path.abspath(args.input_dir)
    metadata_file = os.path.abspath(args.metadata_file)
    output_dir = os.path.abspath(args.output_dir)
    threads = args.threads

    # 检查输入文件和目录是否存在
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"输入目录 {input_dir} 不存在。")
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"元数据文件 {metadata_file} 不存在。")
    os.makedirs(output_dir, exist_ok=True)

    # 检查所需工具
    tools = {
        "cutadapt": check_tool("cutadapt"),
        "vsearch": check_tool("vsearch"),
        "usearch": check_tool("usearch"),
        "seqkit": check_tool("seqkit"),
        "csvtk": check_tool("csvtk")
    }

    # 读取元数据，假设文件无表头，列顺序为：
    # run_id, sample_id, forward_primer, reverse_primer, forward_file, reverse_file
    logging.info("正在加载元数据...")
    metadata = pd.read_csv(metadata_file, sep='\t', header=None)
    metadata.columns = ['run_id', 'sample_id', 'forward_primer', 'reverse_primer', 'forward_file', 'reverse_file']

    # Step 1: 样本解复用（双端模式）
    # 使用 cutadapt 同时处理正向和反向原始数据
    demultiplex_dir = os.path.join(output_dir, "1-demultiplex")
    os.makedirs(demultiplex_dir, exist_ok=True)
    for _, row in metadata.iterrows():
        sample_id = row['sample_id']
        barcode_for = row['forward_primer']
        barcode_rev = row['reverse_primer']
        forward_file = os.path.join(input_dir, row['forward_file'])
        reverse_file = os.path.join(input_dir, row['reverse_file'])
        output_forward = os.path.join(demultiplex_dir, f"{sample_id}.R1.fastq")
        output_reverse = os.path.join(demultiplex_dir, f"{sample_id}.R2.fastq")
        
        if not os.path.exists(forward_file) or not os.path.exists(reverse_file):
            logging.warning(f"跳过 {sample_id} - 正向或反向文件未找到。")
            continue
        
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
        run_command(cmd, f"样本解复用 {sample_id}")

    # Step 2: 合并双端测序文件（对解复用后的数据进行合并）
    merge_dir = os.path.join(output_dir, "2-merge")
    os.makedirs(merge_dir, exist_ok=True)
    for _, row in metadata.iterrows():
        sample_id = row['sample_id']
        demux_forward = os.path.join(demultiplex_dir, f"{sample_id}.R1.fastq")
        demux_reverse = os.path.join(demultiplex_dir, f"{sample_id}.R2.fastq")
        merge_out = os.path.join(merge_dir, f"{sample_id}.merged.fastq")
        
        if not os.path.exists(demux_forward) or not os.path.exists(demux_reverse):
            logging.warning(f"跳过 {sample_id} - 解复用文件未找到。")
            continue
        
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
    os.makedirs(quality_dir, exist_ok=True)
    marker_type = get_valid_input(
        "PCR 产物长度类型: 1. 在一定范围内(default 200~400bp); 2. 固定为一个/多个值: ",
        validator=lambda x: x in ["1", "2"],
        error_msg="请输入 1 或 2。"
    )
    if marker_type == "1":
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
        def validate_lengths(s):
            tokens = s.split()
            return all(token.isdigit() and int(token) > 0 for token in tokens) and len(tokens) > 0

        lengths_input = get_valid_input(
            "请输入目标标记的一个或多个长度值（用空格分隔）: ",
            validator=validate_lengths,
            error_msg="请输入一个或多个正整数，用空格分隔。"
        )
        lengths = lengths_input.split()
        for sample_id in metadata['sample_id']:
            input_file = os.path.join(merge_dir, f"{sample_id}.merged.fastq")
            temp_file = os.path.join(quality_dir, f"{sample_id}.temp.fastq")
            output_file = os.path.join(quality_dir, f"{sample_id}.filtered.fasta")
            if not os.path.exists(input_file):
                logging.warning(f"跳过 {sample_id} - 输入文件 {input_file} 未找到。")
                continue
            with open(temp_file, 'w') as temp_f:
                for length in lengths:
                    cmd = [
                        tools["seqkit"], "seq", "-j", str(threads),
                        "-m", length, "-M", length, input_file
                    ]
                    subprocess.run(cmd, stdout=temp_f, check=True)
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
            os.remove(temp_file)

    # Step 4: 去重复
    dereplicate_dir = os.path.join(output_dir, "4-dereplicate")
    os.makedirs(dereplicate_dir, exist_ok=True)
    all_derep_file = os.path.join(dereplicate_dir, "all_samples_derep.fasta")
    with open(all_derep_file, 'w') as all_f:
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
    os.makedirs(cluster_dir, exist_ok=True)
    cluster_method = get_valid_input(
        "请输入聚类策略: 1. UPARSE3; 2. UNOISE3: ",
        validator=lambda x: x in ["1", "2"],
        error_msg="请输入 1 或 2。"
    )
    otu_file = os.path.join(cluster_dir, "otus.fasta")
    if cluster_method == "1":
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
    os.makedirs(chimera_dir, exist_ok=True)
    nochim_file = os.path.join(chimera_dir, "otus_nochim.fasta")
    chimera_method = get_valid_input(
        "嵌合体检测方法: 1. de novo; 2. 参考数据库: ",
        validator=lambda x: x in ["1", "2"],
        error_msg="请输入 1 或 2。"
    )
    if chimera_method == "1":
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
    os.makedirs(otu_dir, exist_ok=True)
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

    # Step 8: 分类（可选）
    classify = get_valid_input(
        "是否执行 SINTAX 分类注释？(yes/no): ",
        validator=lambda x: x.lower() in ["yes", "no"],
        error_msg="请输入 yes 或 no。"
    )
    if classify.lower() == "yes":
        sintax_dir = os.path.join(output_dir, "8-SINTAX")
        os.makedirs(sintax_dir, exist_ok=True)
        sintax_out = os.path.join(sintax_dir, "otus_sintax.txt")
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

    logging.info("扩增子测序数据处理完成！")
    print("处理完成！输出结果位于:", output_dir)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(f"发生错误: {e}")
        print(f"发生错误: {e}")
