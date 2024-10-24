from alphafold3_pytorch.data import mmcif_parsing
from alphafold3_pytorch.data.data_pipeline import make_mmcif_features, to_mmcif
from Bio.Data import PDBData

def extract_chain_sequences_and_lengths(cif_file_path: str):
    """提取CIF文件中每条链的序列及其长度"""
    
    # 使用 mmcif_parsing 解析 CIF 文件
    mmcif_object = mmcif_parsing.parse_mmcif_object(cif_file_path, cif_file_path.split('/')[-1])

    if mmcif_object is None:
        print(f"无法解析 CIF 文件: {cif_file_path}")
        return

    chain_sequences = {}

    # 遍历结构中的每条链
    for chain in mmcif_object.structure.get_chains():
        chain_id = chain.id
        sequence = []
        
        # 遍历链中的每个残基，提取氨基酸单字母代码
        for residue in chain.get_residues():
            resname = residue.get_resname()
            
            # 通过残基名称获取氨基酸单字母代码
            amino_acid = PDBData.protein_letters_3to1.get(resname, 'X')
            sequence.append(amino_acid)

        # 将序列转化为字符串
        chain_sequence = ''.join(sequence)
        chain_length = len(chain_sequence)
        
        # 保存序列和长度
        chain_sequences[chain_id] = (chain_sequence, chain_length)

    # 输出每条链的序列和长度
    for chain_id, (sequence, length) in chain_sequences.items():
        print(f"Chain {chain_id}: Sequence Length = {length}")
        print(f"Sequence: {sequence}\n")


# 示例用法
# cif_file_path = '/mnt/syc/alphafold3-pytorch/test-folder/data/train/7a4d-assembly1_reconstructed.cif'
# cif_file_path =  "/mnt/syc/alphafold3-pytorch/data/pdb_data/train_mmcifs/a4/7a4d-assembly1.cif"
# extract_chain_sequences_and_lengths(cif_file_path)

# cif_file_path = '/mnt/syc/alphafold3-pytorch/test-folder/data/test/0.cif'
# extract_chain_sequences_and_lengths(cif_file_path)

# cif_file_path = '/mnt/syc/alphafold3-pytorch/test-folder/data/valid/0.cif'
# extract_chain_sequences_and_lengths(cif_file_path)


def check_mmcif_file_content(filepath: str):
    """
    检查 mmCIF 文件中是否包含蛋白质、DNA、RNA 或配体，并打印相应信息。

    :param filepath: mmCIF 文件的路径
    """
    # 解析 mmCIF 文件
    try:
        mmcif_object = mmcif_parsing.parse_mmcif_object(filepath, file_id=filepath)
    except Exception as e:
        print(f"Failed to parse the mmCIF file: {e}")
        return

    contains_protein = False
    contains_dna = False
    contains_rna = False
    contains_ligand = False
    
    # 遍历每个链的化学成分
    for chain_id, chem_comps in mmcif_object.chem_comp_details.items():
        for chem_comp in chem_comps:
            chem_type = chem_comp.type.lower()
            if "peptide" in chem_type:
                contains_protein = True
            elif "dna" in chem_type:
                contains_dna = True
            elif "rna" in chem_type:
                contains_rna = True
            else:
                contains_ligand = True
    
    # 打印检测结果
    if contains_protein:
        print("This file contains protein.")
    if contains_dna:
        print("This file contains DNA.")
    if contains_rna:
        print("This file contains RNA.")
    if contains_ligand:
        print("This file contains ligand.")

    # 如果没有检测到任何成分
    if not (contains_protein or contains_dna or contains_rna or contains_ligand):
        print("No recognized biological components found in the mmCIF file.")

# 使用示例
# check_mmcif_file_content('example.cif')

def extract_chain_sequences_and_lengths(cif_file_path: str):
    """提取CIF文件中每条链的序列、长度及其类型（蛋白质、DNA、RNA或其他）"""
    
    # 使用 mmcif_parsing 解析 CIF 文件
    try:
        mmcif_object = mmcif_parsing.parse_mmcif_object(cif_file_path, cif_file_path.split('/')[-1])
    except Exception as e:
        print(f"无法解析 CIF 文件: {e}")
        return

    if mmcif_object is None:
        print(f"无法解析 CIF 文件: {cif_file_path}")
        return

    chain_sequences = {}

    # 遍历结构中的每条链
    for chain in mmcif_object.structure.get_chains():
        chain_id = chain.id
        sequence = []
        chain_type = "Unknown"  # 初始化链类型为未知
        
        # 获取该链的化学成分类型
        if chain_id in mmcif_object.chem_comp_details:
            chem_comps = mmcif_object.chem_comp_details[chain_id]
            for chem_comp in chem_comps:
                chem_type = chem_comp.type.lower()
                if "peptide" in chem_type:
                    chain_type = "Protein"
                elif "dna" in chem_type:
                    chain_type = "DNA"
                elif "rna" in chem_type:
                    chain_type = "RNA"
                else:
                    chain_type = "Ligand"

        # 遍历链中的每个残基，提取氨基酸单字母代码或核苷酸
        for residue in chain.get_residues():
            resname = residue.get_resname()
            
            if chain_type == "Protein":
                # 通过残基名称获取氨基酸单字母代码
                amino_acid = PDBData.protein_letters_3to1.get(resname, 'X')
                sequence.append(amino_acid)
            elif chain_type in ["DNA", "RNA"]:
                # 获取核苷酸代码
                nucleotide = PDBData.nucleic_letters_3to1.get(resname, 'N')
                sequence.append(nucleotide)
            else:
                # 非蛋白质或核酸的链，暂时跳过
                sequence.append('X')
        
        # 将序列转化为字符串
        chain_sequence = ''.join(sequence)
        chain_length = len(chain_sequence)
        
        # 保存序列、长度和链的类型
        chain_sequences[chain_id] = (chain_sequence, chain_length, chain_type)

    # 输出每条链的序列、长度和类型
    for chain_id, (sequence, length, chain_type) in chain_sequences.items():
        print(f"Chain {chain_id} ({chain_type}): Sequence Length = {length}")
        print(f"Sequence: {sequence}\n")


# cif_file_path = '/mnt/syc/alphafold3-pytorch/test-folder/data/train/6yp7-assembly1.cif'
# cif_file_path = '/mnt/syc/alphafold3-pytorch/test-folder/data/test/7a4d-assembly1_reconstructed.cif'
# extract_chain_sequences_and_lengths(cif_file_path) 

import os 

filepath = "/data/syc/af3_fork/alphafold3-pytorch/data/pdb_data/train_mmcifs/a4/7a4d-assembly1.cif"
file_id = os.path.splitext(os.path.basename(filepath))[0]

mmcif_object = mmcif_parsing.parse_mmcif_object(
    filepath=filepath,
    file_id=file_id,
)
mmcif_feats, assembly = make_mmcif_features(mmcif_object)
cropped_assembly, _, _ = assembly.crop(
    contiguous_weight=0.2,
    spatial_weight=0.4,
    spatial_interface_weight=0.4,
    n_res=150, # 384/2=192
    chain_1="A",
    chain_2="B",
)
mmcif_string = to_mmcif(
    # assembly,
    cropped_assembly,
    file_id=file_id,
    gapless_poly_seq=True,
    insert_alphafold_mmcif_metadata=False,
    # unique_res_atom_names=assembly.unique_res_atom_names,
    unique_res_atom_names=cropped_assembly.unique_res_atom_names,
)
with open("/data/syc/af3_fork/alphafold3-pytorch/test-folder/data/train/7a4d-assembly1_reconstructed.cif", "w") as f:
    f.write(mmcif_string)

print(f"Successfully reconstructed {filepath} after mmCIF featurization.")

print("---------------------before crop-----------------------: ")
cif_file_path =  "/data/syc/af3_fork/alphafold3-pytorch/data/pdb_data/train_mmcifs/a4/7a4d-assembly1.cif"
extract_chain_sequences_and_lengths(cif_file_path)

print("---------------------after crop-----------------------: ")
cif_file_path = '/data/syc/af3_fork/alphafold3-pytorch/test-folder/data/train/7a4d-assembly1_reconstructed.cif'
extract_chain_sequences_and_lengths(cif_file_path)
