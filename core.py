# core.py
import random
import itertools
import pandas as pd
from typing import List, Dict, Literal

BASES = ["A", "T", "C", "G"]

# ------------------ 随机突变 ------------------
def single_mutations(seq: str) -> List[str]:
    return [seq[:i] + b + seq[i+1:]
            for i, orig in enumerate(seq)
            for b in BASES if b != orig]

def double_mutations(seq: str) -> List[str]:
    muts = []
    for i, j in itertools.combinations(range(len(seq)), 2):
        orig_i, orig_j = seq[i], seq[j]
        for b_i in BASES:
            if b_i == orig_i:
                continue
            for b_j in BASES:
                if b_j == orig_j:
                    continue
                muts.append(seq[:i] + b_i + seq[i+1:j] + b_j + seq[j+1:])
    return muts

def triple_mutations_sampled(seq: str, sample_size: int = 44_424) -> List[str]:
    from collections import OrderedDict
    upper = len(list(itertools.combinations(range(len(seq)), 3))) * (3 ** 3)
    sample_size = min(sample_size, upper)
    seen = OrderedDict()
    while len(seen) < sample_size:
        pos = sorted(random.sample(range(len(seq)), 3))
        new_seq = list(seq)
        for p in pos:
            new_seq[p] = random.choice([b for b in BASES if b != seq[p]])
        seen[''.join(new_seq)] = None
    return list(seen.keys())

# ------------------ 随机缺失 ------------------
def del1(seq: str) -> List[str]:
    return [seq[:i] + seq[i+1:] for i in range(len(seq))]

def del2(seq: str) -> List[str]:
    return [seq[:i] + seq[i+1:j] + seq[j+1:]
            for i, j in itertools.combinations(range(len(seq)), 2)]

# ------------------ 主流程 ------------------
def process_sequence(
    seq: str,
    primer_f: str,
    primer_r: str,
    mut_type: Literal["single", "double", "triple", "del1", "del2"],
    *,
    keep_length: bool = True,
    delete_fill_side: Literal["5", "3"] = "5",
    sample_size: int = 44_424
) -> Dict[str, List[str]]:
    """
    返回：
        {
            "sequences":  [全长突变序列列表（已加引物）],
            "core_only":  [仅核心（去掉引物）的突变序列列表],
            "report":     [每条序列的简要描述]
        }
    """
    seq = seq.upper()
    primer_f = primer_f.upper()
    primer_r = primer_r.upper()

    # 1. 根据突变类型生成核心列表
    if mut_type == "single":
        core_list = single_mutations(seq)
    elif mut_type == "double":
        core_list = double_mutations(seq)
    elif mut_type == "triple":
        core_list = triple_mutations_sampled(seq, sample_size)
    elif mut_type == "del1":
        core_list = del1(seq)
    elif mut_type == "del2":
        core_list = del2(seq)
    else:
        raise ValueError("Unsupported mutation type")

    # 2. 缺失时保持长度
    if mut_type in ("del1", "del2") and keep_length:
        n_del = abs(len(core_list[0]) - len(seq))
        fill_seq = "".join(random.choices(BASES, k=n_del))
        if delete_fill_side == "5":
            core_list = [fill_seq + c for c in core_list]
        else:  # 3'
            core_list = [c + fill_seq for c in core_list]

    # 3. 加引物
    full_list = [primer_f + c + primer_r for c in core_list]

    # 4. 生成报告
    report = [f"{mut_type}_{i+1}" for i in range(len(full_list))]

    return {
        "sequences": full_list,
        "core_only": core_list,
        "report": report
    }

# ------------------ 额外便捷 ------------------
def to_dataframe(result: Dict[str, List[str]]) -> pd.DataFrame:
    return pd.DataFrame({
        "sequences": result["sequences"],
        "core": result["core_only"],
        "report": result["report"]
    })
