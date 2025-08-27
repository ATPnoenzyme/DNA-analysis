# streamlit_app.py
import streamlit as st
import itertools
import random
from typing import List, Dict
import pandas as pd
import io

# 从 core.py 里只引用“加引物、补端、包装”的函数
from core import process_sequence, to_dataframe

BASES = ["A", "T", "C", "G"]

# -------------------------------------------------
# 以下函数只在 streamlit 侧使用，core.py 完全不知情
# -------------------------------------------------
def _replace_k_mutations(seq: str, k: int) -> List[str]:
    """k 个碱基同时替换（全部穷举，不抽样）"""
    seq_list = list(seq)
    n = len(seq)
    all_mut = []
    for pos_tuple in itertools.combinations(range(n), k):
        # 每个位置可选的碱基
        choices = [
            [b for b in BASES if b != seq_list[p]]
            for p in pos_tuple
        ]
        for new_bases in itertools.product(*choices):
            tmp = seq_list[:]
            for p, new in zip(pos_tuple, new_bases):
                tmp[p] = new
            all_mut.append("".join(tmp))
    return all_mut


def _delete_k_mutations(seq: str, k: int) -> List[str]:
    """k 个碱基同时缺失（全部穷举）"""
    n = len(seq)
    dels = []
    for pos_tuple in itertools.combinations(range(n), k):
        new_seq = list(seq)
        for p in sorted(pos_tuple, reverse=True):
            del new_seq[p]
        dels.append("".join(new_seq))
    return dels


# -------------------------------------------------
# Streamlit UI
# -------------------------------------------------
st.set_page_config(page_title="适配体随机突变工具", page_icon="🧬")
st.title("🧬 适配体随机突变工具")

seq      = st.text_area("原始适配体序列（仅碱基字母）", "ATGCGTACGTAGCTAGCTAGCTAGCTAGC").replace("\n", "").upper()
primer_f = st.text_input("正向引物", "")
primer_r = st.text_input("反向引物", "")

# ---------- 突变类型 ----------
mut_type = st.selectbox(
    "突变类型",
    ["replace", "delete"],
    format_func=lambda x: {"replace":"碱基替换", "delete":"碱基缺失"}[x]
)

count = st.number_input("突变碱基个数 k", min_value=1, max_value=len(seq) if seq else 1, value=1)

sample_size = st.number_input("生成序列条数（当组合爆炸时随机抽样）", min_value=1, max_value=100_000, value=100)

# ---------- 缺失补端 ----------
if mut_type == "delete":
    delete_fill_side = st.radio("缺失后补端位置", ["5", "3"], format_func=lambda x: f"{x}' 端")
else:
    delete_fill_side = "5"  # 替换类型用不到

# ---------- 运行 ----------
if st.button("生成突变序列"):
    if not seq:
        st.error("请输入原始序列")
        st.stop()

    # 1. 先在 streamlit 侧生成“纯核心”突变序列
    with st.spinner("正在计算…"):
        if mut_type == "replace":
            core_list = _replace_k_mutations(seq, count)
            # 如果穷举结果太多就抽样
            if len(core_list) > sample_size:
                core_list = random.sample(core_list, sample_size)
        else:  # delete
            core_list = _delete_k_mutations(seq, count)
            if len(core_list) > sample_size:
                core_list = random.sample(core_list, sample_size)

    # 2. 把核心列表交给 core.py 做“加引物、补端、包装”
    #    这里 trick：传入一个“假 mut_type”告诉 core.py 不要二次突变
    res = process_sequence(
        seq=seq,
        primer_f=primer_f,
        primer_r=primer_r,
        mut_type="del1",  # 随便填一个不会二次处理的类型
        keep_length=True,
        delete_fill_side=delete_fill_side,
        sample_size=sample_size
    )

    # 3. 用我们算好的 core_list 覆盖 core.py 生成的结果
    res["core_only"] = core_list
    res["sequences"] = [primer_f + c + primer_r for c in core_list]
    res["report"]    = [f"{mut_type}_k{count}_{i+1}" for i in range(len(core_list))]

    df = pd.DataFrame({
        "sequences": res["sequences"],
        "core":      res["core_only"],
        "report":    res["report"]
    })

    st.success(f"✅ 已生成 {len(df)} 条序列！")
    show_mode = st.radio("显示方式", ["含引物全长", "仅核心序列"])
    col = "sequences" if show_mode == "含引物全长" else "core"
    st.dataframe(df[[col, "report"]].rename(columns={col: "序列"}), use_container_width=True)

    # ---------- 下载 ----------
    buffer = io.BytesIO()
    df.to_excel(buffer, index=False, sheet_name="结果")
    buffer.seek(0)
    st.download_button(
        label="⬇️ 下载 Excel",
        data=buffer,
        file_name=f"{mut_type}_k{count}_mutations.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )
