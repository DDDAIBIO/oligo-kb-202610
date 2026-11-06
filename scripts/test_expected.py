# -*- coding: utf-8 -*-
"""
test_expected.py — D16 验收测试（期望值已用独立程序算死，勿改数字）

用法：
  1. 先把 oligo_decorate_skeleton.py 复制为 oligo_decorate.py
  2. 在 oligo_decorate.py 里实现 9 个 TODO 函数
  3. 运行：python test_expected.py
     全部 ok = 通过验收，把完整输出贴进笔记，然后找 AI review

期望值来源：2026-11-06 用独立脚本（不含你的实现）逐条计算并人工核对；
T0 与 D15 手算练习完全一致（gap 内 CpG=2，全串 CG=3）。
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), "oligo_decorate.py")):
    sys.exit("未找到 oligo_decorate.py —— 先把 oligo_decorate_skeleton.py 复制为该文件再实现")

from oligo_decorate import (  # noqa: E402
    parse_notation,
    summarize_residues,
    count_cpg,
    find_gu_motifs,
    find_au_rich,
    max_self_complement,
    gc_content,
    tm_wallace,
    find_hairpins,
    ps_profile,
    g_run_count,
)

T0 = "mC*mGmUmAmCdTdCdGdAdTdCdGdTdAdAmUmGmUmA*mC"          # D15 手算的合成 gapmer
T1 = "eU*eC*eA*eC*eU*eU*eU*eC*eA*eU*eA*eA*eU*eG*eC*eU*eG*eG"  # Nusinersen [文献序列，待核]
T2 = "GmCAmUGmCAmUGmCAmUGmCAmUGmCAmUG"                      # Patisiran 风格 [合成，非真实序列]
T3 = "GGGUGUCCUUCAAGGGAAAA"                                  # 含 GUCCUUCAA [合成]
T4 = "AAUUUAGCAUUUAGCAUUUAG"                                 # 含 AUUUA x3 [合成]
T5 = "UG*mGICAUmGICC"                                        # 含 I（肌苷）[合成]


def prep(seq):
    residues, ps = parse_notation(seq)
    info = summarize_residues(residues)
    # ps_positions 换算：本测试假设你已在 main/summarize 里处理；
    # 这里独立重算一次位置，避免测到解析器而非你的特征函数
    pos, seen = [], 0
    for ch in seq:
        if ch == "*":
            pos.append(seen - 1)
        elif ch.isupper():
            seen += 1
    info["ps_positions"] = pos
    return info


checks = []


def check(name, got, want):
    ok = got == want
    checks.append(ok)
    print(("PASS" if ok else "FAIL"), name, "| got:", got, "| want:", want)


# ---- T0：D15 手算对照（最重要的一个）----
i = prep(T0)
check("T0 len", len(i["bases"]), 20)
check("T0 PS total", ps_profile(i["ps_positions"])["total"], 2)
check("T0 DNA-CpG", count_cpg(i), 2)
check("T0 GC%", round(gc_content(i["bases"]) * 100, 1), 45.0)
check("T0 mod counts: 2'-OMe 个数", sum(1 for m in i["mod"] if m == "m"), 10)
check("T0 motif hits 为空", find_gu_motifs(i["bases"]), {})

# ---- T1：Nusinersen（全 MOE 全 PS，无 DNA 区）----
i = prep(T1)
check("T1 len", len(i["bases"]), 18)
check("T1 PS total", ps_profile(i["ps_positions"])["total"], 17)
check("T1 2'-MOE 个数", sum(1 for m in i["mod"] if m == "e"), 18)
check("T1 DNA-CpG=0（全是 RNA！）", count_cpg(i), 0)
check("T1 U-run>=3 命中 1 处（0-based 4）", find_gu_motifs(i["bases"]).get("U3"), [4])

# ---- T2：Patisiran 风格（只修饰嘧啶，PS=0 —— 回连 C5）----
i = prep(T2)
check("T2 len", len(i["bases"]), 21)
check("T2 PS=0（Patisiran 无硫代！）", ps_profile(i["ps_positions"])["total"], 0)
check("T2 2'-OMe 个数", sum(1 for m in i["mod"] if m == "m"), 10)
check("T2 GC%", round(gc_content(i["bases"]) * 100, 1), 52.4)

# ---- T3：GU-rich 命中 ----
i = prep(T3)
check("T3 GUCCUUCAA @0-based 4", find_gu_motifs(i["bases"]).get("GUCCUUCAA"), [4])

# ---- T4：ARE 命中 3 处 ----
i = prep(T4)
check("T4 AUUUA x3", find_au_rich(i["bases"]), [1, 8, 15])

# ---- T5：含 I 不能崩，I 按 G 计 GC ----
i = prep(T5)
check("T5 len（含 I 照常计数）", len(i["bases"]), 11)
check("T5 PS=1", ps_profile(i["ps_positions"])["total"], 1)
check("T5 2'-OMe=2", sum(1 for m in i["mod"] if m == "m"), 2)
check("T5 GC% (I 按 G)", round(gc_content(i["bases"]) * 100, 1), 72.7)

# ---- 加分项（不算硬验收，但强烈建议）----
# a) 非法输入报错：parse_notation("*mC") / parse_notation("mC*") / parse_notation("") 应 raise ValueError
# b) Tm：修饰序列的 Tm 只要求 > 同序列 Wallace 基础值不可信时的人工判断，v0 不 assert 具体数

print()
if all(checks):
    print("全部 %d 项通过。下一步：跑 CLI 看完整报告 -> 贴进笔记 -> commit -> 找 AI review。" % len(checks))
else:
    print("有 %d 项未过。先看第一条 FAIL。" % (len(checks) - sum(checks)))
    sys.exit(1)
