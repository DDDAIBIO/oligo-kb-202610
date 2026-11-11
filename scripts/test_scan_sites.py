# test_scan_sites.py —— D21 代码日
# 12 条纯函数断言 + 2 条 I/O 冒烟测试
# 期望值全部独立手算; 改数字前先问 (与 test_expected.py 同一规矩)
# 用法: 与 scan_sites_skeleton.py 同目录, 实现后直接 python test_scan_sites.py

import os

import scan_sites_skeleton as ss

PASS = 0
FAIL = 0


def check(name, cond):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"PASS {name}")
    else:
        FAIL += 1
        print(f"FAIL {name}")


# ---------- make_windows ----------

# S1: 60 nt, L=21, step=5 -> 8 个窗口, 首 (0,21), 尾 (35,56)
w = ss.make_windows("A" * 60, 21, 5)
check("S1 len==8, first, last",
      len(w) == 8 and w[0] == (0, 21) and w[-1] == (35, 56))

# S2: L > 序列长 -> 空列表, 不是异常
check("S2 empty", ss.make_windows("A" * 10, 21, 5) == [])

# S3: 60 nt, L=21, step=10 -> 起点只能到 30 (30+21=51<=60), 共 4 个
check("S3 step10 len==4", len(ss.make_windows("A" * 60, 21, 10)) == 4)

# ---------- coverage ----------

TXS = ["ACGUACGUACGG", "TTACGUACGUACAA", "GGGGCCCCUUUU", "ACGUACGUAC"]
# S4: 命中 3/4 -> 0.75
check("S4 0.75", abs(ss.coverage("ACGUACGUAC", TXS) - 0.75) < 1e-9)
# S5: 无命中 -> 0.0
check("S5 0.0", ss.coverage("UUUUUUUUUU", TXS) == 0.0)

# ---------- mode_compartment / compartment_ok ----------

# S6: 模态 -> 区室
check("S6 compartment map",
      ss.mode_compartment("siRNA") == "cytosol_mRNA"
      and ss.mode_compartment("gapmer") == "nucleus_preMRNA"
      and ss.mode_compartment("splice_aso") == "nucleus_preMRNA")

# S7-S8: 区室合规 (D20 §3 表, 5 例)
check("S7 siRNA intron False, CDS True",
      ss.compartment_ok("intron", "siRNA") is False
      and ss.compartment_ok("CDS", "siRNA") is True)
check("S8 splice_aso intron True, CDS False; gapmer intron True",
      ss.compartment_ok("intron", "splice_aso") is True
      and ss.compartment_ok("CDS", "splice_aso") is False
      and ss.compartment_ok("intron", "gapmer") is True)

# ---------- is_seed (窗口内 1-based, 2..8 含两端) ----------

check("S9 seed boundaries",
      ss.is_seed(1, 21) is False and ss.is_seed(2, 21) is True
      and ss.is_seed(8, 21) is True and ss.is_seed(9, 21) is False)

# ---------- window_access ----------

ACC = [0.2] * 5 + [0.8] * 5
# S10: 前 5 均值 0.2, 后 5 均值 0.8
check("S10 access slices",
      abs(ss.window_access(ACC, 0, 5) - 0.2) < 1e-9
      and abs(ss.window_access(ACC, 5, 10) - 0.8) < 1e-9)

# ---------- pop_coverage ----------

# S11: 1 - (0.02+0.01) = 0.97; 空列表 = 1.0; 和超 1 时下限 0
check("S11 pop coverage",
      abs(ss.pop_coverage([(0.02, 10, 11), (0.01, 20, 21)]) - 0.97) < 1e-9
      and ss.pop_coverage([]) == 1.0
      and ss.pop_coverage([(0.6, 0, 1), (0.6, 5, 6)]) == 0.0)

# ---------- 变异重叠 (辅助函数, 也算纯函数) ----------

VAR = [(0.05, 200, 201), (0.01, 210, 211)]
# S12: 窗口 [190, 215) 命中 2 条; 窗口 [0, 21) 命中 0 条
check("S12 overlap",
      len(ss.variants_in_window(VAR, 190, 215)) == 2
      and ss.variants_in_window(VAR, 0, 21) == [])

# ---------- I/O 冒烟测试 (文件不存在则跳过) ----------

if os.path.exists("ttr.fa"):
    txs = ss.load_transcripts(".")
    acc = ss.accessibility("ttr.fa")
    seq = list(txs.values())[0]
    check("I1 acc len == seq len", len(acc) == len(seq))
    check("I2 acc values in [0,1]",
          all(0.0 <= p <= 1.0 for p in acc))
else:
    print("SKIP I1/I2 (ttr.fa 不存在; 下载后重跑)")

if os.path.exists("variants.tsv"):
    vs = ss.load_variants("variants.tsv")
    check("I3 variants parsed", len(vs) > 0
          and all(0.0 <= m <= 1.0 and s2 < e2 for m, s2, e2 in vs))
else:
    print("SKIP I3 (variants.tsv 不存在; 造数据后重跑)")

print(f"\n--- {PASS} PASS / {FAIL} FAIL ---")
if FAIL == 0:
    print("纯函数区全绿 -> 允许进入 B3 (RNAplfold 接入)")
