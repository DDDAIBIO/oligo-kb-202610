"""test_offtarget_scan.py —— D29 工具日②测试（24 断言）
用法：python test_offtarget_scan.py   （与 offtarget_scan_skeleton.py 同目录）
Part A：纯函数（toy 序列，硬编码期望值）  Part B：toy 转录本库 + toy 候选 TSV 端到端
"""
import os
import sys
import tempfile

import offtarget_scan_skeleton as m

PASS = 0
FAIL = 0


def check(name, actual, expect):
    global PASS, FAIL
    if actual == expect:
        PASS += 1
        print(f"  PASS {name}")
    else:
        FAIL += 1
        print(f"  FAIL {name}\n    expected: {expect!r}\n    actual:   {actual!r}")


# ---------------- 固定序列 ----------------
GUIDE = "GCGAUCAUACUGUCGGAUACU"           # 21 nt 合成 guide（RNA）
GUIDE_DNA = "GCGATCATACTGTCGGATACT"
SEED_DNA = "CGATCAT"                      # guide[1:8]（1-based 2–8）
SEED_RC = "ATGATCG"                       # 硬编码：revcomp('CGATCAT')
G_RC = "AGTATCCGACAGTATGATCGC"            # revcomp(GUIDE_DNA)，OFF1 的 15-mer 母串
SENSE_DNA = m.revcomp_dna(GUIDE_DNA) if hasattr(m, "revcomp_dna") else ""

F_T = "A" * 15                            # 中性 filler（与查询最长公共子串 < min_len）
F_A = "A" * 5

OFF1_SEQ = F_T + G_RC[5:20] + F_T         # 15-mer block（i=5），含完整 seed_rc
OFF2_SEQ = "ACGTACGTAC" + SEED_RC + "ACGTACGTAC"   # 只有 seed 命中（7 < 12）
OFF3_SEQ = F_A + GUIDE_DNA + F_A          # passenger 全长 exact
NEG_SEQ = ("TTTTAAAACCCCGGGG" * 2) + ("TTTTAAAACCCCGGGG")   # 干净转录本

TOY_FASTA = "\n".join([
    ">TTR toy on-target",
    F_T + G_RC + F_T,
    ">OFF1 guide 15-mer block",
    OFF1_SEQ,
    ">OFF2 seed only",
    OFF2_SEQ,
    ">OFF3 passenger exact",
    OFF3_SEQ,
    ">NEG clean",
    NEG_SEQ,
])


def revcomp_local(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def main():
    # ============ Part A：纯函数 ============
    print("Part A: 纯函数")

    check("A01 _norm 小写+U->T", m._norm("acgu"), "ACGT")
    check("A02 _norm RNA 长串", m._norm("AUGCUU"), "ATGCTT")
    check("A03 revcomp 基础", m.revcomp_dna("ATGC"), "GCAT")
    check("A04 revcomp 回文", m.revcomp_dna("AATT"), "AATT")
    check("A05 header 基因名", m.gene_from_header(">TTR Homo sapiens chr18"), "TTR")
    check("A06 header 无引号", m.gene_from_header("ENST00000221962.4"), "ENST00000221962.4")

    check("A07 exact 单命中", m.find_exact("ATG", "XXATGXX"), [2])
    check("A08 exact 双命中", m.find_exact("ATG", "ATGATG"), [0, 3])
    check("A09 exact 无命中", m.find_exact("ATG", "GGGGGG"), [])

    check("A10 contig 不重复报 exact",
          m.find_contig(G_RC, F_T + G_RC + F_T, 12), [])
    check("A11 contig 15-mer block",
          m.find_contig(G_RC, OFF1_SEQ, 12), [(5, 15, 15)])
    check("A12 contig min_len 过滤",
          m.find_contig(G_RC, OFF1_SEQ, 16), [])
    check("A13 contig 只报极大块",
          m.find_contig("ACGTAGTTC", "TTTACGTAGCCC", 4), [(0, 3, 6)])

    check("A14 seed rc 命中", m.find_seed(GUIDE, "TTTT" + SEED_RC + "TTT"), [4])
    check("A15 seed 方向陷阱：同向不报", m.find_seed(GUIDE, "TTTT" + SEED_DNA + "TTT"), [])

    # ============ Part B：toy 转录本库端到端 ============
    print("Part B: toy 转录本库")

    tmp = tempfile.mkdtemp(prefix="offtarget_d29_")
    fa = os.path.join(tmp, "toy.fa")
    with open(fa, "w", encoding="utf-8") as fh:
        fh.write(TOY_FASTA + "\n")

    t1 = m.scan_transcript(GUIDE, SENSE_DNA.replace("T", "U"), ">OFF1 x", OFF1_SEQ, 12)
    check("B01 scan_transcript 基因名", t1["gene"], "OFF1")
    check("B02 OFF1 exact=0", len(t1["exact"]), 0)
    check("B03 OFF1 contig 15-mer", t1["contig"], [("guide", 5, 15, 15)])
    check("B04 OFF1 seed 位置", t1["seed"], [23])

    t_ttr = m.scan_transcript(GUIDE, SENSE_DNA.replace("T", "U"),
                              ">TTR x", F_T + G_RC + F_T, 12)
    check("B05 TTR exact（本层不排除）", t_ttr["exact"], [("guide", 15)])

    t3 = m.scan_transcript(GUIDE, SENSE_DNA.replace("T", "U"), ">OFF3 x", OFF3_SEQ, 12)
    check("B06 OFF3 passenger exact", t3["exact"], [("passenger", 5)])

    rows = m.scan_fasta(GUIDE, SENSE_DNA.replace("T", "U"), fa, on_target="TTR", min_len=12)
    check("B07 on-target 排除", [r["gene"] for r in rows], ["OFF3", "OFF1", "OFF2"])
    off1 = next(r for r in rows if r["gene"] == "OFF1")
    check("B08 OFF1 汇总", (off1["n_exact"], off1["max_contig"], off1["n_seed"]), (0, 15, 1))
    off2 = next(r for r in rows if r["gene"] == "OFF2")
    check("B09 OFF2 汇总（仅 seed）", (off2["n_exact"], off2["max_contig"], off2["n_seed"]), (0, 0, 1))
    off3 = next(r for r in rows if r["gene"] == "OFF3")
    check("B10 OFF3 汇总（passenger exact）", (off3["n_exact"], off3["max_contig"], off3["n_seed"]), (1, 0, 0))

    # ============ Part B2：候选 TSV 批量扫描 ============
    print("Part B2: toy 候选 TSV")
    sense1 = revcomp_local(GUIDE_DNA).replace("T", "U")
    g2 = "UUUUAACCCGGGUUUUAACCC"
    sense2 = revcomp_local(g2.replace("U", "T")).replace("T", "U")
    tsv = os.path.join(tmp, "toy_candidates.tsv")
    with open(tsv, "w", encoding="utf-8") as fh:
        fh.write("id\tsense\tguide\n")
        fh.write(f"cand1\t{sense1}\t{GUIDE}\n")
        fh.write(f"cand2\t{sense2}\t{g2}\n")

    out = m.scan_candidates_tsv(tsv, fa, on_target="TTR", min_len=12)
    check("B11 候选数", len(out), 2)
    c1 = out[0]
    check("B12 cand1 脱靶基因数", c1["n_off_genes"], 3)
    check("B13 cand1 exact 合计", c1["n_exact"], 1)
    check("B14 cand1 seed 合计", c1["n_seed"], 2)
    check("B15 cand1 最长连续", c1["max_contig"], 15)
    c2 = out[1]
    check("B16 cand2 干净候选", (c2["n_off_genes"], c2["n_exact"], c2["n_seed"], c2["max_contig"]),
          (0, 0, 0, 0))

    out1 = m.scan_candidates_tsv(tsv, fa, on_target="TTR", min_len=12, limit=1)
    check("B17 limit 参数", len(out1), 1)

    print(f"\nPart A+B: {PASS} PASS / {FAIL} FAIL")
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
