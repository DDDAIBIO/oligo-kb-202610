"""D28 siRNA_design 骨架 -- 候选生成 + 三规则打分 (v0)

依赖: scan_sites_skeleton 的 _read_fasta / make_windows (D21 产出, 禁止重写)

两处声明过的简化(实现时必须保留 🟡 标记, 见 D28 doc §2):
  S3 Reynolds 判据3: 用"中心 7-mer(pos7-13) A/U>=4"代理, 待回原文核内部稳定性算法
  S5 Amarzguioui 条件1: 用"guide 两端各 3 nt A/U 差分>0"代理, 待回原文核

硬约定(D26 拍板):
  三套规则三列并存(uitei/reynolds/amarz), 禁止合成总分
  guide 只存 19 nt 核心, dTdT 突出在合成时加, 不进序列
"""

import sys

try:
    import scan_sites_skeleton as ss
except ImportError:  # 允许与骨架同目录直接运行
    ss = sys.modules[__name__]


# ---------- S0-S1: 基础工具 ----------

def _revcomp_rna(seq):
    """S0: RNA 反向互补 (A<->U, G<->C); 空串返回空串; 小写转大写"""
    #raise NotImplementedError("S0")
     s = seq.upper()
    return "".join({"A": "U", "U": "A", "G": "C", "C": "G"}[c] for c in reversed(s))


def gc_frac(seq):
    """S1: GC 占比 0-1; 空串返回 0.0"""
    #raise NotImplementedError("S1")
    s = seq.upper()
    return (sum(c in "GC" for c in s) / len(s)) if s else 0.0


def _au(seq):
    """内部工具: A/U 计数 (可直接实现, 也可留给实现者)"""
    #raise NotImplementedError("S1")
    s = seq.upper()
    return sum(1 for c in s if c in "AU")


# ---------- S2: Ui-Tei (过滤器, AND 逻辑, guide 坐标系) ----------

def uitei_pass(guide):
    """S2: 四条件 AND -> bool
    (i)   guide[0] in AU
    (ii)  guide[18] in GC
    (iii) guide[0:7] A/U >= 4   (siDirect 口径 >=4/7; 原文摘要为 >=5, 见 D26 §4 漂移)
    (iv)  无 >=10 连续 GC
    """
    #raise NotImplementedError("S2")
    return (guide[0] in "AU" and
            guide[18] in "GC" and
            _au(guide[0:7]) >= 4 and
            "GC" * 10 not in guide)
    """S3 (代理, 🟡): 中心 7-mer sense[6:13] A/U >= 4 -> 1, 否则 0"""
    raise NotImplementedError("S3")


def _c3_low_internal_stability(sense):
    """S3 (代理, 🟡): 中心 7-mer sense[6:13] A/U >= 4 -> 1, 否则 0"""
    return 1 if _au(sense[6:13]) >= 4 else 0


def reynolds_score(sense):
    """S4: 八条累加 0-10, 不设下限(可负)
    1 GC 30-52% -> +1
    2 sense[14:19] 每个 A/U +1 (最多 +5)
    3 S3 代理 -> +1
    4 sense[18]=='A' -> +1
    5 sense[2]=='A'  -> +1
    6 sense[9]=='U'  -> +1
    7 sense[18] in GC -> -1
    8 sense[12]=='G' -> -1
    """
    #raise NotImplementedError("S4")
    s = sense.upper()
    sc = 0
    if 0.30 <= gc_frac(s) <= 0.52:          # 判据1 (S1 的 0-1 口径, 别乘 100)
        sc += 1
    sc += _au(s[14:19])                     # 判据2: 切片只有 5 位, 天然上限 +5
    sc += _c3_low_internal_stability(s)     # 判据3: 代理
    if s[18] == "A":                        # 判据4
        sc += 1
    if s[2] == "A":                         # 判据5
        sc += 1
    if s[9] == "U":                         # 判据6
        sc += 1
    if s[18] in "GC":                       # 判据7: 纯罚分!
        sc -= 1
    if s[12] == "G":                        # 判据8: 纯罚分!
        sc -= 1
    return sc


# ---------- S5-S6: Amarzguioui (分列展示, guide 坐标系) ----------

def _amarz_asymmetry(guide):
    """S5 (代理, 🟡): guide[0:3] 与 guide[16:19] 的 A/U 数差 >0 -> 1, 否则 0"""
    #raise NotImplementedError("S5")
    return 1 if (_au(guide[0:3]) - _au(guide[16:19])) > 0 else 0


def amarzguioui_score(guide):
    """S6: 六条 +/- 体系
    +1 asymmetry(S5); +1 guide[18] in GC; +1 guide[0] in AU; +1 guide[13] in AU;
    -1 guide[18]=='U'; -1 guide[0]=='G'
    """
    #raise NotImplementedError("S6")
    g = guide.upper()
    sc = _amarz_asymmetry(g)            # +1: 不对称
    if g[18] in "GC":                   # +1: pos19 紧
        sc += 1
    if g[0] in "AU":                    # +1: pos1 松
        sc += 1
    if g[13] in "AU":                   # +1: pos14 松
        sc += 1
    if g[18] == "U":                    # -1: pos19 是 U (纯罚分!)
        sc -= 1
    if g[0] == "G":                     # -1: pos1 是 G (纯罚分!)
        sc -= 1
    return sc


# ---------- S7-S9: 候选生成与排序 ----------

def design_candidates(seq, step=1):
    """S7: 全枚举 19-mer; 每窗 -> dict(start, end, sense, guide, uitei, reynolds, amarz)
    sense = 转录本窗口本身; guide = _revcomp_rna(sense)
    """
    #raise NotImplementedError("S7")
    seq = seq.upper()
    rows = []
    for s, e in ss.make_windows(seq, 19, step):
        sense = seq[s:e]
        guide = _revcomp_rna(sense)
        rows.append(dict(start=s, end=e, sense=sense, guide=guide,
                         uitei=uitei_pass(guide),
                         reynolds=reynolds_score(sense),
                         amarz=amarzguioui_score(guide)))
    return rows


def rank_design(fasta_path, top=20):
    """S8: _read_fasta -> 逐转录本 S7 -> 排序键 (-reynolds, -amarz, start)
    返回 [(name, rows)]; 排序注释必须写明"不合成总分(D26 §5)"
    """
    #raise NotImplementedError("S8")
    out = []
    for name, seq in ss._read_fasta(fasta_path).items():
        rows = sorted(design_candidates(seq),
                      key=lambda r: (-r["reynolds"], -r["amarz"], r["start"]))
        # 不合成总分(D26 §5): reynolds 与 amarz 量纲不同(0-10 vs -1~5), pos19 结论相反,
        # 合成即替读者仲裁实验证据; uitei 是 AND 过滤器属性, 不参与排序
        out.append((name, rows))
    return out


def format_table(name, rows, top=20):
    """S9: 表头 rank start end sense guide uitei reynolds amarz
    屏幕打 top 行; 返回 (text, all_rows) -- 全表交 main 写盘 (D21 教训)
    """
    #raise NotImplementedError("S9")
    cols = ["rank", "start", "end", "sense", "guide", "uitei", "reynolds", "amarz"]
    lines = [f"转录本: {name}   候选窗口: {len(rows)}   (三列并存, 不合成总分 - D26 §5)",
             "\t".join(cols)]
    for i, r in enumerate(rows[:top], 1):
        lines.append("\t".join([str(i), str(r["start"]), str(r["end"]),
                                r["sense"], r["guide"],
                                str(r["uitei"]), str(r["reynolds"]), str(r["amarz"])]))
    return "\n".join(lines), rows


# ---------- S10: CLI ----------

def main(argv=None):
    """S10: python siRNA_design.py <fasta> [top]
    写 TTR_siRNA_candidates.tsv (全表); 错误进 stderr; 退出码 0/1/2
    """
    #raise NotImplementedError("S10")
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 1:
        print("用法: python siRNA_design.py <fasta> [top]", file=sys.stderr)
        return 2
    fasta_path = argv[0]
    try:
        top = int(argv[1]) if len(argv) > 1 else 20
    except ValueError:
        print(f"错误: top 必须是整数, 收到 {argv[1]!r}", file=sys.stderr)
        return 2

    try:
        results = rank_design(fasta_path)
    except (ValueError, FileNotFoundError) as e:
        print(f"错误: {e}", file=sys.stderr)
        return 1

    text_all = []
    for name, rows in results:
        text_all.append(format_table(name, rows)[0])

    out = "TTR_siRNA_candidates.tsv"
    with open(out, "w", encoding="utf-8", newline="") as fh:
        cols = ["rank", "transcript", "start", "end", "sense", "guide", "uitei", "reynolds", "amarz"]
        fh.write("\t".join(cols) + "\n")
        rank = 0
        total = 0
        for name, rows in results:
            for r in rows:
                rank += 1
                total += 1
                fh.write("\t".join([str(rank), name, str(r["start"]), str(r["end"]),
                                    r["sense"], r["guide"], str(r["uitei"]),
                                    str(r["reynolds"]), str(r["amarz"])]) + "\n")

    print("\n\n".join(text_all))
    print(f"\n已写出 {out} (全表 {total} 行)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
