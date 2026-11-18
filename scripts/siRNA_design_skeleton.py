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
    raise NotImplementedError("S0")


def gc_frac(seq):
    """S1: GC 占比 0-1; 空串返回 0.0"""
    raise NotImplementedError("S1")


def _au(seq):
    """内部工具: A/U 计数 (可直接实现, 也可留给实现者)"""
    return sum(1 for c in seq.upper() if c in "AU")


# ---------- S2: Ui-Tei (过滤器, AND 逻辑, guide 坐标系) ----------

def uitei_pass(guide):
    """S2: 四条件 AND -> bool
    (i)   guide[0] in AU
    (ii)  guide[18] in GC
    (iii) guide[0:7] A/U >= 4   (siDirect 口径 >=4/7; 原文摘要为 >=5, 见 D26 §4 漂移)
    (iv)  无 >=10 连续 GC
    """
    raise NotImplementedError("S2")


# ---------- S3-S4: Reynolds (排序器, sense 坐标系) ----------

def _c3_low_internal_stability(sense):
    """S3 (代理, 🟡): 中心 7-mer sense[6:13] A/U >= 4 -> 1, 否则 0"""
    raise NotImplementedError("S3")


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
    raise NotImplementedError("S4")


# ---------- S5-S6: Amarzguioui (分列展示, guide 坐标系) ----------

def _amarz_asymmetry(guide):
    """S5 (代理, 🟡): guide[0:3] 与 guide[16:19] 的 A/U 数差 >0 -> 1, 否则 0"""
    raise NotImplementedError("S5")


def amarzguioui_score(guide):
    """S6: 六条 +/- 体系
    +1 asymmetry(S5); +1 guide[18] in GC; +1 guide[0] in AU; +1 guide[13] in AU;
    -1 guide[18]=='U'; -1 guide[0]=='G'
    """
    raise NotImplementedError("S6")


# ---------- S7-S9: 候选生成与排序 ----------

def design_candidates(seq, step=1):
    """S7: 全枚举 19-mer; 每窗 -> dict(start, end, sense, guide, uitei, reynolds, amarz)
    sense = 转录本窗口本身; guide = _revcomp_rna(sense)
    """
    raise NotImplementedError("S7")


def rank_design(fasta_path, top=20):
    """S8: _read_fasta -> 逐转录本 S7 -> 排序键 (-reynolds, -amarz, start)
    返回 [(name, rows)]; 排序注释必须写明"不合成总分(D26 §5)"
    """
    raise NotImplementedError("S8")


def format_table(name, rows, top=20):
    """S9: 表头 rank start end sense guide uitei reynolds amarz
    屏幕打 top 行; 返回 (text, all_rows) -- 全表交 main 写盘 (D21 教训)
    """
    raise NotImplementedError("S9")


# ---------- S10: CLI ----------

def main(argv=None):
    """S10: python siRNA_design.py <fasta> [top]
    写 TTR_siRNA_candidates.tsv (全表); 错误进 stderr; 退出码 0/1/2
    """
    raise NotImplementedError("S10")


if __name__ == "__main__":
    sys.exit(main())
