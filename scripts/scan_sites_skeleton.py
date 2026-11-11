# scan_sites.py 骨架 —— D21 代码日
# 规格来源: 02_技术方案/位点层_靶哪一段_D20.md §7 (10 条判据)
# 用法: 原位替换各 stub (保持签名不变); __main__ 块在文件末尾, 不要往后面追加
# v0 三决定: 变异用转录本坐标 / 覆盖率精确子串 / RNAplfold 走文件接口
# 记号: 窗口一律 0-based 半开 [start, end); is_seed 的 pos 是窗口内 1-based
# 上机时把本文件重命名为 scan_sites.py, 并把 test_scan_sites.py 第 6 行
#   import scan_sites_skeleton as ss  改成  import scan_sites as ss

from pathlib import Path


# ---------- 纯函数区 (test_scan_sites.py 12 条断言的对象, B2 必须全绿) ----------

def make_windows(seq, L, step):
    """候选窗口列表。返回 [(start, end), ...], 0-based 半开, end-start == L。
    L > len(seq) 时返回 [] (空列表, 不是异常——调用方按"无候选位点"处理)。"""
    raise NotImplementedError


def coverage(seq, txs):
    """转录本覆盖率: 精确子串匹配 seq 在 txs 列表中的命中数 / len(txs)。
    txs 为空时返回 0.0。"""
    raise NotImplementedError


def mode_compartment(modality):
    """模态 -> 工作区室。入参取值: 'siRNA'|'gapmer'|'splice_aso'
    -> 'cytosol_mRNA' | 'nucleus_preMRNA'。
    siRNA 只打成熟 mRNA (胞质); gapmer 与剪接 ASO 可达 pre-mRNA (核内)。"""
    raise NotImplementedError


def compartment_ok(region, modality):
    """region 取值: 'CDS'|'3UTR'|'5UTR'|'intron'; modality 同上。
    规则 (D20 §3 表): siRNA 排除 intron; splice_aso 必须 intron;
    gapmer 全部允许。"""
    raise NotImplementedError


def is_seed(pos, L):
    """pos: 窗口内位置, 1-based; L: 窗口长。
    种子区 = 位置 2..8 (1-based, 含两端)。是否适用由调用方按模态决定,
    本函数只管几何判定 (与 L 无关)。"""
    raise NotImplementedError


def window_access(acc, s, e):
    """acc: 全转录本 P(unpaired) 数组 (0-based, 与序列等长);
    返回窗口 [s, e) 的算术均值。空切片返回 0.0。"""
    raise NotImplementedError


def pop_coverage(variants):
    """variants: [(maf, s, e), ...] 已按窗口过滤的变异。
    返回 1 - sum(maf), 下限 0.0 (不要返回负数)。
    独立性假设: 变异不连锁 (v0 接受此近似)。"""
    raise NotImplementedError


def variants_in_window(variants, s, e):
    """过滤: 返回与窗口 [s, e) 有重叠的 [(maf, s2, e2), ...] 子集。
    重叠判定: s2 < e and e2 > s (半开区间标准写法)。"""
    raise NotImplementedError


# ---------- I/O 区 (冒烟测试, 无固定期望值; B2 全绿后才碰这里) ----------

def load_transcripts(fasta_dir):
    """读目录下所有 .fa/.fasta -> dict{标题首词: 去空白后的序列}。
    记得去掉序列内部换行。"""
    raise NotImplementedError


def accessibility(fasta_path, workdir="."):
    """subprocess 调 RNAplfold -W 80 -L 40, 解析生成的 lunp:
    跳过 '#' 开头的行, 每行取 [位置, P_unpaired] (第 2 列)。
    返回 0-based 数组, 长度必须 == 序列长度 (用 assert 检查)。
    pos 是 1-based -> 数组下标 0-based, 换算处写注释 (review 重点)。"""
    raise NotImplementedError


def load_variants(tsv_path):
    """读 TSV (表头 pos/ref/alt/maf) -> [(maf, pos, pos+1), ...]
    pos 转 0-based 半开区间, 供 variants_in_window 使用。"""
    raise NotImplementedError


# ---------- 驱动区 (stretch goal: B4 完成后写) ----------

def rank(fasta_dir, variants_tsv, modality, L=21, step=5):
    """生成排序表 TSV (列见 D21 §7)。v0 简单加和打分, 不加权。"""
    raise NotImplementedError


if __name__ == "__main__":
    print("scan_sites skeleton — 全部函数待实现 (NotImplementedError)")
    print("先跑: python test_scan_sites.py  纯函数区全绿后, 再接 I/O 区")
