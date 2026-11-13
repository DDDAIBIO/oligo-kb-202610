# scan_sites.py 骨架 —— D21 代码日
# 规格来源: 02_技术方案/位点层_靶哪一段_D20.md §7 (10 条判据)
# 用法: 原位替换各 stub (保持签名不变); __main__ 块在文件末尾, 不要往后面追加
# v0 三决定: 变异用转录本坐标 / 覆盖率精确子串 / RNAplfold 走文件接口
# 记号: 窗口一律 0-based 半开 [start, end); is_seed 的 pos 是窗口内 1-based
# 上机时把本文件重命名为 scan_sites.py, 并把 test_scan_sites.py 第 6 行
#   import scan_sites_skeleton as ss  改成  import scan_sites as ss

from pathlib import Path

import sys
ss = sys.modules[__name__]


# ---------- 纯函数区 (test_scan_sites.py 12 条断言的对象, B2 必须全绿) ----------

def make_windows(seq, L, step):
    """候选窗口列表。返回 [(start, end), ...], 0-based 半开, end-start == L。
    L > len(seq) 时返回 [] (空列表, 不是异常——调用方按"无候选位点"处理)。"""
    #raise NotImplementedError
    if L > len(seq):
        return []
    return [(s, s + L) for s in range(0, len(seq) - L + 1, step)]

def coverage(seq, txs):
    """转录本覆盖率: 精确子串匹配 seq 在 txs 列表中的命中数 / len(txs)。
    txs 为空时返回 0.0。"""
    #raise NotImplementedError
    if not seq:
        raise ValueError("空序列: coverage 的 seq 不能为空")
    if not txs:
        return 0.0
    return sum(1 for t in txs if seq in t) / len(txs)

COMPT_MAP = {
    "siRNA": "cytosol_mRNA",
    "gapmer": "nucleus_preMRNA",
    "splice_aso": "nucleus_preMRNA",
}

def mode_compartment(modality):
    """模态 -> 工作区室。入参取值: 'siRNA'|'gapmer'|'splice_aso'
    -> 'cytosol_mRNA' | 'nucleus_preMRNA'。
    siRNA 只打成熟 mRNA (胞质); gapmer 与剪接 ASO 可达 pre-mRNA (核内)。"""
    #raise NotImplementedError
    if modality not in COMPT_MAP:
        raise ValueError(f"未知模态 {modality!r}: v0 仅支持 {sorted(COMPT_MAP)}")
    return COMPT_MAP[modality]

REGIONS = ("CDS", "3UTR", "5UTR", "intron")
def compartment_ok(region, modality):
    """region 取值: 'CDS'|'3UTR'|'5UTR'|'intron'; modality 同上。
    规则 (D20 §3 表): siRNA 排除 intron; splice_aso 必须 intron;
    gapmer 全部允许。"""
    #raise NotImplementedError
    if region not in REGIONS:
        raise ValueError(f"未知区域 {region!r}: v0 仅支持 {REGIONS}")
    if modality == "splice_aso":
        return region == "intron"     # 必须打 intron: 剪接信号只存在于 pre-mRNA
    if modality == "gapmer":
        return True                   # RNase H1 在核内, pre-mRNA / mRNA 都可达
    if modality == "siRNA":
        return region != "intron"     # 胞质里只有成熟 mRNA, 内含子已被切走
    raise ValueError(f"未知模态 {modality!r}: v0 仅支持 siRNA/gapmer/splice_aso")

def is_seed(pos, L):
    """pos: 窗口内位置, 1-based; L: 窗口长。
    种子区 = 位置 2..8 (1-based, 含两端)。是否适用由调用方按模态决定,
    本函数只管几何判定 (与 L 无关)。"""
    #raise NotImplementedError
    if pos < 1 or pos > L:
        raise ValueError(f"位置 {pos} 超出窗口 1..{L}")
    return 2 <= pos <= 8

def window_access(acc, s, e):
    """acc: 全转录本 P(unpaired) 数组 (0-based, 与序列等长);
    返回窗口 [s, e) 的算术均值。空切片返回 0.0。"""
    #raise NotImplementedError
    seg = acc[s:e]
    if not seg:
        return 0.0
    return sum(seg) / len(seg)

def pop_coverage(variants):
    """variants: [(maf, s, e), ...] 已按窗口过滤的变异。
    返回 1 - sum(maf), 下限 0.0 (不要返回负数)。
    独立性假设: 变异不连锁 (v0 接受此近似)。"""
    #raise NotImplementedError
    return max(0.0, 1.0 - sum(m for m, _s, _e in variants))

def variants_in_window(variants, s, e):
    """过滤: 返回与窗口 [s, e) 有重叠的 [(maf, s2, e2), ...] 子集。
    重叠判定: s2 < e and e2 > s (半开区间标准写法)。"""
    #raise NotImplementedError
    return [(m, s2, e2) for (m, s2, e2) in variants if s2 < e and e2 > s]

# ---------- I/O 区 (冒烟测试, 无固定期望值; B2 全绿后才碰这里) ----------
import glob, os, subprocess


def _read_fasta(path):
    """单文件 FASTA -> dict{标题首词: 序列}。序列按行拼接去换行。"""
    txs, name, chunks = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None and chunks:      # 多记录: 先落上一条
                    txs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []   # 标题首词
            else:
                chunks.append(line.upper())
    if name is not None and chunks:
        txs[name] = "".join(chunks)
    return txs

def load_transcripts(fasta_dir):
    """读目录下所有 .fa/.fasta -> dict{标题首词: 去空白后的序列}。
    记得去掉序列内部换行。"""
    #raise NotImplementedError
    paths = sorted(glob.glob(os.path.join(fasta_dir, "*.fa"))
                   + glob.glob(os.path.join(fasta_dir, "*.fasta")))
    if not paths:
        raise FileNotFoundError(f"{fasta_dir} 下没有 .fa/.fasta")
    txs = {}
    for p in paths:
        txs.update(_read_fasta(p))
    if not txs:
        raise ValueError(f"{fasta_dir} 下的 FASTA 没有序列")
    return txs

def _parse_lunp(path, seq_len):
    """解析 lunp -> 0-based P(unpaired) 数组。缺位置/越界都显式报错。"""
    by_pos = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):     # 跳过表头/注释
                continue
            cols = line.split()
            if len(cols) < 2:
                continue
            by_pos[int(cols[0])] = float(cols[1])    # 第 1 列 pos(1-based), 第 2 列 P

    if not by_pos:
        raise ValueError(f"{path} 未解析到任何数据行")
    missing = [p for p in range(1, seq_len + 1) if p not in by_pos]
    if missing:
        raise ValueError(f"{path}: {len(missing)}/{seq_len} 位置缺失 "
                         f"(前几个 {missing[:5]}); RNAplfold 只在窗口覆盖到的位置给值")
    if max(by_pos) > seq_len:
        raise ValueError(f"{path}: 最大位置 {max(by_pos)} 超出序列长 {seq_len}")
    # ⭐ pos 是 1-based, 数组下标 0-based: acc[pos-1] = P
    return [by_pos[p] for p in range(1, seq_len + 1)]

def accessibility(fasta_path, workdir="."):
    """subprocess 调 RNAplfold -W 80 -L 40 -u 1 (stdin 喂入), 解析 lunp。"""
    fasta_path = os.path.abspath(fasta_path)
    txs = _read_fasta(fasta_path)
    if len(txs) != 1:
        raise ValueError(f"{fasta_path} 应含且仅含 1 条序列, 实际 {len(txs)} 条")
    seq_name = next(iter(txs))            # ⭐ 产物前缀取自 FASTA 头行首词, 不是文件名
    seq = txs[seq_name]

    lunp = os.path.join(workdir, seq_name + "_lunp")

    with open(fasta_path, "rb") as fh:    # ⭐ stdin 喂入: 2.7.2 会忽略位置参数并死等 stdin
        r = subprocess.run(["RNAplfold", "-W", "80", "-L", "40", "-u", "1"],
                           cwd=workdir, stdin=fh, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"RNAplfold 失败({r.returncode}): {r.stderr.strip()[:400]}")
    if not os.path.exists(lunp):
        raise FileNotFoundError(f"没找到 {lunp}; 实际产物: {sorted(os.listdir(workdir))}; "
                                f"stderr: {r.stderr.strip()[:400]}")

    acc = _parse_lunp(lunp, len(seq))
    assert len(acc) == len(seq), f"acc 长 {len(acc)} != 序列长 {len(seq)}"
    return acc

"""
def accessibility(fasta_path, workdir="."):
    subprocess 调 RNAplfold -W 80 -L 40, 解析生成的 lunp:
    跳过 '#' 开头的行, 每行取 [位置, P_unpaired] (第 2 列)。
    返回 0-based 数组, 长度必须 == 序列长度 (用 assert 检查)。
    pos 是 1-based -> 数组下标 0-based, 换算处写注释 (review 重点)
    #raise NotImplementedError
    fasta_path = os.path.abspath(fasta_path)
    txs = _read_fasta(fasta_path)
    if len(txs) != 1:
        raise ValueError(f"{fasta_path} 应含且仅含 1 条序列, 实际 {len(txs)} 条")
    seq = next(iter(txs.values()))

    base = os.path.splitext(os.path.basename(fasta_path))[0]
    lunp = os.path.join(workdir, base + "_lunp")

    r = subprocess.run(["RNAplfold", "-W", "80", "-L", "40", "-u", "1", fasta_path],
                       cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"RNAplfold 失败({r.returncode}): {r.stderr.strip()[:400]}")
    if not os.path.exists(lunp):
        raise FileNotFoundError(f"没找到 {lunp}; 实际产物: {sorted(os.listdir(workdir))}")

    acc = _parse_lunp(lunp, len(seq))
    assert len(acc) == len(seq), f"acc 长 {len(acc)} != 序列长 {len(seq)}"
    return acc
"""

def load_variants(tsv_path):
    """读 TSV (表头 pos/ref/alt/maf) -> [(maf, pos, pos+1), ...]
    pos 转 0-based 半开区间, 供 variants_in_window 使用。"""
    #raise NotImplementedError
    out = []
    with open(tsv_path, encoding="utf-8") as fh:
        header = None
        ip = im = None
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            cols = [c.strip() for c in line.split("\t")]
            if header is None:
                header = [c.lower() for c in cols]
                need = {"pos", "maf"}
                if not need <= set(header):
                    raise ValueError(f"表头缺列 {need - set(header)}; 实际 {header}")
                ip, im = header.index("pos"), header.index("maf")
                continue
            pos, maf = int(cols[ip]), float(cols[im])
            if pos < 1:
                raise ValueError(f"pos 应 >= 1 (1-based), 得到 {pos}")
            if not 0.0 <= maf <= 1.0:
                raise ValueError(f"MAF 越界 {maf} @pos{pos}")
            out.append((maf, pos - 1, pos))     # 1-based -> 0-based 半开 [pos-1, pos)
    return out


# ---------- 驱动区 (stretch goal: B4 完成后写) ----------
# ---------- 驱动区 ----------

def _score_dim_coverage(cov):
    if cov < 0.5: return 0
    if cov < 0.9: return 1
    return 2


def _score_dim_snp(win_variants, seed_variants):
    """种子区 MAF>=5% -> 0; 窗内 MAF>=1% -> 1; 其余 -> 2
    v0.1: 原来写的是"窗内有变异就扣 1 分", 于是 MAF 0.001 与 MAF 0.40 同罚,
          MAF 大小在非种子区完全不起作用。改为按 MAF 分档 (1% 以下视为不影响疗效)。"""
    if any(m >= 0.05 for m, _s, _e in seed_variants):
        return 0
    if any(m >= 0.01 for m, _s, _e in win_variants):
        return 1
    return 2


def _score_dim_seq(seq, gc):
    """同碱基连跑>=5 或 GC<25%/>75% -> 0; GC<30%/70% -> 1; 否则 2"""
    run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        if run >= 5:
            return 0
    if gc < 0.25 or gc > 0.75:
        return 0
    if gc < 0.30 or gc > 0.70:
        return 1
    return 2


def rank(fasta_dir, variants_tsv, modality, L=21, step=5):
    """生成排序表 TSV (列见 D21 §7)。v0 简单加和打分, 不加权。"""
    if modality == "splice_aso":
        raise ValueError("v0 的 FASTA 是成熟 mRNA (无内含子), splice_aso 无候选位点; "
                         "需要 pre-mRNA + 注释, 排 v1")

    txs = ss.load_transcripts(fasta_dir)
    if len(txs) != 1:
        raise ValueError(f"v0 只支持单条转录本, 实际 {len(txs)} 条")
    name, seq = next(iter(txs.items()))

    variants = ss.load_variants(variants_tsv)
    if variants:                                 # v0.1: 坐标体系校验, 防"TSV 对旧 fasta 造的"静默跑空
        worst = max(e for _m, _s, e in variants)
        if worst > len(seq):
            raise ValueError(
                f"{variants_tsv} 最大 pos={worst} 超出 {name} 长度 {len(seq)}: "
                f"坐标不匹配 (基因组坐标? 别的转录本/异构体?), 拒绝静默跑空")
    fa_files = sorted(p for p in os.listdir(fasta_dir) if p.endswith((".fa", ".fasta")))
    if not fa_files:
        raise FileNotFoundError(f"{fasta_dir} 下没有 .fa/.fasta")
    fa = os.path.join(fasta_dir, fa_files[0])
    acc = ss.accessibility(fa)

    windows = ss.make_windows(seq, L, step)
    if not windows:
        raise ValueError(f"序列 {len(seq)} nt 短于窗口 {L}")

    means = sorted(ss.window_access(acc, s, e) for s, e in windows)
    n = len(means)
    p30, p70 = means[int(0.30 * (n - 1))], means[int(0.70 * (n - 1))]
    tx_list = list(txs.values())

    rows = []
    for s, e in windows:
        wseq = seq[s:e]
        gc = sum(b in "GC" for b in wseq) / L
        mean_acc = ss.window_access(acc, s, e)
        win_vars = ss.variants_in_window(variants, s, e)
        seed_vars = ss.variants_in_window(variants, s + 1, s + 8)
        pop_cov = ss.pop_coverage(win_vars)

        dim_cov = _score_dim_coverage(ss.coverage(wseq, tx_list))
        dim_comp = 2    # v0 无注释: 整条按 CDS 处理 (已知简化)
        if mean_acc < p30:
            dim_acc = 0
        elif mean_acc > p70:
            dim_acc = 2
        else:
            dim_acc = 1
        dim_snp = _score_dim_snp(win_vars, seed_vars)
        dim_seq = _score_dim_seq(wseq, gc)

        rows.append({
            "start": s, "end": e, "seq": wseq,
            "gc": f"{gc * 100:.1f}", "mean_access": f"{mean_acc:.3f}",
            "seed_hits": len(seed_vars), "pop_cov": f"{pop_cov:.3f}",
            "snp": ",".join(f"{s2 + 1}:{m:.3f}" for m, s2, _e2 in win_vars) or "-",
            "acc_f": mean_acc,      # tie-breaker 用的原始浮点, 不进 TSV
            "score": dim_cov + dim_comp + dim_acc + dim_snp + dim_seq,
        })

    # v0.1: 同分时按可及性降序 —— 原来只看 start, 10 个满分窗口的 rank 其实是位置顺序
    rows.sort(key=lambda r: (-r["score"], -r["acc_f"], r["start"]))
    return name, rows


import csv

COLS = ["rank", "start", "end", "seq", "gc", "mean_access",
        "seed_hits", "pop_cov", "snp", "score"]


def format_table(name, rows, top=None):
    """rows -> (打印用字符串, TSV 行列表)。不碰文件系统, 可测。
    top=None 取全表 (写盘用), top=N 只取前 N 行 (屏幕摘要用)。"""
    lines = [f"转录本: {name}   候选窗口: {len(rows)}   (v0 五维简单加和, 区室维恒 2)"]
    lines.append("\t".join(COLS))
    tsv = []
    for i, r in enumerate(rows if top is None else rows[:top], 1):
        tsv.append({"rank": i, **{k: r[k] for k in COLS[1:]}})
        lines.append("\t".join(str(tsv[-1][c]) for c in COLS))
    return "\n".join(lines), tsv


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 3:
        print("用法: python scan_sites.py <fasta_dir> <variants.tsv> <modality> [L] [step] [top]",
              file=sys.stderr)
        return 2
    fasta_dir, variants_tsv, modality = argv[0], argv[1], argv[2]
    L = int(argv[3]) if len(argv) > 3 else 21
    step = int(argv[4]) if len(argv) > 4 else 5
    top = int(argv[5]) if len(argv) > 5 else 10

    try:
        name, rows = rank(fasta_dir, variants_tsv, modality, L=L, step=step)
    except (ValueError, FileNotFoundError, RuntimeError, AssertionError) as e:
        print(f"错误: {e}", file=sys.stderr)
        return 1

    text, _ = format_table(name, rows, top=top)        # 屏幕只打摘要
    print(text)
    _, tsv_rows = format_table(name, rows, top=None)   # v0.1: 全表落盘, 否则变异窗口查不到
    out = "TTR_sites_ranked.tsv"
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(tsv_rows)
    n_var = sum(1 for r in rows if r["snp"] != "-")
    print(f"\n已写出 {out} (全表 {len(tsv_rows)} 行, 其中 {n_var} 行进过变异筛查; "
          f"屏幕只打 top {top})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""
def _score_dim_coverage(cov):
    if cov < 0.5: return 0
    if cov < 0.9: return 1
    return 2

def _score_dim_snp(win_variants, seed_variants):
    种子区 MAF>=5% -> 0; 窗内有变异 -> 1; 无 -> 2
    if any(m >= 0.05 for m, _s, _e in seed_variants):
        return 0
    if win_variants:
        return 1
    return 2

def _score_dim_seq(seq, gc):
    同碱基连跑>=5 或 GC<25%/>75% -> 0; GC<30%/70% -> 1; 否则 2
    run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        if run >= 5:
            return 0
    if gc < 0.25 or gc > 0.75:
        return 0
    if gc < 0.30 or gc > 0.70:
        return 1
    return 2

def rank(fasta_dir, variants_tsv, modality, L=21, step=5):
    生成排序表 TSV (列见 D21 §7)。v0 简单加和打分, 不加权。
    #raise NotImplementedError
    if modality == "splice_aso":
        raise ValueError("v0 的 FASTA 是成熟 mRNA (无内含子), splice_aso 无候选位点; "
                         "需要 pre-mRNA + 注释, 排 v1")

    txs = ss.load_transcripts(fasta_dir)
    if len(txs) != 1:
        raise ValueError(f"v0 只支持单条转录本, 实际 {len(txs)} 条")
    name, seq = next(iter(txs.items()))

    variants = ss.load_variants(variants_tsv)
    #fa = os.path.join(fasta_dir, sorted(os.listdir(fasta_dir))[0])
    fa_files = sorted(p for p in os.listdir(fasta_dir) if p.endswith((".fa", ".fasta")))
if not fa_files:
    raise FileNotFoundError(f"{fasta_dir} 下没有 .fa/.fasta")
    fa = os.path.join(fasta_dir, fa_files[0])
    acc = ss.accessibility(fa)

    windows = ss.make_windows(seq, L, step)
    if not windows:
        raise ValueError(f"序列 {len(seq)} nt 短于窗口 {L}")

    # 可及性 30/70 分位阈值 (在全窗口均值上算, 不是全碱基)
    means = sorted(ss.window_access(acc, s, e) for s, e in windows)
    n = len(means)
    p30, p70 = means[int(0.30 * (n - 1))], means[int(0.70 * (n - 1))]
    tx_list = list(txs.values())      # v1 换成多异构体列表, 这里不用改

    rows = []
    for s, e in windows:
        wseq = seq[s:e]
        gc = sum(b in "GC" for b in wseq) / L
        mean_acc = ss.window_access(acc, s, e)
        win_vars = ss.variants_in_window(variants, s, e)
        seed_vars = ss.variants_in_window(variants, s + 1, s + 8)  # 种子 1-based 2..8 -> 0-based [s+1, s+8)
        pop_cov = ss.pop_coverage(win_vars)

        dim_cov = _score_dim_coverage(ss.coverage(wseq, tx_list))   # 逐窗口算, v1 直接复用
        dim_comp = 2    # v0 无注释: 整条按 CDS 处理 (已知简化, 见下)
        if mean_acc < p30:    dim_acc = 0
        elif mean_acc > p70:  dim_acc = 2
        else:                 dim_acc = 1
        dim_snp = _score_dim_snp(win_vars, seed_vars)
        dim_seq = _score_dim_seq(wseq, gc)

        rows.append({
            "start": s, "end": e, "seq": wseq,
            "gc": f"{gc * 100:.1f}", "mean_access": f"{mean_acc:.3f}",
            "seed_hits": len(seed_vars), "pop_cov": f"{pop_cov:.3f}",
            "score": dim_cov + dim_comp + dim_acc + dim_snp + dim_seq,
        })

    rows.sort(key=lambda r: (-r["score"], r["start"]))   # 分数降序, 平手按位置(保证可复现)
    return name, rows

import csv
import sys

COLS = ["rank", "start", "end", "seq", "gc", "mean_access", "seed_hits", "pop_cov", "score"]


def format_table(name, rows, top=10):
    rows -> (打印用字符串, TSV 行列表)。不碰文件系统, 可测
    lines = [f"转录本: {name}   候选窗口: {len(rows)}   (v0 五维简单加和, 区室维恒 2)"]
    lines.append("\t".join(COLS))
    tsv = []
    for i, r in enumerate(rows[:top], 1):
        tsv.append({"rank": i, **{k: r[k] for k in COLS[1:]}})
        lines.append("\t".join(str(tsv[-1][c]) for c in COLS))
    return "\n".join(lines), tsv


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 3:
        print("用法: python scan_sites.py <fasta_dir> <variants.tsv> <modality> [L] [step]",
              file=sys.stderr)
        return 2
    fasta_dir, variants_tsv, modality = argv[0], argv[1], argv[2]
    L = int(argv[3]) if len(argv) > 3 else 21
    step = int(argv[4]) if len(argv) > 4 else 5

    try:
        name, rows = rank(fasta_dir, variants_tsv, modality, L=L, step=step)
    except (ValueError, FileNotFoundError, RuntimeError, AssertionError) as e:
        print(f"错误: {e}", file=sys.stderr)
        return 1

    text, tsv_rows = format_table(name, rows)
    print(text)
    out = "TTR_sites_ranked.tsv"
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(tsv_rows)
    print(f"\n已写出 {out} (top {len(tsv_rows)}/{len(rows)})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#if __name__ == "__main__":
#    print("scan_sites skeleton — 全部函数待实现 (NotImplementedError)")
#    print("先跑: python test_scan_sites.py  纯函数区全绿后, 再接 I/O 区")
"""
