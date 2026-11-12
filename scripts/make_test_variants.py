# make_test_variants.py —— 造 v0 的 variants.tsv (D21 B1, 5 行假数据)
# 用法: python make_test_variants.py ttr.fa [pos1 pos2 ...]
#   默认: 从序列均匀取 5 个位置。MAF 按预设梯度, 刻意覆盖打分分支:
#   >=5% (触发 SNP 维 0 分) / 1-5% (1 分) / <1% (无常见变异, 2 分)
# 注意: pos 是 1-based 转录本坐标; ref/alt 从 FASTA 实际碱基取, alt 换成另一种。
# 这是假数据 —— 真 SNP 走 gnomAD, 那是 v1 (见 D20/D21)。

import sys

MAF_PRESETS = [0.40, 0.05, 0.02, 0.001, 0.003]   # 梯度覆盖三个打分分支
BASES = "ACGU"


def main(argv):
    if len(argv) < 2:
        print("用法: python make_test_variants.py <fasta> [pos1 pos2 ...]", file=sys.stderr)
        return 2
    fasta_path = argv[1]

    seq_lines = []
    for line in open(fasta_path, encoding="utf-8"):
        line = line.strip()
        if line and not line.startswith(">"):
            seq_lines.append(line.upper())
    seq = "".join(seq_lines)
    if not seq:
        raise ValueError(f"{fasta_path} 没有序列")

    if len(argv) > 2:
        positions = [int(p) for p in argv[2:]]
    else:  # 默认: 均匀取 5 个位置 (避开两端)
        span = len(seq) - 40
        positions = [20 + i * span // 5 + 1 for i in range(5)]
    positions = sorted(set(positions))

    if len(positions) != len(MAF_PRESETS):
        print(f"提示: {len(positions)} 个位置 vs {len(MAF_PRESETS)} 个 MAF 预设, "
              f"多出的位置循环用预设", file=sys.stderr)

    rows = []
    for i, pos in enumerate(positions):
        if not 1 <= pos <= len(seq):
            raise ValueError(f"pos {pos} 超出序列长 {len(seq)} (1-based)")
        ref = seq[pos - 1]                       # 1-based -> 0-based 取实际碱基
        alt = next(b for b in BASES if b != ref)
        maf = MAF_PRESETS[i % len(MAF_PRESETS)]
        rows.append((pos, ref, alt, maf))

    out = "variants.tsv"
    with open(out, "w", encoding="utf-8", newline="") as fh:
        fh.write("pos\tref\talt\tmaf\n")
        for pos, ref, alt, maf in rows:
            fh.write(f"{pos}\t{ref}\t{alt}\t{maf}\n")
    print(f"已写出 {out} (基于 {fasta_path}, {len(seq)} nt):", file=sys.stderr)
    for pos, ref, alt, maf in rows:
        print(f"  pos={pos:>6}  {ref}>{alt}  maf={maf}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
