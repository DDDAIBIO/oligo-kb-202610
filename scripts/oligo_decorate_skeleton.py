# -*- coding: utf-8 -*-
"""
oligo_decorate_skeleton.py — D16 骨架（2026-11-07 实操日）

用法：
  1. 把本文件复制一份命名为 oligo_decorate.py（在同一目录），你只改副本
  2. 逐个实现 9 个 TODO 函数（函数签名 / 返回值 = D15 冻结的合同，不许改）
  3. 跑测试：python test_expected.py
  4. 跑 CLI：python oligo_decorate.py --seq "mC*mGmUmAmCdTdCdGdAdTdCdGdTdAdAmUmGmUmA*mC"

分工约定（重要，别越界）：
  - parse_notation / summarize_residues ：已提供，直接用（解析是管道，不是今天的练习）
  - 9 个特征函数 ：你来写（这是 D15 §3 合同的实现）
  - report / main ：格式参考已给，输出字段不能少，排版随意

记号规则（v2，已冻结，解析器按此实现）：
  - 大写碱基字母：A C G U T I（I = 2'-O-甲基肌苷，plozasiran 约定，见 D6/D7）
  - 无前缀的碱基 = RNA（如 "G" "A" "U"）
  - "d" 前缀 = DNA（脱氧），如 dC dT —— gapmer 的 gap 写成 dC*dG*dT...
  - 小写前缀 = 2' 修饰 RNA：m = 2'-OMe，f = 2'-F，e = 2'-MOE，l = LNA
  - "*" = 两个残基之间是硫代磷酸（PS）；没有 * = 正常磷酸二酯（PO）
  - 非法输入（空串 / 非法字符 / 以 * 开头结尾 / ** 相连 / 以小写结尾）必须 raise ValueError
"""

import argparse
import sys

# ----------------------------------------------------------------------
# 解析器（已提供 —— 不要改，有 bug 告诉我，别自己悄悄修）
# ----------------------------------------------------------------------

VALID_PREFIX = {"m": "2'-OMe", "f": "2'-F", "e": "2'-MOE", "l": "LNA", "d": "DNA"}
VALID_BASE = set("ACGUTI")


def parse_notation(s: str):
    """把记号字符串解析成 (residues, ps_links)。

    residues : list[str]，每个元素如 "mC" "dA" "I" "G"
    ps_links : int，PS 连接数

    >>> parse_notation("mC*mG*dTdC")
    (['mC', 'mG', 'dT', 'dC'], 1)
    """
    if not s:
        raise ValueError("空序列")
    if s.startswith("*") or s.endswith("*") or "**" in s:
        raise ValueError("PS 记号 * 必须位于两个残基之间")
    residues = []
    ps = 0
    buf = ""
    for ch in s:
        if ch == "*":
            if not residues:
                raise ValueError("位置 %d 的 * 前没有残基" % len(residues))
            ps += 1
            buf = ""
        elif ch.isupper():
            if buf and buf[-1] not in VALID_PREFIX:
                raise ValueError("非法前缀: %r" % buf)
            if ch not in VALID_BASE:
                raise ValueError("非法碱基: %r" % ch)
            residues.append(buf + ch)
            buf = ""
        elif ch.islower():
            buf += ch
        else:
            raise ValueError("非法字符: %r" % ch)
    if buf:
        raise ValueError("以未完成的前缀结尾: %r" % buf)
    if not residues:
        raise ValueError("没有解析到任何残基")
    return residues, ps


def summarize_residues(residues):
    """已提供：把 residues 拆成三条平行信息（特征函数的公共输入）。

    返回 dict：
      bases      : list[str]  纯碱基序列（"mC"->"C"，"dT"->"T"，"I"->"I"）
      is_dna     : list[bool] 是否 DNA（以 d 开头）
      mod        : list[str|None] 2' 修饰（"m"/"f"/"e"/"l"/None）
      ps_positions: list[int]  PS 连接所在的下标 i（连接第 i 与第 i+1 个残基）
    """
    bases, is_dna, mod = [], [], []
    for r in residues:
        core = r
        if core.startswith("d"):
            is_dna.append(True)
            core = core[1:]
        else:
            is_dna.append(False)
        if core and core[0] in "mfel":
            mod.append(core[0])
            core = core[1:]
        else:
            mod.append(None)
        bases.append(core)
    return {
        "bases": bases,
        "is_dna": is_dna,
        "mod": mod,
        "ps_positions": [],  # TODO(main): 在 main 里把 ps_links 换算成位置填进来
    }


# ----------------------------------------------------------------------
# 下面 9 个特征函数 = 你今天的工作（签名与返回值 = D15 合同，不许改）
# 约定：bases 用 list[str]（不用 str），便于将来扩展修饰碱基；I 按 G 计入 GC
# ----------------------------------------------------------------------

def count_cpg(info: dict) -> int:
    """特征 1：DNA 区 CpG 计数 → TLR9 风险。

    只统计相邻两位都是 DNA（is_dna 两个 True）且碱基为 C,G 的位置数。
    RNA 区的 CG 不算（TLR9 认未甲基化 CpG DNA）。返回 int。
    """
    raise NotImplementedError  # TODO 1


def find_gu_motifs(bases: list) -> dict:
    """特征 2：GU/U-rich 基序 → TLR7/8 风险。

    v0 最小表（写死在函数里，注明"最小表，非完备"）：
      "GUCCUUCAA"、"UGUGU"、连续 U >= 3（键名固定用 "U3"）
    返回 {motif: [起始下标(0-based), ...]}，没有命中的 motif 不放进 dict。
    """
    raise NotImplementedError  # TODO 2


def find_au_rich(bases: list) -> list:
    """特征 3：AU-rich 元件 → 炎症/不稳定。

    搜 "AUUUA"，返回命中起始下标列表（0-based）。允许与特征 2 重叠，各自报告。
    """
    raise NotImplementedError  # TODO 3


def max_self_complement(bases: list) -> int:
    """特征 5：最长自互补区（粗筛双链倾向，非折叠预测）。

    v0：找最长的连续子串，使其与其反向互补链相同（回文式自互补），
    窗口长度限定 4–8。返回最长长度；找不到返回 0。
    提示：先写 reverse_complement（A<->U, C<->G，I 按 G、T 按 A 配对）。
    """
    raise NotImplementedError  # TODO 5


def gc_content(bases: list) -> float:
    """特征 6a：GC 含量（比例 0–1）。

    规则（写进注释）：C + G + I 都计入分子（I 按其配对行为归为 G 类）。
    返回 float，空序列 raise ValueError。
    """
    raise NotImplementedError  # TODO 6a


def tm_wallace(bases: list) -> float:
    """特征 6b：基础 Tm（占位公式，注明"占位"）。

    Wallace 规则（仅适合 <14 nt 的粗估）：Tm = 2*(A+U 数) + 4*(G+C+I 数)
    修饰 ΔTm 修正是 D11 欠账，v0 明确不做 —— 在报告里输出基础值并标注。
    """
    raise NotImplementedError  # TODO 6b


def find_hairpins(bases: list) -> list:
    """特征 7：简单发夹粗筛。

    v0：找 (stem 长度 4–8, loop 固定 4) 的模式，即位置 i..i+k-1 与
    j..j+k-1（j = i+k+4）互为反向互补。返回 [(i, k), ...]。
    """
    raise NotImplementedError  # TODO 7


def ps_profile(ps_positions: list) -> dict:
    """特征 8：PS 计数与分布。

    返回 {"total": int, "positions": [...]}。位置在 main 里已换算好，本函数只整理。
    """
    raise NotImplementedError  # TODO 8


def g_run_count(bases: list, k: int = 4) -> int:
    """特征 9：连续 >= k 个 G 的 run 数（G4 倾向粗筛）。返回 int。"""
    raise NotImplementedError  # TODO 9


# ----------------------------------------------------------------------
# 报告与 CLI（字段不能少；排版随意）
# ----------------------------------------------------------------------

def report(seq: str):
    residues, ps = parse_notation(seq)
    info = summarize_residues(residues)
    info["ps_positions"] = list(range(ps))  # TODO(main): 改成真实位置——
    # 提示：parse_notation 现在只返回 PS 总数。两条路任选：
    #   a) 改 parse_notation 让它同时返回位置（推荐，改provided要经过我）；
    #   b) 在 main 里重扫一遍 seq 数 * 的位置。
    # 这是一个设计决定，你选哪条路、为什么，写进注释。

    print("=== oligo_decorate 设计体检报告 ===")
    print("输入记号 :", seq)
    print("长度     :", len(residues), "nt")
    print("PS 数    :", ps)
    print("GC%      :", "%.1f" % (gc_content(info["bases"]) * 100))
    print("DNA-CpG  :", count_cpg(info), "(仅 DNA 区，TLR9)")
    print("GU/U 基序:", find_gu_motifs(info["bases"]), "(TLR7/8)")
    print("AUUUA    :", find_au_rich(info["bases"]))
    print("最长自互补:", max_self_complement(info["bases"]), "nt（粗筛）")
    print("基础 Tm  :", "%.1f" % tm_wallace(info["bases"]), "°C（Wallace 占位，未做修饰修正）")
    print("发夹粗筛 :", find_hairpins(info["bases"]))
    print("PS 分布  :", ps_profile(info["ps_positions"]))
    print("G-run(>=4):", g_run_count(info["bases"]))
    # TODO: 修饰统计一行（m/f/e/l 各多少、d 多少）——字段必须有，实现方式你定


def main():
    ap = argparse.ArgumentParser(description="寡核苷酸设计体检（v0，D16 练习）")
    ap.add_argument("--seq", required=True, help="记号序列，见文件头说明")
    args = ap.parse_args()
    try:
        report(args.seq)
    except ValueError as e:
        print("输入错误:", e, file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
