"""D28 test_siRNA_design -- Part A 纯函数断言 + Part B 冒烟

期望值全部经过参考实现预核验 (2026-11-18), 包括两个手算坑:
  - Reynolds 判据7/8 是纯罚分项, 无 G/C19/G13 不加分 -> patisiran sense = 6 (压线)
  - Amarzguioui ⑤⑥ 是罚分不是反向奖分 -> patisiran guide = 2 (<cutoff 3, 真实分歧!)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import siRNA_design_skeleton as sd

PASSED = FAILED = 0


def check(name, got, want):
    global PASSED, FAILED
    ok = got == want
    PASSED, FAILED = PASSED + ok, FAILED + (not ok)
    print(f"{'PASS' if ok else 'FAIL'} {name}: got {got!r}, want {want!r}")


# ---- 期望锚点 ----
SENSE = "GUAACCAAGAGUAUUCCAU"      # patisiran sense 19 nt (双源核验)
GUIDE = "AUGGAAUACUCUUGGUUAC"      # patisiran guide = revcomp(SENSE)


def t_revcomp():
    check("T0 revcomp 基础", sd._revcomp_rna("AUCG"), "CGAU")
    check("T0b revcomp 空串", sd._revcomp_rna(""), "")
    check("T0c revcomp 小写", sd._revcomp_rna("aucg"), "CGAU")


def t_gc():
    check("T1 gc_frac", round(sd.gc_frac("AAACCC"), 4), 0.5)
    check("T1b gc_frac 空串", sd.gc_frac(""), 0.0)


def t_uitei():
    check("T2 uitei 全过 (A*17+GC)", sd.uitei_pass("A" * 17 + "GC"), True)
    check("T2b uitei fail (i) pos1=G", sd.uitei_pass("G" + "A" * 16 + "C"), False)
    check("T2c uitei fail (iii) AU<4/7", sd.uitei_pass("AACCCCCAAAAAAAAACCC"), False)
    check("T2d uitei fail (iv) GC run=10", sd.uitei_pass("AAAUUAG" + "C" * 10 + "AC"), False)
    check("T2e uitei patisiran guide", sd.uitei_pass(GUIDE), True)


def t_reynolds():
    # 罚分项语义: 有 G/C19 / G13 才 -1, 没有不加分
    check("T3 patisiran sense = 6 (压线)", sd.reynolds_score(SENSE), 6)
    check("T3b all-G = -2 (可负, 不 clamp)",
          sd.reynolds_score("G" * 19), -2)


def t_amarz():
    check("T4 patisiran guide = 2 (<cutoff 3, 真实分歧!)", sd.amarzguioui_score(GUIDE), 2)
    check("T4b asym>0 例 = 3", sd.amarzguioui_score("AUUUU" + "C" * 11 + "GCC"), 3)
    # pos1=G: 失去 +1 且罚 -1 -> 相对 T4b 净 -2
    check("T4c pos1=G = 1", sd.amarzguioui_score("GUUUU" + "C" * 11 + "GCC"), 1)


def t_design_e2e():
    """T5: 方向断言 -- patisiran 埋进测试转录本"""
    mini = "ACGUACGUAC" + SENSE + "UGCAUGCAUG"
    rows = sd.design_candidates(mini)
    hit = [r for r in rows if r["sense"] == SENSE]
    check("T5 靶位点被枚举到", len(hit), 1)
    if hit:
        r = hit[0]
        check("T5b guide = patisiran antisense (方向!)", r["guide"], GUIDE)
        check("T5c 该窗 uitei", r["uitei"], True)
        check("T5d 该窗 reynolds", r["reynolds"], 6)
        check("T5e 该窗 amarz", r["amarz"], 2)
        check("T5f start", r["start"], 10)
    # 反向情形: 若实现把 guide 写成窗口本身, 下面必 FAIL
    check("T5g guide != sense (链方向哨兵)", r["guide"] == r["sense"], False)


def t_sort():
    mini = "ACGUACGUAC" + SENSE + "UGCAUGCAUG"
    rows = sd.design_candidates(mini)
    keys = [(-r["reynolds"], -r["amarz"], r["start"]) for r in rows]
    check("T6 design_candidates 未排序(保持窗口顺序)", keys == sorted(keys), False)


# ---- Part A 跑批 ----

def part_a():
    for t in (t_revcomp, t_gc, t_uitei, t_reynolds, t_amarz, t_design_e2e, t_sort):
        try:
            t()
        except NotImplementedError as e:
            global FAILED
            FAILED += 1
            print(f"SKIP-FAIL {t.__name__}: {e}")
    print(f"\nPart A: {PASSED} PASS / {FAILED} FAIL")
    return FAILED == 0


if __name__ == "__main__":
    ok = part_a()
    sys.exit(0 if ok else 1)
