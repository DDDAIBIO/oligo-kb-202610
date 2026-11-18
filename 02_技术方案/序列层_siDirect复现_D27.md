# D27 ｜ 序列层（三）：siDirect 2.0 在线复现 + 与 D21 排序对齐（对比作业）

> 日期：2026-11-17/18（W4-D27）｜ 预算 45 min ｜ 前置：D25 Reynolds + D26 三规则  
> **今天不学新规则，做一次校准**：拿 Ui-Tei 组的官方实现（siDirect 2.0）当"标准答案"，回答两个问题——① 它给的 top10 和我手动按规则推的是否一致？② 规则排序（siDirect）和可及性排序（D21 `scan_sites`）的重合度如何？后者是"可及性要不要进序列层打分"（第六维度伏笔）的第一手证据。

---

## 0. siDirect 2.0 是什么、为什么信它

- Ui-Tei/Naito 组官方在线工具（`http://sidirect2.rnai.jp/`），D26 三套规则的**权威实现**——之前发现的"≥5 vs ≥4"漂移就发生在它身上，今天顺带验证它用的是哪个口径。
- 三规则打分（Ui-Tei / Amarzguioui / Reynolds 全给）+ **seed Tm 过滤**：默认 seed-target 双链 Tm < 21.5℃（这是 Ui-Tei 2008 的 B 类脱靶防线——**D22 seed 富集、D26 量纲讨论的实体化**）。
- 序列长度限制 🟡（二手文献称 600 bp，回官网确认）——所以今天只喂 200 nt 片段。

## 1. 输入片段（已备好，直接复制）

ttr.fa 全长 616 nt，取 **1-based 240–439**（覆盖 D21 发现的高可及性岛 260–331），GC 53.5%：

```text
TCTGGAGAGCTGCATGGGCTCACAACTGAGGAGGAATTTGTAGAAGGGATATACAAAGTGGAAATAGACACCAAATCTTACTGGAAGGCACTTGGCATCTCCCCATTCCATGAGCATGCAGAGGTGGTATTCACAGCCAACGACTCCGGCCCCCGCCGCTACACCATTGCCGCCCTGCTGAGCCCCTACTCCTATTCCAC
```

## 2. 操作记录模板（跑完填进 §2 表）

**输入参数（逐项记录，别凭记忆）**：

| 参数           | 我用的值      | 官网默认      | 备注                        |
| ------------ | --------- | --------- | ------------------------- |
| seed Tm 上限   | 21.5      | 21.5℃     | 调没调要记录                    |
| GC 范围        | 0-100%    | 30–64% 🟡 | 与 Reynolds 30–52% 不同，注意口径 |
| Target range | start-end | —         | 片段全长                      |
| 规则选项         | 三规则全选     | 三规则全选     | —                         |

**输出 top10（每条一行）**：

| # | 片段内位置 | → 全局 | guide 21 nt（5′→3′） | Ui-Tei 徽标 | seed Tm guide / passenger | 脱靶命中 guide +1/+2/+3 | GC | 备注 |
|---|---|---|---|---|---|---|---|---|
| 1 | 46-68 | 285-307 | `UAUUUCCACUUUGUAUAUCCC`（**反向互补推得**，页面目视确认一次） | **U**（=通过） | 20.1 / 8.2 ℃ | 0 / 2 / 52 | ok | passenger 侧 2/6/61 |
| 2 | | | | | | | | |
| … | | | | | | | | |

### 2.1 读表说明（截图逐列对齐，2026-11-18）

siDirect 2.0 输出共 14 列，按用途分四组：

| 组 | 列 | 含义 | 坑 |
|---|---|---|---|
| 定位 | target position | 输入序列上的靶位点范围（1-based，**含 2 nt 突出**，故长 23 nt） | +239 换算全局；长度 23 不是 21 |
| 序列 | target sequence (21 nt target + 2 nt overhang) | 靶 mRNA 上的 23 nt | **这不是 guide**——它和 guide 反向互补 |
| 序列 | RNA oligo sequences：21 nt guide / 21 nt passenger | 实际合成的两条 RNA（含 U，5′→3′） | guide 与 target 前 21 nt 反向互补（T→U），可当场校验转写对不对 |
| 功能 | functional siRNA selection：**Ui-Tei 徽标** | 通过该规则显示徽标（U = Ui-Tei） | 只显示 Ui-Tei 徽标；Amarzguioui/Reynolds 若要打分得自己算（D28）或看页面其它选项 |
| 稳定性 | seed-duplex stability (Tm)：guide / passenger | seed–靶双链 Tm，双阈值默认 <21.5℃ | **两个值**：guide 20.1（接近上限）、passenger 8.2——不对称在数值上可见 |
| 特异性 | minimum number of mismatches against any off-targets（guide / passenger） | 对任一非靶转录本的最小错配数 | 各为 2 → 满足"≥2 错配"要求 |
| 特异性 | **number of off-target hits with indicated mismatches (strand)** | 0/1/2/3 错配下各有多少个脱靶命中 | ⭐ **D29 `offtarget_scan.py` 的现成对照**：guide 3 错配仍有 **52** 个命中，passenger 侧 61 个 |
| 质检 | target position / GC content | ok / 不过 | 布尔标记，别漏读 |

**一个待核的读表细节 🟡**：guide 列 `0(+) 1` = 0 错配命中 1 个。这条 0 错配命中**大概就是靶 mRNA 本身**（on-target 被计入 off-target 统计），需回页面 [detail] 链接确认——**如果不确认就照抄，将来会把"1 个完美脱靶"这件事虚报出去**（D21 虚警的同型风险）。

### 2.2 ⚠️ C12 更正（我的判断错误）

上一轮我把你填在 Ui-Tei 列的 "U" 判为"疑为列名误抄"——**错了，那是页面上的 Ui-Tei 通过徽标**，你的记录方式是对的。同时你也确实抄错了另外一处：**guide 列填的是 target sequence（23 nt）不是 guide（21 nt）**，两者反向互补，用上面的校验规则一查就出来。**这条记进 C 系列：判断别人记录对不对之前，先看原始界面。**

**坐标换算（本日最易错的坑，同 D25 链方向坑同型）**：

- siDirect 输出位置是**输入片段内 1-based** → 全局位置 = 片段内位置 + **239**。
- guide（反义链）对应 transcript 的**同一段反向**——所以"siRNA 位置"和"D21 窗口 [start, end)"对齐时，用窗口区间比对，不逐位比。
- 记录时直接换算成全局坐标填最后一列，后面全用全局坐标。

## 3. 对比作业：siDirect top10 vs D21 全表

**前置检查**：本地 `TTR_sites_ranked.tsv` 还是旧版 **11 行（top10）**——需要先用 v0.1 版重跑出**全表 120 行**（服务器同步后 `python scan_sites.py . variants.tsv siRNA`，v0.1 起写盘全表）。

**三个指标（v0，足够回答"重合度如何"）**：

1. **满分命中数**：siDirect top10 里有几条落在 D21 score=10 的窗口集合？
2. **可及性 top20 交集**：D21 按可及性排序的前 20 个窗口，与 siDirect top10 的窗口重叠几条？
3. **分歧样本**：挑 2 条"规则通过但可及性排名 >50"的窗口 + 2 条"可及性满分但 Reynolds <6"的窗口，逐条对比——**这两类分歧就是第六维度讨论的素材**。

**预登记预期（先写下你的预测再跑！）**：我预测交集偏低（可及性和序列规则测的是不同性质），你的预测＿＿＿。预测-结果差本身就是发现。

## 4. 结果写回（§4 待填）

- 重合度结果与三个指标数字：
- 最大的意外：
- 对第六维度的初步判断（可及性该做过滤、加权、还是独立列）：
- siDirect 实际用的 Ui-Tei iii 口径（≥5 还是 ≥4）：

## 5. 自检 6 题

1. siDirect 的 seed Tm 过滤对应我们知识库里哪一类脱靶？为什么 Tm 低反而好？
2. 片段内位置换算全局坐标的公式？为什么用区间比对而不是逐位比对？
3. siDirect 的 GC 默认范围和 Reynolds 判据 1 的范围差在哪？这算冲突吗？
4. 如果 siDirect top10 与可及性 top20 完全不重叠，你能立刻下结论"可及性没用"吗？为什么？
5. 三套规则在 siDirect 里都有，为什么它还需要第四道 seed Tm 过滤？（提示：三规则管药效，不管什么？）
6. 今天哪一步最容易产生静默错误？你打算怎么防？
7. 页面上 guide 的 `0(+)` 命中数 = 1，这个"完美匹配的脱靶"最可能是什么？如果不核实就写进报告，会造成什么后果？

## 6. 事实边界

| 断言                                | 状态                 |
| --------------------------------- | ------------------ |
| siDirect 2.0 = Ui-Tei/Naito 组官方工具 | ✅ 多源一致             |
| seed Tm 默认 <21.5℃（Ui-Tei 2008）    | ✅ 多源一致             |
| 序列限制 600 bp                       | 🟡 二手（专利文本），官网确认   |
| GC 默认 30–64%                      | 🟡 单源，官网确认         |
| Ui-Tei iii 实际口径                   | 🟡 今天跑完当场记录（§4 待填） |

## 7. 三档预算

- **A（45 min）**：跑 siDirect + 填 §2 表 + 自检 6 题。
- **B（+20 min）**：完成 §3 三个指标 + §4 写回。
- **C（+并入 D28）**：把 §3 的对齐逻辑写成脚本函数（D28 顺带实现），不再手数。

---

## 复测区（自检批改后留题）

（待填）
