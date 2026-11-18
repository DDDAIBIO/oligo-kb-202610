# D28 ｜ 工具日①（3h）：siRNA_design.py 前半——候选生成 + 三规则打分

> 日期：2026-11-19（W4-D28）｜ 预算 **3 h**（D21 实测值）｜ 前置：D25/26/27  
> **第一次吃自己代码的红利**：`scan_sites` 的 `_read_fasta` / `make_windows` 直接复用，不重写。今天的产出是 D25–D26 三天纸面规则的**代码化**，也是 M1 预演的最大素材。

---

## 0. 任务分解（Part A / Part B，与 D21 同构）

- **Part A（约 2 h）**：实现骨架 11 个 stub，`test_siRNA_design.py` 断言全绿（含 ⭐Patisiran 方向断言）。
- **Part B（约 1 h）**：服务器真跑 `ttr.fa` → `TTR_siRNA_candidates.tsv` 落库 + 自查。

## 1. 复用清单（禁止重写）

```python
from scan_sites_skeleton import ss   # 或 import scan_sites_skeleton as ss
ss._read_fasta(path)      # FASTA 读取（含描述行处理）
ss.make_windows(seq, 19, 1)   # 窗口枚举：设计用 L=19、step=1（全枚举，不像 D21 粗筛用 step=5）
```

## 2. 函数规格（S0–S10，骨架同序）

| #   | 函数                                  | 规格                                                                                                      | 对应知识              |
| --- | ----------------------------------- | ------------------------------------------------------------------------------------------------------- | ----------------- |
| S0  | `_revcomp_rna(seq)`                 | 反向互补，T→U；空串返回空串                                                                                         | 链方向的根基            |
| S1  | `gc_frac(seq)`                      | GC 占比 0–1；空串 0.0                                                                                        | Reynolds 判据 1     |
| S2  | `uitei_pass(guide)`                 | 四条件 AND（**guide 坐标系**）：pos1∈{A,U}；pos19∈{G,C}；pos1–7 A/U **≥4**（siDirect 口径，D26 §4 漂移声明）；无 ≥10 连续 GC    | Ui-Tei = 过滤器      |
| S3  | `_c3_low_internal_stability(sense)` | **代理实现（🟡 待回原文）**：中心 7-mer（pos7–13）A/U ≥4 → 1 分                                                         | Reynolds 判据 3     |
| S4  | `reynolds_score(sense)`             | 八条累加 0–10（判据 3 用 S3 代理）；**不设下限 0**，可负（all-G 例 = −2）                                                     | Reynolds = 排序器    |
| S5  | `_amarz_asymmetry(guide)`           | **代理实现（🟡）**：guide 两端各 3 nt 的 A/U 差分 >0 → 1                                                             | Amarzguioui 条件 1  |
| S6  | `amarzguioui_score(guide)`          | 六条 ± 体系，cutoff ≥3                                                                                       | 分列展示              |
| S7  | `design_candidates(seq, step=1)`    | 全枚举 19-mer → 每窗产出 `{start, end, sense, guide, uitei, reynolds, amarz}`；**sense = 转录本窗口本身，guide = 反向互补** | 候选生成              |
| S8  | `rank_design(fasta_path, top=20)`   | `_read_fasta` → 逐转录本 S7 → 排序键 `(-reynolds, -amarz, start)`                                              | **不合成总分**（D26 拍板） |
| S9  | `format_table(name, rows, top)`     | 表头：`rank start end sense guide uitei reynolds amarz`；全表写盘、屏幕只打 top                                      | D21 教训：写盘必须全表     |
| S10 | `main(argv)`                        | `python siRNA_design.py <fasta> [top]`；输出 `TTR_siRNA_candidates.tsv`；错误进 stderr、退出码 0/1/2               | D21 约定沿用          |

**两条硬约定**：① guide 19 nt 核心，**dTdT 突出在合成时加，脚本不管**（表头注明）；② 三个分数**三列并存**，禁止任何合并。

## 3. ⭐ 方向断言（今天最重要的测试，来源已双源核验）

Patisiran（Onpattro，首个上市 siRNA，TTR 3′UTR）：

```text
sense     : GUAACCAAGAGUAUUCCAU + dTdT   （= TTR mRNA 靶位点，19 nt）
guide     : AUGGAAUACUCUUGGUUAC + dTdT   （= sense 的反向互补）
来源：PMC9330994（药化综述）与药学学报 2023 表格一致 ✅
```

**断言 T8**：把 `GUAACCAAGAGUAUUCCAU` 埋进测试转录本 → `design_candidates` 必须产出 sense = 该 19-mer、guide = `AUGGAAUACUCUUGGUUAC`。  
**如果链方向写反**（guide = 窗口本身），T8 立刻红——这就是 D20 SNP 链归属坑的疫苗。

**顺手的人肉验证**（写测试前先算，别让测试替你发现 bug）：

> ⚠️ **预核验抓出的手算坑（写进 C 系列的教训）**：判据 7/8 是**纯罚分项**——有 G/C19 / G13 才 −1，**没有不加分**。"没有惩罚"≠"+1 奖分"，这是计分表里最容易混的语义。

Patisiran sense 逐条：判据1 +1（GC 36.8%）／判据2 +3（s[14:19]=UCCAU → 3 个 A/U）／判据3 代理 +1（中心 7-mer AAGAGUA = 5 A/U）／判据4 0（s[18]=U）／判据5 +1（A3）／判据6 0（s[9]=A）／判据7、8 无罚 → **合计 6 分，正好压线 ≥6**。真实上市药恰好落在 cutoff 上。

**另一个预核验发现（不修数据，只记录）**：patisiran guide 在 Amarzguioui 六条 ± 体系下只得 **2 分，低于 cutoff 3**——**三套规则对一个真实上市药都会给出"过/不过"相反的结论**。这就是 D26"分列报告"拍板的直接实证：任何"三合一总分"都会把这个分歧抹掉。

## 4. 与 D26 决策的衔接（写进 S8 注释）

- Ui-Tei = 布尔**过滤列**（AND 逻辑天然过滤）
- Reynolds = **主排序键**（连续分天然排序）
- Amarzguioui = **分列展示**（± 体系，不当排序键第一优先）
- 表尾打印一句：`三套规则不合成总分（D26 §5）`

## 5. Part B 步骤（服务器）

```bash
cd ~/project/AOC/scripts/D28        # 新目录
# 同步骨架/测试/依赖 scan_sites_skeleton.py 后：
python test_siRNA_design.py          # 应全绿
python siRNA_design.py ../D21/ttr.fa
# 预期：屏幕 top20、TTR_siRNA_candidates.tsv 全表 ~598 行（616-19+1）
```

**自查三问**（B 段跑完当场答）：top20 里 Ui-Tei 通过几条？Reynolds 最高几分？patisiran 靶位点（若在 ttr.fa 内）排第几？

## 6. review 主题预告（Part A 绿了之后喊我）

1. S3/S5 两个代理实现是否需要在注释里挂 🟡（我在骨架里已挂，看你实现时会不会丢）
2. 排序键稳定性（同分 start 升序）
3. 空串/越界等静默边界（S0/S1）

## 7. 自检 5 题（Part A 全绿后答）

1. 为什么 guide 列存 19 nt 而表头注明 dTdT？把 overhang 存进序列会有什么后果？不知道
2. S3 的代理实现和原文的差别是什么？为什么现在允许代理、D28 之后不允许？ 不知道
3. Reynolds 打分对 all-G 的 19-mer 会给几分？为什么允许负分比 clamp 到 0 更诚实？ 不知道
4. 如果把排序键改成 `(uitei, -reynolds, ...)`（先过滤再排），输出会少什么信息？不知道
5. T8 断言防的是哪一类错误？它防不了哪一类？（提示：防链方向，不防打分错误）不知道

## 8. 三档预算

- **A（1.5 h）**：Part A 前半（S0–S6 纯函数 + T0–T7 绿）
- **B（3 h，今日标准）**：Part A 全绿 + Part B 跑通 + 三问

## 9. Part B 结果与发现（11-18 实跑，TTR 全表 598 候选）

**过程记录**：Part A 首跑 22/23（T2d 边界差一：`mx<=10` → 改 `mx<10`，(iv) 原文口径"无 ≥10 连续 GC"）→ 23/23。Part B 首跑 **KeyError: 'T'**——S0 纯 RNA 查表正确拒收 DNA 字母（ttr.fa 以 DNA 字母存 RNA，FASTA 惯例）；修复 = `rank_design` 进口处 `seq.upper().replace("T","U")` 单点归一化，S0 防御强度不变。重跑成功：候选 598（=⌊(616−19)/1⌋+1 ✓）、TSV 599 行、三数自洽。

**五条发现**（本地参考实现交叉分析）：

1. **⭐ Patisiran 真实位点（start 517）= 全表第 68/598 名**（uitei ✓ / reynolds 6 压线 / amarz 2 不过 cutoff）。真实上市药在自己转录本里排不进前 11%——"规则是先验不是判决"（D25）的定量实证；药选它还有可及性、化学、递送各层的功劳。
2. **top1（521–540，reynolds 9）是 patisiran 位点右移 4 nt 的邻居窗**。规则认为移 4 个更好，真实药物选了 6 分位——谁对要靠可及性/保守性/脱靶数据仲裁，即第六维度存在的意义。
3. **三规则收敛率 ~6.2%**：uitei 过 88（14.7%）、reynolds≥6 共 89、amarz≥3 共 143，三者交集 37。规则间不一致是常态——D26"分列报告"的定量版。
4. **D27 siDirect top1 闭环**：其 antisense 前 19 nt = `UAUUUCCACUUUGUAUAUC` = 我们表 start 286 行的 guide（逐字吻合，坐标口径核对通过）。注意 start 286 不在我们的 top20（start 285 是 rank 14）——siDirect 叠加 seed Tm 21.5℃ 过滤，规则细节差异实例。
5. **top20 与 D21 可及性岛（260–331）只交 4/20**：规则 top 群集在 500–612（A/U 富的 3′ 尾区），可及性岛在 260–331——**序列规则与可及性各自选出不同区域**，这是第六维度整合（D22 §5 规格回流）最真实的一手素材。待办：用服务器 `TTR_sites_ranked.tsv` 查 top20 各窗 mean_access。

- **C（3 h + 顺延）**：把 D27 §3 的对齐逻辑并入（对比作业代码化，可滚入 D29）
