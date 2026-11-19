# 代码日 D29（工具日②）：`offtarget_scan.py` v0 落地

> 预算：**3 h**（工具日标准）｜ 日期：周五 11-20 ｜ 服务器目录：`~/project/AOC/scripts/D29/`
> 关联：D22 §5 规格（本文的契约来源）、W4 实施细则 D29 行、D21 `scan_sites.py`（镜像对象）、D28 `siRNA_design.py`（候选来源）。
> 交付判定：**测试 32/32 全绿 + toy 库实跑通过 → 地图 v1 ③ 行 🟡→✅**（按四值语义改）。

---

## 1. 今天解决的一句话

D21 的问题是"**给定转录本，找候选窗口**"（打哪）；今天问的是反方向：

> **给定药物序列（21-mer guide），在小转录本库里找所有还能打到的地方（打歪到哪）。**

D28 的 `siRNA_design.py` 产出了 598 个候选——但那 120 个窗口 / 598 个候选**没有一格问过"这个 guide 还能配到别的转录本吗"**。今天补上这一格，输出直接回流 D30 整合日的 top10 表作 `dim_offtarget`。

## 2. v0 契约（动手前先读完，这就是验收标准）

| 约定 | 内容 | 出处 |
| --- | --- | --- |
| 字母表 | 统一 DNA 大写 ACGT（含 N）；RNA 输入**单点归一化** U→T（S0） | D28 KeyError 教训复用 |
| 方向 | guide 与靶 RNA **反向平行** → 一律搜 `revcomp(query)` | D20/D28 方向坑再现位，A15 断言专门盯 |
| 三档输出 | `exact`（全长完全匹配）/ `contig`（≥min_len 极大连续 block，默认 12 🟡）/ `seed`（guide 2–8 位 revcomp） | D22 §5 |
| 双链 | exact/contig **guide + passenger 都扫**；seed 只算 guide（v0 口径，passenger seed 留 v1） | D22 §2"两条链都要算" |
| on-target 排除 | `scan_fasta` 层按基因名排除 TTR——打中的不叫脱靶 | W4 D29 行 |
| 排序 | n_exact 降序 → max_contig 降序 → gene 名升序 | D21 可复现纪律 |
| 位置 | 全部 0-based（转录本内偏移） | 工程约定 |
| 搜索方式 | **朴素字符串搜索，不下 Bowtie、不建索引**；只对小转录本集 | W4 风险预案 |
| 阈值口径 | min_len 默认 12（D22 🟡 待核）；W4 写的"≥15-mer"作为报告层过滤（`max_contig ≥ 15` 或 `n_exact > 0`），两档都输出 | 两文档口径差异，参数化处理 |

## 3. 函数清单（9 个桩，逐个实现）

| 桩 | 函数 | 一句话 | 断言 |
| --- | --- | --- | --- |
| S0 | `_norm(seq)` | 大写 + U→T，单点归一化 | A01–A02 |
| S1 | `revcomp_dna(seq)` | 反向互补（N 保留） | A03–A04 |
| S2 | `gene_from_header(header)` | FASTA 头 → 基因名（第一个空白 token） | A05–A06 |
| S3 | `find_exact(query, subject)` | 全长匹配的所有起点 | A07–A09 |
| S4 | `find_contig(query, subject, min_len)` | ≥min_len 的**极大** block [(i, j, L)]；L==len(query) 不报（那是 exact 的活）；L 降序 | A10–A13 |
| S5 | `find_seed(guide, subject)` | guide[1:8] 的 revcomp 全部起点 | A14–A15（方向陷阱） |
| S6 | `scan_transcript(guide, sense, header, seq, min_len)` | 单转录本：双链 exact/contig + guide seed | B01–B06 |
| S7 | `_read_fasta` + `scan_fasta(...)` | 全库扫描、排除 on-target、三档汇总排序 | B07–B10 |
| S8 | `scan_candidates_tsv(...)` | 读 D28 候选 TSV 批量扫 → 回流接口 | B11–B17 |

## 4. 执行流程（沿用 D28 节奏）

1. **上传**：`offtarget_scan_skeleton.py` + `test_offtarget_scan.py` → 服务器 `~/project/AOC/scripts/D29/`，conda activate `aoc_dep`。
2. **Part A 逐函数实现**（建议顺序 S0→S5 纯函数，先 `python test_offtarget_scan.py` 跑到 A15 全 PASS）。
3. **Part B 端到端**：toy 转录本库（5 条合成转录本：on-target / 15-mer 块 / 仅 seed / passenger exact / 干净）+ toy 候选 TSV。**32/32 全绿为止**。
4. **实跑（Part C，toy 全绿后）**：
   - 转录本库 `liver20.fa`：TTR（D21 的 ttr.fa 可直接并入）+ 10–20 条肝脏高表达基因的 **MANE Select** mRNA（建议：ALB, APOB, APOE, PCSK9, LDLR, HBB, TFRC, GUSB, HPRT1, ACTB, GAPDH, KRT18, TTR, SERPINA1, F9, CYP3A4……手动 NCBI/Ensembl 网页下载即可，v0 不写下载脚本）。⚠️ 先看一眼 FASTA 头格式——`gene_from_header` 取第一个 token，如果你的头是 `>ENST…|TTR|…` 形式，header 解析要相应调整（这是实跑第一坑）。
   - 候选：D28 的 `~/project/AOC/scripts/D28/TTR_siRNA_candidates.tsv`（确认列名 `guide`/`sense`/`id`，不同就传参改）。
   - 跑法（REPL 或小脚本都行）：
     ```python
     from offtarget_scan_skeleton import scan_candidates_tsv
     rows = scan_candidates_tsv("TTR_siRNA_candidates.tsv", "liver20.fa",
                                on_target="TTR", min_len=12, limit=20)
     ```
     输出写 `TTR_offtarget_summary.tsv`（id, guide, n_off_genes, n_exact, max_contig, n_seed）。
5. **收尾**：地图 v1 ③ 行 🟡→✅（产物列补实跑数字）；`dim_offtarget` 回流 `scan_sites` 留 D30 整合日做。

## 5. 实跑注意事项（三条坑，全是本周踩过的）

1. **header 解析**：NCBI fasta 头常见 `>gene | acc | title` 或 `>lcl|…`，先 `head -1 liver20.fa` 看格式再定 S2 的取法；toy 断言只保证简单格式。
2. **方向**：如果实跑发现"全库 0 命中但 TTR 明明在库里"，第一嫌疑是把 guide 原向搜了（没 revcomp）——A15 就是防这个的。
3. **seed 命中多是正常现象**：7-mer 随机命中期望 ~6%/kb（D22 §7.2 的 1,200 基因账）——**逐条 seed 命中不构成定罪**，今天的输出只是"计数清单"，B 类的富集判定是 RNA-seq 阶段（D22 S4）的事。报告里写计数，不下结论。

## 6. 完成后的地图动作

`03_SOP/计算切入点地图_v1.md` §2 第 ③ 行：状态 🟡→✅，产物列补：
`offtarget_scan.py v0 落地（32 断言，toy 库 + liver20 实跑）——A 类 exact/contig 双链 + B 类 seed 计数，min_len=12 🟡 待核`。

## 7. <details> 自检（5 题，学完先自测）

<details>
<summary><b>自检清单（点开做，做完打勾）</b></summary>

1. 为什么搜索方向一律用 `revcomp(query)` 而不是直接搜 guide？如果直接搜会漏掉什么还是错报什么？
2. `find_contig` 为什么要求"极大 block"且跳过 L==len(query)？不这样的话输出会怎样重复？
3. exact/contig 为什么双链都扫、seed 为什么只算 guide？（提示：passenger 装载概率 vs seed 机制的归属）
4. toy 库里 OFF2 只有 seed 命中（7-mer）——能不能据此说"OFF2 会被 siRNA 压低"？为什么？这和 D22 §7.2 的 1,200 基因账是什么关系？
5. `scan_fasta` 排除 on-target 用基因名精确匹配——什么时候会排除失败？失败后果是什么（多报还是漏报）？

</details>

---

## 附：本文事实边界

| 陈述 | 状态 |
| --- | --- |
| gapmer/RNase H 切割所需连续互补长度（≥12 一说） | 🟡 D22 已挂待核，v0 参数化不写死 |
| seed = guide 1-based 2–8（7 nt） | ✅ miRNA/siRNA 通行结论（D22/D25 已核） |
| 双链都要算脱靶 | ✅ D22 §2 通行建议 |
| min_len=12、报告层 15 档、双链/单链 seed 口径 | 工程约定，非事实陈述 |
| liver20 基因清单 | 建议清单，实跑时按兴趣调整 |
