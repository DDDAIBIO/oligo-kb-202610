# D12 学习包：GalNAc——一个受体撑起一个赛道

> **日期**：2026-11-03（周二，可顺延）｜ **时长**：1 h ｜ **前置**：D11
> **材料**（不用找资料）：
> - `99_原料/Vutrisiran_Amvuttra_label_要点.md`（主材料，读全文）
> - `99_原料/vutrisiran_label.json`（原文，练习时自己定位句子）
> - `03_SOP/化学修饰速查表_v0.md`（偶联层那一段）
>
> **本页是什么**：练习骨架。**所有 ✏️ 处先自己动手，再展开 `<details>` 对答案。** 对不上就以 label 原文为准，并写进待核清单。

---

## §1 一句话目标

能不看资料说出：**GalNAc 为什么只能打肝细胞、为什么连三份、以及为什么 2025 年三款新药全都复用这个底盘。**

---

## §2 机制骨架：GalNAc → ASGPR（先读要点 §1–§2，再补全这里）

 GalNAc（三触角 N-乙酰半乳糖胺）→ 被肝细胞表面 **ASGPR**（去唾液酸糖蛋白受体）识别 → **受体介导内吞** → 内体酸化后配体解离 → siRNA 进入 RISC。

三个记忆锚点（每个用一句话核对 label 原文）：

1. **组织选择性来自受体分布**：ASGPR 几乎只在肝细胞表达 → 所以 GalNAc 偶联 = **肝靶向**，不是"全身给药后选择性进肝"。
2. **三价不是装饰**：label 原文写的是 *"a ligand containing **three** N-acetylgalactosamine (GalNAc) residues"*。
   ✏️ 多价（三价）带来什么效应？为什么单价不够？

   <details><summary>参考（对完再合上）</summary>

   **多价效应 / avidity**：一个配体同时占用受体上多个结合位点，结合强度远超单价之和（而非简单相加）。单个 GalNAc 对 ASGPR 亲和力只有微摩尔级，三价排列可达纳摩尔级 —— 这是 GalNAc-siRNA 能以 25 mg 级剂量工作的前提。注意区分：**affinity（单个位点的强度） vs avidity（多价协同的总强度）**。

   </details>

3. **代价**：失去 LNP 包封后，保护 RNA 的工作要由化学修饰自己扛 → 对比 Patisiran（PS=0）与 Vutrisiran（S=6），这就是 **C5 更正**（修饰与递送功能重叠、可互相替代）的直接证据。

---

## §3 核心练习（20 min）：自己数一遍分子式

> label 原文（DESCRIPTION 段，已在 json 中核实）：
> *"The molecular formula of vutrisiran sodium is **C530 H672 F9 N171 Na43 O323 P43 S6** with a molecular weight of **17,290 Da**. The molecular formula of the free acid is **C530 H715 F9 N171 O323 P43 S6** with a molecular weight of **16,345 Da**."*

逐项反推（先填，再展开）：

| 元素 | 数值 | ✏️ 你的推断（一句话） | 来源依据（D7 套路） |
|---|---|---|---|
| P | 43 | | 磷原子 ≈ 核苷酸数 |
| S | 6 | | 硫 = PS 骨架数 |
| F | 9 | | 氟 = 2′-F 数 |
| Na | 43 | | ⚠️ 钠不是核苷酸 |

<details><summary>参考答案（先自己写完）</summary>

- **P=43** → 共 43 个磷酸基 ≈ 43 nt（每核苷酸骨架 1 个 P；GalNAc 连接子本身不含 P）。
- **S=6** → 6 个 PS（硫代磷酸）连接。对比：Inclisiran 也是 6 个 PS —— GalNAc-siRNA 的典型"够用就行"取值；对照 Nusinersen 全长 PS 17/17。
- **F=9** → 9 个 2′-F 修饰（分布在 antisense 链中部，与 ESC+ 设计一致；具体位置 label 结构图里数）。
- **Na=43** → 43 个钠离子 = 43 个酸性基团（每核苷酸一个磷酸钠盐），**不是核苷酸数**（D7 最大的坑）。
- **公式自洽检验（必做）**：游离酸 H715 − 钠盐 H672 = **43 = Na 数**。每个钠盐位点失去 1 个 H 换成 Na（H 与 Na 分子量差解释 17,290 − 16,345 ≈ 945 ≈ 43×22）。✅

</details>

### ⚠️ 解不开的点（这正是真实工作的样子）

✏️ 用 P=43 反推链长：43 nt 总量。但常见的二手说法（Inclisiran 的 sense 21 + antisense 23 = 44 nt）对不上；`99_原料/Vutrisiran_Amvuttra_label_要点.md` §3 已记了这个疑问。

**怎么处理**：
1. 回 json 原文找 sense/antisense 或 sequence 相关段落，看 label 是否直接写了链长；
2. 想一个能造成 ±1 差异的化学原因（提示：连接子占位、末端核苷酸缺 5′-P、凸起 bulge……都可能是答案，别硬凑）；
3. 把结论/未解写进 `NOW.md` 待核清单的那一行。

<details><summary>兜底提示（实在卡住 10 分钟再看）</summary>

Label 通常**不给碱基序列表**，只给结构式图。目前确认的事实只有两条：P=43、分子量 17,290 Da。链长问题**允许保持未解**——它在 `NOW.md` 待核清单里挂着（"确认 Vutrisiran 链长"）。本练习的验收点是"知道哪里解不开 + 不硬凑"，不是凑出答案。

</details>

---

## §4 定位练习（10 min）：在原文里找到这三句

用 `99_原料/vutrisiran_label.json`（或任意 JSON 查看器）逐一定位，把 section 名填在括号里：

1. "three N-acetylgalactosamine (GalNAc) residues" →（　　）
2. "25 mg administered by subcutaneous injection once every 3 months" →（　　）
3. 出现 **DLin-MC3-DMA 或 LNP 组分**的句子找得到吗？✏️ 找不到说明了什么？

<details><summary>参考</summary>

1. DESCRIPTION（§11，化学描述段）；2. DOSAGE AND ADMINISTRATION（§2.1）；3. **找不到** —— Vutrisiran 不是 LNP 制剂，说明这份 label 里没有 LNP 组分，与 Patisiran 形成对照（Patisiran 的 DESCRIPTION 里有完整四组分含量）。这是"制剂即成分"最直接的读法训练。

</details>

---

## §5 大图：为什么 2025 三款新药都复用这个底盘

| 药 | 模态 | 靶/适应症 | 偶联 | 共同点 |
|---|---|---|---|---|
| Fitusiran | siRNA | 抗凝血酶（血友病） | GalNAc-肝 | 底盘复用 |
| Donidalorsen | ASO gapmer | 5-10-5（遗传性血管性水肿） | GalNAc-肝 | 底盘复用 |
| Plozasiran | siRNA（序列含 I） | APOC3（家族性乳糜微粒血症） | GalNAc-肝 | 底盘复用 |

✏️ 一句话回答验收问题：**为什么大家都用同一个底盘？**

<details><summary>参考（对完合上）</summary>

因为**底盘（GalNAc 偶联 + 修饰策略 + 制剂）解决的是"分布团"，与序列（"药效团"）解耦**（Khvorova 2017 题眼）。换适应症 = 换序列，不必重做递送开发 → 研发边际成本骤降。这正是"信息药物"的工业意义：**递送平台化、药效定制化**。反例是肝外递送（AOC）——底盘还没平台化，所以每家都要自己趟。

</details>

---

## §6 收尾动作（15 min）

- [ ] 术语表补 1 条：**多价效应（avidity）**——ASGPR 那条已有，别重复添加；补完总条数写到 **29→30**（用实际行数核对，不心算）
- [ ] `NOW.md` 待核清单：Vutrisiran 链长那行，把你 §3 的处理结论补一句
- [ ] D11 欠账检查：`03_SOP/化学修饰速查表_v0.md` 里「ΔTm 量级 + 代价」两列若还没加，**今天只加"代价"列即可**（10 min 保底），ΔTm 留给阶段二热力学计算时一起做
- [ ] commit：`D12: GalNAc拆解 + Vutrisiran分子式反推练习 + 术语表30条`
- [ ] 在 VSCode 点 Sync 推远端

### 验收三问（合上所有资料）

1. GalNAc 为什么**只能**打肝细胞？
2. 为什么是**三价**？单价差在哪？
3. 2025 三款新药为什么共用一个底盘，这个模式的商业逻辑是什么？

> 3 问都能答 → D12 通过。答不出第 2 问 → 回 §2 锚点 2；答不出第 3 问 → 回 §5。

---

## 本节新增术语速记

| 术语 | 一句话 |
|---|---|
| ASGPR | 去唾液酸糖蛋白受体，几乎只表达于肝细胞，识别末端 GalNAc → 肝靶向的生物学基础 |
| 多价效应（avidity） | 多个配体同时结合多个受体位点，总亲和力远超单价之和；三价 GalNAc = 纳摩尔级亲和力 |

> 标签：#GalNAc #ASGPR #siRNA #分子式反推 #W2
