# Eteplirsen / EXONDYS 51 —— label 要点摘抄

> 来源：openFDA `drug/label`，version 15，effective_time **2025-08-13**
> 原文 json：`99_原料/eteplirsen_label.json`
> 拆解日期：2026-11-04（W2 · D14 肝外递送）
> 摘抄原则：只抄 label 原文写明的，label 没写的一律标「⚠️ 待核」不补

---

## 1. 一句话定位

**首个获批的 PMO（磷酸二酰胺吗啉寡聚物）**，2016 年经**加速批准**上市，靶 DMD 基因 **exon 51**，用于可行外显子 51 跳跃的 DMD 患者。

它在 D14 里的作用是**反例坐标**：证明"肌肉递送"要靠 **30 mg/kg/周** 的剂量硬堆，与肝脏靶向的 25 mg/3 月不在一个量级。

---

## 2. 化学（label §11 原文）

| 项 | 内容 |
|---|---|
| 骨架类型 | **PMO**：核糖环换成**六元吗啉环**；连接键是**不带电的磷酸二酰胺**（不是天然带负电的磷酸） |
| 长度 | **30 个亚基**（"contains 30 linked subunits"） |
| 分子式 | `C364 H569 N177 O122 P30` |
| MW | **10,305.7 Da** |
| 修饰 | **无 PS、无 2′-OMe、无 2′-F** —— 中性骨架，全靠骨架本身抗核酸酶 |

> ⭐ 与 ASO/siRNA 的根本差别：**PMO 不带电**。既不能靠静电吸附阳离子载体，也不会被"多阴离子"这个问题困扰；代价是**无法用 GalNAc/LNP 这类依赖电荷或受体的策略**，只能靠剂量。

---

## 3. 制剂（label §11）

- 50 mg/mL；单剂量瓶 100 mg / 500 mg
- 等渗磷酸盐缓冲液，pH **7.5**，渗透压 260–320 mOsm
- 辅料：KCl 0.2 mg/mL、磷酸二氢钾 0.2 mg/mL、NaCl 8 mg/mL、磷酸氢二钠 1.14 mg/mL

---

## 4. 给药方案（label §2.1）

- **30 mg/kg，每周一次（Q1W），静脉输注 35–60 分钟**
- 需稀释：用 0.9% NaCl 稀释到 **100–150 mL** 总量，经 **0.2 µm 在线滤器**给药
- 同一静脉通路不得与其他药物同时输注；出现超敏反应时考虑减慢或中断

### 年负担（70 kg 患者）

| 项 | 计算 | 结果 |
|---|---|---|
| 单次剂量 | 30 × 70 | **2,100 mg = 2.1 g** |
| 年给药次数 | 52 | 52 次 |
| **年总剂量** | 2.1 × 52 | **≈ 109 g** |
| 年输注时间 | 52 × 约 47.5 min | **≈ 41 小时** |

> 对照：年剂量是 **Vutrisiran 的 1092 倍**（100 mg）、**Patisiran 的 300 倍**（364 mg）；年输注时间也是 Patisiran 的 1.8 倍（41 h vs 23 h）。

---

## 5. 机制（label §12.1）

> "Eteplirsen is designed to bind to **exon 51 of dystrophin pre-mRNA**, resulting in exclusion of this exon during mRNA processing... allow for production of an **internally truncated dystrophin** protein."

- 剪接**位阻型**（与 Nusinersen 同类），不是 RNase H 降解型
- 靶点是 **pre-mRNA 的内含子/外显子边界**，作用在**剪接之前**
- 产物是"截短但保留部分功能"的 dystrophin —— 不是恢复全长蛋白

---

## 6. PD：这本 label 最值得记住的一张数（§12.2）

| 试验 | 时间 | 肌肉中 dystrophin 水平（占正常人） |
|---|---|---|
| Study 2 | 治疗 **180 周**（3.5 年） | 平均 **0.93%** |
| Study 3 | 治疗前 | **0.16%** |
| Study 3 | 治疗 48 周后 | **0.44%** |
| Study 3 | 中位增幅 | **0.1%** |

> ⚠️ 这组数字必须带上下文引用：dystrophin 恢复受限**不只因为递送**，还因为 exon skipping 本身的效率、以及"截短蛋白"只能是部分功能。但**109 g/年换来不到 1% 的蛋白恢复**，方向性结论没有争议：**肌肉递送是核酸药里最贵的一段路**。

---

## 7. PK（label §12.3）

- **不被肝微粒体代谢**（各物种、含人均无代谢）
- **肾清除约占给药剂量的 2/3，在 24 小时内完成**
- 推断（标 ⚠️ 为我的推理，非 label 原文）：**中性骨架 → 血浆蛋白结合低 → 肾小球自由滤过 → 快速清除**，这正是"必须每周给药 + 30 mg/kg"的直接原因

---

## 8. 监管状态：加速批准（§1 原文）

> "...approved under **accelerated approval** based on an increase in dystrophin in skeletal muscle... Continued approval... may be contingent upon verification of a clinical benefit in confirmatory trials."

⭐ 与 D8 的 **Pelacarsen 教训**连起来看：**替代终点（surrogate endpoint）≠ 临床终点**。
区别在于：Pelacarsen 是 III 期**没达到**主要终点而失败；Eteplirsen 是**以替代终点换取了上市**，临床获益仍在确证中。两个方向、同一课。

---

## 9. 待核验清单

- [ ] 未找到 `WARNING:` 黑框警示段（本 label 无 boxed warning）—— 需确认是否为 openFDA 提取遗漏
- [ ] 血浆半衰期（label 未在本节给出明确数值）
- [ ] 加速批准的确证性试验（MISSION）最终结果与当前状态
- [ ] 与其他 PMO（golodirsen、casimersen、viltolarsen）的剂量/暴露差异
