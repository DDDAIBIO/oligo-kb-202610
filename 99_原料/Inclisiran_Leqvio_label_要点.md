# Inclisiran（Leqvio）FDA 说明书要点摘录

| 项 | 内容 |
|---|---|
| 来源 | openFDA drug/label API（`search=openfda.brand_name:"LEQVIO"`），原始 JSON = `inclisiran_label.json` |
| 抓取日期 | 2026-10-18（D3） |
| 标签 | `#siRNA #GalNAc #ESC #PCSK9 #Inclisiran` |

---

## 1 · Description（序列与化学）

> LEQVIO contains **inclisiran sodium**, a small interfering RNA (siRNA) directed to **PCSK9** mRNA. Inclisiran contains a covalently linked ligand containing **three N-acetylgalactosamine (GalNAc)** residues to facilitate delivery to hepatocytes. With one exception, the 2′ ribose moieties are present as **2′-F or 2′-OMe**. In addition, **six of the terminal phosphodiester backbones** are present as phosphorothioate linkages.

**分子式**：`C529 H664 F12 N176 Na43 O316 P43 S6`
**分子量**：**17,284.72 g/mol**（钠盐形式）

### 完整序列（5′ → 3′）

```
Sense (passenger) 21 nt，3′ 端接 L96：
Cms-Ums-Am-Gm-Am-Cm-Cf-Um-Gf-Um-dT-Um-Um-Gm-Cm-Um-Um-Um-Um-Gm-Um — L96

Antisense (guide) 23 nt：
Ams-Ams-Gm-Am-Um-Cf-Um-Gf-Gm-Af-Cm-Af-Am-Af-Am-Cf-Gm-Af-Af-Af-Ams-Cfs-Am
```

缩写：`f` = 2′-F；`m` = 2′-OMe；`s` = 硫代磷酸（PS）；`dT` = 2′-deoxy thymidine；`L96` = 三触角 GalNAc

> ⚠️ **关键**：GalNAc（L96）连在 **sense（乘客）链的 3′ 端**，**不是**反义链。
> 出处三处互证：FDA label 原文 "conjugated on the sense strand"；MDPI *Pharmaceuticals* 2022, 15(2):222 "attached ... through its 3′ terminus"；JAHA 2024 "at the 3′ end of the passenger strand"。

---

## 2 · Mechanism of Action

> Inclisiran is a double-stranded siRNA, **conjugated on the sense strand with triantennary GalNAc** to facilitate uptake by hepatocytes. In hepatocytes, inclisiran utilizes the RNA interference mechanism and **directs catalytic breakdown of mRNA for PCSK9**. This **increases LDL-C receptor recycling** and expression on the hepatocyte cell surface, which increases LDL-C uptake and lowers LDL-C levels in the circulation.

逻辑链：siRNA 降解 PCSK9 mRNA → PCSK9 蛋白↓ → LDL 受体不再被降解、回收增加 → 肝细胞表面 LDLR↑ → 血中 LDL-C↓

---

## 3 · Dosage and Administration

- **284 mg** 皮下注射（腹部 / 上臂 / 大腿），由医护人员执行
- **第 1 天 → 第 90 天（3 个月）→ 之后每 6 个月**
- 漏服 < 3 个月：补打并维持原方案；漏服 > 3 个月：重新走"首针 + 3 个月"流程
- LDL-C 可在启动后 **30 天**起随时测，不需卡给药时间

---

## 4 · Pharmacokinetics（**最反直觉的一段**）

| 参数 | 数值 |
|---|---|
| Tmax | 约 **4 小时** |
| Cmax（284 mg） | 509 ng/mL |
| **血浆中检测不到的时间** | 给药后 **24–48 小时** |
| **终末消除半衰期** | **约 9 小时** |
| 血浆蛋白结合率 | 87% |
| 表观分布容积 | 约 500 L |
| 代谢 | 主要被**核酸酶**切成更短的核苷酸；**不是 CYP450 底物**、不是转运体底物 |
| 排泄 | 约 **16%** 经肾清除 |
| 蓄积 | 多次给药**无蓄积** |

> **核心矛盾**：血浆半衰期只有 **9 小时**，24–48 小时就测不到了 —— 可给药间隔是 **6 个月**。
> 结论：**药效持续时间完全来自组织（肝细胞内 RISC 库）的持久性，与血浆 PK 无关。**
> 这是 siRNA 药物 PK/PD 脱耦的教科书案例。

## 5 · Pharmacodynamics

- 单次 284 mg 后 **14 天内**即见 LDL-C 下降
- 30–180 天：LDL-C 平均下降 **38%–51%**
- **Day 180 仍下降约 53%**
- PCSK9 蛋白：Day 120 降 **75%**，Day 180 降 **69%**

---

## 6 · Adverse Reactions

常见（≥3%，3 项安慰剂对照试验、1833 例成人）：**注射部位反应、关节痛（arthralgia）、支气管炎**

> 对比 Nusinersen 的血小板减少 + 肾毒性 —— 模态不同、化学不同、毒性谱完全不同。

---

## 7 · 自己算一遍（D3 的验算练习）

1. **F12 → 12 个 2′-F**：分子式里的 F 原子数直接等于 2′-F 修饰的个数。文献给出的修饰统计是 1×dT + 12×2′-F + 31×2′-OMe = **44** 个核苷 = 21 + 23 ✓ 完全自洽
2. **S6 → 6 个 PS**：sense 5′ 端 2 个 + antisense 5′ 端 2 个 + antisense 3′ 端 2 个 = 6 ✓
3. **Na43 = 21 + 22**：label 结构图分别标注 "Sense Strand 21 Na+" 与 "Antisense Strand 22 Na+"（那是钠离子数，不是核苷酸数 —— 别看错）
4. **P43 = 42 + 1**：骨架磷酸二酯键 (21−1) + (23−1) = 42，再加 L96 连接臂上的 1 个
5. **284 mg 是多少个分子？** 0.284 g ÷ 17284.7 g/mol = 1.64×10⁻⁵ mol → ×6.022×10²³ ≈ **9.9×10¹⁸ 个分子**（约 10¹⁹，记住这个数量级）

---

## 待核验

- [ ] ESC vs ESC+ 的具体差别（Alnylam 平台文献；"ESC+" 特指种子区钝化那一步）
- [ ] 为什么 antisense 是 23 nt 而 sense 是 21 nt（不对称设计的理由）
- [ ] dT 出现在 sense 链第 11 位的作用（降低免疫原性？）
