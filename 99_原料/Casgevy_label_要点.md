# Casgevy（exagamglogene autotemcel）FDA label 要点

> 来源：openFDA `drug/label.json?search=openfda.brand_name:"CASGEVY"`，原始 JSON 存于本目录 `casgevy_label.json`
> 整理日期：2026-10-19（D5）
> 状态：有效

## 基本信息

- 通用名：**exagamglogene autotemcel**（简称 exa-cel）
- 商品名：CASGEVY
- 厂家：**Vertex Pharmaceuticals**
- 类别：**细胞基因治疗**（autologous genome edited HSC-based gene therapy）

## DESCRIPTION（原文要点）

> "CASGEVY is a cellular gene therapy consisting of autologous CD34+ HSCs edited by CRISPR/Cas9-technology at the **erythroid specific enhancer region of the BCL11A gene** to reduce BCL11A expression in erythroid lineage cells, leading to increased fetal hemoglobin (HbF) protein production."

- 制备流程：患者自身 HSC（单采获取）→ 富集 CD34+ → **体外（ex vivo）电转导入 CRISPR/Cas9 核糖核蛋白（RNP）复合物**
- sgRNA 引导 Cas9 在 **BCL11A 基因红系特异性增强子区的关键转录因子结合位点（GATA1）** 产生精确 DNA 双链断裂
- 结果：GATA1 结合被破坏 → BCL11A 表达降低 → γ-globin 表达升高 → HbF 形成增加
- 制剂：冻存悬液，含 **5% DMSO + dextran 40**

## INDICATIONS AND USAGE

患者年龄 **≥ 2 岁**：

1. 伴复发性血管阻塞危机（VOC）的**镰状细胞病（SCD）**
2. **输血依赖型 β-地中海贫血（TDT）**

## MECHANISM OF ACTION（原文要点）

回输后编辑过的 CD34+ 细胞在骨髓植入，分化为红系细胞且 BCL11A 表达降低 → γ-globin 与 HbF 增加。

- **SCD**：HbF 降低细胞内 HbS 浓度 → 阻止红细胞镰变 → 消除 VOC
- **TDT**：γ-globin 改善 α/非 α 珠蛋白失衡 → 减少无效造血与溶血 → 提高总血红蛋白 → 摆脱规律输血

## DOSAGE AND ADMINISTRATION

- 仅供自体使用、仅供静脉使用、单次一次性给药
- 最低推荐剂量：**3 × 10⁶ CD34+ cells/kg**（按体重）
- 给药前 **48 小时至 7 天内必须完成清髓预处理（myeloablative conditioning）**
  - 需考虑癫痫预防、肝静脉阻塞病（VOD）预防
- 融化后 20 分钟内完成输注；**禁止使用在线血液过滤器**；禁止取样/改动/照射

## 两个值得琢磨的点

1. **编辑的是增强子，不是基因** —— 见笔记正文。BCL11A 在 B 细胞、神经系统等其它谱系中必需，全身敲除不可接受；只破坏**红系特异性**增强子，用一次编辑换来**谱系选择性**。
2. **"基因治疗"的门槛在流程不在分子** —— 需要动员 + 单采 + 清髓化疗 + 回输植入，整套流程的毒性和成本远高于分子设计本身。

## 待补

- [ ] 临床疗效数据（CLINICAL STUDIES 章节，本次未提取）
- [ ] 脱靶评估方法与 CMC 要求（FDA 对体内/体外编辑产品的非临床要求）
