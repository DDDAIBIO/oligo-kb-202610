# Spikevax（mRNA-1273）说明书要点

| 项 | 内容 |
|---|---|
| 来源 | DailyMed SPL setid `f96b315c-fa57-4876-a7e5-a9b584d8e6e6`，版本 11，发布 **2026-08-31**（2026-2027 Formula） |
| 原始文件 | `99_原料/spikevax_spl.xml` |
| 抓取日期 | 2026-10-19（D4） |
| 备注 | openFDA drug/label 库**查不到**疫苗（归 CBER 管，不在 CDER label 库），改用 DailyMed |

---

## 11 DESCRIPTION（原文摘录 + 译注）

> SPIKEVAX (COVID-19 Vaccine, mRNA) is a sterile white to off-white injectable suspension for intramuscular use. Each 0.5 mL dose of SPIKEVAX (2026-2027 Formula) contains **50 mcg nucleoside-modified messenger RNA (mRNA)** encoding the **pre-fusion stabilized spike glycoprotein (S)** of the SARS-CoV-2 **JN.1-descendent variant XFG**.

**每 0.5 mL（12 岁及以上）剂量组成：**

| 成分 | 量 | 角色 |
|---|---|---|
| nucleoside-modified mRNA | **50 mcg** | 活性成分（编码 pre-fusion 稳定的 S 蛋白） |
| 总脂质 | **1.01 mg** | LNP：SM-102 + PEG2000-DMG + cholesterol + DSPC |
| tromethamine | 0.25 mg | 缓冲（Tris） |
| tromethamine hydrochloride | 1.2 mg | 缓冲 |
| acetic acid | 0.021 mg | pH 调节 |
| sodium acetate trihydrate | 0.10 mg | pH 调节 |
| sucrose | **43.5 mg** | 冻干保护剂 / 渗透压调节 |

**每 0.25 mL（6 月–11 岁）剂量**：mRNA **25 mcg**，总脂质 **0.50 mg**，蔗糖 21.8 mg（其余等比减半）。

- **不含防腐剂**；预充式注射器的胶塞/活塞不含天然乳胶。

## 3 DOSAGE FORMS AND STRENGTHS

- 6 月–11 岁：单剂 0.25 mL
- 12 岁及以上：单剂 0.5 mL

## 2 DOSAGE AND ADMINISTRATION（要点）

- **仅供肌内注射**
- 6 月–23 月龄：按既往接种史给 1–2 剂（0.25 mL/剂），首剂 month 0，第二剂 month 1
- 2 岁及以上：单剂（0.25 mL 或 0.5 mL）；若既往接种过，需间隔 ≥2 个月
- 使用前确认标签为 **2026-2027 Formula**

## 16 HOW SUPPLIED / STORAGE

- 冻存 **-50 °C ~ -15 °C**
- 解冻后 2–8 °C 可存 **最多 60 天**（或至效期，取先到者）
- 解冻后 8–25 °C 可存 **最多 12 小时**
- **不可再冷冻**；避免阳光与紫外直射

---

## ⚠️ 两处必须注意的"label 没写"的事

### 1. label 只写 "nucleoside-modified"，**没有写 m1Ψ**

具体是 **N1-甲基假尿苷（m1Ψ）完全替代 U**，来自公开文献而非说明书：

- Karikó & Weissman 2005 (*Immunity*)：修饰核苷让 RNA 躲过 TLR3/7/8
- Karikó et al. 2008 (*Mol Ther*)：假尿苷同时**提高翻译**
- 2023 年诺贝尔生理学或医学奖授予这两位

**写法差异的原因**：FDA label 按"成分描述"写大类，不写化学细节。**凡涉及化学细节，必须回到文献或专利核对**——这与 D3 那个 GalNAc 位置的教训是同一类。

### 2. 剂量从 100 mcg 降到了 50 mcg

原始 mRNA-1273（2020 EUA）是 **100 mcg**；2026-2027 Formula 成人剂量为 **50 mcg**。这是真实发生的**剂量下调（down-titration）**——免疫原性可控后降低剂量，是 mRNA 疫苗平台成熟的标志。

**算分子数时两个都要算**（见笔记）。

---

## 需要另外找文献确认的（label 里没有）

- [ ] mRNA 全长（公开值约 3,963 nt，需回原始文献/专利核对）
- [ ] poly(A) 尾长度（约 100 nt，工艺相关，通常不公开）
- [ ] UTR 序列来源（Moderna 自有专利序列）
- [ ] dsRNA 副产物的控制标准（纯度关键质量属性）

> 这四项正是 mRNA 设计里**计算/AI 密度最高**的部分，也是阶段二 P3 项目（mRNA 设计工具）要动手做的。
