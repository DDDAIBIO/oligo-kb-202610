# FAQ：PMO / ASO gapmer / siRNA 三类载荷机制对比（刀、找人、化学、区室）

> 来源：D25–D28 后的追问链（2026-11-18）：siRNA vs ASO 四层对比 → gapmer 结构三层拆解 → 占位型（位阻型）澄清。
> 关联文档：`04_FAQ/FAQ_siRNA作用机制五问_v0.md`（Q5 同源机理）、`04_FAQ/FAQ_启动周四问_v0.md`（gapmer/占位型反转逻辑首记）、`02_技术方案/ASO机制笔记_v0.md`（§91 Nusinersen 为何不能有 gap）、`03_SOP/已上市寡核苷酸药物对比表_v0.md`（#1 / #8 / #10 / #11）、`03_SOP/化学修饰速查表_v0.md` §2.4。
> 适用范围：模态选型（"这个靶点该用哪类载荷"）与脱靶分析流程分叉（D22 SOP 的模态前置判断）。

---

## 0. 一句话总览

> 三类载荷 = "刀、找人方式、形态代价、工作区室"四个自由度的不同组合：**siRNA 自带刀（AGO2 PIWI）+ seed 找人 + 双链 + 限胞质；gapmer 招刀（RNase H1，需 DNA gap）+ 全长配对 + 单链 + 胞质核都行；PMO 没有刀（纯位阻）+ 全长配对 + 单链非天然骨架 + 核内 pre-mRNA。**

选型判断链（第一问永远是区室）：

> 靶在胞质、要高效沉默 → **siRNA**；靶在胞质或核、要降解成熟 mRNA / pre-mRNA → **ASO gapmer**；靶在核内 pre-mRNA、要改剪接（挡住而不切） → **占位型 ASO（2′-MOE 全修饰）或 PMO**。

---

## 1. 总对比表（背这张就够）

| 维度 | siRNA | ASO gapmer（降解型） | 占位型 ASO（2′-MOE 系） | PMO |
| --- | --- | --- | --- | --- |
| **代表药** | Patisiran / Inclisiran（对比表 #11/#12） | Donidalorsen（#8，5-10-5）；Mipomersen、Inotersen（已撤） | Nusinersen / Tofersen（鞘内） | Eteplirsen + 3 个 DMD 同门（#10，裸 PMO 静脉） |
| **靶 RNA** | 胞质成熟 mRNA | 成熟 mRNA 或核内 pre-mRNA 均可 | **核内 pre-mRNA**（如 SMN2 intron 7 ISS-N1） | **核内 pre-mRNA**（如 DMD exon 51） |
| **工作区室** | **胞质**（RISC 在胞质，进不了核） | 胞质 + 核（单链小分子均可进） | 核为主 | 核为主 |
| **刀** | **AGO2 自带 PIWI**（RNase H 样催化中心） | **招募内源 RNase H1**——酶本来就在，ASO 递"刀鞘" | **没有刀**：纯空间位阻，挡住剪接机器/核糖体 | **没有刀**：纯空间位阻 |
| **下刀条件** | 引导链全程配对 → PIWI 在靶 10/11 位切 | DNA gap 与靶 RNA 成 **RNA/DNA 杂合体**，RNase H1 只认杂合体 | 无（目标是"别切"） | 无（目标是"别切"） |
| **找人方式** | **seed（2–8）成核**，7 nt 就能扣住靶，再向 3′ 延伸 | 全长（至少 gap 区）持续配对 | 全长持续配对 | 全长持续配对 |
| **主导脱靶类别** | **B 类 seed 介导**（miRNA 样抑制，完全匹配搜索看不见）→ D22 要 seed 富集 | **A 类序列依赖**（≥15-mer 完全匹配即可抓） | A 类 + 位置类（挡错位点） | 同左 |
| **分子形态** | **双链**：passenger（sense）是负担 | 单链嵌合体：翼 + gap 分区 | 单链均质 | 单链均质 |
| **化学设计** | 2′-F/2′-OMe 高密度 + 少量 PS（Inclisiran 6 个）+ 可挂 GalNAc | 翼 2′-MOE + **中间 DNA gap（功能命门）** + PS 全场或近全场（donidalorsen 15/19） | **全场 2′-MOE + 全长 PS**（Nusinersen 17/17，无 gap） | **六元吗啉环 + 中性磷酸二酰胺骨架**：无核糖、无 PS、无 2′ 修饰、**不带电** |
| **这条化学的代价** | passenger 错误装载 = 打错靶（热力学不对称解决，见五问 Q1） | **PS 非特异蛋白结合毒性**（inotersen 血小板毒性撤市实例） | PS 账同样在；鞘内绕开全身暴露 | **活性低 → 靠剂量取胜**：Eteplirsen 30 mg/kg Q1W，年剂量 ≈109 g（Vutrisiran 的 1092 倍），肌肉 dystrophin 仅 0.44–0.93% |
| **给药现状** | LNP 静脉（Patisiran）/ GalNAc 皮下（Inclisiran） | GalNAc 皮下（donidalorsen）/ 裸皮下（inotersen） | 裸寡核苷酸**鞘内**（CNS 隔室） | 裸 PMO **静脉**（肌肉隔室） |

---

## 2. 逐层拆解

### 2.1 刀的来源：自带、招来、没有——三种机制哲学

- **siRNA 自带刀**：切割酶是 RISC 自带的（AGO2 PIWI），引导链只负责认路。一个 AGO2 可循环切多条 mRNA（催化级联），所以药效可以很强、剂量可以很小。
- **gapmer 招刀**：刀是细胞里本来就有的 RNase H1，口味挑剔——只切 RNA/DNA 杂合双链。所以 gapmer 必须留一段**裸 DNA**（donidalorsen 的 5-10-5 中间 10 个），翼部 2′-MOE 只是护甲和亲和力。**gap 的有无由机制决定**：要招刀就得留 DNA。
- **占位型/PMO 没有刀**：目的恰好相反——**绝不能让靶 RNA 被切**。Nusinersen 占据 ISS-N1 挤走 hnRNP A1 保住 exon 7；只要分子里有一段 DNA 成杂合体，RNase H1 就会把 SMN2 pre-mRNA 切掉，机制直接反转（ASO机制笔记 §91）。所以占位型**全场修饰**（Nusinersen 全 2′-MOE），PMO 更是整体换成非天然骨架。
- **PMO 的极端化**：连"核糖 + 带电磷酸骨架"都整个换掉——吗啉环 + 中性磷酸二酰胺。不带电 → 血浆蛋白结合少 → PS 类毒性谱基本消失；代价是靶亲和力低、活性弱，**只能靠剂量取胜**（年 109 g 对比 siRNA 的毫克级）。

### 2.2 找人方式 → 脱靶分析流程必须分模态（D22 前置判断）

- siRNA 靠 seed 成核，6–7 nt 同源就能扣住靶 → 脱靶面广且**完全匹配搜索看不见**（B 类），必须 RNA-seq seed 富集。
- gapmer / 占位型没有 seed 机制，需要全长（或 gap 区）配对才起效 → 脱靶主体是较长完全匹配的 A 类，BLAST/offtarget_scan（≥15-mer）基本够用。
- **这就是 D22 SOP 里 siRNA 与 gapmer 流程不能照抄同一份的原因。**

### 2.3 区室与靶类型：选型的第一问

- RISC 在胞质 → siRNA 打不了核内 pre-mRNA，**打不了内含子**（D20 结论 1）。
- 单链 ASO/PMO 能进核 → 剪接调变（打 intron 7 的 ISS-N1、打 exon 51）是它们独有的功能类别。
- 反过来，胞质成熟 mRNA 的高效沉默，gapmer 与 siRNA 都能做，此时用"要不要皮下给药（GalNAc）、剂量敏感性、passenger 设计成熟度"等工程因素定夺。

### 2.4 三条"死亡账"对照

- **siRNA 付 passenger 的账**：双链里 sense 链纯是包装，设计不当会被误装进 AGO2 打向 sense 互补的转录本（热力学不对称 + 2′ 修饰压 passenger 装载，见五问 Q1）。
- **gapmer 付 PS 护甲的账**：单链裸 RNA 几分钟降解 → 全身披 PS → 非特异蛋白结合毒性（inotersen 血小板减少撤市；mipomersen 肝毒性撤市——速查表"撤市三种死法"之二）。
- **PMO 付剂量的账**：中性骨架低毒但低效，年 109 g 静脉裸奔，疗效还只有百分之零点几的 dystrophin——**"安全"和"有效"在这里是一对直接的交换**。

---

## 3. 与知识库其他文档的关系（不重复的内容放哪里）

| 内容 | 已在文档 |
| --- | --- |
| siRNA 机理四问 + Inclisiran 皮下 | `FAQ_siRNA作用机制五问_v0.md` Q1–Q4 |
| siRNA vs ASO 机理四层 | 同上 Q5（本文是它的模态扩展版） |
| gapmer/占位型反转逻辑（首次记录） | `FAQ_启动周四问_v0.md`（PS 那一问） |
| Nusinersen 案例拆解 | `ASO机制笔记_v0.md` 案例 1/3 |
| 各药化学账（PS 数、修饰数、剂量） | 对比表 #1/#8/#10/#11、速查表 §2.4 |
| 脱靶两类分报 | `验证层_脱靶分析SOP_D22.md` |

> 本文只保留**跨模态对比**这个视角；单模态深挖一律引上面的文档，避免双份维护。
