# D21 · 代码日：scan_sites.py（学习包）

> 日期：2026-11-14 周六 ｜ 时长：2.5 h ｜ W3 第 4 天
> 前置：`位点层_靶哪一段_D20.md` §7（10 条判据 = 今天的合同）、`D16实操_oligo_decorate.md`
> 标签：#代码日 #位点层 #RNAplfold #scan_sites #D21

---

## 0. 今天的定位：两个任务，一硬一软

**Part A（15 min，先做）**：跑 `test_expected.py`，收尾拖了 5 天的 `oligo_decorate.py` code review。
**Part B（135 min，主体）**：把 D20 §7 的 10 条判据写成 `scan_sites.py`，对一条真实转录本产出**候选位点排序表**。

今天不用懂新概念——判据 D20 都给过了，今天是**把规格变成代码**。

**交付物（下班前应存在的东西）**：
1. `scan_sites.py`（10 个函数，纯函数部分测试全绿）
2. `TTR_sites_ranked.tsv`（一条真实转录本的候选位点排序表）
3. 发给我：断言结果 + 两个设计决定（见 §1）

---

## 1. Part A（15 min）：oligo_decorate 断言收尾

在服务器上跑：

```bash
python test_expected.py
```

然后**发我两样东西**（这是拖了 5 天的 code review 的全部材料）：

1. **21 条断言的结果清单**：PASS 几条、FAIL 几条；FAIL 的把报错行原样贴过来
2. **两个设计决定的答案**：
   - **D-A**：`parse()` 遇到未定义的修饰字母（如 `e` 进来但函数没实现）你是**报错**还是**当普通碱基跳过**？为什么？
   - **D-B**：9 个特征的输出是**每个特征一个数**还是**拼成一个 dict**？调用方（打分器）怎么消费最方便？

**评判标准我先给你**（免得你觉得是主观题）：D-A 报错 > 跳过——静默跳过会把"没实现的特征"伪装成"特征值为 0"，这是和你 D20 学的"SNP 掩蔽静默出错"同一类事故。D-B dict > 单数——特征会增删，dict 是稳定接口。

### Part A 结果（2026-11-12 实测）

**21/21 全部通过**（T0–T5 六个用例）。逐项核对无假绿：

| 用例 | 覆盖到的东西 | 关键数字 |
|---|---|---|
| T0 | 裸串 + `d` 前缀 + `*` + `m` 混合 | len 20 / PS 2 / DNA-CpG **2（全串 CG=3）** / GC 45.0 / 2′-OMe 10 |
| T1 | 全 `e`（2′-MOE）+ 满铺 PS | len 18 / PS **17** / MOE 18 / DNA-CpG 0 / U-run@4 |
| T2 | Patisiran 风格：无 `*` | PS **0** / 2′-OMe 10 / GC 52.4 |
| T3 | 单一基序计数 | GUCCUUCAA @4 |
| T4 | 多重基序 | AUUUA @**[1, 8, 15]** |
| T5 | `I`（2′-O-甲基肌苷） | len 11 / PS 1 / 2′-OMe 2 / GC 72.7（I 按 G 计） |

**这 21 条证明了什么**：解析器（前缀优先级、`*` 归属、`I` 不崩）+ 六个特征的算术（含"DNA-CpG 只数 DNA 区"这个 T0=2 vs 全串=3 的判别点，这是最容易写错的地方，过了说明 `d` 前缀的作用域判断是对的）。

**这 21 条没证明什么**（code review 要补的四个缺口）：
1. **未知修饰字母**（D-A）——测例里没有非法输入
2. **双链输入**（Patisiran 是双链，`/` 或两条链怎么给？报告是否分链）
3. **基序重叠匹配**（`AUUUUA` 里 AUUUA 该报 1 还是 2；T4 的三个位点互不重叠，测不出）
4. **CLI 层**（空输入、超长序列、GC 除零）

### Part A code review 结论（2026-11-12，`parse_notation` 已审）

**两道设计题：都答对了。**
- **D-A**：未知碱基/前缀**抛 `ValueError`**（实测 `XmC` → "非法碱基 'X'"；`qC` → "非法前缀 'q'"）；空串也抛。✅
- **D-B**：`summarize_residues()` 返回 **dict**（`info["bases"]` / `info["mod"]`）。✅

**四档结论**（🔴 会算错数 / 🟡 会静默出错 / 🟢 风格）：

| 级别 | 问题 | 实测证据 | 修法 |
|---|---|---|---|
| 🔴 | `*` 分支清空 `buf` 前不检查 → **未完成前缀被静默吞掉** | `Am*C` → `(['A','C'], 1)`，`m` 消失 | `*` 分支先 `if buf: raise` |
| 🟡 | 前缀只校验 `buf[-1]` → 多字符前缀蒙混过关 | `mmC` → 残基 `"mmC"`；`dmC` → `"dmC"` | 校验**整个** `buf ∈ {"d","m","f","e","l"}` |
| 🟡 | 只返回 PS **数量**，丢位置 | `test_expected.py` 被迫自己重扫串还原位置（重复逻辑） | 拆 `parse_notation_full()` 返回 `ps_positions`，原函数做兼容包装 |
| 🟡 | 不容忍空格/换行 | `"mC mG"` → 非法字符 `' '`（从 PDF 粘贴必踩） | 入口 `s = "".join(s.split())` |
| 🟢 | `*` 分支内 `if not residues` 是死代码（首字符 `*` 已被前置 guard 拦掉） | — | 可删；或保留但注明 |
| 🟢 | 报错位置用 `len(residues)` 而非字符下标 | 报"位置 0"但实际在串中间 | 用 `enumerate` 的 `i` |
| 🟢 | 500 nt 无长度警告 | `'A'*500` 正常返回 500 残基 | `MAX_LEN=50` 上限检查 |

**补丁已回归验证**：修后 T0–T5 六个用例输出与原实现**完全一致**（SAME），21 条断言不受影响；`Am*C`/`mmC`/`dmC` 转为显式报错；空格/换行/首尾空白被容忍。
详见 `05_踩坑记录/20261112_解析器静默降级.md`。

**一句话教训**：**21 条正向断言全绿，三个非法输入就露了底。** 反向测试比正向测试更能说明质量。

---

## 2. Part B 范围冻结：v0 的三个决定（先读再写码）

写码前先冻结范围，不然 2.5 h 必然溢出：

| # | 决定 | v0 | v1（不做） | 理由 |
|---|---|---|---|---|
| 1 | 变异坐标 | **转录本坐标**（TSV：`pos,ref,alt,maf`） | gnomAD 基因组坐标 + 负链翻转 | 把 D20 说的"链翻转静默出错"**挡在模块外**：v0 输入就在转录本坐标系，根本不碰基因组。翻转逻辑等 v1 接 gnomAD 时再写、再单独测试 |
| 2 | 覆盖率 | **精确子串匹配**（ACGU） | Smith-Waterman 容错比对 | 21-mer 全匹配才是 siRNA 的真实场景；容错比对是序列层（③）的事 |
| 3 | 可及性 | **subprocess 调 RNAplfold，解析 lunp 文件** | 装 pybind 绑定 / Cython 化 | 别在环境上花时间；文件接口最稳 |

> 这三个决定本身就是 D18 的 A 级判据在起作用：每个都写得出伪代码，且不偷换概念。

---

## 3. 模块结构（照这个写，别自由发挥）

```
scan_sites.py
├── 纯函数区（不碰 I/O，今天必须全绿）
│   ├── make_windows(seq, L, step) -> list[(start, end)]     # 0-based 半开区间
│   ├── coverage(seq, txs) -> float                          # 精确子串，命中数/总数
│   ├── mode_compartment(modality) -> str                    # 模态 → 工作区室
│   ├── compartment_ok(region, modality) -> bool
│   ├── is_seed(pos, L) -> bool                              # pos 为窗口内 1-based
│   ├── window_access(acc, s, e) -> float                    # P(unpaired) 切片均值
│   └── pop_coverage(variants) -> float                      # 1 − Σ maf（下限 0）
├── I/O 区（今天跑通即可，不做固定期望值测试）
│   ├── load_transcripts(fasta_dir) -> dict{名字: 序列}
│   ├── accessibility(fasta_path) -> list[float]             # 调 RNAplfold 解析 lunp
│   └── load_variants(tsv_path) -> list[(maf, s, e)]         # s,e = 变异在转录本的区间
└── 驱动区
    └── rank(fasta_dir, variants_tsv, modality) -> 打印排序表（stretch goal）
```

**骨架已放在 `scripts/scan_sites_skeleton.py`**，签名和 docstring 都写好了，你原位填实现（**不许加在文件末尾**——D16 的教训，`__main__` 块在文件中部，追加会撞上）。

**测试在 `scripts/test_scan_sites.py`**：12 条纯函数断言 + 2 条 I/O 冒烟测试（文件不存在自动跳过）。期望值全部是我独立算的，**改数字前先问我**（和 test_expected.py 同一规矩）。

---

## 4. 时间块（135 min 怎么切）

| 时段 | 分钟 | 做什么 | 完成标志 |
|---|---|---|---|
| B0 | 15 | 环境：`conda install -c bioconda viennarna`（或 `apt install vienna-rna`）；验证 `RNAplfold --help` | 命令有输出 |
| B1 | 20 | 数据：下载 TTR 转录本 FASTA（§5）；造一个 5 行的变异 TSV（§6） | 两个文件在服务器上 |
| B2 | 30 | 填 7 个纯函数，跑 `test_scan_sites.py` → **12/12 绿** | 全 PASS |
| B3 | 30 | `accessibility()`：跑 RNAplfold → 解析 lunp → 接通 `window_access` | 对 TTR 输出长度 = 序列长度的数组 |
| B4 | 25 | 变异叠加 + `pop_coverage` + 排序打印 | 产出 TSV |
| B5 | 15 | 三条发现写进笔记 + commit + push | 记录在案 |

> B2 全绿之前**不许碰 B3**。纯函数不过就接 I/O，出问题你分不清是解析错还是逻辑错——这和写测试一个道理：先隔离。

---

## 5. 数据：用什么转录本

**用 TTR（转甲状腺素蛋白）**——Patisiran 的靶点，全程呼应。MANE Select 是 `NM_000371.x`（具体版本号以 NCBI 当前 MANE 页为准，🟡 我没核）：

```bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_000371&rettype=fasta&retmode=text" -o ttr.fa
```

注意：这条 FASTA 是**mRNA 全长**（含 UTR），序列本身永远是 5′→3′ 正向——所以 v0 内部不存在链翻转问题（决定 #1 的好处）。

变异 TSV 造 5 行就够（v0 是管道验证，不是真数据）：

```
pos    ref    alt    maf
120    A      G      0.02
135    C      T      0.01
210    G      A      0.05
301    T      C      0.003
450    A      G      0.12
```

---

## 6. RNAplfold 接入细节（B3 的全部内容）

```bash
RNAplfold -W 80 -L 40 < ttr.fa   # 在当前目录生成 lunp 文件
```

解析 `lunp` 的三条要点：

1. **跳过以 `#` 开头的行**（注释头），数据行是：位置（1-based）＋ P(unpaired)（还有其它列，取第 2 列即可）
2. 输出数组长度应等于序列长度——**拿 `len(seq) == len(acc)` 当断言**，这是你接通 I/O 的第一个测试
3. 默认单位是**每碱基**的 P(unpaired)，直接切片取均值就是 `window_access`

三条坑提前说：
- `pos` 是 1-based，你的窗口是 0-based 半开——换算处写一行注释，明天 review 时我会专门看
- FASTA 里可能有换行断行，读入后**先去掉所有空白**再喂给 RNAplfold 和自己的函数
- 大转录本没问题，但**别一次喂多转录本的文件**——lunp 只保留最后一条（v0 单转录本，就绕开了）

---

## 7. 产出表：TTR_sites_ranked.tsv 的列

```
rank  start  end  seq  gc  mean_access  seed_hits  pop_cov  score
```

- `score` = D20 §6 五维的简单加和（v0 不加权——和 oligo_decorate 的全局决定一致）
- 排序按 score 降序；**top 10 里肉眼看两眼**：高分的位点是不是偏向 CDS/3′UTR？有没有 MAF≥5% 的变异混进种子区？——看出来的东西写进三条发现（B5）

---

## 8. 自检 5 问（做完再答）

1. 你的 `make_windows` 对 L > 序列长返回什么？为什么不抛异常？
2. `pop_coverage` 直接 Σ maf 的独立假设在什么情况下会低估覆盖率？
3. lunp 的 pos 是 1-based、你的窗口是 0-based——你在哪一行做的换算？怎么保证没 off-by-one？
4. B2 为什么必须在 B3 之前全绿？
5. Part A 的 D-A 设计题：报错和静默跳过的区别，和 D20 哪个知识点同源？

---

## 9. 明日预告（D22 · 验证层）

RNA-seq 脱靶分析 SOP 骨架：GUIDE-seq/CLEAN-seq 数据长什么样、比对到转录组还是基因组、怎么从 BAM 里数脱靶位点。这是你⑥⑦验证层的第一格——老本行，但要把它写成**能给湿实验同事用的 SOP**。

---

## 附：本文事实边界

| 陈述 | 状态 |
|---|---|
| TTR MANE Select = NM_000371.x | 🟡 版本号待核 NCBI MANE 页 |
| RNAplfold `-W 80 -L 40` 输出 lunp、每碱基 P(unpaired) | ✅ ViennaRNA 文档 |
| efetch 接口拉 FASTA | ✅ NCBI E-utilities 标准用法 |
| conda 装 viennarna | ✅ bioconda 官方渠道 |
| D-A/D-B 评判标准 | ✅ 工程判断，非事实陈述 |
