"""offtarget_scan_skeleton.py —— D29 工具日②骨架：药物序列 → 转录本库脱靶扫描

契约（D22 §5 规格 + W4 细则 D29 行，v0）：
  - 字母表统一为 DNA 大写 ACGT（含 N）；RNA 输入单点归一化 U→T（D28 KeyError 教训的复用）
  - 方向：guide 与靶 RNA 反向平行 → 一律在转录本里搜 revcomp(query)（D20/D28 方向坑再现位）
  - 三档输出（D22 §5）：
      exact  : 全长完全匹配 → A 类铁证
      contig : >=min_len 连续匹配的极大 block（默认 min_len=12，D22 标 🟡 待核；
               与 exact 不重复报——contig 跳过 L==len(query) 的块）
      seed   : guide 1-based 2–8 位（7 nt）的 revcomp 在转录本中的匹配 → B 类候选
               （v0 口径：只算 guide 链的 seed；无 3'UTR 注释，位置标注"未定性"）
  - 双链：A 类（exact/contig）guide 与 passenger 都扫（D22 §2：两条链都能带脱靶）
  - on-target 排除：scan_fasta 层按基因名排除（默认 TTR——打中的不叫脱靶）
  - 排序：n_exact 降序 → max_contig 降序 → gene 名升序（可复现纪律）
  - 位置全部 0-based（转录本内偏移）
  - v0 朴素字符串搜索：不下 Bowtie、不建索引（W4 风险预案）；只对小转录本集
"""
import csv


# ---- S0 归一化 ----------------------------------------------------------
def _norm(seq):
    """大写 + RNA→DNA（U->T）。任何进搜索器的序列都先过这里（单点归一化）。
    >>> _norm("acgu")
    'ACGT'
    """
    #raise NotImplementedError("S0")
    return seq.upper().replace("U", "T")


# ---- S1 反向互补 --------------------------------------------------------
def revcomp_dna(seq):
    """DNA 反向互补。N 保留为 N。
    >>> revcomp_dna('ATGC')
    'GCAT'
    """
    #raise NotImplementedError("S1")
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


# ---- S2 FASTA 头 -> 基因名 ----------------------------------------------
def gene_from_header(header):
    """取 FASTA 头第一个空白分隔 token（去掉开头的 '>'）。
    >>> gene_from_header('>TTR Homo sapiens chr18')
    'TTR'
    >>> gene_from_header('ENST00000221962.4')
    'ENST00000221962.4'
    """
    #raise NotImplementedError("S2")
    return header.strip().lstrip(">").split()[0]



# ---- S3 全长精确匹配 ----------------------------------------------------
def find_exact(query, subject):
    """query 在 subject 中所有出现的起始位置（0-based 升序）。无命中返回 []。"""
    #raise NotImplementedError("S3")
    n = len(query)
    if n == 0 or n > len(subject):
        return []
    return [j for j in range(len(subject) - n + 1) if subject[j:j + n] == query]



# ---- S4 连续匹配（极大 block）------------------------------------------
def find_contig(query, subject, min_len=12):
    """query 与 subject 的所有极大连续匹配 block，返回 [(i, j, L), ...]：
      i = query 内起点, j = subject 内起点, L = 长度，均 0-based。
      极大 = 左端（s[j-1]==q[i-1] 且两端都在界内时不成立）和右端都不能再延伸。
      L == len(query) 的块不报（那是 exact 的活，避免重复）。
      排序：L 降序 -> j 升序 -> i 升序。L < min_len 的块不报。
    """
    #raise NotImplementedError("S4")
    n, m = len(query), len(subject)
    hits = []
    if n == 0 or m == 0:
        return hits
    for i in range(n):                    # query 的每个起点
        j = 0
        while j < m:                      # 扫 subject
            if subject[j] == query[i]:
                k = 0
                while j + k < m and i + k < n and subject[j + k] == query[i + k]:
                    k += 1                # 尽量延伸
                left_ok = (j == 0 or i == 0 or subject[j - 1] != query[i - 1])
                right_ok = (j + k == m or i + k == n or subject[j + k] != query[i + k])
                if k >= min_len and left_ok and right_ok and k < n:
                    hits.append((i, j, k))
                j += max(k, 1)            # 跳过这一块，别重扫
            else:
                j += 1
    return sorted(set(hits), key=lambda t: (-t[2], t[1], t[0]))


# ---- S5 seed 匹配 -------------------------------------------------------
def find_seed(guide, subject):
    """guide 的 1-based 2–8 位（切片 [1:8]）反向互补后在 subject 中的所有起始位置。
    返回 [j, ...] 升序。B 类候选（v0 无 3'UTR 注释，位置未定性）。"""
    #raise NotImplementedError("S5")
    g = _norm(guide)
    seed_rc = revcomp_dna(g[1:8])
    return find_exact(seed_rc, subject)

# ---- S6 单转录本扫描 ----------------------------------------------------
def scan_transcript(guide, sense, header, seq, min_len=12):
    """对一条转录本扫描 guide + passenger 两链。
    返回 dict：
      gene   : gene_from_header(header)
      exact  : [(strand, j), ...]           strand in ('guide','passenger')
      contig : [(strand, i, j, L), ...]
      seed   : [j, ...]                     （只 guide，v0 口径）
    """
    #raise NotImplementedError("S6")
    s = _norm(seq)                                # 转录本归一化，一次就够
    res = {"gene": gene_from_header(header),
           "exact": [],
           "contig": [],
           "seed": find_seed(guide, s)}           # seed 只算 guide（v0 口径）
    for strand, q in (("guide", _norm(guide)), ("passenger", _norm(sense))):
        q_rc = revcomp_dna(q)
        for j in find_exact(q_rc, s):
            res["exact"].append((strand, j))
        for (i, j, L) in find_contig(q_rc, s, min_len):
            res["contig"].append((strand, i, j, L))
    return res


# ---- S7 全库扫描 --------------------------------------------------------
def _read_fasta(path):
    """极简 FASTA 读取，yield (header, seq)。服务器上可换成
    from scan_sites import _read_fasta（复用 D21），接口一致即可。"""
    #raise NotImplementedError("S7a")
    header, buf = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)   # 新头出现 = 上一条记录结束
                header, buf = line, []
            elif line:                            # 跳过空行
                buf.append(line)
    if header is not None:
        yield header, "".join(buf)                # 最后一条没有"下一个头"来触发
    return


def scan_fasta(guide, sense, fasta_path, on_target="TTR", min_len=12):
    """扫整个 fasta 库，排除 on_target 基因；无任何信号（三档全空）的基因不输出。
    返回 [dict, ...]，dict = {gene, n_exact, max_contig, n_seed}，按契约排序。"""
    #raise NotImplementedError("S7b")
    rows = []
    for header, seq in _read_fasta(fasta_path):
        gene = gene_from_header(header)
        if gene.upper() == str(on_target).upper():    # 过滤器 1：on-target 排除
            continue
        t = scan_transcript(guide, sense, header, seq, min_len)
        if not (t["exact"] or t["contig"] or t["seed"]):   # 过滤器 2：无信号不输出
            continue
        rows.append({
            "gene": gene,
            "n_exact": len(t["exact"]),
            "max_contig": max((c[3] for c in t["contig"]), default=0),
            "n_seed": len(t["seed"]),
        })
    rows.sort(key=lambda r: (-r["n_exact"], -r["max_contig"], r["gene"]))
    return rows


# ---- S8 候选批量扫描（回流接口）----------------------------------------
def scan_candidates_tsv(tsv_path, fasta_path, guide_col="guide",
                        sense_col="sense", name_col="id",
                        on_target="TTR", min_len=12, limit=None):
    """读候选 TSV（D28 的 TTR_siRNA_candidates.tsv），逐条（可选 limit）扫库。
    返回 [dict, ...]：{id, guide, n_off_genes, n_exact, max_contig, n_seed}
      n_exact/n_seed 为该候选全库合计，max_contig 为全库最长连续匹配。"""
    #raise NotImplementedError("S8")
    out = []
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if limit is not None and len(out) >= limit:   # 扫描前判断，省算力
                break
            rows = scan_fasta(r[guide_col], r[sense_col], fasta_path,
                              on_target=on_target, min_len=min_len)
            out.append({
                "id": r.get(name_col, ""),
                "guide": r[guide_col],
                "n_off_genes": len(rows),
                "n_exact": sum(x["n_exact"] for x in rows),
                "n_seed": sum(x["n_seed"] for x in rows),
                "max_contig": max((x["max_contig"] for x in rows), default=0),
            })
    return out

