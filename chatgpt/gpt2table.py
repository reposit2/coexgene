import re

infile = "ChatGPT_genefunction.md"
outfile = "ChatGPT_genefunction_table2.tsv"

with open(infile, "r", encoding="utf-8") as f:
    lines = f.readlines()

nlines = len(lines)
output_lines = ["Source\tGeneSetID\tGeneSetName\tGeneList\tn_Genes\tLLM Name\tLLM Analysis\tScore"]

prev2 = ""
prev1 = ""
i = 0

while i < nlines:
    l = lines[i].strip()

    # スコアとデバイスIDの取得
    match_score = re.match(r"\*\*Process\:.+\(([\d\.]+?)\)\*\*", l)
    if match_score:
        score = match_score.group(1)
        drid = ''
        match_drid = re.match(r"\*\*(.+?)\*\*", prev2)
        if match_drid:
            drid = match_drid.group(1)
        i += 1
        continue

    # 機能セクションの検出
    match_func = re.match(r"\*\*(\d+)\. (.+?)\*\*", l)
    if match_func:
        num = match_func.group(1)
        func = match_func.group(2)
        cluster = f"Cluster_{drid}_{num}"

        i += 1
        if i >= nlines:
            break
        l2 = lines[i].strip()

        i += 1
        if i >= nlines:
            break
        l3 = lines[i].strip()

        paragraph_lines = []
        while i + 2 < nlines and not lines[i + 2].startswith("**"):
            paragraph_lines.append(l3)
            i += 2
            l3 = lines[i].strip()
        paragraph_lines.append(l3)
        full_paragraph = " ".join(paragraph_lines)

        # イタリック部（遺伝子候補）の抽出
        italic_matches = re.findall(r'(?<!\*)\*([^\*]+?)\*(?!\*)', full_paragraph)
        gene_set = set()
        for match in italic_matches:
            candidates = match.split(', ')
            for y in candidates:
                gene_match = re.match(r'^([A-Z][A-Z0-9]+)', y)
                if gene_match:
                    gene_set.add(gene_match.group(1))

        genes = sorted(gene_set)
        glist = " ".join(genes)
        gn = len(genes)

        output_lines.append(f"NeST\t{cluster}\t{cluster}\t{glist}\t{gn}\t{func}\t{full_paragraph}\t{score}")

    prev2 = prev1
    prev1 = l
    i += 1

# 出力ファイルへ書き込み
with open(outfile, "w", encoding="utf-8") as f:
    for line in output_lines:
        f.write(line + "\n")

print(f"✅ 出力ファイルを作成しました: {outfile}")

