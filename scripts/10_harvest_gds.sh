#!/usr/bin/env bash
# eutils 多查询并集收割 GEO PE 单细胞数据集 → gds_summary.tsv
set -u
D=/root/pe_signature/data
mkdir -p $D
QUERIES=(
  "preeclampsia+AND+single+cell+RNA-seq"
  "preeclampsia+AND+scRNA-seq"
  "preeclampsia+AND+single-cell"
  "preeclampsia+AND+single+nucleus"
  "preeclampsia+AND+single+cell"
  "preeclampsia+AND+snRNA-seq"
)
> $D/gds_ids_all.txt
for Q in "${QUERIES[@]}"; do
  curl -s --max-time 90 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=${Q}&retmax=500&retmode=json" \
    | /root/miniconda3/bin/python -c "import sys,json; print('\n'.join(json.load(sys.stdin)['esearchresult']['idlist']))" >> $D/gds_ids_all.txt
  sleep 1
done
sort -u $D/gds_ids_all.txt -o $D/gds_ids_all.txt
IDS=$(paste -sd, $D/gds_ids_all.txt)
echo "unique ids: $(wc -l < $D/gds_ids_all.txt)"
curl -s --max-time 120 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=${IDS}&retmode=json" -o $D/gds_summary.json
/root/miniconda3/bin/python - <<'PY'
import json
d = json.load(open('/root/pe_signature/data/gds_summary.json'))['result']
rows = []
for uid in d['uids']:
    r = d[uid]
    acc = r.get('accession', '')
    if not acc.startswith('GSE'):
        continue
    rows.append((acc, r.get('n_samples', ''), r.get('gpl', ''), r.get('title', '')[:110], r.get('summary', '')[:200].replace('\n', ' ')))
rows.sort()
with open('/root/pe_signature/data/gds_summary.tsv', 'w') as f:
    f.write('accession\tn_samples\tplatform\ttitle\tsummary\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
print('GSE rows:', len(rows))
PY
echo HARVEST-DONE
