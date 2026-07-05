## 从 cBioPortal 拉 ClearCode34 基因在 CPTAC ccRCC 的 RNA-seq TPM(做公平外部头对头)
import urllib.request, json, csv
API="https://www.cbioportal.org/api"; PROF="rcc_cptac_gdc_mrna_seq_tpm"; SLIST="rcc_cptac_gdc_all"
ccA=["MAPT","STK32B","FZD1","RGS5","GIPC2","PDGFD","EPAS1","MAOB","CDH5","TCEA3","LEPROTL1","BNIP3L","EHBP1","VCAM1","PHYH","PRKAA2","SLC4A4","ESD","TLR3","NRP1","ST13","ARNT"]
ccB=["SERPINA3","SLC4A3","MOXD1","KCNN4","ROR2","FOXM1","UNG","GALNT10","GALNT4"]
genes=ccA+ccB
def post(url,body):
    req=urllib.request.Request(url,data=json.dumps(body).encode(),
        headers={"Content-Type":"application/json","Accept":"application/json"},method="POST")
    return json.load(urllib.request.urlopen(req,timeout=120))
# symbol -> entrez
gm=post(f"{API}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", genes)
e2s={g["entrezGeneId"]:g["hugoGeneSymbol"] for g in gm}
print("映射到 entrez:", len(e2s), "/", len(genes))
md=post(f"{API}/molecular-profiles/{PROF}/molecular-data/fetch?projection=SUMMARY",
        {"sampleListId":SLIST,"entrezGeneIds":list(e2s)})
expr={}
for r in md:
    s=e2s.get(r["entrezGeneId"])
    if s: expr.setdefault(r["sampleId"],{})[s]=r["value"]
present=sorted({g for e in expr.values() for g in e})
with open("data/raw/extval/cptac_cc34.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["sample"]+present)
    for s,e in expr.items(): w.writerow([s]+[e.get(g,"") for g in present])
print("写出 cptac_cc34.csv | 样本", len(expr), "| 基因", len(present))
print("ccA命中", len([g for g in present if g in ccA]), "/", len(ccA), "| ccB命中", len([g for g in present if g in ccB]),"/",len(ccB))
