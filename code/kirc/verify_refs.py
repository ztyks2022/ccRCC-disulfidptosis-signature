## 用 CrossRef 核实每条参考文献 → 真实 DOI/作者/年/卷页(防止编造)
import urllib.request, urllib.parse, json, time
queries = [
 ("disulfidptosis","Actin cytoskeleton vulnerability disulfide stress disulfidptosis Liu 2023"),
 ("SLC7A11 xCT","Cystine transporter SLC7A11 xCT cancer Koppula 2021"),
 ("TCGA-KIRC","Comprehensive molecular characterization clear cell renal cell carcinoma 2013"),
 ("UCSC Xena","Visualizing interpreting cancer genomics data Xena platform Goldman 2020"),
 ("GSVA/ssGSEA","GSVA gene set variation analysis microarray RNA-seq Hanzelmann 2013"),
 ("glmnet/LASSO-Cox","Regularization paths generalized linear models coordinate descent Friedman 2010"),
 ("ConsensusClusterPlus","ConsensusClusterPlus class discovery visualization Wilkerson 2010"),
 ("ESTIMATE","Inferring tumour purity stromal immune cells expression ESTIMATE Yoshihara 2013"),
 ("WGCNA","WGCNA weighted gene co-expression network analysis R package Langfelder 2008"),
 ("oncoPredict","oncoPredict drug response biomarker discovery Maeser 2021"),
 ("scTenifoldKnk","scTenifoldKnk virtual knockout single cell gene regulatory Osorio"),
 ("Seurat","Integrated analysis multimodal single-cell data Seurat Hao 2021"),
 ("CPTAC ccRCC","Integrated proteogenomic characterization clear cell renal cell carcinoma Clark 2019"),
 ("FinnGen","FinnGen genome-wide association studies 500000 Finnish Kurki 2023"),
]
def cr(q):
    url="https://api.crossref.org/works?"+urllib.parse.urlencode({"query.bibliographic":q,"rows":1})
    try:
        it=json.load(urllib.request.urlopen(url,timeout=30))["message"]["items"]
        return it[0] if it else None
    except Exception as e: return {"_err":str(e)}
def fmt(m):
    if not m or "_err" in (m or {}): return "  [未核到] "+str((m or {}).get("_err",""))
    au=m.get("author",[])
    a1=(au[0].get("family","?")+" "+("".join(w[0] for w in au[0].get("given","").split()) if au[0].get("given") else "")) if au else "?"
    etal=" et al." if len(au)>1 else ""
    yr=(m.get("published-print",m.get("published-online",{})).get("date-parts",[[None]])[0][0])
    jr=(m.get("container-title") or ["?"])[0]
    vol=m.get("volume",""); pg=m.get("page",""); doi=m.get("DOI","")
    ti=(m.get("title") or ["?"])[0]
    return f"  {a1}{etal} ({yr}) {ti}. {jr} {vol}:{pg}. doi:{doi}"
for tag,q in queries:
    print(f"### {tag}"); print(fmt(cr(q))); time.sleep(0.3)
