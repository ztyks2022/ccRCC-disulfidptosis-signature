## 从 cBioPortal 拉 CPTAC ccRCC (rcc_cptac_gdc) RNA-seq TPM(12签名基因)+ OS
import urllib.request, json
API="https://www.cbioportal.org/api"; SID="rcc_cptac_gdc"
PROF="rcc_cptac_gdc_mrna_seq_tpm"; SLIST="rcc_cptac_gdc_all"
ENTREZ={60:"ACTB",23607:"CD2AP",23657:"SLC7A11",2997:"GYS1",10787:"NCKAP1",
        10128:"LRPPRC",23191:"CYFIP1",7094:"TLN1",4628:"MYH10",80224:"NUBPL",
        4719:"NDUFS1",10163:"WASF2"}
def get(url):
    return json.load(urllib.request.urlopen(url, timeout=90))
def post(url, body):
    req=urllib.request.Request(url, data=json.dumps(body).encode(),
        headers={"Content-Type":"application/json","Accept":"application/json"}, method="POST")
    return json.load(urllib.request.urlopen(req, timeout=120))

# 1. sample <-> patient
samps=get(f"{API}/studies/{SID}/samples?projection=SUMMARY&pageSize=100000")
s2p={s["sampleId"]:s["patientId"] for s in samps}
print("样本数:", len(s2p))

# 2. 表达 (12 基因 x 样本)
md=post(f"{API}/molecular-profiles/{PROF}/molecular-data/fetch?projection=SUMMARY",
        {"sampleListId":SLIST,"entrezGeneIds":list(ENTREZ)})
expr={}
for r in md:
    g=ENTREZ.get(r["entrezGeneId"]);
    if g: expr.setdefault(r["sampleId"],{})[g]=r["value"]
print("有表达的样本:", len(expr))

# 3. 临床 OS (per patient)
cd=get(f"{API}/studies/{SID}/clinical-data?clinicalDataType=PATIENT&pageSize=1000000&projection=SUMMARY")
osd={}
for r in cd:
    osd.setdefault(r["patientId"],{})[r["clinicalAttributeId"]]=r["value"]

# 4. 合并写出
genes=list(ENTREZ.values())
rows=[]
for s,e in expr.items():
    p=s2p.get(s); cl=osd.get(p,{})
    osm=cl.get("OS_MONTHS"); oss=cl.get("OS_STATUS")
    rows.append([s]+[e.get(g,"") for g in genes]+[osm,oss])
import csv
with open("data/raw/extval/cptac_ccrcc.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["sample"]+genes+["OS_MONTHS","OS_STATUS"]); w.writerows(rows)
print("写出 data/raw/extval/cptac_ccrcc.csv  行数:", len(rows))
nos=sum(1 for r in rows if r[-1] and r[-2])
print("有OS的样本:", nos, "| OS_STATUS 取值示例:", list({r[-1] for r in rows if r[-1]})[:4])
