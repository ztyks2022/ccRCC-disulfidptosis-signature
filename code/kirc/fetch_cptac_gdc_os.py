## 从 GDC 拿 CPTAC-3 完整生存(死亡时间 + 末次随访),补全 cBioPortal CPTAC 表达队列的 OS
import urllib.request, json, csv, io
body={"filters":{"op":"in","content":{"field":"project.project_id","value":["CPTAC-3"]}},
      "fields":"submitter_id,demographic.vital_status,demographic.days_to_death,diagnoses.days_to_last_follow_up",
      "size":3000,"format":"TSV"}
req=urllib.request.Request("https://api.gdc.cancer.gov/cases",
    data=json.dumps(body).encode(), headers={"Content-Type":"application/json"}, method="POST")
tsv=urllib.request.urlopen(req, timeout=120).read().decode()
rd=list(csv.DictReader(io.StringIO(tsv), delimiter="\t"))
print("GDC CPTAC-3 病例:", len(rd), "| 列:", rd[0].keys() if rd else "无")
def col(d, key):
    for k in d:
        if key in k: return d[k]
    return ""
os_map={}
for r in rd:
    sid=col(r,"submitter_id"); vit=col(r,"vital_status").lower()
    dod=col(r,"days_to_death"); fu=col(r,"days_to_last_follow_up")
    event=1 if "dead" in vit else (0 if "alive" in vit else None)
    t = dod if (dod not in ("","--",None)) else fu
    try: t=float(t)
    except: t=None
    if sid and event is not None and t is not None:
        os_map[sid]={"event":event,"time_days":t}
print("有完整OS的病例:", len(os_map))

# 合并表达
rows=list(csv.DictReader(open("data/raw/extval/cptac_ccrcc.csv")))
out=[]; n_os=0
for r in rows:
    case="-".join(r["sample"].split("-")[:2])   # C3L-00004-01 -> C3L-00004
    o=os_map.get(case)
    if o: r["time_days"]=o["time_days"]; r["event"]=o["event"]; n_os+=1
    else: r["time_days"]=""; r["event"]=""
    out.append(r)
cols=list(rows[0].keys())+["time_days","event"]
with open("data/raw/extval/cptac_ccrcc_os.csv","w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out)
print("写出 cptac_ccrcc_os.csv | 表达样本", len(out), "| 其中有OS", n_os)
