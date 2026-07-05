## 把图编号改成与正文引用顺序一致(顺序引用,期刊要求)。old->new 映射,token 两遍防碰撞。
import re, os, glob
M = {1:1, 2:3, 3:2, 4:4, 5:6, 6:7, 7:8, 8:10, 9:11, 10:12, 11:5, 12:9}  # old -> new

# 1) 正文 + 图注里的 "Fig. N"
f = "manuscript_KIRC_disulfidptosis.md"
s = open(f).read()
s = re.sub(r'Fig\. (\d+)', lambda m: f'Fig. __T{M[int(m.group(1))]}__', s)   # 第一遍打 token
s = re.sub(r'__T(\d+)__', r'\1', s)                                          # 第二遍落地
open(f, "w").write(s)
print("manuscript 重编号完成")

# 2) 结果图文件(v2)同步改名,temp 两遍防碰撞
files = glob.glob("results/Fig[0-9]*_v2.*")
tmp = []
for p in files:
    d, b = os.path.dirname(p), os.path.basename(p)
    mm = re.match(r'Fig(\d+)(\D.*)$', b)
    if not mm: continue
    old, rest = int(mm.group(1)), mm.group(2)
    if M[old] == old: continue
    t = os.path.join(d, f"FigTMP{old}{rest}")
    os.rename(p, t); tmp.append((t, old, rest, d))
for t, old, rest, d in tmp:
    os.rename(t, os.path.join(d, f"Fig{M[old]}{rest}"))
print(f"重命名 {len(tmp)} 个图文件")
