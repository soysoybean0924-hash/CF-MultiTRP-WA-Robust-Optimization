from pathlib import Path
import hashlib,json,shutil
S=Path(r'F:\CF-MultiTRP-WA-Robust-Optimization_staging')
D=Path(r'D:\GitHubProjects\CF-MultiTRP-WA-Robust-Optimization')
assert S.is_dir() and D.is_dir() and not any(D.iterdir())
files=[p for p in S.rglob('*') if p.is_file() and 'matlab_pref' not in p.relative_to(S).parts]
need=sum(p.stat().st_size for p in files)
assert shutil.disk_usage(D.anchor).free>need+1024**3
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
done=0
for src in files:
 rel=src.relative_to(S);dst=D/rel
 dst.parent.mkdir(parents=True,exist_ok=True)
 if dst.exists():raise FileExistsError(dst)
 shutil.copy2(src,dst)
 if dst.stat().st_size!=src.stat().st_size or sha(dst)!=sha(src):
  raise RuntimeError(f'migration verification failed: {rel}')
 done+=1
 if done%50==0:print(f'copied and verified {done}/{len(files)}',flush=True)
report={'source':str(S),'destination':str(D),'files_copied_and_sha256_verified':done,
        'bytes_copied':need,'excluded_environment_directory':'08_一致性验证/matlab_pref',
        'failures':[],'d_free_after':shutil.disk_usage(D.anchor).free}
(D/'MIGRATION_VERIFICATION.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
