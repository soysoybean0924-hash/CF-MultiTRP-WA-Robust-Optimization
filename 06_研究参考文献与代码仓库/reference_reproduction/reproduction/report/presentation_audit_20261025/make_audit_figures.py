"""Read only frozen, complete MATLAB batch files; regenerate the threshold plot."""
from pathlib import Path
import csv, json, sys
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,
                     'pdf.fonttype':42,'ps.fonttype':42})

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1] / 'chasi_paired' / 'results'
SCENE='umi19p4_s27_u10_t2_fc3.5_g17_paired_formal'
BATCHES={'ideal_sumse':220,'nominal_sumse':200,'ideal_pf':180,'nominal_pf':180}
LABELS={'ideal_sumse':'理想 CSI＋sum-SE','nominal_sumse':'非理想 CSI＋sum-SE',
        'ideal_pf':'理想 CSI＋PF','nominal_pf':'非理想 CSI＋PF'}
COLORS=['#174d86','#ca6b30','#287955','#8550a0']
STYLES=['-','--','-.',':']
rows=[]
fig,axes=plt.subplots(2,1,figsize=(10,7.4),layout='constrained',facecolor='white')
for (branch,outer),color,style in zip(BATCHES.items(),COLORS,STYLES):
 stem=f'{SCENE}_{branch}_outer{outer}'
 summary=json.loads((ROOT/f'{stem}_summary.json').read_text(encoding='utf-8'))
 assert summary['evaluation']['outerCount']==outer and summary['evaluation']['stopReason']=='outer_limit'
 mat=loadmat(ROOT/f'{stem}.mat',simplify_cells=True)['diagnostic']
 objective=np.asarray(mat['objective'],float).reshape(-1)
 scale=np.asarray(mat['objectiveScale'],float).reshape(-1)
 diag=np.loadtxt(ROOT/f'{stem}_diagnostics.csv',delimiter=',')
 assert len(objective)==outer+1 and len(scale)==outer+1 and diag.shape==(outer,9)
 x=np.arange(1,outer+1)
 obj=np.abs(np.diff(objective))/scale[:-1]
 cov=diag[:,5]
 assert np.all(obj>0) and np.all(cov>0)
 for a,y in zip(axes,[obj,cov]):
  a.plot(x,y,color=color,linestyle=style,lw=1.5,label=f'{LABELS[branch]}（第{outer}轮）')
 for n,o,c,d in zip(x,obj,cov,diag):
  rows.append([branch,int(n),f'{o:.17g}',f'{c:.17g}',f'{d[1]/d[2]:.17g}',f'{d[3]/d[4]:.17g}',int(d[7])])
for a,thr,ylabel in zip(axes,[1e-5,1e-4],['归一化设计目标变化','发射协方差相对变化']):
 power=int(np.log10(thr))
 a.axhline(thr,color='0.25',ls=':',lw=1.2,label=rf'停止阈值 $10^{{{power}}}$')
 a.set_yscale('log');a.set_ylabel(ylabel);a.set_xlabel('WMMSE 外层迭代次数');a.set_facecolor('white');a.grid(alpha=.15)
 a.legend(frameon=False,fontsize=8,ncol=2)
 a.spines[['top','right']].set_visible(False)
fig.suptitle('WMMSE–ADMM 外层收敛诊断（阶段结果，尚未收敛）')
for ext in ['png','pdf']:
 fig.savefig(HERE/f'figure5_outer_log_thresholds.{ext}',dpi=300,facecolor='white')
with (HERE/'figure5_outer_log_thresholds.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f);w.writerow(['branch','outer','normalizedDesignObjectiveChange','transmitCovarianceChange','primalRatio','dualRatio','innerSolved']);w.writerows(rows)
print('FROZEN_BATCHES',BATCHES,'ROWS',len(rows))
