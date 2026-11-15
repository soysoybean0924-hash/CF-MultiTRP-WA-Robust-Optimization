"""Chinese CDF and paired final metrics from independent evaluation JSON."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager

LABEL={'ideal':'理想 CSI','nominal':'名义 CSI','robust':'误差感知'}


def main():
    p=argparse.ArgumentParser(); p.add_argument('evaluation'); p.add_argument('output_dir')
    p.add_argument('--manifest')
    p.add_argument('--font-file',required=True,help='Chinese TrueType/OpenType font')
    a=p.parse_args(); data=json.loads(Path(a.evaluation).read_text(encoding='utf-8'))
    font_path=Path(a.font_file)
    if not font_path.is_file(): p.error(f'Chinese font not found: {font_path}')
    font=font_manager.FontProperties(fname=str(font_path))
    from matplotlib.ft2font import FT2Font
    cmap=FT2Font(str(font_path)).get_charmap()
    if any(ord(c) not in cmap for c in '每用户真实速率功率误差收敛'):
        p.error(f'font lacks required Chinese glyphs: {font_path}')
    font_manager.fontManager.addfont(str(font_path))
    plt.rcParams['font.family']=font.get_name()
    out=Path(a.output_dir)
    if out.exists(): raise FileExistsError(out)
    out.mkdir(parents=True)
    for goal in ('sumse','pf'):
        fig,ax=plt.subplots(figsize=(7,5)); plotted={}
        for method in LABEL:
            key=f'{method}_{goal}'; row=data.get(key,{})
            if 'final' not in row: continue
            x=np.sort(np.asarray(row['final']['perUserTrueSE'],float))
            plotted[key]=x.tolist()
            if len(x):
                ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',
                        label=f'{LABEL[method]} (5%={np.percentile(x,5):.3g}, 近零={np.mean(x<1e-3):.1%})')
        ax.set(xlabel='每用户真实 SE (bit/s/Hz)',ylabel='经验 CDF',ylim=(0,1),
               title=f'{goal.upper()} 优化目标：真实信道评价')
        ax.grid(True); ax.legend()
        fig.tight_layout()
        for ext in ('png','pdf','svg'): fig.savefig(out/f'真实SE_CDF_{goal}.{ext}')
        plt.close(fig)
        (out/f'真实SE_CDF_{goal}.json').write_text(
            json.dumps(dict(curves=plotted,missing=[f'{m}_{goal}' for m in LABEL
                        if f'{m}_{goal}' not in plotted]),ensure_ascii=False,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(2,1,figsize=(9,8))
    for key,row in data.items():
        history=row.get('history',[])
        xs=[]; primal=[]; dual=[]; ep=[]; ed=[]; powers=[]
        for h in history:
            powers.append(max(h.get('powerByTrp',[float('nan')])))
            trace=h.get('admmTrace',[])
            if trace:
                last=trace[-1]; xs.append(h['outer'])
                primal.append(last.get('primal_residual',float('nan')))
                dual.append(last.get('dual_residual',float('nan')))
                ep.append(last.get('eps_primal',float('nan')))
                ed.append(last.get('eps_dual',float('nan')))
        axes[0].plot(xs,primal,label=key+' primal'); axes[0].plot(xs,dual,label=key+' dual')
        axes[0].plot(xs,ep,linestyle=':',label=key+' tol p')
        axes[0].plot(xs,ed,linestyle=':',label=key+' tol d')
        axes[1].plot([h['outer'] for h in history],powers,label=key)
    axes[0].set_title('ADMM 残差与容差'); axes[0].set_yscale('log'); axes[1].set_title('最大发射功率 (W)')
    for ax in axes: ax.grid(True);ax.legend(fontsize=6)
    fig.tight_layout()
    for ext in ('png','pdf','svg'): fig.savefig(out/f'ADMM与功率.{ext}')
    plt.close(fig)
    if a.manifest:
        manifest=json.loads(Path(a.manifest).read_text(encoding='utf-8'))
        ages=[r['ageMs'] for r in manifest['resources']]
        nmse=[r['sampleNmse'] for r in manifest['resources']]
        theory=[r['theoryNmse'] for r in manifest['resources']]
        fig,ax=plt.subplots(figsize=(7,5)); ax.scatter(ages,nmse,label='样本 NMSE')
        ax.scatter(ages,theory,label='理论协方差/样本真信道能量')
        ax.set(xlabel='SRS 测量年龄 (ms)',ylabel='样本 NMSE',title='年龄与估计误差')
        ax.grid(True);ax.legend();fig.tight_layout()
        for ext in ('png','pdf','svg'): fig.savefig(out/f'年龄_NMSE.{ext}')
        plt.close(fig)
    for goal in ('sumse','pf'):
        fig,axes=plt.subplots(2,2,figsize=(11,8))
        for method in LABEL:
            key=f'{method}_{goal}'; history=data.get(key,{}).get('history',[])
            if not history: continue
            t=[row['outer'] for row in history]
            for ax,field,title in zip(axes.flat,
                ('designedObjective','trueSumSE','outerRelativeChange','covarianceRelativeChange'),
                ('设计目标','真实总 SE','外层目标相对变化','发射协方差相对变化')):
                ax.plot(t,[row.get(field,float('nan')) for row in history],label=LABEL[method]);ax.set_title(title)
        for ax in axes.flat: ax.grid(True); ax.legend()
        axes[1,0].set_yscale('log'); axes[1,1].set_yscale('log')
        fig.tight_layout()
        for ext in ('png','pdf','svg'): fig.savefig(out/f'目标与停止_{goal}.{ext}')
        plt.close(fig)
    keys=[k for k,v in data.items() if 'final' in v]
    for metric in ('trueSumSE','truePF','trueRate5'):
        fig,ax=plt.subplots(figsize=(8,5))
        ax.plot(keys,[data[k]['final'][metric] for k in keys],marker='o',label=metric)
        ax.tick_params(axis='x',rotation=45); ax.legend();ax.grid(True);fig.tight_layout()
        for ext in ('png','pdf','svg'): fig.savefig(out/f'{metric}.{ext}')
        plt.close(fig)
    (out/'原始评价数据.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__': main()
