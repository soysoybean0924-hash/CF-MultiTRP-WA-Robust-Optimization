"""Server-only diagnostic figures from actual layered.py state.json files."""
import argparse
import csv
import json
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

COLORS={'uniform':'#1769aa','variable':'#d95f02','fixed':'#39a852','mapped':'#8e44ad'}
def font():
    path=os.environ.get('SRS_CJK_FONT')
    if not path or not Path(path).is_file():raise RuntimeError('set SRS_CJK_FONT to the uploaded NotoSansSC.otf')
    font_manager.fontManager.addfont(path)
    f=font_manager.FontProperties(fname=path).get_name()
    if any(font_manager.get_font(path).get_char_index(ord(x))==0 for x in '统一变量鲁棒保护功率'):
        raise RuntimeError('font lacks Chinese glyphs')
    plt.rcParams.update({'font.family':f,'axes.unicode_minus':False,'svg.fonttype':'none','pdf.fonttype':42})
def read(path):
    obj=json.loads(Path(path).read_text(encoding='utf-8'))
    if not obj.get('history') or not obj.get('identity'):raise ValueError(f'bad state: {path}')
    return obj
def save(fig,out,stem):
    for ext in ('png','svg','pdf'):fig.savefig(out/f'{stem}.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
def measures(state,threshold):
    ue=np.asarray(state['history'][-1]['perUserTrueSE'],float)
    return dict(p5=float(np.percentile(ue,5)),violation_proxy=float(np.mean(ue<threshold)),
                power=float(np.sum(state['history'][-1]['powerByTrp'])),
                zero=int(np.count_nonzero(ue==0)),status=state['status'],
                strength_beam=state['mechanism']['strengths']['beam'],
                strength_service=state['mechanism']['strengths']['service'],
                strength_power=state['mechanism']['strengths']['power'])
def main():
    p=argparse.ArgumentParser();p.add_argument('--kind',choices=['layer2','layer3','strength','ablation','trigger'],required=True)
    p.add_argument('--states',nargs='+',required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--threshold',type=float,default=0.01);a=p.parse_args();font()
    states=[read(x) for x in a.states];a.out.mkdir(parents=True,exist_ok=False)
    ids=[x['identity'] for x in states]
    if len({x['input_manifest_sha256'] for x in ids})!=1 or len({x['layer1_sha256'] for x in ids})!=1:
        raise ValueError('states use different frozen input or Layer 1 output')
    data=[measures(x,a.threshold) for x in states]
    names=[x['identity']['mode']+('_'+x['identity']['ablation'] if x['identity']['ablation']!='none' else '')+
           ('_'+x['identity']['trigger'] if a.kind=='trigger' else '') for x in states]
    with (a.out/'metrics.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['name',*data[0].keys()]);w.writeheader()
        for n,d in zip(names,data):w.writerow(dict(name=n,**d))
    eligible=all(x['status']=='converged' for x in states)
    for lang in ('zh','en'):
        fig,axes=plt.subplots(1,2,figsize=(11,4.6));colors=[COLORS.get(x['identity']['mode'],'#888888') for x in states]
        axes[0].bar(names,[x['p5'] for x in data],color=colors);axes[0].set_ylim(bottom=0)
        axes[0].set_title('全体 UE 5% SE' if lang=='zh' else 'All-UE 5th percentile SE')
        axes[0].set_ylabel('bit/s/Hz')
        axes[1].bar(names,[x['violation_proxy'] for x in data],color=colors);axes[1].set_ylim(0,1)
        axes[1].set_title('SE 阈值低于率（代理）' if lang=='zh' else 'SE below threshold (proxy)')
        axes[1].set_ylabel('fraction')
        for ax in axes:ax.tick_params(axis='x',labelrotation=20);ax.grid(axis='y',alpha=.2)
        fig.suptitle(('机制诊断｜' if lang=='zh' else 'Mechanism diagnostic | ')+a.kind)
        fig.text(.5,.01,('单场景模型代理；不能称正式体验违约率。' if lang=='zh' else
                         'Single-scene model proxy; not formal experienced-rate violation.'),ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.04,1,.95));save(fig,a.out,f'{a.kind}_diagnostic_{lang}')
        if a.kind=='layer3':
            fig,ax=plt.subplots(figsize=(8,4.5));xx=np.arange(len(states));width=.25
            for j,key in enumerate(('beam','service','power')):
                ax.bar(xx+(j-1)*width,[d['strength_'+key] for d in data],width,label=key)
            ax.set_xticks(xx,names);ax.set_ylabel('dimensionless strength');ax.legend()
            ax.set_title('实际保护强度' if lang=='zh' else 'Applied protection strength')
            fig.tight_layout();save(fig,a.out,f'layer3_strength_{lang}')
    (a.out/'STATUS.json').write_text(json.dumps(dict(kind=a.kind,formal_result=False,
      reason='Independent test seeds and strict experience constraints are absent; diagnostic only',
      all_converged=eligible,states=a.states),indent=2),encoding='utf-8')
    print(json.dumps(dict(kind=a.kind,all_converged=eligible,output=str(a.out),formal_result=False)))
if __name__=='__main__':main()
