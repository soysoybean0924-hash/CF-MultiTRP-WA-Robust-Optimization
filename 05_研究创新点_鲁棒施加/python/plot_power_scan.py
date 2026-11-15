"""Diagnostic power-cap comparison from immutable server state files."""
import argparse,csv,json,os
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    a.out.mkdir(parents=True)
    font=Path(os.environ['SRS_CJK_FONT']);assert font.is_file()
    font_manager.fontManager.addfont(str(font))
    plt.rcParams['font.family']=font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams['axes.unicode_minus']=False
    rows=[]
    for cap in (10,15,20):
        for mode in ('uniform','variable','fixed','mapped'):
            if cap==20:
                branch='full17_stage2_5bc1' if mode in ('uniform','variable') else 'full17_stage3_5bc1_pJ7PNd'
                state=a.results/branch/mode/'state.json'
            else:state=a.results/'full17_remaining_5bc1'/f'power_{cap}w_{mode}'/'state.json'
            d=json.loads(state.read_text(encoding='utf-8'))
            assert d['status']=='converged'
            h=d['history'][-1];m=d['mechanism']
            rows.append(dict(cap_w=cap,mode=mode,iteration=d['iter'],
                effective_cap_w=m['effective_cap_w'],actual_total_tx_power_w=sum(h['powerByTrp']),
                true_weighted_sum_se=h['trueSumSE'],input_manifest_sha256=d['identity']['input_manifest_sha256'],
                source_state=str(state)))
    with (a.out/'metrics.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,rows[0]);w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(figsize=(7,4.5))
    for mode in ('uniform','variable','fixed','mapped'):
        subset=[x for x in rows if x['mode']==mode]
        ax.plot([x['cap_w'] for x in subset],[x['true_weighted_sum_se'] for x in subset],'-o',label=mode)
    ax.set(xlabel='Physical TRP cap (W)',ylabel='Resource-weighted true sum SE (bit/s/Hz)',
        title='Single frozen input: diagnostic power-cap scan')
    ax.grid(alpha=.3);ax.legend();fig.tight_layout();fig.savefig(a.out/'power_cap_diagnostic.png',dpi=180);plt.close(fig)
    (a.out/'STATUS.json').write_text(json.dumps(dict(formal_result=False,
        kind='power_cap',input_count=len(rows),
        reason='Single frozen scene; actual transmit power and effective cap are separate CSV fields',
        sources=[x['source_state'] for x in rows]),indent=2),encoding='utf-8')
if __name__=='__main__':main()
