"""Run a selected project Python branch with only this project's module paths."""
import argparse,runpy,sys
from pathlib import Path
def main():
 p=argparse.ArgumentParser();p.add_argument('branch',choices=['basic','layered'])
 p.add_argument('script',choices=['campaign.py','layered.py','preflight.py','plot_layered.py'])
 a=p.parse_args(sys.argv[1:3]);rest=sys.argv[3:];root=Path(__file__).resolve().parents[1]
 if a.branch=='basic':
  if a.script!='campaign.py':p.error('basic branch runs campaign.py')
  paths=[root/'04_研究方法_WA运行实现'/'python',root/'03_研究场景构造'/'python']
 else:
  paths=[root/'05_研究创新点_鲁棒施加'/'python']
 script=paths[0]/a.script
 if not script.is_file():p.error('script missing from current project')
 sys.path[:0]=[str(x) for x in paths]
 sys.argv=[str(script),*rest]
 runpy.run_path(str(script),run_name='__main__')
if __name__=='__main__':main()

