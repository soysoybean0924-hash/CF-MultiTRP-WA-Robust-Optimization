import re
import xml.etree.ElementTree as ET
from pathlib import Path

root=Path(__file__).resolve().parent
ET.register_namespace('', 'http://www.w3.org/2000/svg')
NS='{http://www.w3.org/2000/svg}'
for html in sorted(root.glob('archify/*/*.html')):
    raw=html.read_text(encoding='utf-8')
    m=re.search(r'<svg\b[^>]*>.*?</svg>',raw,re.S)
    if not m: continue
    svg=m.group(0).replace('<svg ','<svg xmlns="http://www.w3.org/2000/svg" ',1)
    tree=ET.fromstring(svg)
    styles=re.findall(r'<style[^>]*>(.*?)</style>',raw,re.S)
    css=next(s for s in styles if 'THEME VARIABLES' in s)
    light=re.search(r'\[data-theme="light"\]\s*\{(.*?)\}',css,re.S).group(1)
    vals=dict(re.findall(r'(--[a-zA-Z0-9-]+)\s*:\s*([^;]+);',light))
    rules={}
    for cls,body in re.findall(r'^\s*\.([ctam]-[a-zA-Z0-9_-]+)\s*\{([^{}]+)\}',css,re.M):
        decl={k.strip():v.strip() for k,v in re.findall(r'([a-z-]+)\s*:\s*([^;]+);',body)}
        for k,v in list(decl.items()):
            decl[k]=re.sub(r'var\((--[a-zA-Z0-9-]+)\)',lambda m:vals.get(m.group(1),'#777'),v)
        rules.setdefault(cls,{}).update(decl)
    allowed={'fill','stroke','opacity','stroke-width','stroke-dasharray','stroke-linecap','stroke-linejoin','fill-opacity','stroke-opacity'}
    for parent in tree.iter():
        for child in list(parent):
            if child.tag==NS+'pattern' and child.get('id')=='grid': parent.remove(child)
            elif child.tag==NS+'rect' and child.get('fill')=='url(#grid)': parent.remove(child)
        for cls in parent.get('class','').split():
            for k,v in rules.get(cls,{}).items():
                if k in allowed:
                    rgba=re.fullmatch(r'rgba\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*([0-9.]+)\s*\)',v)
                    if rgba and k in ('fill','stroke'):
                        parent.set(k,'#%02x%02x%02x'%tuple(int(rgba.group(i)) for i in (1,2,3)))
                        parent.set(k+'-opacity',rgba.group(4))
                    else:parent.set(k,v)
        if parent.tag in (NS+'text',NS+'tspan'):
            if not parent.get('fill'): parent.set('fill','#111827')
            parent.set('font-family','Microsoft YaHei, Noto Sans CJK SC, Arial, sans-serif')
    vb=[float(x) for x in tree.get('viewBox').split()]
    tree.set('width',str(vb[2]));tree.set('height',str(vb[3]))
    bg=ET.Element(NS+'rect',{'x':str(vb[0]),'y':str(vb[1]),'width':str(vb[2]),'height':str(vb[3]),'fill':'#ffffff'})
    tree.insert(0,bg)
    ET.ElementTree(tree).write(html.with_suffix('.svg'),encoding='utf-8',xml_declaration=True)
