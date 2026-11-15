"""Export standalone white-background SVGs from validated Archify HTML."""
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", NS)
style = {
    "c-grid": {"stroke": "none", "fill": "none"},
    "c-mask": {"fill": "#ffffff", "stroke": "none"},
    "c-lane": {"fill": "#ffffff", "stroke": "#d3dce5"},
    "c-backend": {"fill": "#e9f7f1", "stroke": "#09815d"},
    "c-database": {"fill": "#f3efff", "stroke": "#7651bc"},
    "c-security": {"fill": "#fff0f2", "stroke": "#bb3351"},
    "t-primary": {"fill": "#17212c"},
    "t-muted": {"fill": "#465767"},
    "a-default": {"fill": "none", "stroke": "#56697c"},
    "m-default": {"fill": "#56697c"},
    "m-dashed": {"fill": "#7651bc"},
    "m-emphasis": {"fill": "#09815d"},
    "m-security": {"fill": "#bb3351"},
    "s-backend": {"color": "#09815d", "stroke": "#09815d"},
    "s-database": {"color": "#7651bc", "stroke": "#7651bc"},
    "s-security": {"color": "#bb3351", "stroke": "#bb3351"},
    "semantic-sigil": {"fill": "none", "stroke-width": "1.35"},
}
for name in ("current_execution", "current_algorithm"):
    html = (ROOT / name / f"{name}.html").read_text(encoding="utf-8")
    svg = html[html.index("<svg viewBox="):html.index("</svg>", html.index("<svg viewBox=")) + 6]
    svg = svg.replace("<svg ", '<svg xmlns="http://www.w3.org/2000/svg" ', 1)
    root = ET.fromstring(svg)
    root.set("fill", "#17212c")
    root.set("font-family", "Microsoft YaHei, Noto Sans CJK SC, sans-serif")
    for element in list(root.iter()):
        if element.tag == f"{{{NS}}}style":
            root.remove(element)
            continue
        for class_name in element.get("class", "").split():
            for key, value in style.get(class_name, {}).items():
                element.set(key, value)
    # The white rectangle is part of the export, not the interactive viewer.
    vb = root.get("viewBox").split()
    background = ET.Element(f"{{{NS}}}rect", {"x": "0", "y": "0", "width": vb[2], "height": vb[3], "fill": "#ffffff"})
    root.insert(0, background)
    # The browser's interactive legend is not included in the SVG export.
    # Keep an explicit Chinese legend on the standalone PPT/PDF artwork.
    legend_y = float(vb[3]) - 29
    for x, color, label in ((48, "#09815d", "计算与更新"),
                            (190, "#bb3351", "分支与停止判断"),
                            (370, "#7651bc", "保存与结果")):
        ET.SubElement(root, f"{{{NS}}}rect", {
            "x": str(x), "y": str(legend_y), "width": "13", "height": "13",
            "rx": "3", "fill": "#ffffff", "stroke": color, "stroke-width": "1.4"})
        label_node = ET.SubElement(root, f"{{{NS}}}text", {
            "x": str(x + 19), "y": str(legend_y + 11), "fill": "#465767",
            "font-size": "12", "font-weight": "500"})
        label_node.text = label
    output = ROOT / name / f"{name}.svg"
    ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)
    print(output)
