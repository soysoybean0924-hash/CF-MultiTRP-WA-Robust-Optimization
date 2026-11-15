"""Build a separate five-figure briefing PDF from frozen batches only."""

from pathlib import Path
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph
from PIL import Image

HERE = Path(__file__).resolve().parent
PAIRED = HERE.parents[1] / "chasi_paired" / "figures" / "outer60"
OUT = HERE / "WMMSE_ADMM_汇报修订稿_20261025.pdf"
pdfmetrics.registerFont(TTFont("SimHei", r"C:\Windows\Fonts\simhei.ttf"))
PAGE_W, PAGE_H = landscape(A4)
MARGIN = 38
BODY = ParagraphStyle("body", fontName="SimHei", fontSize=10.2, leading=15.3,
                      textColor="#1F2937")
SMALL = ParagraphStyle("small", parent=BODY, fontSize=9.3, leading=14)

def paragraph(pdf, content, x, y_top, width, style=BODY):
    para = Paragraph(content, style)
    _, height = para.wrap(width, PAGE_H)
    para.drawOn(pdf, x, y_top - height)
    return y_top - height

def title(pdf, number, heading, status):
    pdf.setFillColorRGB(0.09, 0.17, 0.27)
    pdf.setFont("SimHei", 19)
    pdf.drawString(MARGIN, PAGE_H - 43, f"图{number}  {heading}")
    pdf.setFont("SimHei", 9.5)
    pdf.setFillColorRGB(0.35, 0.40, 0.46)
    pdf.drawRightString(PAGE_W - MARGIN, PAGE_H - 40, status)
    pdf.setStrokeColorRGB(0.83, 0.86, 0.89)
    pdf.line(MARGIN, PAGE_H - 53, PAGE_W - MARGIN, PAGE_H - 53)

def footer(pdf, n):
    pdf.setFont("SimHei", 8.5)
    pdf.setFillColorRGB(0.42, 0.46, 0.52)
    pdf.drawString(MARGIN, 20, "单场景受控配对实验；阶段结果，四组尚未满足连续五轮停止判据。")
    pdf.drawRightString(PAGE_W - MARGIN, 20, str(n))
    pdf.showPage()

def picture(pdf, path, x, top, max_w, max_h):
    with Image.open(path) as image:
        width, height = image.size
    scale = min(max_w / width, max_h / height)
    w, h = width * scale, height * scale
    pdf.drawImage(ImageReader(str(path)), x + (max_w - w) / 2, top - h,
                  width=w, height=h, mask="auto")
    return top - h

pdf = canvas.Canvas(str(OUT), pagesize=(PAGE_W, PAGE_H), pageCompression=1)
pdf.setTitle("WMMSE-ADMM 四组配对实验汇报修订稿")

title(pdf, 1, "当前四组 WMMSE–ADMM 修订框架图", "代码核验图；原图见第6页")
picture(pdf, HERE / "archify" / "paired-wmmse-admm" / "paired-wmmse-admm.png",
        MARGIN, PAGE_H - 73, PAGE_W - 2*MARGIN, 265)
y = PAGE_H - 360
y = paragraph(pdf, "<b>实际链路。</b>共同随机初值 X0；理想分支由 Htrue 设计，非理想分支由 Hhat 设计。"
              "外层更新接收矩阵 U、MSE 矩阵 E 与权重 W；固定 U/W 后，内层 ADMM 更新预编码 X、"
              "功率投影变量 Q 和对偶变量 Y。独立评价器只用共同 Htrue 计算真实指标。", MARGIN, y, PAGE_W-2*MARGIN)
y = paragraph(pdf, "<b>边界。</b>当前非理想分支没有把误差协方差 C_rx 输入求解器。"
              "两流、服务集及调度均固定；ADMM 是集中计算，图不表示 SRS 空口观测、动态 Rank 或 TRP 本地分布式求解。",
              MARGIN, y-11, PAGE_W-2*MARGIN)
paragraph(pdf, "<b>停止。</b>内层原始与对偶残差均达标，归一化设计目标变化≤10^-5，"
          "发射协方差相对变化≤10^-4，且三项连续五轮同时满足。20轮仅为保存和轮转节点。",
          MARGIN, y-55, PAGE_W-2*MARGIN)
footer(pdf, 1)

title(pdf, 2, "真实总频谱效率演化", "共同第60轮 · sum-SE 目标")
picture(pdf, PAIRED / "paired_formal_true_sumse.png", MARGIN, PAGE_H-70,
        PAGE_W-2*MARGIN, 415)
paragraph(pdf, "两条曲线分别为理想 CSI＋sum-SE、非理想 CSI＋sum-SE；横轴是 WMMSE 外层轮数，"
          "纵轴是独立评价器在共同 Htrue 上计算的 210 用户真实总 SE（bit/s/Hz）。"
          "共同初值的第0轮实际评价值为 60.103736；第60轮分别为 3163.645841 与 1984.363257。"
          "两组第60轮仅是同轮阶段比较，尚非收敛终点。", MARGIN, 86, PAGE_W-2*MARGIN, SMALL)
footer(pdf, 2)

title(pdf, 3, "真实 PF 效用演化", "共同第60轮 · PF 目标")
picture(pdf, PAIRED / "paired_formal_true_pf.png", MARGIN, PAGE_H-70,
        PAGE_W-2*MARGIN, 415)
paragraph(pdf, "两条曲线分别为理想 CSI＋PF、非理想 CSI＋PF；评价信道仍为共同 Htrue。"
          "PF=Σk ln[(max(rk,0)+10^-6)/(1 bit/s/Hz)]，rk 先按 RB 数跨17组聚合。"
          "第0轮 PF 为 -275.682697；第60轮分别为 559.284267 与 478.546489。"
          "负初值来自低速率用户的对数项，并非负频谱效率。", MARGIN, 86, PAGE_W-2*MARGIN, SMALL)
footer(pdf, 3)

title(pdf, 4, "四种方案的用户 SE 经验 CDF", "共同第60轮 · 单场景")
picture(pdf, PAIRED / "paired_formal_user_cdf.png", MARGIN, PAGE_H-70,
        PAGE_W-2*MARGIN, 415)
paragraph(pdf, "四组各取固定场景第60轮的 210 名 UE，合计四条各210样本的曲线；"
          "没有把迭代轮次或资源组当作独立场景。全部用户均调度，零 SE 用户保留："
          "sum-SE 两支各3名，PF 两支各0名。每用户 SE=Σg(RBg/273)SEkg。"
          "横轴始终为每用户 SE（bit/s/Hz），PF 只是优化目标。", MARGIN, 86, PAGE_W-2*MARGIN, SMALL)
footer(pdf, 4)

title(pdf, 5, "外层停止诊断（对数纵轴）", "各组最新完整批次")
picture(pdf, HERE / "figure5_outer_log_thresholds.png", MARGIN, PAGE_H-70,
        PAGE_W-2*MARGIN, 415)
paragraph(pdf, "理想 CSI＋sum-SE 到220轮，非理想 CSI＋sum-SE 到200轮，两组 PF 到180轮；曲线保留各自实际长度。"
          "上图阈值为归一化设计目标变化 10^-5，下图阈值为发射协方差相对变化 10^-4。"
          "两项末值仍高于阈值；即使内层 ADMM 达标，也须连续五轮三项同时满足才可宣布外层停止。"
          "图2–4为共同第60轮阶段性能，图5为各组最新完整批次诊断，轮数不同。",
          MARGIN, 86, PAGE_W-2*MARGIN, SMALL)
footer(pdf, 5)

title(pdf, 6, "原始 Archify 框架图与适用范围", "桌面原图；非当前四组框架")
picture(pdf, HERE / "archify" / "original_05_wmmse_admm.png", MARGIN,
        PAGE_H-86, PAGE_W-2*MARGIN, 320)
paragraph(pdf, "<b>原图实际嵌入：</b>此页直接使用桌面 05_wmmse_admm.png 的原始字节。"
          "图中 H_est、C_rx、SVD 初始化、ADMM＋CG 与统计误差项更接近茶思主仓库原生 CHASI 链；"
          "它不对应正在运行的四组配对入口。", MARGIN, 157, PAGE_W-2*MARGIN, SMALL)
paragraph(pdf, "<b>关键差异：</b>四组采用共同随机 X0、非理想分支只用 Hhat、没有 C_rx 输入，"
          "内层使用集中式 ADMM 与直接线性求解；停止还要求两项外层变化连续五轮达标。"
          "第1页因此嵌入专为四组入口核验的修订图；各框和箭头的逐项证据见配套审计报告。",
          MARGIN, 98, PAGE_W-2*MARGIN, SMALL)
footer(pdf, 6)

pdf.save()
print(OUT, OUT.stat().st_size)
