# -*- coding: utf-8 -*-
"""45 湿实验方案 docx 生成器"""
from docx import Document
from docx.shared import Pt

doc = Document()
style = doc.styles["Normal"]
style.font.name = "Times New Roman"
style.font.size = Pt(11)

doc.add_heading("湿实验验证方案", level=0)
p = doc.add_paragraph("PE × SIGnature 论文功能验证 | 设计日期：2026-09-29")
p.runs[0].italic = True

doc.add_heading("一、背景与目标", level=1)
doc.add_paragraph(
    '本研究通过单细胞基础模型归因分析，在 7 个公共队列、464,727 个胎盘/蜕膜细胞中发现：'
    '(1) 髓系补体-固有免疫模块（C1QA/C1QC/C5AR1/CTSS/LST1/TYROBP）在胎盘型 PE（早发/重度）中上调，'
    '而在母体型 PE（晚发/混合）中反向；(2) 该两分子面架构与 Admati et al. (Med 2023) 的分子分型独立吻合；'
    '(3) 低表达转录因子 GATA2、NR2F2（EVT 分化与蜕膜化核心调控因子）被归因分析系统性救回。'
    '本湿实验方案旨在通过功能实验验证上述计算发现，将两分子面从计算观察升级为实验证据。'
)

doc.add_heading("二、方案 A：PE 血清诱导实验（首选）", level=1)
doc.add_paragraph("科学逻辑：SIGnature 原文用川崎病血清诱导单核细胞验证 MS1 模块。"
                  "我们用 PE 血清做同样实验，验证补体-髓系模块。")

table = doc.add_table(rows=5, cols=3)
table.style = "Table Grid"
data = [
    ("分组", "血清来源", "样本量"),
    ("EOPE（胎盘型）", "早发型 PE（<34 周，术前采血）", "5–8 例"),
    ("Severe PE", "重度 PE（≥34 周，伴严重特征）", "5–8 例"),
    ("LOPE（母体型）", "晚发型 PE（≥37 周）", "5–8 例"),
    ("正常对照", "同孕周正常妊娠（选择性剖宫产）", "5–8 例"),
]
for i, row in enumerate(data):
    for j, v in enumerate(row):
        table.rows[i].cells[j].text = v

doc.add_paragraph("步骤：THP-1 培养 → 10% PE/对照血清处理 24h → qPCR 检测 C1QA, C1QC, C5AR1, CTSS, LST1, TYROBP（6 基因 panel）")
doc.add_paragraph("预期：EOPE/severe 血清诱导补体基因上调；LOPE/对照不诱导——功能性验证两分子面")
doc.add_paragraph("成本：~¥8,000–12,000 | 周期：~2 个月 | 伦理：血清标本，低风险")

doc.add_heading("三、方案 B：胎盘组织免疫组化", level=1)
doc.add_paragraph("多重免疫荧光：C1q（补体）+ CD68（巨噬细胞）+ CK7（滋养层）；IHC：C5aR1、FLT1")
doc.add_paragraph("分组：EOPE 5 例 + LOPE 5 例 vs 正常对照 5 例（石蜡切片）")
doc.add_paragraph("预期：EOPE 胎盘 C1q 沉积↑ + CD68 浸润↑（共定位验证）")
doc.add_paragraph("成本：~¥5,000–8,000 | 周期：1–1.5 个月")

doc.add_heading("四、方案 C：滋养层 C5a 功能实验", level=1)
doc.add_paragraph("HTR-8/SVneo 划痕/Transwell：对照 / rC5a / rC5a+C5aR拮抗剂（PMX53）")
doc.add_paragraph("读出：迁移率、管状形成、GATA2/NR2F2 qPCR")
doc.add_paragraph("预期：C5a 抑制迁移→拮抗拯救→补体-迁移轴功能证据")
doc.add_paragraph("成本：~¥15,000–20,000 | 周期：2–3 个月")

doc.add_heading("五、方案 D：血清 ELISA 验证队列", level=1)
doc.add_paragraph("EOPE/LOPE/正常各 20 例存留血清，ELISA 检测 C3a、C5a、sC5b-9")
doc.add_paragraph("分析：组间比较 + 与补体模块相关系数做相关")
doc.add_paragraph("成本：~¥10,000–15,000 | 周期：3–4 周（最快启动）")

doc.add_heading("六、推荐组合", level=1)
table2 = doc.add_table(rows=4, cols=5)
table2.style = "Table Grid"
combos = [
    ("预算", "方案", "总成本", "总周期", "Revision 说服力"),
    ("最低", "D（ELISA）", "¥10k", "1 个月", "★★"),
    ("推荐", "A + B（血清诱导 + IHC）", "¥15–20k", "2 个月", "★★★★"),
    ("充裕", "A + B + C", "¥30–40k", "3 个月", "★★★★★"),
]
for i, row in enumerate(combos):
    for j, v in enumerate(row):
        table2.rows[i].cells[j].text = v

doc.add_heading("七、执行优先级", level=1)
for s in [
    "第 1 周：提交 IRB 申请（方案 A + D；可同时申请 B 的胎盘组织）",
    "第 1–2 周：订购 THP-1、qPCR 引物、ELISA kit",
    "第 2–8 周：样本收集",
    "第 8–10 周：方案 A（血清诱导 + qPCR）",
    "第 8–12 周：方案 B（免疫组化，与 A 并行）",
    "第 10–14 周：方案 D（ELISA，可回顾性不等新样本）",
    "方案 C 留给 major revision 二轮或后续课题",
]:
    doc.add_paragraph(s)

doc.add_heading("八、湿实验结果与论文衔接", level=1)
table3 = doc.add_table(rows=5, cols=3)
table3.style = "Table Grid"
links = [
    ("湿实验结果", "论文对应位置", "论文增量"),
    ("A：EOPE/severe 血清诱导补体基因", "R4 亚型段 + Discussion", "两分子面从计算→功能"),
    ("B：C1q 沉积 + CD68 浸润共定位", "R3 共识模块段", "补体模块组织学验证"),
    ("C：C5a 抑制迁移 + GATA2 变化", "R5 TF 捞回段", "机制深度"),
    ("D：血清补体蛋白与模块相关", "R6 bulk 验证段", "临床可及性"),
]
for i, row in enumerate(links):
    for j, v in enumerate(row):
        table3.rows[i].cells[j].text = v

doc.add_heading("九、伦理与样本采集", level=1)
for s in [
    "• 所有方案需广医三院医学伦理委员会审批",
    "• 血清：常规剖宫产术前多采 5 mL 外周血（不额外穿刺），患者知情同意",
    "• 胎盘组织：手术废弃标本（分娩后 30 分钟内取材），患者知情同意",
    "• 样本编码匿名化，临床数据按 CRF 表收集（分娩孕周、血压、尿蛋白、出生体重）",
    "• 排除标准：多胎妊娠、宫内感染、自身免疫病、慢性肾病、既往 CVD",
]:
    doc.add_paragraph(s)

doc.save(r"F:/pe_signature/submission_GB/WetLab_Protocol.docx")
print("saved: WetLab_Protocol.docx")
