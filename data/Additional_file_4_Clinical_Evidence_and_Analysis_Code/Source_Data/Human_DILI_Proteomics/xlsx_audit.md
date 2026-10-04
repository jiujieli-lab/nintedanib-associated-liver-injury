# `fresh_data_dili_source.xlsx` 只读结构审计

审计对象：公开论文 Figure source-data 工作簿 `fresh_data_dili_source.xlsx`。本报告仅依据该工作簿本身及单元格间的一致性检查，不使用任何旧稿结果，也未改动源文件。

## 一、全局结论

- 工作簿含 10 个可见工作表；无隐藏行列、合并单元格、公式、命名区域、Excel Table、批注或超链接。
- 所有看似数值的强度、ALT 和 score 单元格均实际保存为 **Excel 文本**，不能直接按数值运算。提取时须保留原字符串，再经严格解析生成数值列；`ND`、空白和带 `*` 的值不得直接强制转为 0。
- 这是面向作图的个体观测值集合，不是完整的受试者级数据库：没有 subject ID、药物、年龄、性别、病因、采样时间、批次或其他协变量。除 Figure 5 / Supp Fig 7 已在同一行给出的 biomarker–ALT 对以外，不能仅凭 Excel 行号认定跨蛋白或随访配对。
- 组码按来源语境应为 HV（healthy volunteers）、DO（DILI onset）、DF（DILI follow-up）、NDO（non-DILI onset）、NDF（non-DILI follow-up）和 NAFLD；但工作簿内无 codebook，正式分析应在论文方法/数据字典中再次核实。
- `Figure 3`、`Supp Fig 4`、`Supp Fig 5` 大量重复同一批 onset 数据；它们是同一数据的不同作图视图，不是独立验证队列。`Figure 6` 的 HV 30 行也在两个表中重复。

## 二、逐表审计

| Sheet | 实际非空范围与结构 | 观测层级、表头和组列含义 | 可安全提取的统计字段 | 主要歧义/风险 |
|---|---|---|---|---|
| `Figure 2C` | `A1:G203`；12 个由空行分隔的 16 行区块。蛋白：ACO1、ALDOB、ASS1、CES1、CPS1、DMGDH、FBP1、FAH、GSTA1、HPD、OTC、PCK2。每块第 6 行式表头，B:G 为 HV10、DO10、DF10、NDO5、NDF5、NAFLD10。 | Discovery cohort 的单蛋白个体观测值；测量尺度明确写为 `VSN Normalized Relative Intensity`。短组列末端空白仅为宽表填充。 | protein、group、within-group position、VSN intensity、source cell。每蛋白 50 个数值，共 600 个。 | 所有数值为文本；无 subject ID。**ALDOB 与 HPD 的全部 50 个值逐格完全相同**，高度疑似源文件复制/标签问题，二者不可同时作为独立信号，须回查原始质谱文件或作者数据说明。 |
| `Figure 3` | `A1:D1054`；12 个蛋白区块：ACO1、ALDOB、ASS1、CPS1、DMGDH、FAH、GSTA1、HPD、LECT2、OTC、CES1、FBP1。B:D 表头分别声明 HV n=60、DO n=82、NDO n=34。 | Confirmatory cohort onset 的单蛋白个体值；宽表按组纵向堆叠。工作表未直接写单位/是否为 light:heavy relative intensity。 | protein、group、value、detected/ND、within-group position。可作组间分布、效应量、检测率分析。 | DMGDH 的 HV 列有 61 个数值（`B357:B417`），与 n=60 不符；CES1 的 DO 序列内有一个无标签空白 `C936`，而非末端填充。HPD：HV 7 个 ND；OTC：HV 15、DO 1 个 ND；CES1：HV 23、DO 5、NDO 2 个 ND；FBP1：HV 23、DO 6、NDO 2 个 ND。无 ID，不能跨蛋白直接拼接。 |
| `Figure 5` | `A1:D1191`；8 个 148 行区块：GLDH、CK18、ASS1、ACO1、ALDOB、FAH、CPS1、FBP1。每块 142 条记录，HV60 + DO82。 | 同一记录内已配对的 confirmatory biomarker 与 ALT，均为 log2 尺度；GLDH/CK18 是 U/L，候选蛋白标注为 light/heavy relative intensity。 | group、log2_ALT、log2_biomarker、protein；适合组内/总体 Spearman、稳健回归及 DILI-onset 与 HV 的效应比较。 | GLDH 有 2 个 ND；FBP1 有 29 个 ND。CK18 下限 100 U/L，`log2(100)=6.643856...` 的重复值属于左删失/下限值，不宜当普通连续值。FAH/CPS1 区块列顺序为 biomarker–Type–ALT，其他多数区块为 group–ALT–biomarker，必须按表头解析，不能固定列号。无临床协变量。 |
| `Figure 6` | `A1:C91`（E1 为说明）及 `A93:C153`（E93 为说明）两张长表。表 1：90 行；表 2：60 行。 | 字段为 `SigName`（Zone 1/2/3）、`scores`、`Group`。表 1 各 zone×HV/DO/DF 均 10 行；表 2 各 zone×HV 10 行、zone×NDO/NDF 各 5 行。 | zone、score、group；可作每组每 zone 的描述统计和可视化。 | 每位样本看起来贡献 3 个 zone score，且 DO–DF/NDO–NDF可能重复测量，但无 subject ID；90/60 行不能当独立受试者做普通检验。表 2 的 HV 30 行与表 1 HV 30 行完全重复。灰色区域说明是标准误，但表中未提供独立 SE 列。建议仅作描述性/探索性证据，除非获得映射表。 |
| `Supp Fig 1-LECT2` | `A1:G16`；单一 LECT2 discovery 区块；组列与 Figure 2C 相同。 | Discovery cohort 的 LECT2 个体 VSN intensity。 | group、LECT2 intensity、ND indicator、source cell。含 ND 后各组总数与表头 n 一致。 | HV/DO/DF/NAFLD 各 2 个 ND，NDO/NDF 各 1 个 ND；末端空白为 n=5 组的结构填充。该蛋白与 Figure 2C 分表存放，合并时需统一尺度和缺失编码。 |
| `Supp Fig 2` | `A1:E886`；10 个区块：ACO1、ALDOB、ASS1、CPS1、DMGDH、FAH、GSTA1、HPD、LECT2、OTC。B:E 为 DO、DF、NDO、NDF，表头不写 n。 | Confirmatory onset/follow-up 的单蛋白个体值，宽表。理论上可描述四组分布及随访变化。 | protein、group/time-state、value、ND flag、within-column position。 | 无 subject ID，故不能确认 DO–DF 或 NDO–NDF 的逐例配对。观察到的非空数随蛋白变化：常见为 DO82、DF83、NDO34、NDF21；CPS1 为 DO83（含 1 ND）、DF84、NDF20；DMGDH 为 DO83（含 1 ND）、DF84、NDF21；HPD 为 DO83（含 1 ND）、DF84（含 7 ND），且 **`E627` 在应为数值的 NDF 列中错误写成字符串 `NDF`**；OTC 为 DO82（1 ND）、DF83（11 ND）、NDF21（3 ND）。这些计数须与样本清单核对。 |
| `Supp Fig 3-PCK2` | `A1:F88`；PCK2 单一区块；B:F 为 HV、DO、DF、NDO、NDF，表头不写 n。 | Confirmatory PCK2 个体值。 | group、PCK2 value、ND flag、asterisk flag、within-column position。 | 观察数为 HV60（54 数值+6 ND）、DO82、DF83、NDO34、NDF23。`C28="3077.0266*"` 是未解释的星号/极端值，工作簿无脚注；必须保留 flag 并做含/不含该值的敏感性分析。单位/尺度未在 sheet 中说明。 |
| `Supp Fig 4` | `A1:C1042`；12 个区块，蛋白和数值与 Figure 3 对应；B/C 为 HV/DO。 | Figure 3 的 HV–DO 作图子集。 | 仅可作为 Figure 3 提取结果的复核。 | 对 12 个蛋白，HV 和 DO 序列均与 Figure 3 **逐值完全一致**，不是独立数据。DMGDH 多出的第 61 个 HV 值和 CES1 DO 的内部空白也被原样复制。 |
| `Supp Fig 5` | `A1:C1042`；12 个区块，蛋白和数值与 Figure 3 对应；B/C 为 DO/NDO。 | Figure 3 的 DO–NDO 作图子集。 | 仅可作为 Figure 3 提取结果的复核。 | 对 12 个蛋白，DO 和 NDO 序列均与 Figure 3 **逐值完全一致**，不是独立验证。OTC、CES1、FBP1 的 ND 编码与 Figure 3 相同。 |
| `Supp Fig 7` | `A1:D1035`；7 个长表区块：HPD `A1:D147`、OTC `A149:D295`、GSTA1 `A297:D443`、DMGDH `A445:D591`、CES1 `A593:D740`、LECT2 `A741:D887`、PCK2 `A889:D1035`。每块 142 行，HV60+DO82。 | 每行给出 protein、Type(HV/DO)、log2 biomarker 和 log2 ALT；属于同一记录内的成对图源。 | protein、group、log2_biomarker、log2_ALT、ND flag；可做稳健相关/回归。 | CES1 保留 29 个 ND。与原始值顺序核对显示：HPD 的 7 个 ND 和 OTC 的 16 个 ND 在此表被替换为 `-14.287712...`（对应原尺度 0.00005）；PCK2 的 6 个 ND 被替换为 `-2.867752...`（对应 0.137），但 sheet 未说明替代规则。PCK2 数值与 Supp Fig 3 的四位小数值仅存在舍入差。不能将这些区块当成新的队列。 |

## 三、跨表一致性与独立性

1. `Figure 3` 是 confirmatory onset 数据的首选规范来源；`Supp Fig 4` 和 `Supp Fig 5` 仅作 QA，不应再次进入模型或 meta-analysis。
2. `Supp Fig 2` 与 Figure 3 的 onset 列明显重叠，但版本并非全都逐行一致。数值标准化后，ALDOB、ASS1、DMGDH、FAH、GSTA1、HPD、OTC 的 DO/NDO 序列一致；ACO1、CPS1、LECT2 出现插入/位移或额外值。不能把两个 sheet 合并扩充样本量，应回查原始 manifest 后指定唯一版本。
3. `Figure 5` 和 `Supp Fig 7` 提供 biomarker–ALT 同行配对，可用于相关性分析，但仍来自与 Figure 3 相同的 confirmatory 受试者，不构成外部验证。
4. `Figure 6` 两表共用同一 HV 数据；只应计一次。
5. `Figure 2C` 中 ALDOB/HPD 完全重复是当前最严重的源数据一致性问题。在未解决前，发现队列 HPD 结果应标为不可判定或仅作敏感性排除，不可据此增加证据权重。

## 四、推荐的安全提取规范

生成长表时保留以下字段：`source_sheet`, `source_block`, `source_cell`, `cohort`, `protein`, `group`, `within_group_position`, `measurement_scale`, `value_raw`, `value_numeric`, `is_ND`, `is_structural_padding`, `is_flagged`, `qc_note`。

- 数值解析：仅接受标准十进制/科学计数法；`3077.0266*` 解析为 3077.0266 并令 `is_flagged=TRUE`；`ND` 保持左删失标记；空白先区分短组列的结构填充和组内真实空缺。
- 不用 0 代替 ND。主分析报告检测率与检测值分布；连续分析采用左删失方法或至少进行“排除 ND / 预设 LOD 分数替代 / 秩检验”三重敏感性分析。
- 不以 `within_group_position` 冒充 subject ID。只有 Figure 5 / Supp Fig 7 同一行内的 ALT–protein 可以直接配对；DO–DF 配对和跨蛋白多变量矩阵必须先由外部样本清单、共享受试者 ID 或可验证的生物标志物指纹确认。
- 不能把重复 sheet 当成独立复现。证据单位应按实际 cohort 计数：discovery、confirmatory onset/follow-up，而不是按图号计数。
- 当前工作簿足以支持单蛋白组间效应、检测率、稳健 biomarker–ALT 相关和数据驱动候选排序；不带 subject ID 直接训练多蛋白患者级分类器、做配对随访检验或调整混杂因素均不安全。

## 五、主分析前必须冻结的 QC 决策

1. 回查并裁定 Figure 2C 的 ALDOB/HPD 重复。
2. 裁定 DMGDH-HV 的第 61 值（`Figure 3!B417`）和 CES1-DO 的内部空白（`Figure 3!C936`）。
3. 获取/核对 Supp Fig 2 与 Supp Fig 3 的实际组样本数和 DO–DF/NDO–NDF 配对清单。
4. 将 `Supp Fig 2!E627="NDF"` 作为数据错误，不作数值使用，并追溯正确值。
5. 明确 Supp Fig 7 对 HPD、OTC、PCK2 的 ND 替代规则；在未确认前，相关性结果必须带删失敏感性分析。
6. 对 `Supp Fig 3-PCK2!C28` 的星号值做保留/剔除双版本稳健性检查，不自行删除。
