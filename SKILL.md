---
name: crystal-pore-blender
description: 从真实 CIF 计算周期可达孔道并在 Blender 中绘制原子框架、客体与孔道的完整图或纵向剖面图。用于 MOF 等晶体的孔道表面可视化、单孔选择、去客体计算但保留客体显示，以及透明背景可编辑科学图。
---

# CIF 孔道与 Blender 剖面

交付从 CIF 坐标与周期孔隙计算得到的表面、真实原子框架和客体；包括可编辑 `.blend`、透明/白底 PNG、原始与显示网格、参数及必要验证。不得用手画圆管、任意路径或 AI 图像替代计算孔道。

## 工作边界

- 源 CIF 保留不变，计算模型、显示模型分开。水合晶格去水计算与脱水相 CIF 是两种不同模型。
- 在多个 CIF、客体、配色或视角中，优先沿用用户指定的最终版本；不要根据相似文件名另选源文件。
- 不绑定 Zn、Co、TTF 或水。客体选择、节点元素和显示范围依当前结构确定。本次 TTFP 风格是可选预设，不是所有材料的固定要求。
- 仅在关键输入缺失时问用户：哪个 CIF/相、移除哪些客体、无序的哪一个明确模型、哪个互不连通孔道。已能核查的内容先从文件中确认。
- 源占位或无序未经处理时停止计算并请求明确模型；不要自动选 major occupancy、补氢、优化或重定位客体。

## 科学决策与执行

1. **读取源结构并记录索引。** 显式展开对称性，核查晶格、组成、占位和客体。需要修复 CIF 时使用可用的 `clean-cif`；修复输出另存。所有分数坐标统一以 `X = F @ A` 转成同一个行存储晶格。与旧 Blender 合并时核对晶格 Gram 矩阵及基底，不能直接混用不同解析器的笛卡尔坐标。
2. **锁定孔道定义。** 确认探针半径、原子 vdW 半径表、可达/贯通条件和表面类型。用户未另指定且沿用本流程时可取 0.9 Å、Contact Surface、Zeo++ 默认半径。半径表是数值模型的一部分，不能为了视觉效果更换。
3. **孔网和周期表面。** 优先用实际可执行的 Zeo++ 做 `-res`、`-chan` 分析；必要时用周期场 winding 独立检查孔道维度和轴向。Contact Surface 是可达探针中心集合按完整球形探针扩张后的空间边界，不是中心可达空间的 SAS。详见 [方法与边界](references/methods.md)。不能把替代计算称为 Mercury 原生导出。
4. **选择一个实际孔道。** 联合选择客体的整数晶格平移，以连通分量判定是否属于同一孔。不要给各客体分别找“最近的配体”，这会把它们放入不同平行孔。裁剪后核查单个孔分量、侧壁不被边框截断、指定锚点在孔内。沿轴端口可以开放。
5. **补齐真实孔壁。** 按原子到计算孔壁的邻近关系选取周期框架，并补全有限配体及配位节点。分子整体展开和移动，不能逐原子 wrap。无穷共价网络不适用“整分子全部补齐”；明确有限展示边界并保持周期键，或另用相应结构工具。
6. **渲染。** 先用配置确定视角、原子/客体比例、孔壁颜色与透明度。使用下面脚本或与已有 `mof-blender` 场景衔接；无需重复做与孔道无关的拓扑分析。源原子位置保持不变，只调整显示半径、材质和灯光。
7. **真实输出验收。** 保存并重开 `.blend` 核对原子、客体、网格和默认图层；检查实际最终 PNG，保证没有意外裁边、文字、多面体、暗原子或不可见客体。保留原始计算网格；显示平滑不能反馈到孔径等数值结论。

## 脚本入口

先读 [参数与文件接口](references/parameters.md)。依当前 CIF 填写 `job.json`，其中索引基于本脚本实际解析结果，不可沿用另一解析器的索引。

```sh
python scripts/pore_pipeline.py inspect --cif input.cif --out intake
python scripts/pore_pipeline.py compute --job job.json --out calculation
python scripts/pore_pipeline.py display --job job.json --calculation calculation --out display
blender --background --factory-startup --disable-autoexec --python-exit-code 1 \
  --python scripts/render_pore.py -- --data display/display.json \
  --mesh display/pore_mesh.npz --style assets/ttfp-approved-style.json \
  --out render --mode cutaway
python scripts/white_background.py render/figure_transparent.png
```

依赖：几何计算需 NumPy、SciPy、pymatgen、scikit-image；渲染用 Blender 内置 Python/NumPy；白底合成用 Pillow。工具位置需在当前环境定位，不把此项目的绝对路径写入新任务。只有确实缺少依赖且需要执行时才安装。

`display` 可直接复用包含 `contact_level_A`、`lattice_A` 的已验证周期场；同目录需有 `intake.json` 和匹配的 `calculation.json`。不要把其他 CIF 或另一半径模型的缓存拿来拼接。

## 本次认可的视觉预设

加载 [TTFP 配置](assets/ttfp-approved-style.json) 时默认：框架正常元素色、统一细球棍；仅客体突出；单孔侧视沿轴横向展开；孔壁轻微平滑、青蓝透明；无文字、无晶体形貌。水的紧凑空间填充采用 vdW 显示半径的 0.55 倍。完整侧视和半剖图可独立配置客体表示。

纵向半剖是与孔道轴平行的真实平面裁切，只移除朝向相机的半个表面，不封口、不移动客体。孔壁外侧浅青蓝，内侧稍深蓝灰。默认不加切口细边缘。透明度控制着色混合而非折射玻璃，避免孔壁把原子压暗。

仅温和平滑显示网格并固定原有开口，记录最大和平均位移。预设三轮 Taubin `lambda=0.25, mu=-0.26`，本次网格最大位移不到 0.03 Å；换结构后必须重新测量，不能承诺仍是该数值。需要单独导出的科学孔道图，保留未平滑网格。

## 常见纠错

- **看上去是两根孔道：** 检查所有客体锚点是否在同一有限连通分量；联合调整整个孔族的整数平移，不能搬动分子内部坐标。
- **孔道悬在原子外：** 孔隙必须由完整周期计算框架得到；增加真实孔壁邻近配体，不缩放孔道来适应较小原子片段。
- **剖面过暗或客体变色：** 先使用明亮填充光及透明着色；不用折射玻璃或黑色内壁。透明着色的多层叠加效果与单层不相同。
- **水过大：** 改显示 scale，保留 O/H 坐标。空间填充不意味着必须按 100% vdW 半径展示。
- **“只显示孔道”仍留有原子：** 隐藏所有非孔道显示集合，包括配位节点、共价键、客体，不仅检查名为 `atoms` 的集合。
- **沿轴被裁边：** Blender 正交相机 `ortho_scale` 按水平宽度拟合，明确 `sensor_fit=HORIZONTAL`，用图幅比例约束高度。
- **Mercury 高级接口不可用：** 说明实际限制并用已授权的可行方法；不要反复安装大 SDK、逆向接口或把未成功计算写成结果。

具体已验证数据和适用范围见 [本次案例记录](references/verified-case.md)，仅在复现这套图或查错时读取。
