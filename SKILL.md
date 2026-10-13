---
name: crystal-pore-blender
description: 从有序 CIF 计算周期孔隙表面并生成可编辑 Blender 场景。用于 MOF、COF、沸石和分子晶体的孔道可视化、整晶胞或单孔显示、去客体计算但保留客体显示、纵向或横向剖面、透明底和白底科学图，也可批量导出网格。
---

# 晶体孔隙与 Blender

先计算周期几何，再生成显示。保留源 CIF、原始网格和参数，交付用户请求的 PNG、`.blend` 或 OBJ。不要以手画管道或生成图像替代孔隙计算。

## 执行

1. 确定源 CIF 和相态；从文件核查结构，不按相似文件名替换。无序或部分占位需明确的有序模型。保留解析警告，必要时检查原文件。去客体的源晶格与另一结构相的 CIF 分别处理。
2. 运行 `crystal-pore-blender doctor` 检查依赖。若未安装，在本仓库运行 `python -m pip install .`。渲染需要 Blender 4.2+ 或匹配的 `bpy`；仅导出网格用 `--no-render`。不要自动安装多个大软件或把未完成渲染写成成功。
3. 一般整晶胞图用 `run --cif input.cif --out output`。定制任务先 `init --cif input.cif --out job.json`，从生成的索引表配置客体、半径、孔轴和范围，再 `validate --job job.json`。索引只属于本次解析结果。
4. 读取[参数](references/parameters.md)配置任务。默认不删除客体；仅移除用户指定的完整分子。`guest_groups` 控制计算移除，`guest_images` 控制原位显示。按周期键图展开跨胞客体，不逐原子找最近锚点。客体结构变化只允许整数晶格平移。
5. 读取[方法](references/methods.md)确定表面及连通条件。默认 Contact Surface、1.2 Å 探针和 pymatgen 表格半径是起始配置；按当前问题选择，记录实际数值。`all` 包含孤立几何空腔，不能称为外部气体可达。用于数值比较时检查网格收敛及孔网维度。
6. 用整晶胞 `cell` 展示一般孔网；确需单孔时设 `single`，选择锚点、联合客体镜像与裁剪。检查一个有限分量及未裁断的侧壁。二维、三维孔网不自动假定唯一孔轴；使用用户指定方向。无限共价网络用 `network` 有限裁剪，有限分子补全用 `molecules`。
7. 运行请求的视角或剖面。配色与图层通过风格 JSON 调整；渲染已有显示时用 `render`，不重复计算。只运行必要的变体。显示球径和平滑不改变计算半径或数值结果。
8. 检查实际最终 PNG 和重开的 `.blend`。查看坐标误差、网格和透明边界。保证客体可辨、无意外裁边，白底版本来自同一透明图。核查完整图层的可见状态，不能只检查原子集合名。

## 常用命令

```sh
crystal-pore-blender run --cif input.cif --out output
crystal-pore-blender init --cif input.cif --out job.json
crystal-pore-blender run --job job.json --out output/section --mode cutaway --view side
crystal-pore-blender render --display output/display --out output/recolored --style style.json
```

横向剖面使用 `--mode cutaway --section transverse`；侧视、端视和剖切需明确 `pore_axis_uvw` 或 `--axis`，一维唯一孔族可自动推断方向。整晶胞自由视角不需要孔轴。

源和数值模型相同时，`display` 可复用 `compute` 缓存；改变探针、计算半径、客体移除或连通选择需重新算场。保留原始 mesh，平滑只用于显示。Zeo++ 是可选独立孔网分析；只有报告实际成功运行时才声称使用了它。

## 查错

- 多条平行孔：核对单孔锚点、裁剪与所有客体的联合周期镜像。
- 缺少孔壁或跨胞键：检查 `framework_mode`、裁剪、孔壁邻近范围与显示键阈值，不缩放孔隙去适配原子。
- 无可达孔道：检查客体、探针和分辨率；若研究孤立空腔才选 `accessibility: "all"`。
- 孔壁过暗：调透明度、配色与光照，不改变科学半径。
- 缺元素半径：查明有来源的 vdW 值并显式配置；显示颜色可使用默认灰色。

只对关键且无法从文件确认的科学输入向用户澄清。工具路径在当前环境定位。不要执行 SHA256 验证或无关的 smoke test，用任务相关的真实输出检查验收。
