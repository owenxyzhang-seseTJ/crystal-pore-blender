# 参数与文件接口

## 入库和索引

`inspect` 生成 `intake.json`：源 CIF、完整对称展开后的源索引、元素、标签、分数坐标及行存储晶格。分数坐标内部模 1，仅用于周期坐标；显示时对完整分子再展开。不得把 VESTA、Gemmi、pymatgen 或其他旧文件的整数索引直接互换。

`job.json` 的 `cif` 路径相对于 job 文件，其他本机可执行路径为实际定位的路径。计算脚本默认使用 pymatgen 的晶格方向；若要接现有 Blender，设置 `lattice_A` 为该场景的行存储矩阵，脚本检查与源晶格度量一致。

下面是已经确认客体索引后的配置形状。示例的索引、平移与范围来自本次 Zn 水合案例，新结构必须先重新检查。

```json
{
  "cif": "input.cif",
  "guest_groups": [[156,52,53],[157,54,55]],
  "node_elements": ["Zn"],
  "probe_radius_A": 0.9,
  "grid_step_A": 0.14,
  "convergence_step_A": 0.20,
  "pore_axis_uvw": [1,0,0],
  "open_crop_axes": [0],
  "crop_fractional": {"min":[-0.15,0.01,-0.52],"max":[3.15,0.91,0.28]},
  "guest_images": [
    {"group":0,"shift":[0,0,-1]}, {"group":1,"shift":[0,0,0]},
    {"group":0,"shift":[1,0,-1]}, {"group":1,"shift":[1,0,0]},
    {"group":0,"shift":[2,0,-1]}, {"group":1,"shift":[2,0,0]}
  ]
}
```

`guest_groups` 内每组都是完整客体的源索引，第一个原子是孔道锚点，水应以 O 开头。这些组全部在孔道计算中移除；`guest_images` 指定仍需显示的完整客体周期镜像。非选定客体保持原样；若其它客体也应移除，显式列出。

脚本以锚点展开客体，适合尺寸小于晶胞半尺度的离散分子，如本次水。大的跨胞客体应先按实际分子键连通整体展开，并扩展脚本的数据输入，不可把离锚点最近的各原子镜像当成普遍正确的分子重建。

`node_elements` 指框架配位节点；不要一概排除所有金属，因为 Fc 配体内部的 Fe–C 要保留为完整组分。`donor_elements` 默认 N/O，`coordination_cutoff_A` 默认 2.60；新材料依据实际配位验证后设置。其他键阈值通过 `bond_cutoffs_A` 提供，如 `{"C-P":1.95}`。

`vdw_radii_A` 提供额外元素或有来源的覆盖。改变它必须重新计算场，并使用相同半径模型核对原生孔网；它与 `style` 中显示球径是两套参数。

`max_voxels` 默认 12000000；若超限，选适当晶胞或分辨率，不忽略资源限制。`smooth_display` 默认 true，`max_display_smoothing_A` 默认 0.10（明确的失败上限，并非应达到的位移）。

`zeopp_executable` 可设置实际原生程序绝对路径。未设置时只完成周期网格方法，报告 `zeopp.ran=false`；不得声称已使用 Zeo++ 验证。需要独立原生验证时定位/安装合法可执行版本并设置此项。

## 计算产物

`compute` 保存：

- `intake.json`：与缓存对应的源结构。
- `framework_for_pores.cif`：显式删除指定完整客体后的计算模型。
- `periodic_probe_field.npz`：`contact_level_A`、`accessible_centers`、`lattice_A`、`probe_radius_A` 和计算框架坐标。
- `calculation.json`：半径、组分、精度、winding、去客体索引、原生执行状态。
- 启用原生程序时的 `zeopp.res`、`zeopp.chan`、`zeopp.log`。

`display` 复用这些产物，核对源、晶格、客体选择及半径参数。只修改裁剪、客体周期镜像或显示平滑时不用重复算周期场。已有 legacy 场如缺 `intake.json`、`calculation.json`，只有当前源坐标、晶格、半径、去客体列表均经核对，才可补充对应元数据；不要跳过缓存匹配检查。

## 显示产物

`display/display.json` 包含 framework atoms / bonds、完整 guests / bonds、统一晶格、已验证孔轴、探针及半径模型。`pore_mesh.npz` 有原始顶点 `raw_vertices_A`、显示顶点 `vertices_A`、`triangles`；另有原始/显示 OBJ 和 validation。

有限配体选择按孔壁邻近关系补齐，其结果不承诺任何固定配体或原子数。无限共价网络会显式报错，需按相应周期网络的方法制备显示数据。不要把无限 COF 当作一个有限配体无限补齐。

## 渲染设置

风格文件可从 `assets/ttfp-approved-style.json` 复制修改：

- `guest_representation` 为 `spacefill` 或 `ball-stick`；spacefill 使用计算半径表乘 `guest_vdw_scale` 作为显示半径，默认 0.55。ball-stick 使用 `guest_ball_radii_A`，显示原子间源几何键。
- `palette`、`framework_radii_A` 必须包含实际元素，未知元素不随机着色。
- `pore_outer` / `pore_inner` 分别用于外表面和背向面（从腔内看到的内侧）。
- `pore_opacity` 控制完整图；`cutaway_opacity` 控制半剖。0.16 双层叠加近似 0.29 的单层不透明度；不是把孔道科学模型改成两个嵌套壳。
- `section_outline=false` 为本次认可设置。当前 helper 不添加轮廓线。
- `camera_direction`、`camera_up` 给侧视朝向倾向；脚本将方向投影到孔轴法平面，并使孔轴横向。沿孔端视或横向切片需要另写相机/裁切，不要套用纵向半剖。
- `resolution` 和 `samples` 可为必要预览临时降低；最终本次为 2000×1300、128 samples。

每次修改输出到新 render 目录。`figure.blend` 为主图设置，`figure_transparent.png` 为真实渲染；白底脚本输出 `figure_white.png`。`render_validation.json` 记录保存后重开校验与实际 CPU/GPU 后端，`style_used.json` 固化风格。无需渲染所有变体；只做用户请求的视角和确实有用的图层。
