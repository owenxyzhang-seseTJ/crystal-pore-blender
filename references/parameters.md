# 参数与输出

## 工作配置

使用 `crystal-pore-blender init --cif input.cif --out job.json` 生成配置。`cif` 路径相对 job 文件，其余输出路径相对当前目录。源索引来自本次对称展开的 `job.intake.json`。

| 字段 | 默认值 | 含义 |
| --- | --- | --- |
| `radii_model` | `pymatgen` | 元素表格 vdW 半径；`custom` 要求全部显式提供；`legacy` 为兼容旧数值模型的有限表 |
| `vdw_radii_A` | `{}` | 有来源的半径补充或覆盖，Å |
| `probe_radius_A` | `1.2` | 球形探针半径，Å |
| `grid_step_A` | `0.25` | 晶格方向采样步长上限，Å；斜晶胞不等同笛卡尔正交网格 |
| `surface` | `contact` | `contact` 或 `probe-center`；共用一次计算缓存 |
| `accessibility` | `percolating` | 保留周期贯通分量；`all` 包含孤立几何空腔 |
| `guest_groups` | `[]` | 完整客体的源索引列表，每组第一个原子为锚点；从计算中删除 |
| `guest_images` | 每组零平移一次 | 显示的完整客体镜像；设 `[]` 隐藏所有移除客体 |
| `display_mode` | `cell` | `cell` 显示裁剪中的全部孔；`single` 选一个锚点分量 |
| `crop_fractional` | `[0,0,0]` 至 `[1,1,1]` | 分数坐标矩形裁剪，可跨晶胞；实际端点落在网格上 |
| `pore_anchor_fractional` | 无 | 单孔模式无客体时的三维分数坐标锚点 |
| `pore_axis_uvw` | 仅唯一一维孔族推断 | 直接晶格方向，不是晶面法向；侧视/端视/剖切需要 |
| `open_crop_axes` | `[0]` | 单孔裁剪允许开口的网格轴编号；其它边界不能截到孔壁 |
| `framework_mode` | `network` | 有限周期网络显示；`molecules` 补全有限分子片段 |
| `node_elements` | `[]` | 配位节点元素；避免把配体内部金属误当节点 |
| `donor_elements` | `["N","O"]` | 配位键供体元素 |
| `coordination_cutoff_A` | `2.6` | 指定节点与供体的显示键阈值 |
| `bond_cutoffs_A` | `{}` | 元素对显示键覆盖，如 `{"Si-O":1.9}`；否则共价半径和乘 1.15 |
| `wall_extra_A` | `0.45` | 选邻近孔壁原子的额外距离 |
| `smooth_display` | `true` | 保留开口的三轮 Taubin 显示平滑 |
| `max_display_smoothing_A` | `0.15` | 显示顶点位移上限 |
| `convergence_step_A` | 无 | 额外计算较粗网格并比较 Contact 体积分数 |
| `convergence_tolerance` | `0.025` | 粗细 Contact 体积分数允许绝对差；需另核查孔网维度 |
| `max_voxels` | `4000000` | 周期计算体素上限，不是精确内存限额 |
| `max_display_voxels` | `8000000` | 展示裁剪体素上限 |
| `max_display_atoms` | `30000` | 防止过大展示原子集合 |
| `zeopp_executable` | 无 | 可选 Zeo++ 程序路径或 PATH 名称；写入匹配半径表后调用 |

`guest_images` 的格式为 `[{"group":0,"shift":[1,0,0]}]`。`group` 引用当前 `guest_groups` 中的一组，`shift` 必须是三个整数。原子索引不能跨 CIF 使用；批量模板不接受客体索引。

默认保留所有源原子参与计算。去客体时，保持源晶格不变；这是明确的几何模型，不代表框架弛豫后的实验相。未知 vdW 半径会报错，不自动猜值。

## 命令分步执行

```sh
crystal-pore-blender inspect --cif input.cif --out intake
crystal-pore-blender compute --job job.json --out calculation
crystal-pore-blender display --job job.json --calculation calculation --out display
crystal-pore-blender render --display display --out render
```

`run` 串联计算、显示、渲染和白底合成。`--no-render` 停在几何输出。`run` 不覆盖非空目录，除非使用 `--force`。分步命令允许重建各自目录；缓存失配会拒绝。

`--probe`、`--step`、`--surface`、`--accessibility`、`--axis` 在 `run` 和 `batch` 中覆盖 job 的对应字段。`batch --cif a.cif b.cif` 按文件名建立子目录，重复文件名会报错；每项结果及错误记入 `batch.json`，有失败时返回非零退出码。

## 风格与视角

复制 `crystal_pore_blender/styles/default.json` 后用 `--style` 加载。`palette` 覆盖元素颜色；未覆盖的常见非金属使用固定元素色，其余为中性灰。`framework_radii_A` 按元素覆盖显示半径，缺项用 `framework_radius_A`。这些都不影响计算半径。

`guest_representation` 支持 `spacefill` 和 `ball-stick`；前者用数值 vdW 半径乘 `guest_vdw_scale` 作显示球径，后者用 `guest_ball_radii_A` 或 `guest_ball_radius_A`。`guest_palette` 单独覆盖客体颜色。

`pore_outer`、`pore_inner` 是表面正面与背面颜色；`pore_opacity`、`cutaway_opacity` 分别控制完整和剖面透明度。`show_framework`、`show_guests`、`show_pore` 控制整个集合，保存后的视口也使用相同可见状态。

`--view free` 使用 `camera_direction`，`side` 投影到孔轴法平面并使孔轴横向，`end` 沿孔轴观察。`--mode cutaway --section longitudinal` 保留背向相机的半壳；`transverse` 以孔轴为法向。切面默认经过客体锚点均值或裁剪中心，可用风格中的 `section_origin_A` 指定笛卡尔原点。剖面不加盖、不裁原子。

`resolution` 为 `[宽,高]`，`samples` 为 Cycles 采样数。渲染自动尝试可用 GPU，失败回到 CPU。Blender 图层关闭后，相机仅按可见几何拟合。`section_outline` 目前仅接受 false。

## 输出

- `calculation/`：源索引、去客体 CIF、周期中心/Contact 标量场、实际半径、体积分数、winding、维度、收敛和 Zeo++ 执行记录。
- `display/`：原始/显示 OBJ、包含两套顶点的 NPZ、原子与客体 JSON、坐标及平滑位移核查。
- `render/`：`figure.blend`、`figure_transparent.png`、`figure_white.png`、最终风格和保存后重开核查记录；剖面附剖切网格。
- `job_used.json`、`run.json`：运行参数和各阶段完成状态。

计算文件记录本机源文件路径和用户自己的结构数据；自行分享输出时按实际需要选择文件。公开仓库不需要提交个人任务输出。
