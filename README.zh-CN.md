# crystal-pore-blender

把 CIF 里的孔隙变成可编辑的 Blender 场景。

看晶体里的原子不难，想把原子之间的空间画清楚就麻烦一些：孔壁要连续，剖面要能看见内部，客体还得留在正确的位置。这个工具先计算周期孔隙，再生成 Blender 场景。你可以直接出图，也可以打开 `.blend` 继续调整。

适用于有序晶体，包括 MOF、共价网络、沸石和分子晶体。仓库没有附带具体课题的结构或原子索引。示例是人为构造的碳原子孔道，只用于演示。

[English](README.md) · [参数说明](references/parameters.md) · [计算方法](references/methods.md)

![合成示例孔道的纵向剖面](examples/preview.png)

示例的实际渲染图。孔壁来自原子模型的周期计算。

## 安装

需要 Python 3.11 或更高版本。下载或克隆仓库后，在仓库目录运行：

```sh
python -m venv .venv
# macOS / Linux
source .venv/bin/activate
# Windows PowerShell 使用 .venv\Scripts\Activate.ps1
python -m pip install .
crystal-pore-blender doctor
```

生成 PNG 和 `.blend` 还需要 [Blender 4.2 或更高版本](https://www.blender.org/download/)。程序会查找 PATH，以及 macOS、Windows 的常见安装目录。没找到时，用 `--blender "程序路径"` 指定。

也可以运行 `python -m pip install '.[render]'` 安装 `bpy`，直接通过 Python 渲染。它需要与你的 Python 版本、系统匹配的 wheel；找不到时，使用独立的 Blender 程序即可。

项目还没有发布到 PyPI。除了在下载的仓库中 `pip install .`，也能直接安装 GitHub 版本：

```sh
python -m pip install "git+https://github.com/owenxyzhang-seseTJ/crystal-pore-blender.git"
```

## 先跑一张图

安装后直接运行下面的命令，生成合成孔道的纵向剖面，无需另找结构文件：

```sh
crystal-pore-blender demo --out output/demo --mode cutaway --view side
```

示例输入保存在 `output/demo/input`，结果在 `output/demo/scene`。如果终端找不到命令，用 `python -m crystal_pore_blender` 替换 `crystal-pore-blender` 即可。PowerShell 中把多行命令写成一行。

换成自己的结构：

```sh
crystal-pore-blender run --cif structure.cif --out output/structure
```

默认显示一个晶胞的孔隙，不会猜测或删除溶剂。只要网格、不需要渲染时，加 `--no-render`，这样无需安装 Blender。

输出包括透明底和白底 PNG、可编辑 `.blend`、原始及显示 OBJ 网格、周期场缓存，以及记录参数和核查结果的 JSON。原子坐标只作整数晶格平移；球径、颜色、灯光属于显示设置。

## 常用调整

```sh
# 指定探针和网格步长，仅导出几何文件
crystal-pore-blender run --cif structure.cif --out output/fine --probe 0.9 --step 0.15 --no-render

# 绘制探针中心区域的边界
crystal-pore-blender run --cif structure.cif --out output/centers --surface probe-center

# 同时包含孤立空腔
crystal-pore-blender run --cif structure.cif --out output/cavities --accessibility all

# 沿指定晶格方向看孔道端面
crystal-pore-blender run --cif structure.cif --out output/end --axis 1 0 0 --view end

# 复用已有显示网格，换配色重新渲染
crystal-pore-blender render --display output/structure/display --out output/new-view --style my-style.json

# 批量处理；失败的结构及原因记入 batch.json
crystal-pore-blender batch --cif first.cif second.cif --out output/batch --no-render
```

纵剖用 `--mode cutaway --view side`；横剖再加 `--section transverse`。剖切需要确定 `pore_axis_uvw`，或者用 `--axis u v w` 指定。二维、三维孔网没有唯一的孔轴，需要自己选择方向。剖切只裁孔隙表面，保留框架和客体。

配色、透明度、显示球径、分辨率及图层设置放在风格 JSON 中。从 [default.json](crystal_pore_blender/styles/default.json) 复制一份，用 `--style` 加载。常见非金属使用元素色，其他元素默认灰色。`show_framework`、`show_guests`、`show_pore` 可以关闭整个图层。

## 去客体计算，保留客体显示

先为自己的 CIF 生成配置和索引表：

```sh
crystal-pore-blender init --cif structure.cif --out job.json
crystal-pore-blender validate --job job.json
```

`job.intake.json` 记录对称展开后的源索引、标签和坐标。把需要从计算中删除的完整客体写入 `guest_groups`，再通过 `guest_images` 选择要显示的周期镜像。客体按周期键图整体展开，避免跨晶胞分子被拆散。

默认是整晶胞显示。要选单孔，设置 `display_mode: "single"`，填写裁剪范围和 `pore_anchor_fractional`，或用显示客体作为锚点。程序检查锚点是否属于同一个有限连通分量、侧壁是否被裁断。具体格式见[参数说明](references/parameters.md)。

`framework_mode: "network"` 适合有限范围内的周期网络；`"molecules"` 会补全邻近的有限分子片段，遇到无限共价组分时停止。显示键通过距离判断，特定配位方式应核查或修改阈值。

## 如何理解孔壁

默认 Contact Surface 表示球形探针在所选中心区域内运动后占据空间的边界。`probe-center` 表示探针中心可处的位置边界，两者含义不同。默认只选周期贯通分量；`all` 也包含孤立的几何空腔，不能据此认为外部气体一定能进入。

计算记录 winding 向量和孔网维度。网格精度会影响窄窗口的连通性，定量比较前应检查收敛。可选 Zeo++ 分析会写入与当前计算一致的半径表，并记录实际执行结果；本工具的网格是周期场构网结果。

无序、部分占位需先提供明确的有序模型。缺失的 vdW 半径需显式给出。程序不会自动补氢、选占位、移动客体或优化结构。[方法说明](references/methods.md) 给出定义和适用边界。

## 作为 Codex Skill

把仓库放进 Codex 的 skills 目录，重新打开会话后调用 `$crystal-pore-blender`，例如：

> 根据这个 CIF 生成孔道纵向剖面，客体保持原位，给我白底图、可编辑场景和原始网格。

[SKILL.md](SKILL.md) 是给智能体的执行说明。命令行工具可以独立使用。

## 开发与验证

```sh
python -m pip install '.[dev]'
python -m pytest -q
python -m build
```

测试覆盖周期贯通、倾斜晶胞、跨胞客体、有序模型检查、缓存匹配和命令行。GitHub Actions 在 Python 3.11、3.12 上运行。全部测试结构为合成数据。
