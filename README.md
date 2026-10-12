# crystal-pore-blender

A Codex skill for CIF-derived periodic pore surfaces and editable Blender cutaways.

根据真实 CIF 计算连续周期孔道，将原子框架、客体和孔道拼在同一 Blender 场景中。支持单孔侧视和纵向半剖、透明/白底 PNG，以及可编辑 `.blend`。

计算与显示分离：可以在孔道计算中去掉客体，同时在原位置显示客体；不会把去客体的水合框架换成脱水相，也不会用手画管道代替实际孔隙。

## 安装与调用

将本仓库放入 Codex 的 skills 目录：

```sh
git clone https://github.com/owenxyzhang-seseTJ/crystal-pore-blender.git \
  "${CODEX_HOME:-$HOME/.codex}/skills/crystal-pore-blender"
```

在新会话中使用 `$crystal-pore-blender`，例如：

> 根据这个水合 CIF，用 0.9 Å 探针计算去水后的单孔，但渲染保留水。制作细框架、突出客体的纵向半剖图，孔壁内外稍异色，透明背景，不加切口轮廓线。

使用说明：[SKILL.md](SKILL.md)。参数：[references/parameters.md](references/parameters.md)。计算定义：[references/methods.md](references/methods.md)。已执行的验证：[references/verified-case.md](references/verified-case.md)。

## 依赖与执行

Python 环境安装几何/图像依赖：

```sh
python -m pip install -r requirements.txt
```

另外需要 Blender；独立孔网验证使用合法安装的 Zeo++，通过 `job.json` 的 `zeopp_executable` 指定。软件二进制、CIF 和研究数据均未打包在仓库中。

典型执行路径：

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

`job.json` 要根据 `inspect` 的源索引指定完整客体、节点元素、孔轴、联合周期镜像和裁剪范围。不要直接把参考材料的原子索引套到另一个 CIF。

## 方法与风格

孔道由周期 vdW 净空、贯通探针中心集合和完整球形探针扩张得到 Contact Surface，scikit-image 提取三角网格。Zeo++ 验证孔网，不将网格称为 Mercury 或 Zeo++ 原生表面导出。

默认视觉预设为正常元素色的细框架，紧凑空间填充客体（显示 vdW scale 0.55），浅青蓝外壁和蓝灰内壁，透明材质，无文字、形貌多面体和切口轮廓线。配色、尺寸、客体表示和透明度均可修改。完整侧视需要球棍客体时，将风格中的 `guest_representation` 设为 `ball-stick`。

只对显示网格温和平滑，保留原始 mesh；不得从缩放后的球径或显示剖面推断数值孔径。无序结构需显式模型，未知元素需明确 vdW 半径和颜色，无穷共价网络需相应周期展示实现。
