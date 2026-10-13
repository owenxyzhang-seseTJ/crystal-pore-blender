# 周期孔隙的计算方法

## 定义

行存储晶格为 `A`，分数坐标通过 `X = F @ A` 变成 Å。框架原子为 `r_i`，数值 vdW 半径为 `R_i`，周期净空为：

`c(x) = min_i,t (|x - r_i - t @ A| - R_i)`。

球形探针半径为 `p`。几何允许的探针中心集合满足 `c(x) >= p`。默认保留周期贯通的分量；`accessibility: "all"` 保留全部允许中心，包括孤立空腔。

`probe-center` 的等值场是 `c(x) - p`。Contact 场是 `p - distance(x, selected_center_set)`，表示所选中心集合经完整球形探针扩张后的空间边界。Contact 场通过采样中心构建，因此仍有网格离散误差。两种表面都用 marching cubes 构网。

数值半径来自 pymatgen 的 `Element.van_der_waals_radius` 表，显式覆盖优先；最终使用值写入报告。`custom` 模型要求每个实际元素都有显式值。`legacy` 仅为兼容有限旧半径表，不能假定与任意 Zeo++ 安装的默认值一致。显示球径独立配置。

## 周期性和资源

计算在完整周期晶胞中完成，展示再裁剪。镜像范围使用倒易基底给出的晶胞面间距确定，支持倾斜晶格。原子场的查询按块执行，Contact 扩张复用一个中心 KD tree，逐镜像查询，不存储全部中心镜像。

体素数量上限在分配前检查；该上限不是精确 RAM 预算。高分辨率、大晶胞、病态倾斜晶胞仍可能消耗较多资源。粗网格可用于预览，不能据此保证窄窗口连通。

净空远离所需阈值时被截断。缓存能判断允许中心和构建表面，不能从截断场求最大腔径、Di 或 Df。网格体积分数是此半径和探针模型下的几何量，也不是吸附量或 BET 孔容。

## 孔网维度

局部体素按六邻接标签，然后在周期边界合并。合并过程记录整数晶格平移；闭环的非零 winding 才证明贯通。维度是 winding 向量的秩，范围 0 至 3。

`[1,1,0]` 的一维斜向通道不能报告为独立纯 a 轴与 b 轴。报告同时保留 winding 向量、维度，以及其张成空间包含哪些纯晶格轴。粗细网格的体积分数接近仍不能保证拓扑相同，定量任务应核查维度和方向。

整晶胞内可能出现多个有限孔隙片段。`single` 用客体锚点或指定锚点选择一个有限连通分量，并检查禁止开口的边界；`cell` 展示裁剪中全部分量。多个孔族、2D/3D 孔网没有自动唯一轴，需要明确显示方向。

## 客体和框架

源 CIF 先显式展开对称性，禁用分数坐标的理想化舍入，并记录解析警告。无序及部分占位需明确的有序模型；不自动挑选主占位、补氢或弛豫结构。

客体通过周期键图展开，保证各原子的差异只是整数晶格平移。分组需对应有限连通分子；断开的组或带 winding 的组会报错。键按共价半径或显式阈值识别，特定配位键需要检查。

周期框架默认使用有限网络范围，保留范围内真实周期原子与键。`molecules` 模式补全靠近孔壁的有限组分，不允许把无限共价网络当作一个可补全分子。显示键仅用于绘图，不用于判定电子结构或化学机制。

## 平滑和剖面

显示网格作三轮 Taubin 正负步平滑，原始开放边界固定；原始网格另存。报告实际顶点位移，不承诺不同结构具有相同误差。

剖面按三角形与平面的真实交点裁剪。纵向平面法向垂直显示轴，横向平面法向平行显示轴。切面不加盖。剖切和平滑只影响显示，不能反馈为新的科学孔径或数值孔容。

## 可选 Zeo++

指定 `zeopp_executable` 后，程序生成与当前模型完全一致的两列 `.rad` 半径文件，通过 `-r` 传给 Zeo++，然后运行 `-ha -res -chan`。只有返回成功且产生本次非空输出，才记录原生分析成功。未指定程序时，报告 `zeopp.ran=false`。

本工具的表面仍是周期网格方法的产物。Zeo++ 提供独立孔网分析，不能把网格称为 Zeo++ 或 Mercury 原生表面导出。当前自动检查执行和输出文件；孔径、通道数量及方向应查看原生报告并与网格结果比较。

## 来源

- [pymatgen CIF parser](https://pymatgen.org/pymatgen.io.html)
- [pymatgen element properties](https://pymatgen.org/pymatgen.core.html)
- [scikit-image marching cubes](https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.marching_cubes)
- [Zeo++ methods and publications](https://www.zeoplusplus.org/docs.html)
- [Zeo++ command reference, including radius files](https://github.com/richardjgowers/zeoplusplus/blob/master/documentation/README)
