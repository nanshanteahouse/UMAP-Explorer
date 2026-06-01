# Interactive UMAP Explorer — 交互式UMAP可视化平台

基于 Dash (Plotly) 构建的单细胞RNA测序数据交互式UMAP可视化Web应用。覆盖11个GEO视网膜发育数据集，提供三种探索模式：单数据集分析、跨数据集对比、以及基于Harmony整合的共享UMAP空间。

## 页面功能

| 页面 | 路径 | 功能说明 |
|------|------|----------|
| **单数据集浏览** | `/` | 逐个数据集探索——支持按细胞类型、基因表达、QC指标、Leiden聚类着色。附带小提琴图展示表达分布。 |
| **数据集对比** | `/comparison` | 两个独立数据集的并列UMAP视图，支持联动着色和基因检索。 |
| **整合视图** | `/integrated` | 所有11个数据集投影到Harmony校正后的共享UMAP空间，可跨物种、跨发育阶段观察全局细胞类型关系。 |

## 项目架构

```
app.py                          # Dash入口——初始化、全局布局、页面注册、WSGI应用导出
├── pages/
│   ├── single_dataset.py       # 主页 (/)，回调逻辑、控件面板、基因检索、小提琴图
│   ├── comparison.py           # /comparison，双UMAP并列视图、独立数据集选择器
│   └── integrated.py           # /integrated，Harmony整合多数据集视图
├── components/
│   ├── data_loader.py          # 基于PyArrow内存映射的Parquet加载器，模块级缓存
│   ├── umap_figure.py          # 核心UMAP渲染——Scattergl (2D) / Scatter3d (3D)，自适应降采样
│   ├── controls.py             # 侧边栏——数据集下拉框、着色方式选择、基因检索触发
│   ├── gene_selector.py        # 基因搜索/自动补全组件
│   ├── info_panel.py           # 数据集元信息和文献引用展示
│   └── violin_plot.py          # 基因表达分布小提琴图
├── preprocessing/
│   ├── build_data_cache.py     # 从04_clustered.h5ad提取Parquet缓存
│   ├── build_integrated.py     # Harmony整合 → 共享UMAP空间
│   ├── extract_gene_metadata.py # 构建gene_index.json用于自动补全
│   └── validate_cache.py       # 缓存验证：NaN检查、数值范围、列完整性
├── assets/
│   └── style.css               # 单一CSS样式表（Dash自动加载）
├── data/
│   ├── dataset_registry.json   # 数据集元信息注册表（已入库）
│   ├── gene_index.json         # 基因搜索索引（已入库）
│   └── *.parquet / *.h5ad      # 缓存数据（gitignore，由预处理脚本生成）
├── Dockerfile                  # 生产镜像——python:3.11-slim，gunicorn 4 workers，端口8050
├── nginx.conf                  # 反向代理——WebSocket升级、50M请求体、静态资源缓存
├── start.sh                    # 本地开发启动脚本
└── requirements.txt
```

## 快速开始（开发环境）

**前置条件**：Python 3.11+、虚拟环境、`data/` 目录中已构建好缓存数据。

```bash
# 1. 克隆仓库并创建虚拟环境
python3 -m venv ../.venv
source ../.venv/bin/activate
pip install -r requirements.txt

# 2. 构建数据缓存（需要外部目录中的04_clustered.h5ad文件）
python preprocessing/build_data_cache.py --data-dir /path/to/h5ad_files
python preprocessing/extract_gene_metadata.py
python preprocessing/validate_cache.py

# 3. （可选）构建整合视图
python preprocessing/build_integrated.py

# 4. 启动应用
./start.sh                # 默认端口 8050
./start.sh -p 8051        # 自定义端口
```

浏览器打开 **http://0.0.0.0:8050** 即可访问。

## 数据流水线

应用运行时不会直接读取原始 `.h5ad` 文件。预处理流水线会将其转换为基于 PyArrow 内存映射的列式 Parquet 文件，实现快速数据访问。

```
04_clustered.h5ad (外部目录)
    │
    ├──[build_data_cache.py]──→ data/{dataset_id}/umap_metadata.parquet
    │                           data/{dataset_id}/gene_expression.parquet
    │
    ├──[extract_gene_metadata.py]──→ data/gene_index.json
    │
    └──[build_integrated.py]──→ data/integrated/integrated_umap.parquet
                                data/integrated/integrated_info.json

validate_cache.py ← 缓存构建完成后运行，检查数据完整性
```

### 预处理脚本一览

| 脚本 | 输入 | 输出 | 说明 |
|------|------|------|------|
| `build_data_cache.py` | 每个数据集的 `04_clustered.h5ad` | `data/{id}/*.parquet` | 使用 `--data-dir` 指定外部h5ad目录 |
| `extract_gene_metadata.py` | 已缓存的parquet文件 | `data/gene_index.json` | 为前端基因自动补全提供索引 |
| `build_integrated.py` | 所有已缓存的parquet文件 | `data/integrated/` | 基于Harmony-PyTorch的批次效应校正 |
| `validate_cache.py` | 已缓存的parquet文件 | 验证报告 | NaN、数值范围、列覆盖度检查 |

## 生产环境部署

```bash
# 构建镜像
docker build -t vis-website .

# 启动容器（挂载数据卷）
docker run -p 8050:8050 -v /path/to/data:/app/data vis-website
```

完整生产栈（Nginx前置）：

```bash
docker run -d --name vis-app -v /path/to/data:/app/data vis-website
# 然后使用提供的 nginx.conf 配置反向代理
```

### Docker镜像详情

- **基础镜像**：`python:3.11-slim`（scanpy/numba兼容性要求）
- **运行时**：Gunicorn，4个工作进程，120秒超时（适配大规模UMAP渲染）
- **用户**：非root用户 `appuser`（UID 1000）
- **健康检查**：每30秒通过HTTP探针对8050端口进行检查
- **数据**：将 `/app/data` 挂载为卷，数据文件不打包进镜像

## 核心设计决策

- **不使用 Plotly Express**——所有图表均通过 `plotly.graph_objects.Figure` 和 `make_subplots` 构建，以获得对渲染的完全控制。
- **内存映射 Parquet**——元数据通过 `pyarrow.parquet.read_table` 以内存映射方式加载；基因表达列按需按基因查询加载，避免一次性加载完整表达矩阵。
- **自适应降采样**——当数据集超过50万个细胞时，UMAP渲染器会进行均匀降采样以保证交互流畅性。
- **客户端基因搜索**——`gene_index.json` 每个会话只加载一次，实现即时自动补全；实际表达值由服务端按需查询返回。
- **Scattergl / Scatter3d**——2D UMAP 使用GPU加速的WebGL渲染；3D视图使用Plotly的Scatter3d轨迹类型。
- **模块级缓存**——已加载的DataFrame缓存在模块作用域中，避免跨回调的重复读取。

## 开发注意事项

- **无 Node.js / npm / 前端构建工具**——纯Python项目，单一 `style.css` 通过Dash静态资源机制自动提供服务。
- **暂无测试**——仓库中目前没有自动化测试。`validate_cache.py` 提供了数据完整性校验。
- **虚拟环境**——约定位于仓库上级目录的 `../.venv`。
- **回调作用域**——回调函数定义在页面模块中，而非组件模块中。仅使用标准Dash回调（无模式匹配回调）。
- **.omo/** 和 **.sisyphus/** 目录为AI助手工作区，已被gitignore。

## 文件大小说明

- `data/` 目录包含GB级别的Parquet文件，已被gitignore排除。
- `gene_index.json` 可能很大（数千基因 × 多数据集），基因表达查询在回调中由服务端过滤。
- `dataset_registry.json` 保存元数据（物种、发育阶段、细胞数、列模式），需与缓存数据保持同步。

## 依赖项

| 包名 | 用途 |
|------|------|
| `dash>=4.1.0` | Web框架，支持多页面路由 |
| `plotly>=6.7.0` | 交互式UMAP图表（Scattergl、Scatter3d） |
| `gunicorn>=26.0.0` | 生产级WSGI服务器 |
| `pandas>=2.3.0` | DataFrame数据操作 |
| `numpy>=2.4.0` | 数值数组计算 |
| `pyarrow>=24.0.0` | 内存映射Parquet I/O |
| `scipy>=1.17.0` | 稀疏矩阵支持 |
| `anndata>=0.12.0` | 单细胞数据容器（整合视图） |
| `scanpy>=1.12.0` | 单细胞分析工具（预处理） |
| `harmony-pytorch>=0.1.8` | 批次效应校正（整合UMAP） |

## 数据集

项目覆盖11个GEO单细胞RNA测序数据集，涵盖人类和小鼠视网膜发育的多个时间点。各数据集通过 `data/dataset_registry.json` 注册，包含以下元信息：

- 物种（Human / Mouse）
- 发育阶段
- 细胞数量
- 可用标注列（细胞类型、Leiden聚类、QC指标等）
- 高变基因列表
- 文献引用

如需添加新数据集，需按数据流水线构建Parquet缓存，并在 `dataset_registry.json` 中添加对应条目。
