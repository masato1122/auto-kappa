# MLIPS Strict Relaxation

这个文档介绍如何使用新添加的 `MLIPSStrictRelaxation` 类来进行基于机器学习原子间势(MLIPS)的严格结构弛豫。

## 特性

- 支持多种MLIPS计算器 (eSEN, MACE等)
- 使用 Birch-Murnaghan 状态方程进行体积-能量拟合
- 自动优化结构并找到最优体积
- 兼容2D和3D材料
- 完全独立于原有的 `StrictRelaxation` 类

## 支持的MLIPS

1. **eSEN** (fairchem)：需要安装 `fairchem-core`
2. **MACE**：需要安装 `mace-torch`
3. 其他ASE兼容的MLIPS计算器

## 安装依赖

```bash
# 安装 fairchem (for eSEN)
pip install fairchem-core

# 或者安装 MACE
pip install mace-torch

# 或者两者都安装
pip install fairchem-core mace-torch
```

## 基本使用

### 1. 导入必要的模块

```python
from pymatgen.io.vasp import Poscar
from auto_kappa.vasp.relax import MLIPSStrictRelaxation
```

### 2. 设置MLIPS计算器

#### 使用eSEN:
```python
from fairchem.core.models.model_registry import model_name_to_local_file
from fairchem.core.common.relaxation.ase_utils import OCPCalculator

checkpoint_path = model_name_to_local_file('eSCN-S2EF-force-only')
calc = OCPCalculator(checkpoint_path=checkpoint_path)
calc_type = "esen"
```

#### 使用MACE:
```python
from mace.calculators import mace_mp

calc = mace_mp(model="medium", dispersion=False, default_dtype="float32")
calc_type = "mace"
```

### 3. 加载结构并运行弛豫

```python
# 加载结构
struct_init = Poscar.from_file("POSCAR").structure

# 创建弛豫实例
relax = MLIPSStrictRelaxation(
    struct_init, 
    calc, 
    calc_type=calc_type,
    outdir="./mlips_volume_relaxation"
)

# 运行弛豫
Vs, Es = relax.with_different_volumes(
    initial_strain_range=[-0.03, 0.05],  # -3% to +5% volume strain
    nstrains=15,                         # 15 strain points
    fmax=0.01,                          # 0.01 eV/Å force convergence
    maxstep=0.2,                        # 0.2 Å maximum step
    max_steps=500                       # 500 max optimization steps
)
```

### 4. 分析结果

```python
# 绘制 Birch-Murnaghan 拟合结果
relax.plot_bm(figname='bm_fit_mlips.png')

# 打印结果摘要
relax.print_results()

# 获取最优结构
optimal_structure = relax.get_optimal_structure()
optimal_structure.to(filename="POSCAR.optimal")

# 获取拟合误差
mae = relax.get_fitting_error()
print(f"Mean Absolute Error: {mae:.5f} eV")
```

## 完整示例

运行提供的示例脚本：

```bash
cd /path/to/auto-kappa
python examples/mlips_strict_relaxation_example.py
```

确保当前目录下有 `POSCAR` 文件。

## 输出文件

运行后会生成以下文件：

```
mlips_volume_relaxation/
├── 1/                      # 第一个应变点的计算
│   ├── CONTCAR            # 弛豫后的结构
│   ├── energy.dat         # 总能量
│   ├── strain.yaml        # 应变信息
│   └── opt.log           # 优化日志
├── 2/                      # 第二个应变点...
├── ...
├── volume_energy.csv      # 所有体积-能量数据
├── result.yaml           # 拟合结果摘要
├── POSCAR.init          # 初始结构
├── POSCAR.opt           # 最优结构
└── fig_bm_mlips.png     # Birch-Murnaghan拟合图
```

## 参数说明

### 初始化参数
- `initial_structure`: 初始结构 (pymatgen Structure 或 ASE Atoms)
- `mlips_calculator`: MLIPS计算器对象
- `calc_type`: 计算器类型 ("esen", "mace", 等)
- `dim`: 维度 (2 for 2D, 3 for 3D)
- `outdir`: 输出目录

### 弛豫参数
- `initial_strain_range`: 初始应变范围 [最小, 最大]
- `nstrains`: 应变点数量
- `fmax`: 力收敛标准 (eV/Å)
- `maxstep`: 最大步长 (Å)
- `max_steps`: 最大优化步数

## 与VASP StrictRelaxation的区别

| 特性 | StrictRelaxation (VASP) | MLIPSStrictRelaxation |
|------|-------------------------|----------------------|
| 计算器 | VASP | MLIPS (eSEN, MACE, 等) |
| 速度 | 慢 | 快 |
| 精度 | 高 (DFT精度) | 中等 (MLIPS精度) |
| 计算成本 | 高 | 低 |
| 适用性 | 所有材料 | 训练集覆盖的材料 |

## 注意事项

1. **MLIPS精度**：MLIPS的精度取决于训练数据，对于未见过的化学环境可能不准确
2. **能量零点**：不同MLIPS模型的能量零点不同，只关注相对能量差
3. **模型选择**：选择合适的MLIPS模型很重要，建议先在小体系上测试
4. **收敛标准**：MLIPS计算通常可以使用更严格的收敛标准

## 故障排除

1. **计算器导入错误**：确保安装了相应的MLIPS包
2. **内存不足**：大体系可能需要调整批处理大小
3. **优化不收敛**：调整 `fmax`, `maxstep`, `max_steps` 参数
4. **拟合质量差**：增加应变点数量或调整应变范围

## 支持

如有问题请联系开发者或查看auto-kappa文档。