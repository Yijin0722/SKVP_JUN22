# 第一阶段：相同双原子的交换对称化说明

这个文件夹是从原代码复制出来以后修改的版本，用来实现相同双原子分子体系的第一阶段 monomer exchange symmetry。

父目录里的原始源文件没有被修改。原始代码快照在：

```text
../exchange_phase1_original
```

修改后的代码在：

```text
../exchange_phase1_modified
```

## 这一阶段的目标

第一阶段只验证一件事：把原来的 ordered channel basis 投影到 monomer exchange symmetric basis 以后，S matrix 的结果是否和原来一致。

这一阶段暂时不处理 cross section，也不把 cross section 的输出逻辑改成对称化版本。也就是说，这一步的目标是：

```text
先保证新的 channel basis 和 S matrix 求解是对的
```

而不是一步到位把所有输出和截面计算都改完。

## 新 basis 的定义

原来的 channel 是 ordered form：

```text
|a,b>
```

其中

```text
a = (j1,k1)
b = (j2,k2)
```

对于两个相同的 diatom，交换两个 monomer 后，`|a,b>` 和 `|b,a>` 描述的是同一个物理交换对。因此新的 symmetric channel 写成：

```text
|a,b>_sym = (|a,b> + |b,a>) / sqrt(2)
```

如果 `a == b`，也就是交换以后还是自己，那么这个 symmetric channel 只有一个分量：

```text
|a,a>_sym = |a,a>
```

系数是 `1`，不是 `1/sqrt(2)`。

## 这一版代码实际做了什么

这一版还没有直接在小 basis 里面从头构造全部矩阵。它做的是：

```text
1. 按照原代码生成完整 ordered quant_mat
2. 按照 monomer exchange symmetry 建立 symmetric channel map
3. 用原代码路径构造 ordered basis 下的 M_* 矩阵
4. 把 ordered M_* 矩阵投影到 symmetric basis
5. 分别求解：
   ordered basis 的 S matrix
   symmetric basis 的 S matrix
6. 比较两个结果是否一致
```

所以这一阶段的意义是“验证数学变换正确”。它还不是最终的省计算量版本，因为 `mat_M`、`mat_M0`、`mat_M00`、`mat_M10` 仍然先按旧的大 ordered basis 构造了一遍。

真正减少 `M_*` 构造成本，需要下一阶段直接在 symmetric basis 中构造这些矩阵。

## 修改过的文件

### `module_skvp_AtomDiatom.f90`

新增了 exchange-symmetric basis 需要的全局数组和变量：

```text
nsym
n_sym_open
sym_quant_mat
sym_old_idx
sym_coeff
sym_ncomp
old_to_sym
sym_open_idx
old_open_pos_to_sym_open
```

这些变量负责记录：

```text
新的 symmetric channel 有多少个
每个 symmetric channel 对应旧 basis 里的哪一个或哪两个 ordered channels
每个分量的系数是多少
open channel 在新旧 basis 之间怎么对应
```

同时新增了投影后的矩阵和 S matrix：

```text
mat_M_sym
mat_M0_sym
mat_M00_sym
mat_M10_sym
Smat_sym
```

### `skvp_AtomDiatom_sym_phase1.f90`

这是这一阶段真正的新主程序。主要新增了几部分：

```text
build_exchange_symmetric_basis
```

根据原来的 `quant_mat(4,ncf)` 构造 exchange-symmetric channel basis。

```text
build_exchange_symmetric_open_basis
```

在 `solve_target_levels` 得到 `open_idx` 之后，构造 symmetric open-channel map。

```text
project_matrices_to_exchange_symmetric
```

把原来 ordered basis 下的 `M_*` 矩阵投影到 symmetric basis。

```text
SolveSMatrixGeneric
```

把原 `PhaseShift` 里面求解 S matrix 的线性代数过程抽出来，使 ordered basis 和 symmetric basis 可以共用同一个求解流程。

```text
write_exchange_phase1_outputs
```

输出这一阶段的 channel map 和 S matrix 对比结果。

### `Makefile`

新增了两个目标：

```text
make baseline
make sym_phase1
```

其中：

```text
baseline
```

编译原始主程序副本 `skvp_AtomDiatom_baseline.f90`。

```text
sym_phase1
```

编译新的第一阶段程序 `skvp_AtomDiatom_sym_phase1.f90`。

## 怎么运行

在 `exchange_phase1_modified` 文件夹里：

```sh
make baseline
./baseline > baseline_run.log 2>&1
```

然后运行修改版：

```sh
make sym_phase1
./sym_phase1 > sym_phase1_run.log 2>&1
```

当前测试使用的是这个文件夹里的：

```text
input.nml
```

也就是从当前 working tree 复制过来的那一组参数。

## 当前测试结果

当前 input 下得到：

```text
old ncf      = 121
new nsym     = 66
old n_open   = 121
new n_open   = 66
projected N  = 4488
max |dP|     = 5.3290705182E-15
```

其中：

```text
ncf
```

是旧 ordered basis 的 channel 数。

```text
nsym
```

是 exchange-symmetric basis 的 channel 数。

这一组参数下，channel 数从 `121` 降到 `66`。

S matrix 概率对比结果：

```text
sum_ordered_projected = 9.9014618385E-01
sum_sym               = 9.9014618385E-01
max_abs_diff          = 5.3290705182E-15
```

这个差值是数值舍入误差量级，说明 symmetric basis 的 S matrix 求解和 ordered basis 投影后的结果一致。

## 重要解释：为什么不能直接逐行比较旧的 `proba_all.dat`

旧代码里的 `proba_all.dat` 是 ordered channel basis 的概率：

```text
|S(old_i -> old_f)|^2
```

但是新的 symmetric final channel 可能是两个 ordered final channels 的线性组合：

```text
|f>_sym = c1 |old_f1> + c2 |old_f2>
```

所以正确比较方式不是逐行比较旧的 `proba_all.dat`，而是先比较振幅：

```text
S(old_i -> f_sym)
  = c1 S(old_i -> old_f1) + c2 S(old_i -> old_f2)
```

然后再取概率：

```text
|S(old_i -> f_sym)|^2
```

这一版代码输出的 `phase1_smatrix_compare.dat` 做的就是这个比较。

## 输出文件说明

```text
phase1_sym_channel_map.dat
```

记录每个 symmetric channel 对应旧 ordered basis 中的哪一个或哪两个 channel，以及对应系数。

```text
phase1_sym_open_map.dat
```

记录 symmetric open channel 和旧 open channel 的对应关系。

```text
phase1_smatrix_compare.dat
```

逐个 symmetric open channel 比较：

```text
ordered S matrix 投影后的概率
symmetric basis 直接求解得到的概率
二者差值
```

```text
phase1_smatrix_summary.dat
```

给出这一轮测试的总览，包括 `ncf`、`nsym`、`n_open`、`n_sym_open` 和最大概率差。

```text
baseline_run.log
```

原始 baseline 版本运行日志。

```text
sym_phase1_run.log
```

修改后的 phase-1 版本运行日志。

## 这一阶段还没有做的事

这一阶段还没有：

```text
1. 直接在 symmetric basis 中构造 M_* 矩阵
2. 把 CrossSection 改成 exchange-symmetric basis 的版本
3. 重新定义最终物理 cross section 的输出格式
4. 处理 antisymmetric exchange sector
```

现在只做 symmetric sector，因为当前目标是先解决相同双原子体系中最直接需要的对称化 S matrix 验证。

下一阶段如果要真正减少计算量，核心要改的是矩阵构造部分，也就是不要先生成完整 ordered `M_*`，而是在 `BAM_theta`、`M_V`、`mat_M`、`mat_M0`、`mat_M00`、`mat_M10` 的构造过程中直接进入 symmetric basis。
