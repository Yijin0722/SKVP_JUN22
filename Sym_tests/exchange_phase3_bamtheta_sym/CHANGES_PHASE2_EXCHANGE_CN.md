# 第二阶段：直接构造交换对称 basis 下的 M 矩阵

这个文件夹是从第一阶段的 `exchange_phase1_modified` 复制出来的第二阶段版本。

第一阶段做的是：

```text
先构造完整 ordered-basis M_* 矩阵
再投影到 exchange-symmetric basis
```

第二阶段新增的是：

```text
直接在 exchange-symmetric basis 里构造 M_* 矩阵
```

这样后面真正用于 S matrix 求解的矩阵尺寸可以从旧的：

```text
dim_x * ncf
```

变成新的：

```text
dim_x * nsym
```

## 当前实现状态

当前 phase2 仍然保留两条路径：

```text
1. projected reference path
   old ordered M_* -> project -> S_projected

2. direct symmetric path
   directly build symmetric M_* -> S_direct
```

保留第一条路径是为了验证第二条路径是否正确。等 direct path 稳定后，可以把 projected reference path 关掉，那时才是实际省时间和省内存的生产版本。

## 这一阶段改了哪些文件

### `Makefile`

新增目标：

```sh
make sym_phase2
```

默认目标也改成：

```text
sym_phase2 baseline
```

也就是说直接 `make` 会优先编译 phase2 版本。

### `skvp_AtomDiatom_sym_phase2.f90`

这是 phase2 的新主程序，由 `skvp_AtomDiatom_sym_phase1.f90` 复制而来。

主循环里现在的逻辑是：

```text
CALL basic_aux_mat_calcul

! reference path
CALL potential_mat_calcul
CALL make_scatt_mat
CALL project_matrices_to_exchange_symmetric
CALL SolveSMatrixGeneric(... projected symmetric M_* ...)

! direct path
CALL potential_mat_calcul_sym_direct
CALL make_scatt_mat_sym_direct
CALL SolveSMatrixGeneric(... direct symmetric M_* ...)

CALL write_exchange_phase2_outputs
```

新增的重要子程序：

```text
make_scatt_mat_sym_direct
```

它把 kinetic、threshold、rotational 和 `Wdd` 项直接累加到：

```text
mat_M_sym
mat_M0_sym
mat_M00_sym
mat_M10_sym
```

而不是先构造旧的：

```text
mat_M
mat_M0
mat_M00
mat_M10
```

再投影。

```text
cleanup_phase2_matrix_storage
```

用于在 projected reference path 完成后释放旧的大矩阵和中间 potential arrays，然后再运行 direct path。

```text
write_exchange_phase2_outputs
```

输出 `S_projected` 和 `S_direct` 的比较。

### `sub_potential_aux_mat_cacul.f90`

新增：

```text
potential_mat_calcul_sym_direct
```

它直接构造 symmetric basis 下的 potential block：

```text
mat_M_sym
mat_M0_sym
mat_M00_sym
mat_M10_sym
```

这里用的数学关系还是：

```text
M_sym(s,t) = sum_a sum_b c(a,s) c(b,t) M_old(a,b)
```

区别是 phase2 不再真的生成完整的 old `M_V/M0_V/M00_V/M10_V`，而是在循环里直接把 old component pair 的贡献加到 symmetric matrix 里。

## 当前还没有优化的部分

这一阶段仍然会构造：

```text
BAM_theta(n_pot, ncf, ncf)
```

所以 angular potential coupling 的这一层还没有降到：

```text
BAM_theta_sym(n_pot, nsym, nsym)
```

也就是说，phase2 已经把最终进入 S-matrix solver 的 `M_*` 矩阵改小了，但还没有把 `BAM_theta` 这一层也改小。

下一步如果继续优化，目标应该是：

```text
直接构造 BAM_theta_sym
```

或者在 potential matrix loop 中完全不保存 old `BAM_theta`，而是即时计算并累加到 symmetric matrix。

## 当前测试结果

使用当前文件夹里的 `input.nml`，运行：

```sh
make sym_phase2
./sym_phase2 > sym_phase2_run.log 2>&1
```

得到：

```text
ncf        = 121
nsym       = 66
n_open     = 121
n_sym_open = 66
N_sym      = 4488
```

direct path 和 projected reference path 的比较：

```text
sum_projected = 9.9014618385E-01
sum_direct    = 9.9014618385E-01
max_abs_dP    = 1.3322676296E-15
max_abs_dS    = 1.0547118734E-14
```

这里：

```text
max_abs_dP
```

是从初始 `(0,0,0,0)` 到所有 symmetric open final channels 的概率最大差。

```text
max_abs_dS
```

是整个 symmetric open-channel S matrix 的复数矩阵元最大差。

这两个差值都在数值舍入误差范围内，说明 direct symmetric M_* 构造和 phase1 的 projected reference 是一致的。

## 输出文件

```text
phase2_sym_channel_map.dat
```

记录 symmetric channel 和 ordered channel components 的对应关系。

```text
phase2_sym_open_map.dat
```

记录 symmetric open channel 和 old open channel components 的对应关系。

```text
phase2_smatrix_compare.dat
```

逐个 symmetric final open channel 比较：

```text
P_projected
P_direct
abs_diff
```

```text
phase2_smatrix_summary.dat
```

记录本次测试的总览，包括 `ncf`、`nsym`、`n_open`、`n_sym_open`、`N_sym`、`max_abs_dP` 和 `max_abs_dS`。

```text
sym_phase2_run.log
```

phase2 运行日志。

## 和 phase1 的区别

phase1 验证的是：

```text
old ordered basis S matrix
投影到 symmetric basis 后
和 symmetric projected M_* 求解是否一致
```

phase2 验证的是：

```text
直接构造 symmetric M_*
和 old ordered M_* 投影后是否一致
```

所以 phase2 已经进入了真正减少 `M_*` 矩阵尺寸的阶段。

当前为了验证，程序还同时跑 reference path。等确认更多 input 也稳定后，可以做一个 `sym_phase2_direct_only` 或 runtime flag，把 reference path 关掉。
