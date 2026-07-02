# Direct-Only 交换对称 S 矩阵版本

这个文件夹是从 `exchange_phase2_modified` 复制出来的 direct-only 版本。

它把 phase2 里的两件事合成了一步：

```text
1. 不再跑 projected reference path
2. 直接把 symmetric open channels 当作 S matrix 的真实 channel
```

也就是说，这个版本的 S matrix 是：

```text
Smat_sym(1:n_sym_open, 1:n_sym_open)
```

而不是旧代码里的：

```text
Smat(1:n_open, 1:n_open)
```

## 和 phase2 验证版的区别

phase2 验证版做的是：

```text
old ordered M_* -> project -> S_projected
direct symmetric M_* -> S_direct
compare S_projected and S_direct
```

direct-only 版本只做：

```text
direct symmetric M_* -> Smat_sym
```

所以它不会再构造 old ordered `mat_M/mat_M0/mat_M00/mat_M10`，也不会再调用：

```text
project_matrices_to_exchange_symmetric
```

## 主要文件

### `skvp_AtomDiatom_sym_direct_only.f90`

新的 direct-only 主程序。

主流程是：

```text
CALL build_exchange_symmetric_basis
CALL solve_target_levels
CALL build_exchange_symmetric_open_basis
CALL basic_aux_mat_calcul
CALL potential_mat_calcul_sym_direct
CALL make_scatt_mat_sym_direct
CALL SolveSMatrixGeneric(..., Smat_sym, ...)
CALL write_direct_symmetric_outputs(Smat_sym)
```

这里 `Smat_sym` 的行列就是 symmetric open channel。

### `sub_potential_aux_mat_cacul.f90`

保留 phase2 新增的：

```text
potential_mat_calcul_sym_direct
```

这个子程序直接构造：

```text
mat_M_sym
mat_M0_sym
mat_M00_sym
mat_M10_sym
```

### `Makefile`

新增并默认使用：

```sh
make sym_direct_only
```

直接运行：

```sh
./sym_direct_only > sym_direct_only_run.log 2>&1
```

## S matrix channel 的定义

旧代码里，S matrix channel 由：

```text
open_idx(i)
quant_mat(:, open_idx(i))
```

定义。

direct-only 版本里，S matrix channel 由：

```text
sym_open_idx(i)
sym_quant_mat(:, sym_open_idx(i))
```

定义。

所以：

```text
Smat_sym(i,j)
```

表示：

```text
第 i 个 symmetric open channel -> 第 j 个 symmetric open channel
```

每个 symmetric channel 具体由哪些 old ordered channels 组成，可以看：

```text
direct_sym_open_map.dat
direct_sym_channel_map.dat
```

## 输出文件

### `proba_sym.dat`

类似旧代码的 `proba.dat`。

它写的是从初始 symmetric channel `(0,0,0,0)` 出发，到所有 symmetric open final channels 的概率：

```text
E_eV  P(sym_in -> sym_1)  P(sym_in -> sym_2) ...
```

初始 channel 不是硬编码成第 1 行，而是通过 old `(0,0,0,0)` 找到对应的 `sym_in_open`。

### `proba_sym_all.dat`

类似旧代码的 `proba_all.dat`，但 channel label 换成 symmetric channel。

每一行格式是：

```text
E_eV
Jtot
incoming_sym_open
incoming_sym_index
incoming_rep_j1 incoming_rep_k1 incoming_rep_j2 incoming_rep_k2
final_sym_open
final_sym_index
final_rep_j1 final_rep_k1 final_rep_j2 final_rep_k2
probability
```

其中 `rep` 是这个 symmetric channel 的代表 ordered label。完整组成要看：

```text
direct_sym_open_map.dat
```

例如一个 symmetric channel 可能对应：

```text
(0,0,2,0) 和 (2,0,0,0)
```

系数都是：

```text
1/sqrt(2)
```

### `direct_smatrix_summary.dat`

记录 direct-only 运行的基本信息：

```text
ncf
nsym
n_open
n_sym_open
N_sym
sym_in_open
prob_sum
```

## 当前测试结果

使用当前 `input.nml`，direct-only 运行结果为：

```text
ncf        = 121
nsym       = 66
n_open     = 121
n_sym_open = 66
N_sym      = 4488
sym_in     = 1
prob_sum   = 9.9014618385E-01
```

运行日志：

```text
sym_direct_only_run.log
```

本次运行时间约：

```text
0.611 min
```

作为对比，phase2 验证版因为同时跑 reference path 和 direct path，本组 input 下约为：

```text
1.207 min
```

## 仍然没有做的事

这个版本仍然没有做 cross section。

它也仍然保留：

```text
BAM_theta(n_pot,ncf,ncf)
```

所以后续如果继续优化，下一步可以考虑把 angular potential coupling 也改成 symmetric basis：

```text
BAM_theta_sym(n_pot,nsym,nsym)
```

或者直接在 potential matrix loop 中即时计算 angular coupling，不再保存 old `BAM_theta`。
