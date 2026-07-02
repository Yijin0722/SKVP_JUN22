# Phase4：BAM_theta 与 open basic auxiliary 的对称化

本文件夹从 `exchange_phase3_bamtheta_sym` 复制而来，目标是在 direct-only 的相同双原子求解路径里继续减少旧 ordered-channel 矩阵。

本阶段只做两件事：

1. `BAM_theta_sym` 直接在 exchange-symmetric basis 中生成，不再先分配 `BAM_theta(n_pot,ncf,ncf)`。
2. `BAM_x*` 和 `BAM_xx*` 的 direct-only 求解路径改用 `n_sym_open` 尺寸的矩阵。

`BAM_r0 / BAM_r00 / BAM_r10` 本阶段暂时不改，留到下一步处理。

## 修改的文件

- `module_skvp_AtomDiatom.f90`
  - 新增 `sym_kvec(:)`。
  - 新增 `BAM_x1_sym / BAM_x3_sym / BAM_x4_sym`。
  - 新增 `BAM_xx1_sym / BAM_xx3_sym / BAM_xx4_sym`。
  - 新增 `BAM_xx1b_sym / BAM_xx3b_sym / BAM_xx4b_sym`。

- `sub_basic_aux_mat_calcul.f90`
  - 保留旧 `basic_aux_mat_calcul` 给 baseline/旧路径使用。
  - 新增 `basic_aux_mat_calcul_sym_direct` 给 direct-only 主程序使用。
  - 新 routine 仍生成旧 `kvec/open_idx` 所需的 target 信息，但 direct scattering 使用 `sym_kvec` 和 `n_sym_open` 尺寸的 open auxiliary 矩阵。
  - 新增输出 `basic_sym_summary.dat`，记录 `BAM_x*` 和 `BAM_xx*` 的尺寸缩减。

- `sub_potential_aux_mat_cacul.f90`
  - `potential_mat_calcul_sym_direct` 改为调用 `build_BAM_theta_sym_direct`。
  - 新增 `theta_ordered_value(ipot,n,n_prime)`，用于按需计算 ordered angular 元素。
  - 新增 `build_BAM_theta_sym_direct`，直接生成 `BAM_theta_sym(n_pot,nsym,nsym)`。
  - 旧 `build_BAM_thetas` 和 `build_BAM_theta_sym_from_old` 保留作旧路径/对照，不再被 direct-only 主路径调用。

- `skvp_AtomDiatom_sym_direct_only.f90`
  - 主流程改为调用 `basic_aux_mat_calcul_sym_direct`。
  - `build_exchange_symmetric_open_basis` 现在生成 `sym_kvec(n_sym_open)`。
  - `make_scatt_mat_sym_direct` 的 closed-open / open-open kinetic 与 asymptotic blocks 改用 `BAM_x*_sym` 和 `BAM_xx*_sym`。
  - 清理逻辑新增 sym auxiliary 矩阵与 `sym_kvec` 的 deallocate。

## 当前尺寸变化

使用当前 `input.nml` 测试：

- `ncf = 121`
- `nsym = 66`
- `n_open = 121`
- `n_sym_open = 66`
- `dim_x*nsym = 4488`

`bamtheta_sym_summary.dat`：

```text
old BAM_theta entries = 439230
sym BAM_theta entries = 130680
ratio = 0.29752066116
mode = direct
```

`basic_sym_summary.dat`：

```text
BAM_x*  old entries = 25410
BAM_x*  sym entries = 13860
ratio_x = 0.54545454545

BAM_xx* old entries = 87846
BAM_xx* sym entries = 26136
ratio_xx = 0.29752066116
```

## 验证结果

编译命令：

```bash
make sym_direct_only
```

运行命令：

```bash
./sym_direct_only > sym_direct_only_run.log 2>&1
```

运行结果：

- `direct_smatrix_summary.dat`：

```text
Energy = 0.49987339963 eV
Jtot = 0
ncf = 121
nsym = 66
n_open = 121
n_sym_open = 66
N_sym = 4488
prob_sum = 9.9014618385E-01
```

与 `exchange_phase3_bamtheta_sym` 对比：

- `direct_sym_channel_map.dat`：完全一致。
- `direct_sym_open_map.dat`：完全一致。
- `direct_smatrix_summary.dat`：完全一致。
- `proba_sym.dat` 最大绝对差：`0.0000000000000000E+00`。
- `proba_sym_all.dat` 最大绝对差：`0.0000000000000000E+00`。

## 还没有改的部分

为了把风险分开，本阶段没有改：

- `BAM_r0(n_pot,pb_nbr,n_open)`
- `BAM_r00(n_pot,n_open,n_open)`
- `BAM_r10(n_pot,n_open,n_open)`

因此 direct-only 路径现在仍然会保留旧 `open_idx/kvec`，因为 `BAM_r0/r00/r10` 还需要旧 open-channel 编号。下一阶段的目标就是把这三组 radial potential open-block 也换成：

- `BAM_r0_sym(n_pot,pb_nbr,n_sym_open)`
- `BAM_r00_sym(n_pot,n_sym_open,n_sym_open)`
- `BAM_r10_sym(n_pot,n_sym_open,n_sym_open)`

完成后，direct-only 主路径才可以进一步摆脱旧 `open_idx(n_open)` 对矩阵构造的依赖。
