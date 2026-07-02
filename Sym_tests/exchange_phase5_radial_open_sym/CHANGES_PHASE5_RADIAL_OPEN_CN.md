# Phase5：radial potential open-block 的对称化

本文件夹从 `exchange_phase4_theta_open_sym` 复制而来。本阶段把 direct-only 路径中的三组 radial potential open-block 改成 `n_sym_open` 尺寸：

- `BAM_r0`
- `BAM_r00`
- `BAM_r10`

注意：本阶段没有新增 `BAM_r0_sym` 这类新名字，而是让 direct-only 路径中的原变量 `BAM_r0 / BAM_r00 / BAM_r10` 直接按 symmetric open channel 尺寸分配。

## 核心想法

对同一个 exchange-symmetric open channel，其两个 ordered component 只是单体交换：

```text
(j1,k1,j2,k2) <-> (j2,k2,j1,k1)
```

这两个 component 的 rotational threshold 相同，因此 radial wave number `k` 相同。radial potential integral 只依赖 radial coordinate、PES radial coefficient `A_cache` 和这个 open-channel `k`，不依赖 ordered component 的标签顺序。

因此 open radial block 可以先在 `n_sym_open` 上计算：

```text
BAM_r0(ipot,i,os)
BAM_r00(ipot,os,ot)
BAM_r10(ipot,os,ot)
```

然后 angular 部分直接乘已经对称化的：

```text
BAM_theta_sym(ipot,s,t)
```

这样不需要再对 `old_open_a / old_open_b` 做 component sum。

## 修改的文件

- `sub_potential_aux_mat_cacul.f90`
  - 新增 `build_BAM_rs_sym_direct`。
  - direct-only 路径的 `potential_mat_calcul_sym_direct` 改为调用 `build_BAM_rs_sym_direct`。
  - `BAM_r0` direct 分配尺寸从：

    ```fortran
    BAM_r0(n_pot,pb_nbr,n_open)
    ```

    改为：

    ```fortran
    BAM_r0(n_pot,pb_nbr,n_sym_open)
    ```

  - `BAM_r00 / BAM_r10` direct 分配尺寸从：

    ```fortran
    BAM_r00(n_pot,n_open,n_open)
    BAM_r10(n_pot,n_open,n_open)
    ```

    改为：

    ```fortran
    BAM_r00(n_pot,n_sym_open,n_sym_open)
    BAM_r10(n_pot,n_sym_open,n_sym_open)
    ```

  - closed-open potential block 现在使用：

    ```fortran
    BAM_r0(ipot,ir+1,os) * BAM_theta_sym(ipot,s,t)
    ```

  - open-open potential block 现在使用：

    ```fortran
    BAM_r00(ipot,os,ot) * BAM_theta_sym(ipot,s,t)
    BAM_r10(ipot,os,ot) * BAM_theta_sym(ipot,s,t)
    ```

  - 新增输出 `radial_sym_summary.dat`。

旧的 `build_BAM_rs` 保留给旧 `potential_mat_calcul` / baseline 路径使用；direct-only 主路径不再调用它。

## 当前尺寸变化

使用当前 `input.nml` 测试：

- `n_open = 121`
- `n_sym_open = 66`
- `n_pot = 30`
- `pb_nbr = 70`

`radial_sym_summary.dat`：

```text
BAM_r0:
old entries = 254100
sym entries = 138600
ratio = 0.54545454545

BAM_r00 + BAM_r10:
old entries = 878460
sym entries = 261360
ratio = 0.29752066116
```

结合前一阶段：

- `BAM_theta` 已经从 `ncf*ncf` 变成 `nsym*nsym`
- `BAM_x*` 已经从 `n_open` 变成 `n_sym_open`
- `BAM_xx*` 已经从 `n_open*n_open` 变成 `n_sym_open*n_sym_open`
- 本阶段 `BAM_r0/r00/r10` 也已经从 `n_open` 体系变成 `n_sym_open` 体系

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

与 `exchange_phase4_theta_open_sym` 对比：

- `direct_sym_channel_map.dat`：完全一致。
- `direct_sym_open_map.dat`：完全一致。
- `direct_smatrix_summary.dat` 最大绝对差：`0.0000000000000000E+00`。
- `proba_sym.dat` 最大绝对差：`0.0000000000000000E+00`。
- `proba_sym_all.dat` 最大绝对差：`0.0000000000000000E+00`。

本阶段运行时间：

```text
execution time: 0.214 min
```

上一阶段 phase4 同一 input 约为：

```text
execution time: 0.591 min
```

这个时间下降主要来自 radial open-block 从 `n_open` 缩到 `n_sym_open`。

## 仍然保留的旧结构

当前 direct-only 主路径中的大矩阵已经基本都是 symmetric 尺寸。不过以下 bookkeeping 仍然保留：

- `quant_mat(4,ncf)`
- `open_idx(n_open)`
- `kvec(n_open)`
- `old_to_sym(ncf)`
- `old_open_pos_to_sym_open(n_open)`

这些现在主要用于：

- 构造 symmetric basis 的 map
- 输出 channel map
- 找初始态 `(0,0,0,0)` 对应的 symmetric open row
- 保留和旧 ordered basis 的可验证对应关系

也就是说，本阶段已经解决了 `BAM_r0/r00/r10` 的矩阵尺寸问题；剩下如果还要继续瘦身，就是把这些 bookkeeping 也改成不依赖完整 ordered `ncf/n_open` 的生成方式。
