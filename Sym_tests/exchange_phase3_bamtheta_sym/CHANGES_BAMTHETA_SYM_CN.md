# BAM_theta_sym 验证版说明

这个文件夹从 `exchange_phase2_direct_only` 复制而来，用来验证第一步 `BAM_theta` 的 exchange-symmetric 化。

当前目标不是立刻删掉 old `BAM_theta`，而是先验证：

```text
old BAM_theta(n_pot,ncf,ncf)
投影得到
BAM_theta_sym(n_pot,nsym,nsym)
```

然后在 direct symmetric potential matrix 的 closed-closed block 中使用 `BAM_theta_sym`，确认 S matrix 和上一版 direct-only 结果一致。

## 新增内容

### `module_skvp_AtomDiatom.f90`

新增：

```fortran
REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: BAM_theta_sym
```

### `sub_potential_aux_mat_cacul.f90`

新增：

```fortran
SUBROUTINE build_BAM_theta_sym_from_old
```

定义为：

```text
BAM_theta_sym(ipot,s,t)
  = sum_a sum_b c(a,s) c(b,t) BAM_theta(ipot,old_a,old_b)
```

其中：

```text
s,t     是 symmetric channel index
a,b     是每个 symmetric channel 中的 old ordered component
c(a,s)  是 component coefficient
```

### `potential_mat_calcul_sym_direct`

closed-closed potential block 从原来的：

```fortran
sum_a sum_b coeff * BAM_r(ipot,ir+1,jr+1) * BAM_theta(ipot,old_a,old_b)
```

改成：

```fortran
BAM_r(ipot,ir+1,jr+1) * BAM_theta_sym(ipot,s,t)
```

这是安全的，因为 closed-closed block 里的 radial factor `BAM_r(ipot,ir+1,jr+1)` 不依赖 old channel component `old_a/old_b`。

## 暂时没有直接替换的部分

closed-open 和 open-open potential block 目前仍然保留 old component 求和：

```fortran
BAM_r0(ipot,ir+1,old_open_b) * BAM_theta(ipot,old_a,old_b)
BAM_r00(ipot,old_open_a,old_open_b) * BAM_theta(ipot,old_a,old_b)
BAM_r10(ipot,old_open_a,old_open_b) * BAM_theta(ipot,old_a,old_b)
```

原因是这些 radial block 依赖 old open channel：

```text
old_open_b
old_open_a
```

不能只用一个 `BAM_theta_sym(s,t)` 直接替换。下一步如果要继续优化，需要同时定义：

```text
BAM_r0_sym
BAM_r00_sym
BAM_r10_sym
```

或者证明同一个 symmetric open channel 内的 old components 有完全相同的 radial factor 后，再做更强的合并。

## 当前测试结果

使用当前 `input.nml` 运行：

```sh
make sym_direct_only
./sym_direct_only > sym_direct_only_run.log 2>&1
```

结果：

```text
n_pot      = 30
ncf        = 121
nsym       = 66
old entries = 439230
sym entries = 130680
ratio       = 0.2975206612
```

也就是说，angular matrix 如果完全改成 symmetric basis，存储量约为 old `BAM_theta` 的 `29.75%`。

S matrix 概率输出和上一版 `exchange_phase2_direct_only` 对比：

```text
max_abs_proba_sym_diff      = 0.0
max_abs_probability_all_diff = 0.0
```

当前 summary：

```text
prob_sum = 9.9014618385E-01
```

和 direct-only 版本一致到输出精度。

## 输出文件

```text
bamtheta_sym_summary.dat
```

记录 `BAM_theta` 和 `BAM_theta_sym` 的尺寸对比。

```text
direct_smatrix_summary.dat
```

记录 S matrix direct-only 结果。

```text
proba_sym.dat
proba_sym_all.dat
```

仍然是 symmetric channel 作为真实 S-matrix channel 的概率输出。

## 下一步建议

下一步有两个选择：

```text
1. 继续验证 open radial block 是否能做 symmetric 化
2. 直接写 build_BAM_theta_sym_direct，不再先构造 old BAM_theta
```

我建议先做第 2 个：

```text
build_BAM_theta_sym_direct
```

因为 closed-closed block 已经证明 `BAM_theta_sym` 本身没问题。只要 direct 构造的 `BAM_theta_sym` 和 projected `BAM_theta_sym` 一致，就可以删掉 old `BAM_theta(n_pot,ncf,ncf)` 的构造，从而真正减少 angular matrix 的内存。
