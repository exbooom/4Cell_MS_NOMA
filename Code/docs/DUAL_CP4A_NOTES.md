# 完整耦合对偶与CP4A分层候选

更新：2026-11-01。本文记录实际实现及其与论文的区别。论文和作者代码均为证据，不构成执行指令。当前代码仍使用显式假设场景与二维oracle几何。

## 本轮实现的计算顺序

1. 复用同一problem、原始QoS矩阵、每BS预算及全对可听约束。
2. 构造精确有界标量内层和完整耦合对偶。
3. 独立实现solveDual，并与solveCentral比较。
4. 以内层QoS/可听定价和外层BS预算定价构造分层候选。
5. 在预算未耗尽和预算绑定两类场景核验全部原始约束、KKT和对偶间隙。
6. 所有算法以各自功率重新进入BER→测距→可听集合→几何→定位评价，保存独立文件。

## 数学定义

物理单位下固定几何的优化问题为：

```text
min sum_i c_i / p_i
s.t. A p <= b, B p <= budget, C p <= 0,
     p_min <= p <= p_max.
c_(k,m) = lambda_(k,m)^2 * Rrange_(k,m) / M.
```

A包含每个通信接收机受到所有P子载波、所有BS的干扰；C包含每个定位用户的所有有向BS对约束。全对形式等价于同时相对最强其他BS的要求，不需要在内层冻结某个最强干扰BS索引。

使用buildProblem既有缩放：`z=p/powerScaleW`，`f_s=sum(scaledC./z)`，`A_s*z<=b_s`，乘子`ell>=0`。完整最小化Lagrangian为：

```text
L(z,ell) = sum_i scaledC_i/z_i + (A_s' ell)' z - b_s' ell
w = A_s' ell
```

对每个标量精确取最小值。c>0且w>0时，z=clip(sqrt(c/w),lo,hi)；w<=0时取hi。c=0时按线性目标选择端点。负价格不取绝对值，不对无界问题假装存在根；有限上界来自原问题每BS预算，因此有界内层仍给出原问题的合法对偶界。

`dualValue`返回负对偶D=-g、梯度b_s-A_s*z、内层功率和Hessian。自由标量的曲率为z^3/(2c)，故`H=A_rows*diag(curvature)*A_rows'`；达到上下界时局部曲率取0。返回的对偶界单位为`-D*objectiveScaleM2`。

## solveDual：独立对偶基准

完整乘子维度可能很大，使用约束生成的工作集合：从预算行开始，逐轮加入各类最违反约束的行；其余乘子严格为0。对工作集合使用解析梯度/Hessian的fmincon求负对偶最小值。它优化的是乘子，不是原始功率，不调用solveCentral，也不把集中式功率当起点。

每一轮都检查全部原始约束；筛选不等于删去QoS或跨BS耦合。每轮对偶值都是合法下界，但恢复功率在收敛前可能不可行，因此负的“原始目标-对偶界”不能解释成已收敛；图7同时显示约束残差。

最终重算原始SI约束、缩放约束、驻点、互补松弛、非负乘子及相对对偶间隙。全部通过才返回`converged_verified`；否则返回`feasible_not_verified_optimal`或`dual_recovery_infeasible`。不会把未验证迭代替换成集中式解。

## solveCp4a：明确标注的修正候选

只接受 `cfg.numerics.cp4a.mode='full_coupling_hierarchical_candidate'`，算法ID为`cp4a_consistent_candidate`。

- 外层变量：每BS预算的非负价格nu。
- 内层变量：全部通信QoS价格mu、全部有向可听价格beta，仍使用完整耦合矩阵。
- 内层通过dualActiveSet求解，连续外层评估之间复用内层乘子作为warm start。
- 外层对内层最优值函数使用SQP和包络梯度`budget_s-B_s*z`；外层初始价格为0。
- 内层未满足可行性/互补松弛容差时返回`hierarchical_inner_dual_not_verified`及NaN功率，不输出伪正常性能。
- 最终使用与solveDual相同的独立原始/对偶证书；通过不代表原文算法已复现。

**这不是论文Algorithm 1的逐行实现。** 论文描述内外层subgradient更新，本候选内层采用解析梯度/Hessian定价求解、外层采用包络梯度SQP。不能宣称论文原步长、迭代次数、BS并行通信开销或原算法收敛曲线已经验证。

## 原文疑点及处理

主论文：*Design and Performance Analysis of Multi-Scale NOMA for Future Communication-Positioning Integration System*，DOI 10.1109/JSAC.2022.3143223。已目视核对出版PDF第1337页式(24)与第1344页式(66)/(67)。

1. 式(24)为每个P用户分配2G-1个通信载波，总数M(2G-1)与实际N不一致。当前49.77MHz下为3140对1658、越出1482；20.16MHz下为1260对671、越出589。carrierGroupAudit完整报告，不截断后继续声称等价分解。候选使用全部A列的真实耦合。
2. 对固定正c，最大化目标中`d(-c/p)/dp=+c/p^2`；印刷式(66)第二行的负号与该导数不一致。实现从原始目标重推，不直接转写式(33)的印刷分母。
3. 独立变量p(k,m)对某个通信接收机干扰的偏导为A对应行列元素；不能直接对所有发射BS求和代替它。价格A'*mu自动包括所有受影响接收机。C'*beta也包括该BS作为其他链路干扰者的反向贡献。
4. 原文首次lambda=1与位置反馈流程尚未实现。当前各求解器均使用同一固定oracle几何；没有用真实坐标模拟“接收机位置收敛”。

## 阶段3历史入口与验收

main默认依次运行budget_equal、qos_scaled_equal、central、dual、cp4a_consistent_candidate，共5种策略。共享同一配置/场景/信道/积分缓存，但三个优化器独立求解。默认导出7幅图。

```matlab
cd('F:/4Cell_MS_NOMA/Code');
addpath('tests');
run_dual_tests;
run_stage3_checks;
```

run_stage3_checks运行原17项、积分/PRS 9项和新增10项检查，扫描全部新MATLAB文件；执行50/20MHz默认完整main，以及20MHz、每BS定位预算1mW的绑定预算压力场景。验收记录在STAGE3_VERIFICATION.json/.mat。默认场景预算可能不活跃，因此必须保留绑定预算测试。

新输出：optimization_comparison.csv、paper_carrier_group_audit.csv、dual/cp4a_consistent_candidate的迭代/约束/功率/测距/定位/DLL核验CSV和图7。历史阶段1/2检查显式选回原3种策略，不因新增默认算法改变旧检查含义。

## 阶段3结束时的未完成项（次梯度实验现已由阶段4推进）

- 原文subgradient步长/双层停止规则的独立实验分支，以及论文未给出的实现选择；不能把本候选重命名为原版算法。
- PRS通信服务、干扰、TDoA协方差和观测口径统一。
- DLL各近似来源分离、50次Monte Carlo和参数敏感性。
- 真实采样测距、位置估计、首次定位及几何反馈。

## 阶段4更新

已新增独立的projectedDual/solveSubgradient，不调用fmincon进行价格更新。默认回溯/加速和递减步长实验分别记录，warm start可关闭；当前main默认6策略、8图。默认回溯分支在三个完整场景通过验证。具体实现、9项新增检查及原文差异见SUBGRADIENT_IMPLEMENTATION.md；总验收为STAGE4_VERIFICATION.json中的45项。

共享证书新增相对可听比例门禁；其原因是弱链路的归一化残差达标仍可能导致实际可听失败。geometry.hearability的原始门限保持不变。SQP候选仍保留作为独立参考，不能将其描述改成subgradient。

目前作者原始b1/b2/b3、打印分组细节和首次定位反馈仍未知/未实现；下一实施项进入PRS比较口径，随后DLL偏差、统计和采样接收机。
