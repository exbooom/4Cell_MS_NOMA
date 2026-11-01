# 分层投影价格／次梯度实验分支

更新：2026-11-01。代码位于 `F:/4Cell_MS_NOMA/Code`。本阶段完成计划中独立的次梯度更新实验；作者未公开的步长和细节仍按工程选择记录，不宣称恢复了作者程序。

## 实际调用链

```text
main
  -> experiment.runSingle
    -> optim.solveSubgradient
       -> optim.projectedDual (内层通信QoS和全对可听价格)
       -> 投影更新外层BS预算价格
       -> optim.dualSolution / diagnostics 独立证书
    -> BER -> 两/三项测距及DLL积分 -> 可听集合 -> 几何 -> 定位
  -> 原始数据、逐次更新CSV、8幅图与源码SHA256
```

算法ID为 `cp4a_subgradient_candidate`。新增乘子更新中不调用fmincon、solveCentral、solveDual、solveCp4a或dualActiveSet，也不借用参考求解器的乘子／功率作为起点。所有算法仍共享linprog可行性预检；“无价格优化器”不表示整个工程没有Optimization Toolbox依赖。

保留 `central`、`dual`、`cp4a_consistent_candidate` 三条独立参考路径，不通过回退覆盖本分支失败的输出。

## 投影更新与坐标

沿用阶段3负对偶D=-g。固定外层价格nu后，对内层的通信价格mu及有向可听价格beta，定义正的行坐标尺度：

```text
s_i = 1 / max(norm(A_s(i,:),2), 1e-12)
ell_i = s_i * x_i
Q = diag(s) * A_s
grad_x D = diag(s) * (b_s - A_s*z)
x_next = max(0, x - alpha * grad_x D)
```

因此在原缩放价格坐标中等价于 `ell_next=max(0,ell+alpha*s.^2.*(A_s*z-b_s))`。违反约束提高价格，宽松约束降低价格，投影保证非负。这是完整耦合矩阵下的修正实验，不使用论文式(24)越界分组或式(66)/(67)有疑点的打印导数。

内层功率依然来自有界标量Lagrangian的精确最小值，不是手动功率修补。全部受影响通信接收机和反向可听贡献通过A_s'×ell进入价格。

外层更新为 `nu_next=max(0,nu-alphaOuter*(budget_s-B_s*z))`。内层价格在外层求值间可warm start，默认启用；设warmStartInner=false时重新从零开始，测试覆盖两种方式。外层回溯产生的被拒绝试探仍完整保存。

## 两种步长模式

`cfg.numerics.subgradient.stepRule`：

- `backtracking`（默认）：内层采用投影梯度、可选FISTA动量及梯度重启，通过二次上界检查步长；外层对内层最优值使用投影梯度和Armijo回溯。默认accelerate=true，不等同于原文朴素subgradient。
- `diminishing`：内外层分别使用初始步长除以 `(t+1)^eta`，eta默认0.75且限定0.5<eta<=1；该模式不使用动量。它是显式的工程步长实验，不把有限迭代上限内的未收敛结果称为成功。

默认内层初始步长1e-3、外层0.1、回溯缩小0.5、增长1.15；这些值是本工程配置，不是作者b1/b2/b3。行归一化同时改变每行的有效步长，不能直接与论文原乘子坐标下的步长比较。

### 稳定的内层步长检查

直接相减O(1)对偶值会丢失接近最优点的微小改进。本实现使用精确的标量对偶Bregman余项，避免靠绝对的1e-14容差放行不合适的步长：

```text
w_next = fixedPrice + Q' * x_next
R = -sum((z_next-z_y) .* (w_next-c./(z_next.*z_y)))
accept if R <= norm(x_next-y)^2 / (2*alpha)
```

使用极小的相对浮点容差；动量通过梯度方向重启，避免用舍入后的对偶值差判断重启。该处理仅稳定计算，不改变原始功率目标或约束。

## 停止条件与可听性边界

内层必须同时满足：

1. 最大归一化可变约束违例<=2e-8。
2. `sum(abs(ell.*residual))/scaledRecoveredObjective<=2e-7`，避免正负互补项抵消。
3. 逐有向可听约束的相对违例<=1e-8，且不大于配置comparisonTolerance的十分之一。

最终还必须通过共享dualSolution的原始可行性、KKT、非负乘子、相对对偶间隙和相对可听性检查。共享diagnostics新增 `maxRelativeHearabilityViolation`：

```text
max(0, (C*p) ./ max(max(C,0)*p, realmin))
```

该比例与实际接收功率可听判定对应。归一化残差足够小并不保证弱链路相对误差足够小；本阶段完整回归曾发现1名用户因此只听到3个BS，已在根因处补充停止/证书门限，未放宽hearability中的comparisonTolerance。新增单元检查专门验证这一边界。

满足数值容差不等于精确零违例；很小的负原始-对偶间隙可能伴随容差内的正违例，仍完整记录而不归零。

## 状态、轨迹和接口

`projectedDual(problem,opts,variableRows,fixed,warm)` 返回与阶段3对偶state兼容的数据，另有converged/reason、iterations/evaluations、coordinateScale、innerRelativeComplementarity及maxRelativeHearabilityViolation。opts为cfg.numerics；variableRows默认全部约束，分层求解时仅为QoS与可听行。固定价格和warm向量为全部约束行长度。

`solveSubgradient(problem,opts,feasibility)` 返回通用solution、innerCalls/innerIterations、innerTrace、outerPriceHistory、stepConfiguration、usedPriceOptimizer=false、usedCentralFallback=false、literalPaperAlgorithm=false。第三参数可省略；调用方提供时复用可行性预检。

成功状态为converged_verified。未收敛明确标记subgradient_inner_iteration_limit、subgradient_inner_line_search_failed、subgradient_outer_iteration_limit或subgradient_outer_line_search_failed。保留最后已接受的功率、证书与失败轨迹；原问题不可行时沿用原不可行状态和NaN功率。不能把失败状态下的条件精度当成功算法结果。

导出新增：

- `subgradient_inner_calls.csv`：每次内层求值、预算价格、更新数、停止原因及残差，包含外层拒绝的试探。
- `subgradient_inner_history.csv`：按innerCall索引拼接的内层迭代采样记录，默认每100次更新采样，保留最终点。
- `subgradient_outer_prices.csv`：只记录已接受外层迭代的预算价格。
- `cp4a_subgradient_candidate_solver_history.csv`及通用功率/约束/测距/定位/DLL CSV。
- 图8：外层间隙/违例、预算价格、每次内层更新数和最终BS功率。

优化比较CSV新增相对可听违例、usesFmincon和stepRule，完整配置及未采样终态保存在run.mat。cfg/run版本为4.0；main默认6策略、8图。旧stage3验证显式固定原5策略。

## 复跑

```matlab
cd('F:/4Cell_MS_NOMA/Code');
addpath('tests');
run_subgradient_tests;  % 9项
run_stage4_checks;      % 共45项及3个完整场景
```

只运行新分支及集中式对照：

```matlab
cfg=msnoma.config.base;
cfg.experiment.algorithmIds={'central','cp4a_subgradient_candidate'};
cfg.io.showFigures=false;
main(cfg);
```

朴素递减步长与每轮重新初始化的实验：

```matlab
cfg=msnoma.config.base;
cfg.experiment.algorithmIds={'central','cp4a_subgradient_candidate'};
cfg.numerics.subgradient.stepRule='diminishing';
cfg.numerics.subgradient.accelerate=false;
cfg.numerics.subgradient.warmStartInner=false;
cfg.io.showFigures=false;
main(cfg); % 检查真实status；不能假定在配置迭代上限内达到验收容差
```

## 与原文仍不同的内容

固定oracle几何、无首次定位反馈；完整耦合代替打印分组；重新推导的正确导数；显式但非作者公开的步长和数值尺度；默认动量/回溯/warm start。原文逐行和作者参数的一致性仍未确认，但计划中的独立一阶更新实验链已建立并可与参考解比较。

下一项工作按交接计划进入PRS通信服务／干扰／TDoA与有限观测口径统一，随后做DLL近似来源分离、Monte Carlo，最后再接采样接收机与估计几何。
