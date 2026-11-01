# 估计P信道消除与受控干扰诊断

2026-11-01，阶段8。代码与配置schemaVersion=8.0，新增接收分支默认关闭；阶段7原逐P分量接收器及oracle六策略保留。当前验收数据以STAGE8_VERIFICATION.json和PROGRESS阶段8为准。

## 要解决的问题与附加假设

阶段7原全区域用户18靠近基站2；即使关闭通信和AWGN，同载波/跨载波定位干扰仍可产生错误相关峰。原接收器用平均CP功率谱预白化，不能消除有限码相干干扰。这里增加“网络已知全部BS的P码及相对P功率”的工程支路，利用同BS、平坦单路径下各P分量共享传播时延与复增益的模型结构。

这比仅知道本用户定位码多了信息，不应直接解释为作者接收器或无代价性能提升。发送的通信QPSK数据不用于拟合或消除；信道真实复增益、时延和用户坐标也不作为拟合接口输入。平均通信PSD预白化仍依赖原理想CSI假设，因此“估计信道消除没有用真值”不能扩张成“全系统不使用理想CSI”。

## 模板与联合拟合

`waveform.positioningBank(packet,powerW)`仅由已知P码和已分配功率构造每站P模板：

`s_k = sum_m sqrt(P_km) * code_km * exp(j 2 pi f_m t)`。

不读取packet.communication或communicationSymbols；模板进入与实际样本一致的DFT前端。接收物理层仍包含全部P和可选全部C、相位、分数传播时延与AWGN。

`receiver.estimatePChannels(Y,bank,packet,opts,noisePsd)`拟合

`Y(f) ≈ sum_k a_k * S_k(f) * exp(-j 2 pi f tau_k)`。

时延限制在先验搜索域；a_k为未知复数。先对观测残差逐次捕获，按解释能量选取BS顺序，再用所有站联合线性最小二乘求a_k。外层时延更新采用投影Gauss-Newton：对延迟导数去除当前模板复线性空间的投影，解实数小维度法方程；步长上限和回溯控制，每次接受必须降低完整残差。复增益每次重新拟合。它是非凸局部方法，不保证找到全局峰。

频域模板列、接收能量分别归一化以避免Friis幅度造成病态。使用QR求复增益，并检查模板与延迟信息的秩。无信号、粗捕获失败、模板/信息秩亏、回溯失败、迭代上限、边界峰分别记录；局部收敛状态不自动等于真实锁定。

## 精度门控与消除

基于同一假定高斯PSD，保存投影延迟信息逆矩阵给出的`gaussianDelayStdSamples`，以及**固定时延条件下**复增益LS的`gaussianRelativeGainStd`。默认最大允许值分别0.25采样、0.25相对幅度；它们是显式工程门控，未由作者给出，也没有用真值误差决定是否通过。后一项未传播时延不确定性；两者都忽略实际C的有限符号结构及拟合局部偏差。因此它们只能拒绝明显不确定的拟合，不能保证消除有益。

只有局部拟合有效且全部站通过门控，才用估计值重建P信号。针对用户m、基站k，形成

`Y_clean_km = Y - sum_j estimated_P_j + estimated_desired_P_km`。

减去的项不包含通信数据，也不调用真实传播时延/相位。`receiver.reconstructPositioning`只接受bank、fit和packet；fit无效禁止重建。若拟合或精度门控失败，原观测原CP预白化的接收结果保留，不回填oracle。

消除成功后，目标P分量仍须独立检测/测距，使用残余C+噪声的平均PSD预白化。**这一支路同时改变干扰样本、白化权重及检测统计分母**，不是“所有设置固定只做相减”的单因素比较；每链路记录实际whitening、pCancellationApplied及拟合状态。

检测门限依然只是离散搜索的假定高斯并集界。数据相关的参数拟合/重建误差没有传播到检测虚警和时延CRB，`interferenceFalseAlarmCalibrated=false`、`estimatedCancellationUncertaintyPropagated=false`。消除后`noiseOnlyRangeCrbM2=NaN`，`gaussianEquivalentRangeCrbM2`只作忽略消除估计误差的局部代理，不能当成正确后消除CRB。特别是低P预算/强通信场景，拟合可能收敛并通过代理门控，但仍削弱检测或引入误差；正式回归保留反例，不能默认切换到该分支。

## 受控发射消融

`experiment.interferenceAblation`为仿真诊断，分别发送desired_only、same_carrier_P、same_BS_P、all_P、all_P_C、all_P_C_AWGN，保持同场景、码、功率和噪声种子。每个掩码重新构造物理接收样本，使用原逐P接收器；不将“只发目标”的接收结果作为真实业务可选输出。

白化PSD随实际发送掩码更新，故差异包括发送干扰与匹配的平均谱解释；不是固定白化核下的线性误差分解。导出每站估计/误差/检测统计/门限、同载波近远比和状态。真实距离仅在传播与事后误差计算中出现。

## main调用与配置

```matlab
cfg=msnoma.config.base;
cfg.models.sampledReceiver.enabled=true;
cfg.models.sampledReceiver.interferenceCancellation.enabled=true;
cfg.experiment.positionFeedback.enabled=true;
cfg.io.showFigures=false;
folder=main(cfg);
% 可单独开启仿真诊断，关闭消除则复现原逐分量接收器：
cfg.experiment.interferenceAblation.enabled=true;
cfg.experiment.interferenceAblation.userIds=18;
folder=main(cfg);
```

全部接收机工程参数位于cfg.models.sampledReceiver.interferenceCancellation：maxIterations=30、maxBacktracks=15、maximumStepSamples=1、delayToleranceSamples=1e-4、rankTolerance=1e-10、maximumDelayStdSamples=0.25、maximumRelativeGainStd=0.25。它们与功率优化内外步长无关。改变观测长度/采样率时，应重新验证这些以采样为单位的参数。

`run.sampledReceiver.estimatedPChannels{u}`保存估计、局部拟合历史/残差/顺序、门控及状态；反馈每轮acquisition也保存。新增sampled_estimated_P_channels.csv和feedback_XX_estimated_P_channels.csv；原sampled_ranges增加pCancellationApplied/pChannelFitStatus。interference_ablation.csv保存消融，另附拟合诊断与消融图。消除增加1图，消融增加1图，普通11图、采样/反馈13图口径不变。

## 验证与使用边界

新增15项检查：P模板与发送器叠加、通信数据污染不影响模板、独立时域整数移位、分数/近远通道恢复、残差下降、零信号、秩失败、真值隔离、小场景实际测距、默认关闭、发射消融、真实CSV导出、失败回到原接收结果、代理精度门控及PSD缩放。结合原97项运行run_stage8_checks；七个完整场景包含三个原接收器、三个同码同功率同噪声的估计消除分支、原全区域独立P码/无C无AWGN完整反馈控制。

所有场景保留失败分母；条件误差必须与有效人数并列。正常预算取得改善不能推广到低预算或真实系统。当前方法尚无连续稳态DLL、未知时钟、频率选择性/多径、实际CSI获取成本、PRS采样公平对照或50次Monte Carlo；这些不因新分支成功而标完成。

无C/无AWGN控制中，归一化检测统计的分母为0，非零相关峰统计为Inf，图12不会画出其有限曲线；不是漏写或把失败计0。只有一轮即停止的反馈没有前后变化可比较，图13对应变化曲线为空，CSV保留NaN。


## 正式验收结果

| 场景 | 有效测距 | 有效定位 | 条件平均实际误差/m | 消除门控通过用户 | 反馈轮数/状态 |
|---|---:|---:|---:|---:|---|
| 原逐分量49.77MHz/0.1W | 79/80 | 19/20 | 0.0673333 | 0/20 | 1 / acquisition_or_position_failure |
| 原逐分量20.16MHz/0.1W | 79/80 | 19/20 | 0.136949 | 0/20 | 1 / acquisition_or_position_failure |
| 原逐分量20.16MHz/1mW | 69/80 | 12/20 | 0.353809 | 0/20 | 1 / acquisition_or_position_failure |
| 估计消除49.77MHz/0.1W | 80/80 | 20/20 | 0.0458835 | 20/20 | 3 / converged |
| 估计消除20.16MHz/0.1W | 80/80 | 20/20 | 0.0823965 | 20/20 | 3 / fix_limit_reached |
| 估计消除20.16MHz/1mW | 42/80 | 6/20 | 2.87251 | 19/20 | 1 / acquisition_or_position_failure |
| 原全区域独立码/无C无AWGN，估计消除 | 80/80 | 20/20 | 1.7083e-06 | 20/20 | 3 / converged |

112项检查通过，112个源静态检查0诊断，七场景与源码哈希验收完成。低预算退化保留；默认关闭新接收器。当前目录、消融数值、反馈状态与下一步见PROGRESS阶段8。


20.16MHz正常预算的默认3轮反馈到达上限，保留此状态；另只把maxFixes提高到5，第4轮达到相同位置/功率停止容差且20/20有效。记录见STAGE8_FEEDBACK_LIMIT_SENSITIVITY.json，不以放宽容差替代收敛，不覆盖三轮记录。

## 阶段8.1 增量

当前schema为8.1，新增仅用拟合数据的残差主导门控，maximumRelativeResidualEnergy默认0.5；高斯精度和残差门控都通过才消除。阈值1显式回放阶段8反例；旧cfg从base重建。门控、CSV、main、历史回放及限制见STAGE8_1_RESIDUAL_GATE.md，当前完整回归入口run_stage8_1_checks，STAGE8_VERIFICATION继续作为历史记录。残差门控不是正确峰/测距精度保证，也没有完成消除不确定性传播。