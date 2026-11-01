# 4Cell MS-NOMA MATLAB 实现

当前状态（2026-11-01）：已跑通4基站、20定位用户的单场景理论计算路径，包含场景、跨基站信道、通信干扰/QoS、测距方差、几何、可听性与功率优化。阶段5/6接入相关TDoA、PRS服务审计、DLL八组合方差差异和C接收滤波解析核验；阶段7新增有限采样波形、实际捕获/测距、位置估计、首次lambda=1与估计几何反馈，以及可配置初始步长敏感性实验。阶段8接入可选全部P码库联合时延/复增益估计、估计干扰消除及发射消融，正常预算改善与低预算退化均保留。oracle基准单独保留，捕获或定位失败明确停止反馈。采用显式工程假设，不声称复现作者CP4A或论文全部数值。

**接手入口：`PROGRESS.md`最新阶段8.1。** 接口在 `docs/MODULE_INTERFACES.md`；当前验证见 `docs/STAGE8_1_VERIFICATION.json`，原阶段8验收保留，阶段1–7记录保留为历史源码快照。采样与反馈说明见 `docs/SAMPLED_RECEIVER_AND_FEEDBACK.md`，可选估计P消除见 `docs/ESTIMATED_P_CANCELLATION.md`；DLL解析核验、PRS比较及对偶差异分别见对应专题文档。`docs/ARCHITECTURE_PLAN.md`是原设计快照。

## 运行

在 MATLAB R2024b 命令窗口执行：

```matlab
cd('F:/4Cell_MS_NOMA/Code');
main
```

默认名义50MHz，按整数G得到有效49.77MHz。main显示六种功率策略、PRS服务和DLL分解/解析接收端审计，默认打开11个图窗，并保存MAT、CSV、PNG与清单到 `F:/4Cell_MS_NOMA/Results/run_*`。

修改带宽、可听性假设或关闭图窗：

```matlab
cfg = msnoma.config.base('smoke_assumed');
cfg.paper.requestedBandwidthHz = 20e6;  % 有效20.16MHz
cfg.models.hearability.ratio = 0.1;     % 显式复合q，不是作者已知Omega
cfg.io.showFigures = false;
runDir = main(cfg);
```

噪声的实际计算输入是 `cfg.assumed.noisePsdWPerHz`；默认值为k_B×290K。单独修改温度或噪声系数字段不会自动重算此PSD，改变噪声时须同时记录一致的PSD和来源。配置的原文值与工程假设分别位于paper、assumed和models。

检查模型与完整运行：

```matlab
addpath('tests');
run_tests;          % 17项数学检查
run_precision_tests; % 9项DLL积分与PRS检查，其中LTE图样检查需要LTE Toolbox
run_dual_tests;     % 10项对偶/分层数学和独立求解检查
run_subgradient_tests; % 9项投影/次梯度更新及物理证书检查
run_prs_comparison_tests; % 18项TDoA、双向PSD、通信服务、网格零点及失败检查
run_bias_tests;     % 17项DLL单因素/闭合、滤波独立积分、导出和失败检查
run_sampled_tests;  % 17项有限波形、真值隔离、packet复用、测距/位置/反馈和步长证书检查
run_cancellation_tests; % 15项P码库、联合估计、消除/失败/精度门控、消融与导出
run_stage8_checks;  % 112项检查、7个完整场景、接收器配对与源码哈希核验
run_stage7_checks;  % 原97项检查、5个场景，阶段8分支默认关闭
run_stage6_checks;  % 历史80项检查、3个完整场景、11图、阶段7模块默认关闭
run_stage5_checks;  % 历史63项检查、3个场景、9图，显式关闭阶段6附加模块
run_stage4_checks;  % 历史45项检查、3个完整场景、8图导出，显式关闭新PRS比较
run_stage3_checks;  % 历史五策略36项检查、静态检查、两带宽+绑定预算压力main、7图导出
run_stage2_checks;  % 历史三策略26项检查、静态检查、两带宽main、积分/能量/前端敏感性与6图导出
run_stage1_checks;  % 历史阶段1检查，主动关闭DLL积分/PRS，含重放和功率下限敏感性
```

PowerShell批处理：

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); c=msnoma.config.base; c.io.showFigures=false; main(c);"
```

需要Optimization Toolbox（linprog、fmincon）；使用MATLAB自带图形、表格和Java运行时。默认主计算路径不依赖LTE Toolbox、CVX或原始信号生成器；LTE单符号参考模式及完整阶段2测试中的图样检查需要LTE Toolbox。只需将Code加入路径，不要对整个工作区使用genpath。

## 六种策略

- `budget_equal`：每BS预算均分，可能违反通信QoS，完整保留其违例。
- `qos_scaled_equal`：在均分结构下用共同尺度满足全部通信QoS，不平衡近远效应。
- `central`：求解完整耦合的凸模型，包含QoS、BS预算、所有基站对的相对可听性和正功率定义域。
- `dual`：独立完整耦合对偶，精确有界标量内层、工作集合乘子定价，最终检查所有约束和对偶界。
- `cp4a_consistent_candidate`：完整耦合分层修正候选，内层QoS/可听价格、外层BS预算价格；不是原文subgradient逐行实现。
- `cp4a_subgradient_candidate`：独立内外层投影价格更新，提供回溯/递减步长实验，不调用fmincon更新乘子；原六策略报告使用固定oracle，另有可选首次定位与估计几何反馈路径。

目标函数为平均位置**方差**；summary同时保存条件平均位置**标准差**、有效人数和覆盖率。不同策略的有效用户集合可能不同，不能只比较条件均值判断优劣。无条件CDF终点为覆盖率。

## 输出

每次运行保存完整 `run.mat`（cfg、num、scene、ch、模型、功率、乘子、残差、测距与定位结果）、config.json、manifest.json（源码SHA256、MATLAB/工具箱版本、模型假设）、summary.csv和定位坐标。

每种策略还有功率、用户、测距及全部约束的CSV；central、dual和分层候选另有各自迭代轨迹。约束CSV的signed residual≤0表示可行，同时给出SI单位与归一化值。结果目录不会覆盖上次运行。

## 当前边界

当前使用二维oracle几何、独立测距误差、无未知接收机时钟状态、自由空间功率增益和假设的q。优化采用论文式(10)两项近似；式(9)三项近似另作求解后评价。相对可听不包含绝对检测门限。

DLL有限PSD数值积分和PRS能量明确的AWGN下界已实现；完整耦合对偶、分层修正候选及独立投影/次梯度实验分支已验证。PRS的通信服务、双向PSD与TDoA比较已实现，但默认严格服务合同不匹配，有限观测与稳态DLL仍不等价，不能称为完全公平对照。有限采样捕获已实现，连续稳态DLL、PRS采样对照及50次Monte Carlo尚未实现。作者原始步长及逐行一致性仍未确认；工程初始步长允许自定并以证书和敏感性检查验收。`MS_NOMA`中的原作者单小区演示保留独立入口。公开论文在工作区References目录。

## 阶段2：积分核验与PRS

默认抽查用户[1,10,20]的全部4条BS链路，积分缓存由同次运行的全部策略共用。采用固定全带前端 `global_frontend` 和精确正弦鉴别器，返回噪声、通信和定位干扰分量及网格收敛报告；两项/三项近似相对误差分别保存。积分只作后评价，优化目标保持式(10)。

```matlab
cfg = msnoma.config.base;
cfg.models.rangingIntegral.userIds = 1:cfg.paper.M; % 扩大检查，耗时增加
cfg.models.rangingIntegral.frontendMode = 'local_symmetric'; % 另一前端解释
cfg.io.showFigures = false;
main(cfg);
% 只运行原理论主干时，可关闭新增模块：
cfg.models.rangingIntegral.enabled = false;
cfg.models.prs.enabled = false;
cfg.models.dllBias.enabled = false;
```

PRS默认是广义OFDM广播信号，匹配central每个BS的总平均定位功率；在20ms周期内显式分配符号、峰值功率与BS顺序时隙。其他发送在PRS时隙静默。50MHz推广不是标准LTE。该基准采用已知相位AWGN延时CRB和共享OLS几何，没有匹配通信QoS/检测覆盖，也不等同于MS稳态DLL估计器；目前不能称为论文PRS曲线复现。

新增输出：`*_dll_validation.csv`、`dll_quadrature.csv`、`prs_broadcast_budget.csv`、`prs_positions.csv`、`prs_ranging.csv`与图5/6。`run_stage2_checks`还输出 `dll_frontend_sensitivity.csv/.mat`。本次两项近似最大方差误差约13.5%，完整数值与有效结果目录见PROGRESS.md。

## 阶段3：对偶与分层候选

阶段3为5策略、7图；当前默认6策略、11图。`optimization_comparison.csv`保存原始目标、对偶界、间隙、相对central目标差、原始残差/KKT和耗时；`paper_carrier_group_audit.csv`明确展示论文式(24)分组与实际载波数的不一致。两种新求解器均不调用central回退。

```matlab
cfg = msnoma.config.base;
cfg.experiment.algorithmIds = {'central','dual','cp4a_consistent_candidate'};
cfg.io.showFigures = false;
main(cfg);
```

cfg/run当前版本为8.1；旧cfg请重新从base建立后覆盖参数。`cfg.numerics.dual`控制对偶容差、工作集合和迭代上限；`cfg.numerics.cp4a`控制明确标注的候选模式和内外层容差。未通过原始可行性/KKT/对偶间隙的输出不会标为converged_verified。

固定oracle几何下，修正候选内层使用解析梯度/Hessian的fmincon、外层使用包络梯度SQP；另有阶段4独立投影/次梯度实验及阶段7实际采样反馈路径。作者Algorithm 1逐行一致性仍未确认。详细公式与差异在docs/DUAL_CP4A_NOTES.md，最新场景及数值见PROGRESS.md。

## 阶段4：投影/次梯度实验

`docs/SUBGRADIENT_IMPLEMENTATION.md`记录完整公式、坐标/步长、失败状态、导出接口和原文差异。默认回溯加速分支在20/50MHz与绑定预算场景均通过证书；diminishing为可复跑实验，必须检查真实status，不能假定有限迭代内收敛。

新配置位于cfg.numerics.subgradient。只选集中式和新分支：

```matlab
cfg=msnoma.config.base;
cfg.experiment.algorithmIds={'central','cp4a_subgradient_candidate'};
cfg.io.showFigures=false;
main(cfg);
```

新增subgradient_inner_calls.csv、subgradient_inner_history.csv、subgradient_outer_prices.csv及图8。共享证书新增相对可听检查，避免弱链路“归一化残差达标但物理可听失败”。未更改实际可听性门限。参考解和投影失败轨迹分别保留，不作集中式回退。

## 阶段5：PRS的TDoA、双向干扰和服务审计

`run.prsBenchmark`继续保存旧AWGN/OLS示例。新 `run.prsComparison` 在相同oracle掩码下给出MS/PRS的相关TDoA传播，并分别审计全局通信静默与其他BS通信继续两种PRS时隙策略。默认参考最小TOA方差BS，支持固定参考；缺站、秩失败及不可行能量计划保留具体原因。

双向有限OFDM PSD使用现有gCP/gPC的正确接收机维度。PRS端采用对角高斯频率噪声近似，C端逐PRS槽核算BER；原模型未计的C-to-C干扰仍未计。活动C功率保持原值，通信占时、实际平均能量、维持原能量所需峰值、最坏活动BER和尝试速率分开导出。默认严格满尝试速率合同会因PRS占时失败，不会自动降速放宽结果。

```matlab
cfg=msnoma.config.base;
cfg.models.prs.comparison.tdoa.referencePolicy='fixed';
cfg.models.prs.comparison.tdoa.referenceBs=1;
cfg.io.showFigures=false;
main(cfg);
% 只关闭阶段5附加比较，保留旧PRS基准：
cfg.models.prs.comparison.enabled=false;
```

新输出为 `prs_comparison_summary.csv`、两种策略各自的service/communication/tdoa/ranging CSV及图9。定位广播能量、几何口径、通信服务、通信能量和观测是否匹配各有单独标志。`observationMatched=false`、`fullyMatched=false`；PRS的几何有效不能解释为检测覆盖率。详细公式、限制和下一步在 `docs/PRS_COMPARISON_IMPLEMENTATION.md`。

## 阶段6：DLL独立方差分解与接收滤波解析审计

同一central功率、信道与环路下，交叉组合精确正弦/线性化、有限C谱瓣/同一有限载波的窄线谱、全带/局部前端共8种设置。保存全部12条单因素配对边，另给出6步可闭合路径，将其余渐近近似、N与N+1、定位干扰省略单独列出。三个单因素百分比受上下文影响，不能直接相加。这里的“偏差”指预测方差差异，不是实际时延估计的经验均值偏差。

```matlab
cfg=msnoma.config.base;
cfg.models.dllBias.userIds=[1,10,20]; % 默认抽查；扩大到1:M会增加积分计算
cfg.models.dllBias.referenceAlgorithm='central';
cfg.io.showFigures=false;
main(cfg);
% 只关闭阶段6，保留阶段5的9图口径：
cfg.models.dllBias.enabled=false;
cfg.models.prs.comparison.receiverAuditEnabled=false;
```

DLL独立研究沿用rangingIntegral数值设置，普通积分模块禁用时仍可运行。缺少指定参考会返回reference_solution_unavailable，不换算法。新增6个dll_bias CSV和图10；完整核、收敛历史和矩阵在run.mat。

PRS接收端另计算矩形符号PSD经矩形匹配滤波的全频积分，与旧载波中心PSD取样并列，新增2个prs_receiver_filter CSV和图11。该支路假设平稳独立符号和平均相对时序；不是相干同步OFDM解，也没有生成样本、估计时延或位置。旧点模型、PRS测距界和oracle优化保留，fullyMatched继续false。公式、接口边界和数值在 `docs/DLL_BIAS_AND_RECEIVER_AUDIT.md` 与PROGRESS.md。

## 文献版本核验：arXiv 早期稿

新增 `F:/4Cell_MS_NOMA/1910.04350v1.pdf`（2019-10-10，30页）作为交叉验证来源；最终参数仍以JSAC 2022为准。比较见 `docs/ARXIV_VERSION_COMPARISON.md`，来源登记见 `docs/SOURCE_LEDGER.json`。关键差异包括PCJPA/CP4A、式(10)三项factor与两项variance、0.8W/0.1W示例，以及PRS参考由旧[35]变为最终[28]。单小区目录的Joint channel/time delay论文对应旧[35]，不能替代当前最终版PRS参考。来源登记当时仅修改文档；阶段7的采样链路是后续独立工程实现。

## 阶段7：有限采样链路、估计几何反馈与步长试验

```matlab
cfg=msnoma.config.base;
cfg.models.sampledReceiver.enabled=true;
cfg.experiment.positionFeedback.enabled=true;
cfg.experiment.stepSensitivity.enabled=true;
cfg.experiment.stepSensitivity.multipliers=[0.01,1,100];
cfg.io.showFigures=false;
folder=main(cfg);
```

新增模块默认关闭，原main仍为6策略/11图；采样、反馈各增加一图，步长试验再增加一图。默认采用作者L1CP码表、4倍带宽采样、请求1023码片并补齐完整通信符号。实际观测长度、码片/符号、理想CSI假设和随机种子均保存；可显式选择独立随机BPSK模式，不能把受控模式的结果当作者码场景的覆盖率。

`waveform -> receiver.estimateDelay -> geometry.estimateRanges -> experiment.positionFeedback`已接入main。估计器不读取真实时延和用户位置；首次全部lambda=1，后续lambda只取上一轮实际估计坐标。完整反馈要求全部用户及全部BS测距有效；失败停止并保留记录。实测位置误差和理论oracle标准差分开导出；有限包捕获不等同0.2Hz稳态DLL，PRS完全匹配标志保持false。

工程步长由cfg.numerics.subgradient配置；敏感性实验同时缩放内外初始步长，保留回溯后的实际接受值、原始可行性、KKT和对偶证书。它们是缩放坐标下的工程参数，不声称恢复作者b1/b2/b3。

新增CSV为sampled_ranges/positions/packet_audit、feedback_iterations/geometry_powers/ranges/positions及step_sensitivity；完整对象在run.mat。方法、接口与失败状态见 `docs/SAMPLED_RECEIVER_AND_FEEDBACK.md`，验收和后续优先级见PROGRESS阶段7。

## 阶段8：可选估计P信道消除与干扰诊断

```matlab
cfg=msnoma.config.base;
cfg.models.sampledReceiver.enabled=true;
cfg.models.sampledReceiver.interferenceCancellation.enabled=true;
cfg.experiment.positionFeedback.enabled=true;
cfg.io.showFigures=false;
folder=main(cfg);
% 仿真发射消融，可单独开启：
cfg.experiment.interferenceAblation.enabled=true;
cfg.experiment.interferenceAblation.userIds=18;
folder=main(cfg);
```

新接收器额外假设全部P码及相对功率已知，通过同BS全部P分量联合估计时延/复增益，只用估计值消除P干扰；通信数据与真实时延不用于消除。局部拟合/高斯代理精度门控失败时保留原接收器。门控未标定真实多址干扰，正常预算改善不能保证低预算改善，默认仍关闭；不能为追求覆盖而自动挑选更好的一组误差。

成功消除后使用C+噪声预白化，白化权重及检测分母也随之改变，明确作为另一接收模式。样本相关的消除估计误差没有传播至CRB/虚警；PRS完全公平和稳态DLL标志维持false。新增P信道估计CSV及诊断图；消除与发射消融分别增加1图。细节、参数与反例见 `docs/ESTIMATED_P_CANCELLATION.md` 和PROGRESS阶段8。

原区域/正常预算下，49.77MHz参考20/20有效且3轮反馈收敛；20.16MHz参考20/20有效，默认3轮到上限，另仅把maxFixes改5即在第4轮达到原停止容差。低预算反例仍为6/20，未默认推广。112项检查及七个配对/控制场景见STAGE8_VERIFICATION.json；额外轮数敏感性见STAGE8_FEEDBACK_LIMIT_SENSITIVITY.json，未放宽容差。


阶段8.1已验收：残差主导门控默认上限0.5，低预算回退原接收器，正常预算改善保留；117项检查及六次完整main通过。门控为保守工程使用条件，不是精度保证；消除开关仍默认关闭。详见docs/STAGE8_1_RESIDUAL_GATE.md与PROGRESS末尾。
