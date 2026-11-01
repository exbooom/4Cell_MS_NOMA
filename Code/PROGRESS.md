# 多小区复现进度与交接

当前阶段8.1已验收：新增只使用拟合数据的残差主导门控，低预算反例回退原接收器，正常预算改善保留；117项检查通过，115个MATLAB源静态0诊断，六次完整main及独立哈希/图面核验。新接收器仍默认关闭，门控不是精度保证。接手见末尾阶段8.1、docs/STAGE8_1_RESIDUAL_GATE.md和STAGE8_1_VERIFICATION.json；阶段1–8记录保留为历史源码快照。

更新日期：2026-11-01（Asia/Shanghai）。本文件是后续接手入口，持续维护的原件为 `F:/4Cell_MS_NOMA/Code/PROGRESS.md`。

## 当前任务

用户已授权分步编写场景、跨基站信道、测距、几何、可听性和优化模块，并由main调用。阶段1建立单场景理论主干；阶段2接入有限带宽DLL与能量PRS基准；阶段3/4完成完整耦合对偶、分层修正及独立投影候选；阶段5完成TDoA/PRS服务审计；阶段6接通独立DLL方差近似分解与平均时序接收滤波解析审计；阶段7接入实际采样测距、位置估计、首次lambda=1与估计几何反馈；阶段8加入可选估计P消除和受控发射消融。接手以末尾阶段8记录为准。用户允许数值步长自定，工程选择及敏感性单独验收；不把有限包捕获当成稳态DLL，也不把修正候选称为作者原版CP4A。

## 阶段状态

| 阶段 | 状态 | 验收 |
|---|---|---|
| 资料与架构 | 完成 | docs/ARCHITECTURE_PLAN.md、SOURCE_LEDGER.json |
| A 配置、场景、信道 | 已实现并验证 | 固定局部随机流、四角基站、P/C两类用户、Friis功率增益 |
| B 通信、测距、几何、可听性 | 已实现并验证 | BER正反算、两项/三项测距方差、二维QR几何、相对可听 |
| C 约束、可行性、集中式优化 | 已实现并验证 | 全耦合A、预算和全对可听，LP预检，fmincon，KKT与手算例 |
| D main、输出与交叉检查 | 已实现并验证 | 20/50MHz名义分支均跑通4×20；MAT/CSV/4图/源码哈希 |
| E 有限带宽DLL积分 | 已实现并验证 | 有限C载波PSD、正弦鉴别器、分项积分、网格收敛与前端敏感性 |
| F 显式能量PRS基准 | 已实现并验证公式与预算 | [28]式(5)/(16)、广播能量、峰值功率和分时调度；尚非论文公平曲线 |
| G 完整耦合对偶、分层候选 | 已实现并验证 | 独立对偶界、原始可行性/KKT、20/50MHz和绑定预算压力场景 |
| H 独立投影/次梯度实验 | 已实现并验证默认分支 | 内外投影更新、回溯/递减步长、完整证书和相对可听门限 |
| I PRS可比性审计 | 已实现，真实全公平未达成 | 相关TDoA、两种时隙策略、双向PSD、活动BER、时长/能量/峰值与观测标志 |
| J DLL方差分解、解析C接收滤波 | 已实现并验证 | 八组合、12配对、6步闭合；矩形PSD积分与平均时序假设审计 |
| K 有限采样、实际测距/位置与首次反馈 | 已实现并验收 | 原区域失败保留；中心区域实际接收20/20、3轮反馈收敛；oracle独立保留 |
| L 工程步长敏感性 | 已实现并验收 | 0.01/1/100倍率，回溯后实际接受步长和完整证书 |
| M 全P码库拟合与估计干扰消除 | 已实现并验收可选分支 | 正常预算改善；低预算反例保留，未默认启用；联合拟合不使用真值/通信数据 |
| N 干扰发射消融与精度门控 | 已实现并验收 | 六类发送掩码，独立估计/误差，门控未等于真实干扰校准 |
| 后续：低预算稳健接收、稳态DLL/PRS与统计 | 继续推进 | 有限捕获不等于稳态跟踪；不以条件均值替代覆盖率 |

## 接手规则

1. 先读本文件和README，再读docs/ARCHITECTURE_PLAN.md；后者是设计快照，具体已实现接口以源码为准。
2. 原始作者代码保留在MS_NOMA。本轮新代码统一进入+msnoma与main.m，不混入单小区工程。
3. 每次运行保存完整配置、坐标、模型、求解状态、残差与性能。论文未公开参数标为assumed。
4. MATLAB使用D:/DevTools/MATLAB2024b/bin/matlab.exe；项目Python如需使用，遵守用户指定Anaconda解释器。
5. 后续开始工作时，先检查本文件末尾的验证记录与下一步，避免把计划模块当成已完成。

## 本轮模型决策

- 首轮采用整数G回算有效带宽，20/50MHz分别约20.16/49.77MHz。
- 二维已知高度、无接收机钟差状态；真实坐标计算oracle几何。
- 自由空间Friis功率增益；通信用户在最近基站分区内抽样；不额外添加快衰落。
- 默认标准QPSK每比特能量分支；另保留印刷BER及其配套反演。
- q为显式假设的复合可听性比例；不是已确认的作者Omega。
- 优化使用两项测距近似、完整跨基站QoS和全对可听约束。
- 失效定位为NaN，覆盖率与条件平均分别报告；仅相对可听不等于真实检测成功。

## 验证记录

本轮验证命令：

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); addpath('tests'); run_stage1_checks;"
```

- 17项数学检查通过，详见tests/run_tests.m。覆盖载波、随机流隔离、用户服务区、信道接收点、非对称多维干扰、两种BER反演、负QoS阈值、测距标度/单位/第三项、几何对称及秩亏、可听性、目标导数、手算最优解、q>1预检、等功率策略和CDF分母。
- 新MATLAB文件的Code Analyzer诊断为0。
- 两个真实规模场景均通过独立QoS、全对可听、目标重算、KKT与文件导出核验。
- 将数值功率下限从1e-10 W缩小到1e-12 W，在相同20MHz场景重放集中式求解，目标相对变化为2.6201e-14；实际解的最小功率大于原下限100倍。
- 验证摘要：docs/STAGE1_VERIFICATION.json；更详细数据在相邻MAT文件和各次run.mat。

### 本轮有效运行目录

| 分支 | 结果目录 |
|---|---|
| 名义50MHz，有效49.77MHz，G=79，N=1658 | `F:/4Cell_MS_NOMA/Results/run_20261031_231420_dd29de9c` |
| 名义20MHz，有效20.16MHz，G=32，N=671 | `F:/4Cell_MS_NOMA/Results/run_20261031_231432_80ec02e3` |

此前230924/230937目录是初轮试跑，最后一轮有效结果以上表及STAGE1_VERIFICATION.json为准，不删除试跑数据。

### 数值结果（仅为显式假设下的单场景）

| 指标 | 49.77MHz集中式 | 20.16MHz集中式 |
|---|---:|---:|
| 状态 | converged_verified | converged_verified |
| 有效定位用户 | 20/20 | 20/20 |
| 全部通信链路最大BER | 0.01 | 0.01 |
| 最大正归一化约束违例 | 0 | 0 |
| 缩放KKT驻点残差 | 7.845e-10 | 2.687e-9 |
| 最大互补松弛残差 | 2.2217e-8 | 2.4279e-8 |
| 平均位置方差（m²） | 0.01758 | 0.27014 |
| 有效用户平均位置标准差（m） | 0.063887 | 0.24185 |

相同q下两个等功率策略均有16/20用户满足相对可听规则；budget_equal最大BER约0.058，违反QoS；qos_scaled_equal最大BER为0.01。集中式覆盖更全但有效用户集合不同，不能将其条件平均精度与基准条件均值直接作公平比较。本轮数值不用于宣称论文表II已复现。

## 已实现文件的职责

- main.m：主入口，调用runSingle，默认6种策略、8幅图；可接收用户修改的cfg。
- +config：base、validate；严格检查已支持的几何、时钟与分布模式。
- +model：numerology、sincpi、spectralKernels。
- +scenario：samplePositions、make、fromSnapshot、freeSpaceGain、channels。
- +comm：makeBerModel、interferenceMatrix、interference、ber、qosLimit。
- +ranging：equivalentCommunicationGain、factorTwoTerm、varianceTwoTerm、varianceThreeTerm、integralKernels、varianceIntegral。
- +prs：allocation、powerBudget、delayCrb、positionBenchmark。
- +geometry：weights、hearability、positionAccuracy。
- +optim：buildProblem、objective、checkFeasibility、solveCentral、equalPower、diagnostics，以及阶段3新增dualValue、dualActiveSet、dualSolution、solveDual、carrierGroupAudit、solveCp4a，以及阶段4新增projectedDual、solveSubgradient。
- +experiment：runSingle、evaluate。
- +report：aggregate、unconditionalCdf、plotRun、plotPrecision、plotOptimization、saveRun、exportConstraints、exportPrecision、exportOptimization、exportSubgradient、plotSubgradient、sourceHashes。
- tests：run_tests、run_stage1_checks、run_precision_tests、run_stage2_checks、run_dual_tests、run_stage3_checks、run_subgradient_tests、run_stage4_checks。

实际签名、形状和返回状态见docs/MODULE_INTERFACES.md。注意report.plotRun是实际名称，原架构草案中的report.plot只是规划名称。

## 下一步优先级

1. **PRS可比性（当前下一项）**：独立投影/次梯度链已通过阶段4，优化实验主干可作为参考。接下来统一PRS通信服务、干扰、TDoA协方差及有限观测口径；作者实际步长仍未知，不再用等待作者作为阻塞。
2. **优化稳健性扩展**：继续保留集中式、完整对偶、SQP分层和投影分层四条独立路径；默认回溯分支已验收，递减步长仅作显式实验。作者b1/b2/b3、首次定位和原始分组仍不能确认，不按论文数值反调参数。
3. **DLL近似偏差**：分离正弦线性化、有限载波PSD、有限带宽积分及前端频率坐标各自贡献；默认只抽查3名用户，可设置userIds=1:M扩大检查。当前积分仅作求解后评价，优化仍使用式(10)。
4. **统计与敏感性**：复用runSingle增加50次Monte Carlo；q、噪声、C用户分布和带宽模式分别记录，统一随机场景；保存求解失败与覆盖不足的不同分母。
5. **采样接收机与估计几何**：当前source=oracle。首次λ=1、接收机输出位置反馈尚无实现；不要宣称“真实位置迭代收敛”。

## 重要接手细节

- cfg.models.hearability.ratio是q的唯一计算输入；不从L1C码表推断作者Omega。
- cfg.assumed.noisePsdWPerHz是实际N0；温度与噪声系数字段是默认假设说明，改变它们不会自动改PSD。
- cfg.paper.positioningBudgetW当前为每BS相同标量；不同BS预算尚未暴露为配置向量。
- 数组顺序：P为K×M，p=P(:)，变量号k+(m-1)K；通信行号k+(n-1)K。
- gPC为K×N×K×M，gCP为M×K×N。接收点不同，不能互换。
- 时间方差到距离方差仅在ranging层乘c0²。优化目标为mean(Psi²)，报表均值为mean(Psi(valid))。
- rank不足、听到BS不足、通信QoS失败与优化不可行是不同状态；run.mat保存全部输入，CSV无效定位写NaN。
- 当前实现只支持integer_nearest带宽、二维oracle几何、无未知时钟状态、指定用户分布；配置选择其他模式会明确报错。
- sourceHashes对每个新.m文件按原始字节做SHA256。后续改源码应重跑验证，生成新结果目录；旧结果仍由其manifest关联当时源码哈希。

## 阶段2历史验证与交接（2026-11-01）

验证命令：

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); addpath('tests'); run_stage2_checks;"
```

17项原有数学检查和9项新增检查均通过，Code Analyzer诊断0。新增检查包括独立自适应积分对照、功率/噪声标度、频率坐标与对称性、小D极限、PRS导频归一化、独立物理能量频率矩CRB、峰值/联合调度不可行、LTE Toolbox单符号图样和广播能量不随用户数重复计算。完整检查使用LTE Toolbox验证6RB单符号12导频；默认广义OFDM路径不依赖该工具箱。

两带宽均跑通main，集中式解仍通过QoS、全对可听和KKT核验；MAT、CSV、默认6幅PNG导出完成。验证摘要见 `docs/STAGE2_VERIFICATION.json/.mat`。

| 指标 | 49.77MHz | 20.16MHz |
|---|---:|---:|
| DLL选中用户/链路 | [1,10,20] / 12 | [1,10,20] / 12 |
| 积分原始分项最大网格相对变化 | 6.0815e-10 | 3.7066e-9 |
| 两项近似相对积分的最大方差误差绝对值 | 13.5488% | 13.4515% |
| 三项近似相对积分的最大方差误差绝对值 | 13.4323% | 13.3407% |
| global/local前端最大方差差异绝对值 | 11.7001% | 11.7038% |
| PRS四BS合计时隙占20ms周期比例 | 17.3333% | 18.6667% |
| 假设PRS条件平均位置界（m） | 0.000126495 | 0.000297561 |

两/三项误差定义为 `approxVariance/integralVariance-1`，本次抽查均为负，说明近似低估方差；前端差异为 `globalVariance/localVariance-1` 的最大绝对值。网格变化仅表示数值积分收敛，不能用来消除模型假设或近似误差。

本阶段有效目录：

- 49.77MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_000324_2304afbf`
- 20.16MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_000341_1c0cd0cd`

`run_20261101_000013_d354497e`是绘图修正前的试跑，导出不完整；保留作追溯，不作为验收结果。阶段1目录继续保留历史意义；后续默认入口及数据结构以阶段2为准。

### DLL频率与计算范围

默认 `cfg.models.rangingIntegral.frontendMode='global_frontend'`：在未中心化载波坐标中，以 `(N+1)*deltaFc/2` 为固定接收前端中心，对第m个P子载波减去 `m*deltaFp` 后积分。这是明确的工程解释，尚非作者确认的实现。`local_symmetric`则对每个目标P载波使用±Bfe/2；两种解释都保存，未让RF载频3.5GHz进入鉴别器。当前默认选择是阶段2实现决策，覆盖最初架构草案中的建议。

积分计算式(6)/(51)—(55)的有限C载波PSD与正弦鉴别器；P干扰按式(8)只计其他BS同一P子载波。先用x=f_local*Tp无量纲化，8→16→32→64个采样/通信载波间隔逐次细化，记录每次分项积分并要求收敛。默认用户[1,ceil(M/2),M]共享缓存供所有功率策略评价。输出K×选中用户，不能把抽查结果当作全部用户的积分精度。

优化器仍采用式(10)两项凸目标，积分和三项模型均为后评价。选中用户的积分位置精度在 `run.algorithms(i).evaluation.integral.positionStdM`；几何重新使用该策略的相对可听集合。CSV导出逐链路分量，完整缓存与几何在run.mat。

### PRS能量、调度与边界

依据主论文参考[28] *Achievable Localization Accuracy of the Positioning Reference Signal of 3GPP LTE*，DOI 10.1109/ICL-GNSS.2012.6253127，实现式(5)和式(16)。OneCell目录中的 *Joint channel and time delay estimation...* 是另一篇论文，不能混用。

默认15kHz、comb6、去DC的 `generalized_ofdm`，权重平方和=Nc；49.77MHz推广不属于LTE标准。另提供LTE Toolbox生成的单符号参考模式，不将重复通用符号称为完整LTE PRS occasion。采用已知相位、无CP、导频独占和信息独立累加。

每个BS的PRS平均功率等于所选参考策略（默认central）所有P用户功率之和，广播能量只计算一次；E=Pavg×20ms。默认在1.1W峰值内选足够符号数，Pinstant=E/实际观测时长，不做能量裁剪。四BS依次占用互不重叠时隙，期间其他发送静默，要求总时隙≤测量周期。预算、起始时间和实际符号数均导出。不满足功率/时间预算则返回不可行状态，不伪造位置结果。

位置传播复用主论文无未知时钟的二维OLS几何，尚未实现[28]原始TDoA相关协方差。当前PRS没有匹配通信QoS，也没有测量检测门限/覆盖率；`geometryValid`不能解释为可检测覆盖。PRS是有限观测AWGN下界，MS为稳态DLL，两者估计器、带宽利用和服务假设不同。本次PRS界显著更小，不能据此声称算法优劣或论文PRS曲线已复现。

### 阶段2结束时的下一动作（已由阶段3推进）

以下为当时的接手建议：先读本节及阶段2验证摘要，再从 `+msnoma/+optim/buildProblem.m`、`objective.m`、`solveCentral.m` 开始实现对偶基准；当时尚无solveDual/solveCp4a；现已实现，见阶段3。保留集中式解作独立比较，任何新算法先通过相同原始约束与目标重算，再扩展Monte Carlo。每次修改MATLAB源码后新建运行目录，核对manifest的SHA256，不覆盖历史结果。

最终归档核对：两个阶段2有效目录各自manifest中的49个MATLAB源文件SHA256均与当前源码逐一匹配；每个目录含6幅PNG。阶段2新增图5/6已目视检查。文档更新未修改已验证的MATLAB源码。

## 阶段3历史验证与交接（2026-11-01）

按优化依赖链实现：完整问题→有界标量内层→完整耦合对偶→与集中式对照→内层QoS/可听、外层预算的分层候选→BER/测距/几何重新评价→图表与证书。

新增函数均在 `+msnoma/+optim`，由runSingle和main实际调用，无集中式回退。cfg/run版本升至3.0；旧配置应从新base重建后覆盖自定义参数。默认algorithmIds为budget_equal、qos_scaled_equal、central、dual、cp4a_consistent_candidate；两种新增求解器都共享同一个原问题，分别求解价格与恢复功率。

### 算法范围

- solveDual完整使用A'×mu及C'×beta，包括跨接收机干扰与反向可听贡献；约束生成仅减少当前非零乘子集合，最终仍检查全部原始行。
- 标量内层精确处理有限上下界，负价格时取上界，不使用sqrt(abs(price))。对偶界、相对间隙、原始可行性、KKT、乘子均通过后才标converged_verified。
- solveCp4a只支持mode=full_coupling_hierarchical_candidate：内层完整QoS/可听定价，外层BS预算价格。内层用解析梯度/Hessian的fmincon，外层用包络梯度SQP；这与论文Algorithm 1的内外subgradient不同。literalPaperAlgorithm=false，固定oracle几何，没有首次定位或位置反馈。
- 式(24)的打印载波分组在49.77MHz越出1482个、20.16MHz越出589个载波；仅导出审计，不用截断分组求解。式(66)/(67)的导数从实际原问题重新推导。
- 具体推导、原文差异和未实现项见 `docs/DUAL_CP4A_NOTES.md`。当前候选是可校验的分层数学实现，不是作者实现的确定恢复。

### 验证

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); addpath('tests'); run_stage3_checks;"
```

36项检查通过：原17项、积分/PRS 9项、新增对偶/分层10项。新增内容包括有界标量、负价格、零系数、梯度/Hessian有限差分、弱对偶、手算预算最优解、非对称完整耦合、工作集约束恢复、不可行状态、错误模式拒绝及载波分组审计。59个MATLAB源文件的Code Analyzer诊断0。

三个完整main均通过：4BS/20定位用户的50/20MHz默认分支，加20MHz/每BS定位预算1mW的压力分支。两种新求解器都通过原始约束、独立目标、KKT、相对对偶间隙、全部通信QoS与20/20相对可听定位检查。PRS和DLL缓存继续有效。每个场景输出5种策略、7幅PNG和逐项CSV。

| 场景 | dual相对central目标差绝对值 | 分层候选相对central目标差绝对值 | 分层相对对偶间隙 | 外层迭代/内层求值 |
|---|---:|---:|---:|---:|
| 49.77MHz，0.1W/BS | 1.8572e-7 | 1.8572e-7 | 1.0686e-10 | 1 / 1 |
| 20.16MHz，0.1W/BS | 1.4074e-8 | 1.4019e-8 | 3.8084e-10 | 1 / 1 |
| 20.16MHz，0.001W/BS | 5.1706e-9 | 4.9434e-9 | 4.1670e-10 | 19 / 42 |

默认场景预算未耗尽，外层最优价格为0；不能据此宣称完整外层更新已被充分验证。压力分支确实触发预算绑定与19次外层更新，其最大归一化违例1.5340e-10，小于1e-7验收阈值；两个默认分支最大正归一化违例为0。轻微低于集中式目标是数值求解精度差，不是改变模型获得收益。

当前有效结果：

- 49.77MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_103739_44ea0314`
- 20.16MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_103801_2a4b0adb`
- 20.16MHz绑定预算：`F:/4Cell_MS_NOMA/Results/run_20261101_103811_41e5c04c`

验证摘要：`docs/STAGE3_VERIFICATION.json/.mat`；新增结果`optimization_comparison.csv`、`paper_carrier_group_audit.csv`、dual/cp4a_consistent_candidate的solver_history/constraints/powers/ranging/users/dll_validation文件与图7。约束生成早期功率可能不可行；不能只凭对偶间隙小停止，图7同时显示原始约束违例。

阶段1/2的历史验证脚本显式恢复原三策略，避免默认增加算法影响历史导出检查。旧目录的manifest继续关联当时源码，不应与当前源码相同作为历史结果的硬性要求；本阶段目录才与当前源码核对。

### 阶段3结束时的接手顺序（第一项已由阶段4推进）

1. 当时计划优先决定并实现独立的原文subgradient实验分支，保留本轮已验证的完整耦合/分层基准；不沿用有疑点的打印分组和错误导数。
2. PRS统一通信服务/干扰、TDoA相关协方差及有限观测比较口径；已有广播能量预算可继续复用。
3. 分离DLL正弦近似、PSD近似和前端解释的各自偏差。
4. 在统一场景和明确失败分母下增加50次Monte Carlo、q/噪声/带宽敏感性。
5. 最后接采样接收机、估计位置、首次lambda=1与几何反馈。保持oracle与估计几何结果分开。

最终补充核验：三个阶段3有效目录中各59个MATLAB源文件的SHA256均与当前源码逐字节匹配，每个目录7幅PNG；默认拓扑图及默认/压力场景优化图已目视检查。核对记录在docs/STAGE3_SOURCE_HASH_VERIFICATION.json。另将maxPasses强制设为1，独立确认dual返回dual_recovery_infeasible，分层返回hierarchical_inner_dual_not_verified及NaN功率；未调用集中式回退。这两项为补充检查，不计入36项套件数量。

## 阶段4验证与交接（2026-11-01，历史记录）

按计划完成最高优先级的独立一阶价格更新实验，保留阶段3所有参考求解器。新ID `cp4a_subgradient_candidate` 实现内层QoS/可听价格投影、外层BS预算价格投影；乘子更新不调用fmincon或其他参考价格优化器，不调用集中式回退。可行性预检仍复用linprog，不能称整个工程无优化工具箱依赖。

cfg/run版本4.0，默认6策略、8图。新增 `optim.projectedDual`、`optim.solveSubgradient`、`report.exportSubgradient`、`report.plotSubgradient`，已实际接入main。实现推导与接口详见 `docs/SUBGRADIENT_IMPLEMENTATION.md`。

### 方法及边界

- 默认backtracking分支：行归一化的内层投影梯度、FISTA动量/梯度重启和稳定二次上界回溯；外层投影梯度/Armijo回溯。保持完整A/C耦合。
- diminishing分支：内外初始步长/(t+1)^eta，默认eta=0.75，无动量；测试验证实际一步更新、重放确定性和有限迭代失败状态，未宣称该默认递减步长在真实规模内达到验收精度。
- 这些步长、归一化、动量和warm start均为显式工程选择，不是已知作者b1/b2/b3。不使用越界分组或打印导数，几何仍固定oracle，无首次定位反馈。因此仍标candidate而不是原版CP4A。
- 内层二次上界使用稳定的标量Bregman余项，避免相减O(1)对偶值导致近最优阶段浮点抵消。
- 内层停止同时看归一化可行性、绝对互补项总和/目标、实际相对可听违例；外层停止必须通过完整原始/对偶证书，不用乘子相邻变化小代替收敛。

### 本轮修复的物理边界

第一轮真实规模回归发现：归一化残差约8.3e-9虽小于原1e-7门限，弱链路相对接收功率仍可能偏离，使1名用户只听到3个BS。没有放宽geometry.hearability的门限；新增diagnostics.maxRelativeHearabilityViolation、内层更严格的相对可听停止条件，并在共享dualSolution证书中加入同一物理比例检查。新增单元检查专门隔离验证该门禁。

旧的 `run_20261101_105215_6fbb2630` 是捕获问题的试跑，不能作为验收。105403/105430/105449目录为中间回归；当前源码验收以本节有效目录为准，不删除历史证据。

### 验证

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); addpath('tests'); run_stage4_checks;"
```

45项全部通过：17原数学、9积分/PRS、10对偶/分层、9新增投影检查。新增包括实际投影一步、无动量梯度、手算预算最优解、完整耦合与每轮重置、递减重放、有限迭代失败、非法模式拒绝、无价格优化器调用及相对可听证书边界。65个MATLAB源文件Code Analyzer诊断0。

三个完整main都通过，均为4BS/20定位用户：50MHz、20MHz，以及20MHz/1mW每BS预算压力场景。四条优化路径均通过独立目标重算、通信QoS、全对可听、KKT、对偶间隙和20/20有效定位；DLL与PRS后评价继续有效。各目录含6策略数据和8幅PNG。

以下为新投影分支最终结果：

| 场景 | 相对central目标差 | 相对对偶间隙 | 最大归一化违例 | 最大相对可听违例 | 外层更新/内层调用 |
|---|---:|---:|---:|---:|---:|
| 49.77MHz，0.1W/BS | 1.42456e-08 | 1.99981e-07 | 5.08572e-09 | 0 | 0 / 1 |
| 20.16MHz，0.1W/BS | 1.85894e-07 | 1.99988e-07 | 2.24435e-10 | 1.08028e-14 | 0 / 1 |
| 20.16MHz，0.001W/BS | 7.77639e-07 | 7.82854e-07 | 0 | 0 | 61 / 76 |

有效结果目录：

- 49.77MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_105550_e4c148e5`
- 20.16MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_105616_ab5e4b45`
- 20.16MHz绑定预算：`F:/4Cell_MS_NOMA/Results/run_20261101_105634_98a5eb0f`

摘要在 `docs/STAGE4_VERIFICATION.json/.mat`。新增subgradient_inner_calls/inner_history/outer_prices三类CSV；内层调用含外层被拒绝试探，outer_prices只含已接受更新。optimization_comparison.csv新增usesFmincon、stepRule及相对可听残差；完整终态/配置在run.mat。图8显示外层证书、预算价格、内层更新数和最终功率。

历史stage3脚本显式选原5策略，阶段1/2继续3策略；默认算法增加不会改变历史检查口径。旧配置请从新base建立再覆盖参数。

### 下一步（依赖顺序）

1. PRS公平口径：通信服务/资源成本、其他BS干扰、参考[28]TDoA相关协方差、有限观测与MS稳态DLL的比较定义。
2. DLL偏差分离：正弦线性化、有限PSD、前端频率解释分别量化。
3. Monte Carlo和敏感性：统一随机场景、保留求解失败及不同覆盖分母，先小样本回归后50次。
4. 采样接收机、测距、估计位置，再接首次lambda=1及几何反馈；oracle基准继续保留。

最终归档核验：三个阶段4有效目录各65个MATLAB源文件SHA256均与当前源码逐字节匹配，每目录8幅PNG。拓扑图和默认/压力场景投影图已目视检查；记录见docs/STAGE4_SOURCE_HASH_VERIFICATION.json。验收后仅更新交接文档，未修改MATLAB源码。

## 阶段5验证与交接（2026-11-01，历史记录）

按依赖顺序完成PRS比较层。cfg/run版本5.0，默认仍为6策略，新增图9后共9图。所有新增计算均在Code/+msnoma的模块/函数中，由experiment.runSingle及main调用；原作者MS_NOMA未改动。详细公式和边界见 `docs/PRS_COMPARISON_IMPLEMENTATION.md`，接口见 `docs/MODULE_INTERFACES.md`末节。

### 已实现

1. `geometry.tdoaAccuracy`：按参考[28]式(20)–(23)构造共享参考误差产生的非对角协方差；Cholesky白化、QR解算、二维至少3站、参考选择/掩码/零距离/秩和缺站失败。默认最小TOA方差参考，支持固定站且不静默换站。BS同步、公共接收机偏移被差分消除；局部冻结协方差，没有空间导数。MS与PRS另按相同oracle掩码/相关TDoA规则后评价，并保留MS相对可听掩码结果；MS仍是稳态方差代理。
2. `prs.interferenceKernels`、`delayCrbInterference`：有限OFDM PSD的C→PRS和PRS→C两个方向，正确区分gCP/gPC维度，导频分块和共用缓存；顺序PRS没有同时其他PRS。已知相位的对角高斯频率噪声近似，零干扰还原旧CRB，非完整有色噪声似然。
3. `prs.serviceAudit`：全局通信静默和其他BS通信继续两种时隙策略；实际活动C功率保持原值，分开核算机会时长、实际平均C能量、等能量所需功率/峰值、活动槽均值/最坏BER、尝试比特率和预期正确比特率。静默槽不记为BER=0；无通信机会时活动BER为NaN。默认严格满尝试速率合同，不自动降速、补功率或截断。
4. `prs.comparison`、独立CSV/图/manifest集成：保留旧 `run.prsBenchmark` 的AWGN/OLS示例；新 `run.prsComparison` 分别记录定位广播能量、几何、通信服务、通信能量及观测匹配。没有PRS调度/通信补偿联合优化，只审计当前明确时隙方案。

### 不应遗漏的数值边界

- 定位广播能量与TDoA几何口径可以匹配，但PRS有限观测已知信号CRB和0.2Hz MS稳态DLL不是同一估计器；没有从20ms推导等效环路记忆。`observationMatched=false`、`fullyMatched=false`，仍不能称为论文公平PRS曲线复现。
- 几何有效人数不是检测覆盖率，`detectionCoverageEvaluated=false`。主优化仍为无时钟OLS，新增TDoA没有改写优化问题。
- C接收机目前沿用载波中心PSD取样，不是实际匹配滤波带积分；矩形符号的对齐网格有采样零点。49.77MHz默认网格出现此现象，加上erfc下溢产生BER=0，不能解释为真实PRS完全无干扰。专门增加零点/偏移网格、人工非零PSD的BER失败检查，并将限制写入kernel/service/manifest和图9。C-to-C干扰仍与原主模型一样未计入。
- 显示图的1e-16 BER下限不改动导出原值。补偿峰值可行只是审计结论，未实施补偿或重算补偿后的干扰；实际C能量匹配标志仍为false。

### 验证

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); addpath('tests'); run_stage5_checks;"
```

63项全部通过：原数学17、DLL/PRS 9、对偶/分层10、投影9、新增PRS比较18。新增检查还覆盖参考不变性、未知时钟Schur等价、双向PSD显式求和、分块、干扰单调性、活动槽加权、严格服务/能量/峰值/零机会失败、不可用参考/计划及实际CSV导出。74个MATLAB源文件Code Analyzer诊断0。

三个完整main均为4BS/20用户，分别有效49.77MHz、20.16MHz及20.16MHz/每BS定位预算1mW压力分支。四条优化路径均保持20/20有效定位、通信QoS、独立目标、原始约束/KKT和对偶证书；每目录6策略、9图。PRS两种策略各20/20 oracle几何有效，不能作为检测覆盖率。

| 场景 | PRS总占时 | 全静默C机会比 | 全静默保持C能量所需活动功率 | 其他BS继续C的最坏活动BER |
|---|---:|---:|---:|---:|
| 49.77MHz，0.1W/BS | 17.3333% | 82.6667% | 1.20968W，超过1.1W上限 | 0，受PSD网格零点/数值下溢限制 |
| 20.16MHz，0.1W/BS | 18.6667% | 81.3333% | 1.22951W，超过1.1W上限 | 0.347262，超过0.01门限 |
| 20.16MHz，0.001W/BS | 1.3333% | 98.6667% | 1.01351W，峰值可行但未实施 | 0.227646，超过0.01门限 |

三组两种PRS策略均因正占时无法满足同符号时钟下的原满尝试速率合同；20MHz继续其他BS通信还违反活动BER。不能只选其较低的定位CRB而忽略通信损失。详细逐用户精度、逐载波BER和逐BS成本在CSV/run.mat，数值汇总在 `docs/STAGE5_VERIFICATION.json/.mat`。

当前源码有效结果：

- 49.77MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_111410_fa412f7e`
- 20.16MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_111441_7527b46c`
- 20.16MHz绑定预算：`F:/4Cell_MS_NOMA/Results/run_20261101_111500_a8706e86`

110938目录是捕获CSV汇总初始化错误的未完成试跑，不作为验收。111105/111133/111152为61项中间回归目录；其后增加网格零点/BER失败检查和图面修正。最终以上三个63项目录为准。汇总初始化已修复，并在单元检查中加入真实导出回归；不删除历史证据。

新增输出：`prs_comparison_summary.csv`及global_muting/other_bs_communication两种策略各自的service、communication、tdoa、ranging CSV（共9个新CSV）。run.mat保留协方差、差分矩阵、PSD核和每槽BER；manifest保留独立匹配标志和PSD取样限制。历史阶段2/3/4脚本显式关闭comparison，继续各自6/7/8图口径。旧cfg从新base创建后覆盖，不猜缺失字段。

本次验证已逐字节确认三个最终manifest的74个MATLAB源SHA256与当前源码一致，每目录9幅PNG；三组图9及默认场景全部9图已目视检查，零BER显示与多行标题已修正。独立归档记录为 `docs/STAGE5_SOURCE_HASH_VERIFICATION.json`。验收后仅整理文档和归档，不再修改MATLAB源码。

### 后续接手优先级

1. **DLL及接收端近似偏差分离**：固定场景/功率，分别改变正弦线性化、有限PSD和前端解释；同时独立核验C载波PSD点取样与实际匹配滤波带积分/符号对齐差异，保留当前点取样基线。不要直接用零点网格的BER结论推广到接收机。
2. **统计与敏感性**：先小样本回归，再50次Monte Carlo；统一场景、随机种子、失败状态和分母，加入q/噪声/带宽以及PRS网格/时序假设敏感性。未建立有限观测等价前继续标部分匹配。
3. **采样链路**：波形、实际测距、估计位置，再接首次lambda=1和几何反馈；oracle结果与估计几何结果分开。补偿通信功率或改变PRS调度时必须重算双向干扰/BER，不只改单个通过标志。
4. 作者原始步长/未公开参数仍无法确定；保留candidate/assumed与来源，不用“调到论文曲线”替代可验证模型。

## 阶段6验证与交接（2026-11-01，历史快照）

本轮按最新确认的优先级实现DLL独立方差差异量化，并核验PRS干扰在C接收滤波前后的解析差异。cfg/run版本6.0，默认6策略、11图。新增函数均位于Code/+msnoma，经main调用；原作者MS_NOMA未改动。公式与边界集中在 `docs/DLL_BIAS_AND_RECEIVER_AUDIT.md`。

### 已实现与定义

- `ranging.communicationPsd`、扩展 `integralKernels`、`biasKernels`、`asymptoticVariance`、`biasDecomposition`：固定同一参考功率/信道/环路，正弦精确/线性、有限C谱瓣/同有限载波的窄线谱、global/local前端共8种组合。复用已有完全同设置基线核。不是八次重新优化。
- 基线单独翻转三个因素，并导出其余因素四种上下文下的全部12条配对边。另给出六步可闭合路径：正弦→C谱宽→前端→其余A0/噪声/定位项渐近→有限N到名义B的N+1→去定位干扰。未将其余渐近误差混归到PSD，也未把不同分母的百分比直接相加。
- 默认central及用户[1,10,20]（12链路）；核/参考可显式配置。噪声、通信、定位、总方差分别保留，固定oracle OLS传播所选用户位置精度。这里的“偏差”是预测方差差异，不是实际时延估计的均值偏差。前端是解释敏感性，不宣布global或local为已确认的真实接收器。
- `receiver.rectangularPsdOverlap`、`prs.receiverAudit`：用两三角相关函数乘积的解析傅里叶积分计算矩形发送PSD与C矩形滤波重叠；与旧载波中心点取样结果并列保存，活动机会、C/PRS功率及其他假设固定。模型是平稳独立符号/平均相对时序/高斯等效干扰，不是相干对齐OFDM或原作者已确认的PRS时序。
- experiment.dllBias、CSV/图/manifest已接入runSingle/main。旧优化、旧点取样PRS策略和测距界保留；滤波分支不重算PRS测距、不把有限观测标为匹配。

### 验证和数值

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); addpath('tests'); run_stage6_checks;"
```

80项全部通过：原63项加DLL/接收滤波17项。新增覆盖PSD显式逐项求和、线谱还原式(58)、有限N/名义B差异、八组合和缓存、12配对、逐链路闭合/分量隔离、小间隔极限、中间载波前端对称、零通信、独立频域积分、重叠特殊值/对称缩放/长窗极限、真实CSV导出和禁用/缺参考失败。87个MATLAB源文件Code Analyzer诊断0。

三个完整main均为4BS/20用户、6策略/11图。四条优化路径继续通过20/20有效定位、通信QoS、独立目标重算、原始可行性、KKT及对偶证书，未引入回退。每次偏差研究八种核均收敛，基线与已有central积分逐值一致。

下表是**默认选中12条链路**、从exact/global/finite基线单独改变一个因素的最大绝对相对方差差，百分比针对方差，不是标准差：

| 场景 | 正弦线性化 | 有限C谱瓣→窄线谱 | global→local前端 | 最大积分加密变化 |
|---|---:|---:|---:|---:|
| 49.77MHz，0.1W/BS | 13.4917% | 0.0100002% | 10.4746% | 6.0815e-10 |
| 20.16MHz，0.1W/BS | 13.4886% | 0.0271785% | 10.4775% | 3.7066e-9 |
| 20.16MHz，0.001W/BS | 13.4837% | 0.0272056% | 10.4776% | 3.7066e-9 |

该场景的线性化与切局部前端都会降低预测方差；有限谱瓣转线谱略微升高。先线性化后再切前端，其差异会远小于基线单独切前端：有明显上下文交互，不能把13.49%和10.47%相加。程序保存12条边和有顺序的闭合链以避免此误读。三个场景的六步闭合残差在双精度下为0；最终旧两项式相对积分基线的最大绝对差分别13.5488%、13.4515%、13.3381%。

PRS其他BS继续通信策略、相同功率和时隙下，接收模型审计如下：

| 场景 | 旧C载波点PSD最坏活动BER | 新平均时序矩形滤波BER | 新活动BER是否≤0.01 |
|---|---:|---:|---|
| 49.77MHz，0.1W/BS | 0（采样零点/下溢） | 0.201402 | 否 |
| 20.16MHz，0.1W/BS | 0.347262 | 0.273655 | 否 |
| 20.16MHz，0.001W/BS | 0.227646 | 0.126292 | 否 |

滤波解释使49.77MHz原采样零点不再表现为无干扰，但20MHz峰值PSD被平均后BER降低；不能当成普适单方向校正。新模型仍未满足默认严格通信服务与统一观测合同，也不是实测解码结果。

### 当前有效目录和输出

- 49.77MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_113009_3dc53e95`
- 20.16MHz：`F:/4Cell_MS_NOMA/Results/run_20261101_113040_88d74fa4`
- 20.16MHz绑定预算：`F:/4Cell_MS_NOMA/Results/run_20261101_113104_6d9b54b4`

112721目录为1.2MHz/3用户的开发烟测，不作为完整规模验收。最新验证是 `docs/STAGE6_VERIFICATION.json/.mat`。三最终目录manifest的87个MATLAB源SHA256均与当前源码逐字节匹配；源码验收后只整理文档/归档。

新增8类CSV：dll_bias_variants、paired_contrasts、chain、positions、kernels、effect_summary，及prs_receiver_filter_links/service。图10展示受控偏差、六步归因、相同功率下的方差模型与加密检查；图11展示点/滤波PSD、活动BER和期望正确比特率。run.mat保存全部八组合、积分历史和接收模型假设；manifest记录状态及无采样接收机标志。旧阶段1–5脚本关闭阶段6附加模块，保留原4/6/7/8/9图数。

### 明确尚未实现

1. 多小区采样波形接收、实际时延估计、检测门限、失锁统计及从测距值估计位置。
2. 首次lambda=1、从估计位置重算几何及功率/位置迭代反馈。
3. PRS有限观测与MS环路稳态的统一观测接收机或有依据的等价定义。fullyMatched/observationMatched仍false。
4. 50次Monte Carlo及完整敏感性。现阶段三个确定性场景不能代替统计复现。

下一步优先接**有限采样波形→接收滤波/相关→测距→位置估计→首次lambda=1/几何反馈**，oracle单独保留。先用无噪声已知时延和独立几何例验证，再接噪声/干扰；采样率、长度、序列、相位/时钟及相对时序均须显式记录。PRS解析滤波只作此链路的核验基线，不能作为“接收机已实现”的状态依据。

阶段6最终归档：重新逐字节核验三组有效结果各87个MATLAB源文件SHA256全部一致，每组11幅PNG；默认全部11图及三组新增图10/11已目视检查，标题和图例可读、无裁切。独立记录见 `docs/STAGE6_SOURCE_HASH_VERIFICATION.json`。验收后仅更新进程文档、核验记录与交接副本，未修改已验证的MATLAB源码。

## 新增文献核验：arXiv 1910.04350v1（2026-11-01）

已将用户新增的30页早期稿作为补充验证来源登记；首页及arXiv官方版本页确认2019-10-10、PCJPA名称和作者信息。最终复现目标仍是JSAC 2022，本次未修改MATLAB源码或参数。详细对照见 `docs/ARXIV_VERSION_COMPARISON.md/.json`。

- **PRS来源已解释**：早期稿第17/30页引用[35]《Joint channel and time delay estimation for LTE positioning reference signals》，正是OneCell目录原有PRS论文；最终版[28]换成《Achievable Localization Accuracy...》。这是不同论文，不能混用公式。
- **模型边界**：早期式(10)为含其他站定位功率的三项ranging-factor，最终式(10)为省略定位干扰的两项方差；早期式(40)/最终式(33)、早期式(28)/最终式(24)仍保留原印刷符号/分组疑点，修正候选与独立证书不取消。
- **参数不覆盖**：早期图示0.8W及beta_th=8e-3不替换最终版0.1W/BS、Xi_th=1e-2；两版PRS表值分别为8.44/3.17m和8.11/3.28m（20/50MHz），旧表数值不作当前拟合目标。
- **算法线索**：早期稿补充三类次梯度更新和他站功率初始化；数值步长、完整采样测距/位置反馈仍未恢复。首次lambda=1仍以最终版明确分支为依据。

证据抽取和官方页面快照位于 `References/version_comparison/`；SOURCE_LEDGER新增来源并补全早期[35]身份关系。本次不增加算法验收项，既有80项/三场景验证状态不变。新来源及现有源文件哈希均已核对，记录见 `docs/ARXIV_SOURCE_REGISTRATION_VERIFICATION.json`。下一步仍是采样波形、接收、测距、位置估计和首次几何反馈。


## 阶段7验证与交接（2026-11-01，历史快照）

用户允许数值步长自定，因此未等待作者未公开的b1/b2/b3。已按链路实现：有限四站MS波形 -> 跨BS传播叠加 -> 真实样本捕获/测距 -> 二维位置估计 -> 首次lambda=1功率求解与估计几何反馈；oracle六策略独立保留。所有函数在Code/+msnoma，main统一调用，未改作者MS_NOMA文件。cfg/run schemaVersion=7.0；新重实验默认关闭，默认仍11图，开启采样/反馈为13图，另启步长为14图。

### 接口与工程定义

- waveform.makePacket/component/reference/transmit/propagate/communicationPsd：有限矩形BPSK定位、QPSK通信、完整OFDM符号、已知平坦CSI、单路径、同步站/已知接收钟；记录码、相位、采样、前端、名义/有限包功率和种子。默认作者32条L1CP码跨载波复用，独立随机码为另一个显式模式。
- receiver.estimateDelay：不知道真实时延/位置/复增益；先平均CP谱高斯等效预白化和相关捕获，再分数时延细化。门限只在假定高斯模型下解释，未标定真实多址干扰虚警。捕获状态不保证锁定正确。连续稳态DLL尚未实现。
- geometry.estimateRanges：只输入BS坐标、测距与有效标志；代数初值加阻尼Gauss-Newton。至少3非共线站，残差、秩、域外和求解失败分开；不裁剪域外结果，无效实测误差NaN，条件均值附完整覆盖分母。
- experiment.positionFeedback：首次lambda全部1，冻结几何进行带完整约束/证书的功率求解；接真实采样接收，仅全部用户/全部BS测距有效才更新下一轮几何。失败停止，不回填真值、不删除用户。默认共同码/通信数据与噪声；切换commonRandomNumbers仅切换每轮噪声种子。CSI固定，几何更新来自上一轮估计坐标。
- experiment.stepSensitivity：对缩放坐标下内外初始步长乘0.01/1/100，允许回溯。authorStepsKnown=false，验收依据原始约束、相对可听、KKT与对偶界。CSV向后检索最后实际接受步长，最终warm start零更新不再误标步长0。
- 修复sampledReceiver第7参数packet复用条件，新增独立检查保证main预建包不会被重复生成；保存实际chips/通信符号与种子以重放。

方法和完整接口见docs/SAMPLED_RECEIVER_AND_FEEDBACK.md、docs/MODULE_INTERFACES.md。传播真值只供物理信道及事后误差评价，不输入时延/位置估计器。采样前端是有限离散DFT带限模型，尚未验证连续射频抗混叠等价。

### 五个完整场景的结果

| 场景 | 有效测距 | 有效定位 | 条件平均实际位置误差/m | 反馈轮数 | 反馈状态 |
|---|---:|---:|---:|---:|---|
| 49.77MHz，0.1W/BS，作者码+C+P+AWGN | 79/80 | 19/20 | 0.0673333 | 1 | acquisition_or_position_failure |
| 20.16MHz，0.1W/BS，作者码+C+P+AWGN | 79/80 | 19/20 | 0.136949 | 1 | acquisition_or_position_failure |
| 20.16MHz，0.001W/BS，作者码+C+P+AWGN | 69/80 | 12/20 | 0.353809 | 1 | acquisition_or_position_failure |
| 49.77MHz，全区域，独立码/P，无C无AWGN | 79/80 | 19/20 | 0.0703923 | 1 | acquisition_or_position_failure |
| 49.77MHz，中心区域，独立码/P，无C无AWGN | 80/80 | 20/20 | 0.0439273 | 3 | converged |

表中误差仅对应采样参考central有效用户，是实际估计坐标的误差；不是主summary的oracle预测标准差，也不把失败计0。反馈每轮功率/测距/位置、最终覆盖和条件均值在对应feedback CSV与STAGE7_VERIFICATION.json中独立保存。第一轮停止时没有前后两次结果，位置/功率变化保持NaN，图13对应两个子图无曲线，不是漏写或零变化。中心控制末轮20/20，条件平均误差0.0439282m；不将该中心子区域结果冒充全区域性能。

**关键发现**：原0–200m区域即使关闭C和AWGN，其他P仍可能导致失败；近BS用户18的若干测距误捕获/未检测使位置残差超限。保持原门限、种子和失败记录。另以[50,150]×[50,150]m中心区域、无C无AWGN但保留所有P分量验证完整多用户反馈；它达到3轮、全部20用户和80链路有效并满足位置/功率双停止准则。两者是明确不同的场景，不用后者替代前者。平均PSD预白化不能保证消除有限码相干干扰，接下来应研究抗近远捕获/取消干扰。

额外开发诊断保留在Results/stage7_diagnostics/stage7_long_observation.mat及long_observation_user18.json：对原区域用户18把独立码请求长度改4095（实际4108），仍出现错误峰/位置失败；仅增加观测不保证单调改善。这是单用户诊断，不是完整统计敏感性，也不写成作者参数。正式五场景保持请求1023码片：49.77MHz实际1027码片/13个C符号/0.433333ms，20.16MHz实际1024码片/32个C符号/1.066667ms。

### 工程步长结果

20.16MHz/每BS 1mW绑定预算压力场景：

| 初始倍率 | 内层初始步长 | 外层初始步长 | 最后实际内层步长 | 最后实际外层步长 | 总内层迭代 | 目标/m² | 相对对偶间隙 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 1e-05 | 0.001 | 0.00115805 | 0.0393875 | 17304 | 3.13696 | 1.39051e-07 |
| 1 | 0.001 | 0.1 | 0.00152088 | 0.0307715 | 22414 | 3.13696 | 7.82854e-07 |
| 100 | 0.1 | 10 | 0.115 | 0.0316137 | 17772 | 3.13696 | 9.15493e-07 |

三组都通过完整证书，相对目标跨度7.7644e-07，小于1e-6。步长太小可能增加工作量，太大由回溯收缩；本次三个倍率均可用不能证明任意倍率/任意场景收敛。继续保留默认1倍，不按论文图拟合。实际耗时、总迭代、调用数、残差及完整轨迹在CSV/run.mat；缩放坐标参数不能直接等同作者b1/b2/b3。

### 验证与当前源码匹配目录

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); addpath('tests'); run_stage7_checks;"
```

97项全部通过：原80项加有限采样/接口/反馈17项。新检查覆盖独立时域整数移位、复幅度分数时延、有限PSD显式和、完整C符号谱零点、局部随机流、作者PRN复用、纯噪声/零信号、独立测距几何、缺站/秩亏/域外、真值隔离、packet复用、实际小场景反馈/强制失败及步长证书。104个MATLAB源文件Code Analyzer诊断0。

五个完整main均保留4BS/20用户、六策略；四个优化器继续通过通信QoS、目标独立重算、原始可行性/相对可听、KKT/对偶证书和相对central目标差<1e-6。每次后续lambda均与上一轮实际估计坐标重新计算的lambda逐值一致。分别检查80测距/20位置行、图数、步长和源哈希。当前源码对应目录：

- 49.77MHz，0.1W/BS，作者码+C+P+AWGN：`F:/4Cell_MS_NOMA/Results/run_20261101_122149_c2fd95ca`
- 20.16MHz，0.1W/BS，作者码+C+P+AWGN：`F:/4Cell_MS_NOMA/Results/run_20261101_122227_ee48d93b`
- 20.16MHz，0.001W/BS，作者码+C+P+AWGN：`F:/4Cell_MS_NOMA/Results/run_20261101_122317_9135a33a`
- 49.77MHz，全区域，独立码/P，无C无AWGN：`F:/4Cell_MS_NOMA/Results/run_20261101_122420_925f1ff5`
- 49.77MHz，中心区域，独立码/P，无C无AWGN：`F:/4Cell_MS_NOMA/Results/run_20261101_122505_b591a683`

五目录全部104个源SHA256已用独立Python逐字节复核，每幅PNG能解码；默认全部13图和五组新增图12/13及压力场景图14已目视检查。记录：docs/STAGE7_VERIFICATION.json/.mat、docs/STAGE7_SOURCE_HASH_VERIFICATION.json。源码冻结后仅更新文档/归档。

### 开发记录与下一步

- 120028至121329目录为接收/白化/完整通信符号和相关插值逐步调试记录；不能视为当前源码最终结果。121551/121631/121723为首次完整回归已通过的前三组，121829为原区域无C无AWGN但有P干扰的失败控制；其最初“无C无噪声应全成功”的验收假设被实际结果否定，未改检测/位置门限。随后修复packet复用、实际步长导出和图例，并新增中心区域完整反馈正例；以上均由当前五场景重新回归。
- **优先1：有限码干扰与接收稳健性。** 固定种子/全区域，分离同载波P、跨载波P、通信干扰、白化失配与错误峰，研究可验证的多峰/联合捕获或干扰抵消；不能用已知真值时延取消信号。先局部诊断，再重跑全体用户和失败分母。
- **优先2：采样率/前端/观测敏感性及小样本统计。** 维持作者码和独立码两条明确模式，记录近远距离、码复用、时序、检测/误捕获；再扩50次Monte Carlo。条件误差必须与覆盖并列。
- **优先3：稳态DLL与PRS同观测接收。** 有限捕获不等于0.2Hz稳态跟踪；需显式环路启动/收敛/失锁和PRS采样估计及通信成本。PRS observationMatched/fullyMatched维持false，不能因为MS新增采样而改写。
- **持续边界：** 已知接收机时钟与CSI、二维单路径、q/Omega未知、修正候选非原作者CP4A、单场景非论文统计复现；作者原始步长未知继续记录但不阻塞工程推进。


## 阶段8验证与交接（2026-11-01，当前入口）

本轮继续解决阶段7的有限码/近远干扰失败，实现了**额外已知全部P码及相对功率**的联合拟合/估计P消除分支和受控发射消融。schemaVersion=8.0；所有模块在Code/+msnoma，经main统一调用。原六策略oracle、逐P分量接收器、原门限/种子/全区域继续保留；新接收器默认关闭。原作者MS_NOMA文件未改。

### 实现及接口

- waveform.positioningBank：仅已知P码/功率构造每BS全部P模板；通信数据不参与。
- receiver.estimatePChannels：逐次相关初值、QR拟合未知复增益、投影Gauss-Newton细化共同BS时延、回溯检查实际残差下降。接口无真实时延、增益、场景或用户坐标。多址非凸局部拟合可能收敛到错误解，converged不等于正确锁定。
- receiver.reconstructPositioning：只按有效估计值重建P干扰；物理层真值仅用于仿真传播与事后误差。
- 代理精度门控：延迟std≤0.25采样、固定时延条件下相对复增益std≤0.25，全部BS通过才消除。此门控使用假定高斯PSD信息，不读取真实误差，不保证强C/低预算场景有益。拟合或门控失败时原接收器逐值保留。
- sampledReceiver/positionFeedback：消除后加回目标P，仍用真实接收样本执行目标检测/测距；改用残余C+噪声预白化，因此同时改变权重和检测分母，明确作为另一接收模式。没有用真实时延取消信号或填补失败。首次lambda=1、下一轮lambda取上轮实际估计位置、全部用户/全部BS有效要求保持。
- experiment.interferenceAblation：desired_only、same_carrier_P、same_BS_P、all_P、all_P_C、all_P_C_AWGN六类**发送掩码仿真**，使用原接收器，不当作业务可选的oracle抑制。PSD随实际发送掩码匹配变化，不能当成固定核的唯一因素归因。
- 新增P信道CSV、消除应用状态、拟合/门控诊断图及消融图；完整历史在run.mat。docs/ESTIMATED_P_CANCELLATION.md、MODULE_INTERFACES.md和README记录参数/接口/复跑方法。

### 七个完整场景与配对结果

| 场景 | 有效测距 | 有效定位 | 条件平均实际误差/m | 消除门控通过用户 | 反馈轮数/状态 |
|---|---:|---:|---:|---:|---|
| 原逐分量49.77MHz/0.1W | 79/80 | 19/20 | 0.0673333 | 0/20 | 1 / acquisition_or_position_failure |
| 原逐分量20.16MHz/0.1W | 79/80 | 19/20 | 0.136949 | 0/20 | 1 / acquisition_or_position_failure |
| 原逐分量20.16MHz/1mW | 69/80 | 12/20 | 0.353809 | 0/20 | 1 / acquisition_or_position_failure |
| 估计消除49.77MHz/0.1W | 80/80 | 20/20 | 0.0458835 | 20/20 | 3 / converged |
| 估计消除20.16MHz/0.1W | 80/80 | 20/20 | 0.0823965 | 20/20 | 3 / fix_limit_reached |
| 估计消除20.16MHz/1mW | 42/80 | 6/20 | 2.87251 | 19/20 | 1 / acquisition_or_position_failure |
| 原全区域独立码/无C无AWGN，估计消除 | 80/80 | 20/20 | 1.7083e-06 | 20/20 | 3 / converged |

前三组原接收器与第四至第六组在相同场景/功率/全部码/通信符号/噪声种子下配对；程序逐值核验，不通过门控的用户也逐值核验仍为原接收器结果。误差是采样参考central的**有效用户条件实际坐标误差**，不是oracle预测标准差；失败保持NaN。反馈实际覆盖/误差单独保存在STAGE8_VERIFICATION及feedback CSV，不用参考central的覆盖代替反馈覆盖。

正常预算恢复缺失用户，完整反馈状态见表；原全区域无C无AWGN控制保留所有P分量，验证估计消除与位置反馈确实可以在整个原区域执行。**低预算分支依然明显退化，不能默认推广新接收器。** 局部拟合/高斯精度门控不足以识别有限通信干扰、估计偏差和数据相关的消除误差。保留反例，没有根据真实位置误差逐用户挑选更好接收器，也没有把低预算结果藏起来。

### 用户18的六类发射诊断

49.77MHz/作者码/同一central功率，误差仅事后评价；未检测为NaN（表内“—”）：

| 发射掩码 | BS1误差/m | BS2误差/m | BS3误差/m | BS4误差/m |
|---|---:|---:|---:|---:|
| desired_only | 3.10329e-06 | -2.44053e-07 | 3.08174e-06 | 4.09006e-06 |
| same_carrier_P | -0.0219247 | -0.00182487 | 0.0195387 | 0.00807031 |
| same_BS_P | 0.103462 | 0.119199 | 0.0901046 | 0.119757 |
| all_P | 0.0251149 | 0.0838706 | 55.0887 | -21.6612 |
| all_P_C | — | 0.0533226 | -194.114 | -148.133 |
| all_P_C_AWGN | — | 0.0533473 | -194.113 | -148.131 |

只有目标、同载波P或同BS P时误差小；全部P发送后BS3/4出现大错误峰，再加C后错误更大且BS1未检出。说明不能仅用同载波近远比解释此样本；尚不能由六个非线性捕获结果唯一分解各干扰贡献。关键是多分量/多BS有限波形叠加及平均谱接收失配，下一步应在固定接收核/已知标注的合成场景继续分析。

### 验收与当前源码目录

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); addpath('tests'); run_stage8_checks;"
```

112项全部通过：原97项加新15项；112个MATLAB源文件Code Analyzer诊断0。新增覆盖独立整数时域移位、分数近远合成估计、残差下降、通信数据污染隔离、零信号/秩/失败、实际采样小场景、真实导出、失败保留原结果、代理精度门控和PSD方差缩放。七个完整4BS/20用户main、四个优化器通信QoS/可行性/相对可听/KKT/对偶证书、同随机数配对及真实估计几何重算均通过。原全区域无C无AWGN控制要求全部80测距/20用户有效且反馈达到双停止准则；未放宽门限。

- 原逐分量49.77MHz/0.1W：`F:/4Cell_MS_NOMA/Results/run_20261101_125140_8ca850bb`
- 原逐分量20.16MHz/0.1W：`F:/4Cell_MS_NOMA/Results/run_20261101_125228_a750b5aa`
- 原逐分量20.16MHz/1mW：`F:/4Cell_MS_NOMA/Results/run_20261101_125306_14b5e24b`
- 估计消除49.77MHz/0.1W：`F:/4Cell_MS_NOMA/Results/run_20261101_125515_c67a79f2`
- 估计消除20.16MHz/0.1W：`F:/4Cell_MS_NOMA/Results/run_20261101_125638_e655714d`
- 估计消除20.16MHz/1mW：`F:/4Cell_MS_NOMA/Results/run_20261101_125728_4a1d7ea5`
- 原全区域独立码/无C无AWGN，估计消除：`F:/4Cell_MS_NOMA/Results/run_20261101_125930_36c3eb87`

每组图数13或14；默认仍六策略/11图。消除、消融各追加一图。所有七目录共112个源SHA256用独立Python逐字节核验一致，全部PNG可解码；默认全部图和七组新增接收/反馈/消除/消融图已目视检查。详见docs/STAGE8_VERIFICATION.json/.mat及STAGE8_SOURCE_HASH_VERIFICATION.json。源码冻结后仅整理文档/归档。

### 明确保留的限制与接手优先级

1. **低预算/强C稳健性优先。** 当前高斯精度门控不是充分条件。需要在不看真值误差的前提下，研究保留样本的拟合验证、检测校准、消除误差传播、多峰假设或更适当的有限通信统计模型；不能按事后较低误差自动切换模式。已比较保留CP白化与切换C白化，开发反例保留在Results/stage8_diagnostics；两者都没有根治低预算。
2. **观测/前端/采样和小样本统计。** 固定全区域、全部用户与种子，先检查码长、采样率、码模式和模型失配，再扩50次Monte Carlo；同时保存错误峰、检测失败、位置残差/域外和完整分母。
3. **稳态DLL与PRS同观测链路。** 目前仍为有限捕获；没有连续环路启动/稳态/失锁、PRS采样接收与统一通信服务成本，PRS fullyMatched/observationMatched继续false。
4. **新接收器附加知识/误差界。** 全部P码和相对功率是新增工程信息；通信数据没有使用，平均PSD仍用理想CSI。消除后的估计误差未传播到虚警/CRB，noiseOnlyRangeCrbM2=NaN，局部GaussianEquivalent界不是正确完整后消除CRB。
5. **目标和优化边界不变。** 已知接收机时钟、平坦单路径、q/Omega未确认、修正候选非作者逐行CP4A；作者数值步长未知继续记录，工程步长可自定且已有证书/敏感性。当前结果不等同论文统计复现。


### 阶段8补充：20.16MHz反馈轮数敏感性

七场景中的20.16MHz/正常预算/估计P消除分支，在默认3轮上限时为fix_limit_reached：第3轮最大位置变化0.000919261m已低于0.1m容差，但相对功率变化0.00106168仍略高于1e-3，因此没有误标收敛。

另仅把cfg.experiment.positionFeedback.maxFixes从3改为5，保留原随机数、接收器、精度门控和全部容差，第4轮达到双停止准则：最大位置变化1.25467e-05m，相对功率变化1.48609e-05，20/20有效，末轮平均实际位置误差0.0822002m。首次lambda=1和后续估计几何逐值重算均核验。它是增加计算轮数的明确敏感性实验，不覆盖原三轮未收敛记录；默认maxFixes仍3。

补充目录：`F:/4Cell_MS_NOMA/Results/run_20261101_130239_5d09fd17`。记录：docs/STAGE8_FEEDBACK_LIMIT_SENSITIVITY.json。14图已检查，112个源SHA256与七组正式回归及当前源码一致；仅更新文档/归档。当前共七个正式场景加一个额外轮数对照，单元检查仍为112项，未把额外main重复计作新单元测试。


## 阶段8.1：残差主导时拒绝估计P消除（已验收）

本轮在接收链路新增独立函数`receiver.cancellationGate(fit,opts)`，先检查阶段8高斯代理精度，再要求加权拟合残差/接收能量不超过可配置上限（默认0.5）。该上限定义P拟合至少解释一半加权能量的保守工作区，是显式工程选择，不是论文参数或正确峰/精度保证。不输入真实时延、增益、坐标、测距误差或通信数据来决定消除；平均PSD仍沿用理想CSI假设。未通过则原样使用原接收器，所有失败和分母保留。

新增门控保存precisionGatePassed/residualGatePassed/maximumRelativeResidualEnergy，CSV/manifest/诊断图同步。上限1回放旧精度门控，旧反例不覆盖；新接收器仍默认关闭。schema升至8.1，旧cfg从base重建；原作者代码未改。接口及工程边界见docs/STAGE8_1_RESIDUAL_GATE.md。下表预算均为每BS总P预算，实际误差为有效用户条件均值；失败不计0。

| 场景 | 有效测距 | 有效定位 | 条件平均实际误差/m | 消除通过用户 | 反馈状态/轮数 |
|---|---:|---:|---:|---:|---|
| 20.16MHz/1mW原接收 | 69/80 | 12/20 | 0.3538092 | 0/20 | acquisition_or_position_failure / 1 |
| 20.16MHz/1mW旧消除，上限1 | 42/80 | 6/20 | 2.872514 | 19/20 | acquisition_or_position_failure / 1 |
| 20.16MHz/1mW新增门控，上限0.5 | 69/80 | 12/20 | 0.3538092 | 0/20 | acquisition_or_position_failure / 1 |
| 49.77MHz/0.1W新增门控 | 80/80 | 20/20 | 0.04588352 | 20/20 | converged / 3 |
| 20.16MHz/0.1W新增门控 | 80/80 | 20/20 | 0.08239653 | 20/20 | fix_limit_reached / 3 |
| 全区域独立P码/无C无AWGN新增门控 | 80/80 | 20/20 | 1.708296e-06 | 20/20 | converged / 3 |

低预算新增门控在此固定场景拒绝全部用户的消除：19名通过精度门控的用户因残差被拒，另一名未通过原精度门控。测距和定位逐值等于原接收器，首次lambda=1反馈测距也逐值一致。因此避免了本次旧消除的6/20退化，保留原12/20覆盖；**并没有解决剩余8人的定位失败，也没有保证其他场景始终改善**。低预算0.25/0.5/0.75三组上限重跑同样观测，全部为零消除且与默认门控测距逐值一致；它只说明此单场景在该区间不敏感，不能推广为通用最优参数。

两个正常带宽及无C无AWGN控制保留20/20定位；与历史8.0对应成功分支的测距/定位CSV逐值一致。49.77MHz正常预算及无C控制3轮反馈收敛；20.16MHz默认3轮仍为fix_limit_reached，不放宽容差或误标收敛。此前增至5轮并在第4轮收敛的历史敏感性记录仍保留。

验收命令：

```powershell
& 'D:\DevTools\MATLAB2024b\bin\matlab.exe' -batch "cd('F:/4Cell_MS_NOMA/Code'); addpath('tests'); run_stage8_1_checks;"
```

原112项与新增5项共117项通过，115个MATLAB源Code Analyzer为0诊断。六次完整main均检查四条优化路径的原始约束、相对可听、KKT/对偶证书，以及实际估计几何反馈；独立Python SHA256逐字节核验六组manifest全部与当前源码一致。全部PNG可解码，六组新增接收/反馈/消除图已目视检查。核验文件：docs/STAGE8_1_VERIFICATION.json/.mat、STAGE8_1_SOURCE_HASH_VERIFICATION.json。历史阶段8验收保留，不冒充当前源码哈希。

本轮结果目录：

- 20.16MHz/1mW原接收：`F:\4Cell_MS_NOMA\Results\run_20261101_134333_625aa67b`
- 20.16MHz/1mW旧消除，上限1：`F:\4Cell_MS_NOMA\Results\run_20261101_134441_0dec3af6`
- 20.16MHz/1mW新增门控，上限0.5：`F:\4Cell_MS_NOMA\Results\run_20261101_134541_e2ed74f5`
- 49.77MHz/0.1W新增门控：`F:\4Cell_MS_NOMA\Results\run_20261101_134845_6612bf7d`
- 20.16MHz/0.1W新增门控：`F:\4Cell_MS_NOMA\Results\run_20261101_135026_7f3fccf2`
- 全区域独立P码/无C无AWGN新增门控：`F:\4Cell_MS_NOMA\Results\run_20261101_135250_e52d705b`

接手优先级：先扩展观测码长/码模式/前端和随机种子敏感性，保留错误峰/未检出/位置残差和全分母，再研究独立留出验证、消除误差传播及多峰估计；当前残差门控不等于这些模块已经完成。连续稳态DLL、PRS同观测采样公平比较、50次Monte Carlo仍待实现；作者b1/b2/b3/Omega未知不作为工程实现阻塞。
