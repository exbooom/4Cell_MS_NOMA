# 当前实现接口（阶段4）

更新：2026-11-01。此文件记录已存在的接口；未实现的规划函数见PROGRESS.md。以下函数包前缀均为 `msnoma.`。

## 输入结构与约定

- `cfg=config.base('smoke_assumed')`：原文参数在paper，工程假设在assumed，模型分支在models，数值策略在numerics；io控制输出与图窗。`config.validate(cfg)`按已实现分支检查。
- `num=model.numerology(cfg)`：K/M/N/G、TcS/TpS、deltaFcHz/deltaFpHz、requestedBandwidthHz/effectiveBandwidthHz、BHz/BfeHz、PcW/c0Mps。BHz为有效带宽。
- `powerW`始终为K×M；所有信道字段为功率增益。展平顺序为MATLAB列优先。
- 核心计算函数不修改全局随机流；experiment、report和main负责调度及IO。

## 场景和跨基站信道

```matlab
scene = scenario.make(cfg,num,stream)        % stream可省略，使用局部固定种子
scene = scenario.fromSnapshot(snapshot,num)
[gain,distanceM] = scenario.freeSpaceGain(txPosM,rxPosM,prop)
ch = scenario.channels(scene,num,prop)
```

scene字段：bsPosM[K×2]、pUserPosM[M×2]、cUserPosM[K×N×2]、seed、domainM。定位用户均匀抽样，通信用户在最近BS服务区内独立抽样，排除最小距离。

freeSpaceGain输出是接收机×发射机。ch字段：gCServing[K×N]、gPC[K×N×K×M]、gP[K×M]、gCP[M×K×N]及距离。最低距离违规会报错，不悄悄截断增益。

## 通信

```matlab
kernels = model.spectralKernels(num)
berModel = comm.makeBerModel(cfg.models.ber)
[A,rowMap] = comm.interferenceMatrix(ch,num,kernels)
I = comm.interference(powerW,A,num)
BER = comm.ber(gServing,PcW,TcS,N0WPerHz,I,berModel)
qos = comm.qosLimit(gServing,PcW,TcS,N0WPerHz,threshold,berModel)
```

W为N×M，A为KN×KM，I与BER为K×N。qos包含limitWPerHz[K×N]、cleanBer、cleanFeasible、threshold和model。负限值保留，表示无定位干扰时也不能满足QoS。

BER具名分支：qpsk_bit_energy和paper_literal。前者采用Eb=gPcTc/2、分母N0+I；后者按印刷式(2)，Gamma/gamma系数明确在配置中。正向、逆向必须用同一个berModel。

## 测距

```matlab
gCeq = ranging.equivalentCommunicationGain(ch,num)
factor = ranging.factorTwoTerm(ch,num,N0WPerHz,dll)
two = ranging.varianceTwoTerm(powerW,factor)
three = ranging.varianceThreeTerm(powerW,ch,num,factor)
```

gCeq为M×K，对应式(62)有限载波求和。factor保存RtauS2W、RrangeM2W及噪声/通信分量，均为K×M。two/three保存delayVarianceS2、rangeVarianceM2、rangeStdM。three另有positioningDelayVarianceS2，只用于求解后评价，不把随其他功率变化的第三项伪装成固定R。

## 几何与可听性

```matlab
geom = geometry.weights(scene,bsMask,geometrySpec)
hear = geometry.hearability(powerW,ch,hearSpec)
position = geometry.positionAccuracy(rangeVarianceM2,geom)
```

bsMask是K×M逻辑矩阵。geom包括lambda[K×M]、各用户H/bsIds、rank、conditionNumber、valid。QR求解、SVD判秩，rank不足不产生伪精度。

hear包含mask[K×M]、receivedPowerW、ratioToStrongestOther、count[M×1]与enoughBs。只用相对功率条件，零功率不算可听。position包含positionVarianceM2、positionStdM、covarianceM2与valid；experiment.evaluate再结合minBs判定最终有效性。

## 优化

```matlab
problem = optim.buildProblem(A,qos,factor,geom,ch,cfg,num)
feasibility = optim.checkFeasibility(problem,cfg.numerics)
[f,gradient,H] = optim.objective(z,problem)
solution = optim.solveCentral(problem,cfg.numerics,feasibility) % 第三参数可省略
solution = optim.equalPower(problem,policy)
diagnostics = optim.diagnostics(problem,solution)
```

problem包含原始Aqos/bqos、Bbudget/budgetW、Chear、c、定义域下限/上限、原始Aineq/bineq、scaledA/scaledB/scaledC、功率/目标/行缩放、变量/约束映射。正功率原始目标为sum(c./p)。z=p/powerScaleW。

feasibility状态：noise_alone_violates_qos、all_pair_hearability_impossible、linear_constraints_infeasible、linear_solver_residual_failed或feasible_positive_domain_with_numeric_lower_bound。成功提供scaledPoint与powerW；失败保留原因。

solution保存powerW、objectiveM2、exitflag、multipliers、history、solverOutput、diagnostics。集中式仅在退出、原始可行性、驻点和互补松弛全部通过时标converged_verified；其他状态明确区别。无可行解时powerW为NaN，实验层不制造性能输出。

diagnostics保存分组SI残差、scaledResidual、maxNormalizedViolation、KKT驻点/互补残差和最小乘子。基准策略没有乘子，其KKT指标为NaN。等功率基准会保留通信和可听性违例。

## 实验和报告

```matlab
run = experiment.runSingle(cfg,algorithmIds,snapshot) % 后两项可省略
evaluation = experiment.evaluate(powerW,context)
summary = report.aggregate(run)
[x,F] = report.unconditionalCdf(positionStdM,valid)
figures = report.plotRun(run,cfg.io)
runDir = report.saveRun(run,outputRoot)
runDir = main(cfg)                         % cfg可省略，打开桌面图窗
```

run.algorithms每项包含id、solution、evaluation；各算法共享坐标、信道、A、R和初始全BS几何。evaluate从原功率重算BER、两/三项测距、可听集合及重新计算后的几何。

summary为table：algorithm/status、nTotal/nValid、coverage、conditionalMeanM/conditionalMedianM、objectiveMeanM2、maxBer、communicationQosSatisfied及残差。CDF分母始终为全部用户数，无效样本不填0。

saveRun保存完整run.mat、配置JSON、源码SHA256和工具箱版本、摘要、坐标及各策略CSV。plotRun绘图后由main导出PNG；showFigures=false时导出完关闭隐藏图窗。

## 最小复用示例

```matlab
cfg = msnoma.config.base;
cfg.io.showFigures = false;
run = msnoma.experiment.runSingle(cfg,{'qos_scaled_equal','central'},[]);
replay = msnoma.experiment.runSingle(cfg,{'central'},run.scene);
```

重放场景必须与当前K/M/N一致。改变带宽导致N改变时，不可把原通信用户张量直接塞入新维度；定位用户的默认局部随机流可保持相同，而通信用户数量与抽样过程应另行记录。

## 有限带宽DLL积分（阶段2）

```matlab
kernel = ranging.integralKernels(ch,num,dll,opts)
exact = ranging.varianceIntegral(powerW,ch,kernel,N0WPerHz)
```

opts使用 `cfg.models.rangingIntegral`，包括enabled、userIds、frontendMode、discriminator、initialSamplesPerCarrier、maxRefinements、relativeTolerance、requireConvergence和chunkSize。直接调用integralKernels须提供非空userIds；runSingle把空列表自动扩为unique([1,ceil(M/2),M])。frontendMode支持global_frontend/local_symmetric，discriminator支持exact_sine/small_spacing。

kernel保存userIds[1×U]、A0PerS/noiseIntegral/communicationIntegralS/positioningIntegralS[1×U]、localBoundsHz[U×2]、relativeChange/samplesPerCarrier/converged[1×U]、refinementHistory{U}及列名historyColumns。可跨功率策略复用，但不能跨信道、numerology、DLL或前端设置复用。

exact的delayVarianceS2、rangeVarianceM2、rangeStdM以及noise/communication/positioningDelayVarianceS2均为K×U，列顺序严格对应userIds；零接收功率的总方差为Inf。默认requireConvergence=true，未收敛则明确报错。run.algorithms(i).evaluation.integral另有twoTermRelativeVarianceError/threeTermRelativeVarianceError[K×U]和选中用户positionStdM[U×1]，误差定义为approx/integral-1。未选用户没有精确积分结果。

run.integralKernel保存共享缓存，cfg.models.rangingIntegral.enabled=false时为空，evaluation不含integral字段。积分为求解后评价，不进入优化目标。

## PRS（阶段2）

```matlab
allocation = prs.allocation(spec,bandwidthHz)
budget = prs.powerBudget(averagePowerW,measurementPeriodS,peakLimitW,allocation,spec)
crb = prs.delayCrb(allocation,cn0Hz,symbolCount,c0Mps)
benchmark = prs.positionBenchmark(scene,ch,num,averagePowerW,noisePsdWPerHz,spec,geometrySpec)
```

spec使用cfg.models.prs。allocation支持generalized_ofdm及lte_single_symbol_reference（需要LTE Toolbox）。返回activeK/weightsSquared[Na×1]、numCarriers/numPilots、symbolPeriodS、occupiedBandwidthHz、isStandardLtePattern、referenceSymbol；sum(weightsSquared)=Nc。LTE图样模式只表示单符号参考，不代表完整标准occasion。

averagePowerW为K×1、每BS所有P用户的定位功率总和。budget返回averagePowerW、energyPerMeasurementJ、symbolCount、instantaneousPowerW、observationTimeS、dutyCycle、slotStartS、feasible[K×1]与jointDutyCycle等。fit_peak_limit只支持通用符号；single_reference_symbol可能因能量/峰值不兼容而不可行。jointDutyCycle>1使全部BS不可行，不通过截断能量修复。预算是广播一次，不乘M。

cn0Hz为K×M，symbolCount为可广播的K×1。crb返回delayVarianceS2/rangeVarianceM2/rangeStdM、fisherInformationPerS2[K×M]、snrLinear、frequencyMoment、phaseModel；按[28]式(5)已知相位模型，符号信息累加。

benchmark包含allocation、budget、crb、geometry、position、referencePowerPolicy、status、detectionCoverageEvaluated=false。不可行时position=[]，无crb；可行状态computed_assumed_AWGN_broadcast_bound，几何有效性不代表检测覆盖。position使用主路径无时钟OLS传播，不是[28]的TDoA位置CRB。budget.communicationQosMatched=false。

run.prsBenchmark默认以central功率为参考；由spec.referenceAlgorithm选择，缺失时status=reference_solution_unavailable，不改用其他策略。关闭时status=disabled；资源不可行时status=observation_energy_or_peak_limit_infeasible。

## 新增报告接口与版本

```matlab
report.exportPrecision(run,runDir)
figures = report.plotPrecision(run,cfg.io)
```

阶段2时run/config的schemaVersion为2.0；阶段4为4.0，当前为5.0（见末节），旧cfg请从新base重新创建后覆盖自定义字段。阶段2main拼接原4图与新2图；阶段3另有优化图7，阶段4追加投影图8，阶段5追加PRS审计图9，关闭模块或缺失可用结果时相应附图不创建。

exportPrecision保存逐链路积分对照、积分收敛、每BS广播能量/符号数/时隙起点、PRS几何结果和链路CRB。manifest保存DLL模式/抽查用户/收敛与PRS状态；完整原始结构和数值历史保存在run.mat。run_stage2_checks另行导出global/local前端对照，不属于每次main的常规输出。

## 完整耦合对偶与分层候选（阶段3）

```matlab
[negativeDual,gradient,z,H] = optim.dualValue(lambda,problem,rows,fixed)
state = optim.dualActiveSet(problem,opts,variableRows,fixed,warm)
solution = optim.dualSolution(problem,opts,state,algorithm)
solution = optim.solveDual(problem,opts,feasibility)
solution = optim.solveCp4a(problem,opts,feasibility)
audit = optim.carrierGroupAudit(num)
report.exportOptimization(run,runDir)
figures = report.plotOptimization(run,cfg.io)
```

opts为cfg.numerics，需包含dual和cp4a新结构。solveDual/solveCp4a第三参数可省略，会独立LP预检。buildProblem新增G数值用于载波分组审计，所有原始缩放定义不变。

dualValue默认rows=全部行、fixed=0；lambda是对应rows的非负乘子。fixed为全行长度向量，只能在非variableRows上给固定价格；当前rows覆盖其同名坐标。返回负对偶标量、rows维梯度、全变量z及rows×rows半正定Hessian。标量内层精确处理负价格与上下界。

dualActiveSet仅优化variableRows，其余价格取fixed；warm为全行长度历史价格，工作集合只取variableRows中正warm价格和可变预算行。state.lambda为全部行长度，z为KM×1；另有negativeDual/scaledObjective/scaledGap、projectedGradient、innerComplementarity、maxVariableViolation、activeRows、passes/iterations、exitflag/output和history/historyColumns。maxVariableViolation只针对内层可变价格对应行，最终全部约束在dualSolution重算。

solveDual输出沿用central公共字段；另有dual完整state、dualBoundM2、dualityGapM2、relativeDualityGap、elapsedS、usedCentralFallback=false。multipliers.ineqlin/lower/upper均为既有缩放问题的乘子，可直接调用optim.diagnostics。仅原始可行、KKT/互补、非负乘子和相对间隙都通过才标converged_verified。

solveCp4a只接受cp4a.mode=full_coupling_hierarchical_candidate。输出algorithm=cp4a_consistent_candidate、literalPaperAlgorithm=false、innerMethod/outerMethod、geometryMode、innerEvaluations/innerIterations和groupAudit。外层预算价格梯度为scaledB_budget-scaledA_budget*z，内层固定该预算价格求完整QoS/可听对偶。内层warm start与central无关。内层检查失败返回hierarchical_inner_dual_not_verified、NaN功率和错误说明。其他未收敛状态与dualSolution一致，禁止静默回退。

carrierGroupAudit输出printedWidth/printedFirst/printedLast、printedCarrierCount、actualCarrierCount、outOfBandCount、uncoveredCount、isExactPartition及solverUsesPrintedGroups=false。仅用作原文审计，不能用其截断集合替代完整干扰矩阵。

runSingle新增ID dual/cp4a_consistent_candidate；阶段3默认5种策略；当前另含cp4a_subgradient_candidate共6种。阶段3的run/config schemaVersion=3.0，run新增optimizationLabel。report.saveRun保存完整独立证书并导出optimization_comparison.csv、paper_carrier_group_audit.csv；每种新算法的迭代历史、全部约束和评价CSV沿用通用接口。阶段3追加优化图，共7图。历史stage1/2验收显式使用旧3策略。

接口说明中的“候选”不能省略：内层解析梯度/Hessian内点法及外层SQP不等于论文的subgradient Algorithm 1，几何仍是固定oracle。

## 独立投影价格与物理证书（阶段4）

```matlab
state = optim.projectedDual(problem,opts,variableRows,fixed,warm)
solution = optim.solveSubgradient(problem,opts,feasibility)
report.exportSubgradient(run,runDir)
figures = report.plotSubgradient(run,cfg.io)
```

opts为cfg.numerics，需包含subgradient结构。projectedDual可省略后三参数，默认优化所有约束价格；分层时variableRows仅含QoS/可听行，fixed为全行长度非负预算价格，warm为全行长度历史价格。省略使用零向量；固定价格在variableRows上必须为0。rowScale/列顺序不变。

state沿用z/lambda/negativeDual/scaledObjective/scaledGap/exitflag/output/history/historyColumns，新增converged、reason、iterations/evaluations、coordinateScale、innerRelativeComplementarity和maxRelativeHearabilityViolation。inner history列为迭代、对偶界、恢复目标、违例、互补、坐标投影梯度、接受步长、求值数。默认每100步采样和最终点。

solveSubgradient返回algorithm=cp4a_subgradient_candidate、通用解和证书，以及innerCalls/innerIterations、innerTrace{call}（budgetPrice与完整state）、outerPriceHistory[(outer+1)×K]、stepConfiguration、usedPriceOptimizer=false、usedCentralFallback=false、literalPaperAlgorithm=false。feasibility第三参数可省略，仍有共享LP预检。它不调用fmincon更新乘子。

成功同时要求内层停止门禁和共享原始/对偶证书。失败状态带subgradient_前缀和具体inner/outer原因；保留最后已接受解和未验证轨迹。原问题不可行沿用旧状态和NaN功率。runSingle会保留有限功率的物理评价及明确失败status，统计成功算法时不得丢弃该状态信息。

buildProblem新增hearabilityComparisonTolerance；diagnostics新增maxRelativeHearabilityViolation，定义为max(0,Cp/max(C_positive*p,realmin))。dualSolution将该字段纳入converged_verified门禁，与实际可听比例一致。内层投影停止还使用更严格的相对可听容差。该变化对所有对偶候选共享，未修改geometry.hearability门限。

stepRule=backtracking时支持FISTA/梯度重启和内层稳定二次界、外层Armijo；stepRule=diminishing时不使用动量，步长指数在(0.5,1]。warmStartInner控制是否每次内层从零开始。作者原步长未知，不能省略candidate名称。

阶段4的cfg/run版本4.0、6策略、8图；stage3验收固定旧5策略。report.exportOptimization新增usesFmincon、stepRule和相对可听违例。逐内层调用包含被拒绝外层试探，outer_prices和主solver_history只含已接受外层迭代。

## PRS相关TDoA、双向干扰及服务（阶段5，当前接口）

```matlab
position = geometry.tdoaAccuracy(scene,rangeVarianceM2,spec,bsMask)
kernel = prs.interferenceKernels(ch,num,allocation,spec)
crb = prs.delayCrbInterference(allocation,gP,budget,noisePsdWPerHz,kernel,c0Mps)
service = prs.serviceAudit(ch,num,budget,kernel,noisePsdWPerHz,berModel,berThreshold,spec,policy,referenceBer)
comparison = prs.comparison(scene,ch,num,benchmark,reference,cfg,berModel)
report.exportPrsComparison(run,runDir)
figures = report.plotPrsComparison(run,cfg.io)
```

新增配置 `cfg.models.prs.comparison`：enabled、chunkSize、minimumAttemptedRateRatio（默认1）、communicationPeakLimitW（默认初始1.1W）及tdoa子结构。tdoa含referencePolicy=`minimum_variance`或`fixed`、referenceBs、rankTolerance。旧cfg从base重建，不自动猜缺失字段。修改paper定位预算不会自动覆盖PRS/C峰值，需显式保持一致或记录不同限制。

tdoaAccuracy输入rangeVarianceM2[K×M]，bsMask为逻辑[K×M]（省略为全true）；spec取comparison.tdoa。输出positionStdM/positionVarianceM2/valid/rank/referenceBs[M×1]、failureReason[M×1 cell]、usableBsMask[K×M]及每用户cell的covarianceM2[2×2]、rangeDifferenceCovarianceM2[(J-1)×(J-1)]、differenceMatrix[(J-1)×J]、H[(J-1)×2]、bsIds[J×1]。正方差独立TOA、BS同步、oracle位置冻结协方差，公共接收机时钟经差分消除。二维至少3站。固定参考缺失不会换站，失效结果NaN；minimum_variance自动在可用站中选最小TOA方差（不是额外接收功率检测器）。

interferenceKernels的spec取comparison。返回communicationPsdWPerHz[K×M×Na]（当前实际Pc下，其他BS C对目标PRS）、prsPsdPerWAtCommunication[K×N×K]（每W广播PRS对C接收机）、pilotFrequencyHz/basebandFrequencyHz/communicationFrequencyHz、模型及干扰标志。使用gCP[M×K×N]，gPC[K×N×K×M]取相同频率平坦通道的最后维第一项；spec.chunkSize控制导频分块。广播导频已按sum(weightsSquared)=Nc归一化。不是全有色协方差似然，未计C-to-C。

delayCrbInterference输入gP[K×M]、可行budget（K×1符号数/瞬时功率）、正N0及kernel，输出delayVarianceS2/rangeVarianceM2/rangeStdM/fisherInformationPerS2[K×M]。已知相位、对角高斯导频噪声近似；零干扰还原原delayCrb，保留model/fullColoredNoiseLikelihood=false。无同时其他PRS，未伪加其干扰。

serviceAudit的policy仅支持global_muting和other_bs_communication，berModel必须显式qpsk_bit_energy；referenceBer[K×N]来自同一次参考算法。输出communicationOpportunityFraction/activeCommunicationPowerW/averageCommunicationPowerW/requiredPowerForEqualCommunicationEnergyW[K×1]，energyCompensationPeakFeasible/operatingPowerFeasible/communicationEnergyMatched/attemptedRateContractMet/activeSlotBerQosMet/serviceContractMet/energyAndServiceMatched[K×1 logical]。活动C功率固定num.N*num.PcW；补偿功率仅作审计、不实施、不截断。

service还含meanActiveBer/worstActiveBer/referenceBer/expectedCorrectBitRatePerCarrierBps/referenceExpectedCorrectBitRatePerCarrierBps[K×N]、berDuringPrsSlot[K×N×K]（静默站为NaN）、berNoPrs[K×N]、attemptedBitRatePerCarrierBps[K×1]及failureReason[K×1 cell]。尝试比特率为2/Tc×机会，预期正确比特率按活动槽加权，不是packet goodput；无活动机会的BER为NaN。

run.prsComparison为独立后评价；禁用/参考缺失返回disabled/reference_solution_unavailable且无policies；PRS计划不可行返回prs_observation_schedule_infeasible与空policies；成功为computed_partial_match_service_and_observation_audit，仍不保证服务合同可行。含kernel、msTdoaCommon、msTdoaHearable、commonBsMask、policies(1:2)。每policy有id/crb/tdoa/service。定位广播能量和oracle几何匹配标志为true；observationMatched/fullyMatched/detectionCoverageEvaluated保持false。MS TDoA来自独立误差的稳态DLL方差代理，不是有限观测CRB。

阶段5的cfg/run版本5.0、6策略/9图。main调用comparison，saveRun调用新CSV导出并记录manifest.prsComparisonStatus/prsMatching，图9由plotPrsComparison产生。完整矩阵和状态保存在run.mat。阶段2/3/4历史脚本显式关闭comparison保持旧图数；run_stage5_checks执行63项单元检查、三个完整main及源码哈希核对。不得把新几何有效人数当检测覆盖率，或以观察期相同声称估计器等价。

kernel.communicationReceiverModel明确C端为载波中心PSD取样，不是匹配滤波带积分；gridCaveat记录对齐网格的矩形符号PSD零点。service.receiverCaveat和manifest.prsMatching保留同一限制。49.77MHz默认BER=0可能来自此零点和erfc下溢，不能推断实际接收机无干扰；图9的绘图下限不改CSV。

## DLL独立方差分解与解析滤波（阶段6历史接口）

```matlab
unitPsd = ranging.communicationPsd(localFrequencyHz,userId,ch,num,chunkSize)
bundle = ranging.biasKernels(ch,num,dll,opts,existingKernel)
finiteN = ranging.asymptoticVariance(powerW,ch,num,N0,dll)
result = ranging.biasDecomposition(powerW,ch,num,dll,N0,bundle,factor,geom)
study = experiment.dllBias(run)
overlap = receiver.rectangularPsdOverlap(deltaFrequencyHz,transmitSymbolS,receiveSymbolS)
audit = prs.receiverAudit(ch,num,allocation,budget,N0,berModel,berThreshold,spec,pointKernel,referenceBer)
report.exportDllBias(run,runDir)
report.exportReceiverAudit(run,runDir)
figures = report.plotDllBias(run,cfg.io)
figures = report.plotReceiverAudit(run,cfg.io)
```

cfg/run版本6.0、默认6策略/11图。新增models.dllBias(enabled,userIds,referenceAlgorithm)，默认true/自动[1,ceil(M/2),M]/central。comparison.receiverAuditEnabled默认true。旧配置从base重建。积分数值设置沿用rangingIntegral；dllBias独立开关不随普通积分附图关闭。

communicationPsd输入任意实有限形状的局部频率数组，userId为单个P用户；输出同形状PSD（秒，即每单位每C载波功率的1/Hz），已汇总全部BS/有限N载波，未乘Pc。chunkSize默认256。gCP为M×K×N，不使用C接收机张量。

integralKernels增加可选opts.communicationModel=finite_carriers（缺省）或spectral_lines。后者取同一有限载波列表的单位面积零谱宽极限，边界线谱半权重；不是无限平坦PSD。kernel增加communicationModel与quadratureOptions；旧finite模型的归一化和主优化不变。

biasKernels输入opts沿用rangingIntegral，必须显式给非空userIds。existingKernel可省略，全部quadratureOptions严格相同才复用。返回variants(1:8)，每项id、flags[linearizedSine,spectralLines,localFrontend]、kernel、reused；另有userIds、factorOrder、baselineIndex=1、allConverged。八组合是共同功率下的设置，不是八次功率优化。

asymptoticVariance输出全K×M的噪声/通信/定位/总时延方差及距离方差，保留附录B的有限N精确加权和。nominalBCommunicationMultiplier=(N+1)/N，与旧varianceThreeTerm通信项的差异独立列出。它是解析渐近模型，没有数值积分收敛含义。

biasDecomposition对给定正powerW、相同ch、N0、dll和oracle geom评价八组合。selected=L时delay/noise/communication/positioningVarianceS2、rangeVarianceM2、relativeToBaseline为K×L×8；positionStdM为L×8。pairedContrasts(1:12)记录factor、fromIndex、toIndex、absoluteDifferenceS2/relativeToFrom/logVarianceRatio[K×L]。mainEffectRelativeToBaseline为K×L×3，定义为基线单独翻转某因素。

chainVarianceS2[K×L×7]及六个chainStepIds形成有序闭合路径；chainDifferenceS2、chainRelativeToBaseline、chainRelativeToPrevious[K×L×6]保留不同分母，不混用可加性。totalRelativeBias/closureResidualS2[K×L]与maximumRelativeClosureResidual核对最终原式(10)。前三因素之外，另外保存其余积分渐近、N→N+1和去定位干扰三步。偏差是预测方差差异，不是采样时延均值偏差；前端是解释敏感性；归因顺序不唯一。

experiment.dllBias在主优化和PRS后评价后运行。run.dllBias含status、referenceAlgorithm、bundle、result及optimizerChanged=false/sampledReceiverImplemented=false。禁用status=disabled；指定参考缺失/评价失败为reference_solution_unavailable；非正功率为reference_power_outside_positive_domain。不更换参考算法、不重求功率、不修改oracle权重。

rectangularPsdOverlap返回与频差数组同形状的非负单位PSD积分（秒），Ttx/Trx为正标量。零频差采用稳定矩级数，一般频差用三角相关乘积的闭式；对应平稳独立符号与平均相对时序。此函数无样本输入输出，也不执行相关峰、测距或定位估计。

prs.receiverAudit的spec取comparison，使用allocation、budget、pointKernel和referenceBer[K×N]；返回pointUnitPsdWPerHzPerW/filteredUnitPsdWPerHzPerW[K×N×K]、pointService/filteredService。两者复用相同活动时长、实际C/PRS功率和QPSK模型，仅C接收端干扰核不同。status=computed_analytical_average_timing_filter_audit，sampledReceiverImplemented/finiteObservationMatched/prsRangingRecomputed均false；assumptions标注独立平稳符号/平均时序，不代表相干对齐OFDM或已知PRS序列的真实接收器。没有改旧PRS端测距CRB。

新结果位于run.dllBias、run.prsComparison.receiverAudit，默认追加图10/11和8个CSV。manifest记录dllBiasStatus/dllBias及receiverAudit。历史stage1–5脚本显式关闭阶段6附加模块以保留原图数。run_stage6_checks执行80项检查、两带宽加绑定预算的三个main和87源文件哈希核对。完整数学边界见DLL_BIAS_AND_RECEIVER_AUDIT.md。

## 有限采样与实际估计反馈（阶段7接口，阶段8兼容扩展）

阶段7 cfg/run schemaVersion为7.0，当前阶段8为8.0。新增实验默认关闭，不改变原六策略/11图的默认路径；采样与反馈各追加一图，步长试验再追加一图。旧阶段文档的schema和状态属于对应历史快照。

```matlab
opts = config.sampledDefaults(cfg)
packet = waveform.makePacket(num,opts)
signal = waveform.component(packet,k,m)
[spectrum,signal] = waveform.reference(packet,k,m)
tx = waveform.transmit(packet,powerW,includeCommunication)
receivedSpectrum = waveform.propagate(tx,packet,gain,delayS,phaseRad,noisePsd,seed)
psd = waveform.communicationPsd(frequencyHz,num)
estimate = receiver.estimateDelay(receivedSpectrum,templateSpectrum,packet,opts,noisePsd)
fix = geometry.estimateRanges(bsPosM,rangesM,valid,opts)
study = experiment.sampledReceiver(scene,ch,num,powerW,opts,noisePsd,packet)
feedback = experiment.positionFeedback(context,cfg,feasibility,packet)
study = experiment.stepSensitivity(problem,cfg,feasibility)
[ranges,positions] = report.sampledTables(study)
report.exportSampled(run,folder)
figures = report.plotSampled(run,cfg.io)
```

以上均为msnoma包下接口。powerW[K×M]、gain/delayS/phaseRad[K×1]；packet包含chips[K×M×L]、communicationSymbols[K×N×Nc]、communication[K×activeSamples]，其中L将请求码片数向上补齐到G的整数倍。sampleRateHz=oversampling×B，samplesPerChip=(M+1)×oversampling；FFT记录至少覆盖两倍活动包及延迟保护。

component返回单个BS/P用户单位幅度的活动包[1×activeSamples]。reference另返回已投影到前端的FFT[1×fftLength]。transmit按sqrt(powerW)叠加全部定位分量及可选通信，返回spectrum[K×fftLength]和名义/有限观测功率审计。propagate是唯一读取实际时延的物理层，返回接收FFT[1×fftLength]，不返回真值标签。communicationPsd输入输出形状相同，返回有限C载波总功率PSD[W/Hz]，包含Pc、尚未乘信道增益。

estimateDelay输入两个FFT[1×fftLength]及波形元数据，opts还须含speedMps。无位置/信道/真值时延参数。输出delayS/rangeM、detected/valid、status、coarseLagSamples/fractionalLagSamples、complexAmplitude、normalizedPeak/threshold/atSearchBoundary及局部noiseOnlyRangeCrbM2。检测先按coarseOversampling插值选峰，再连续分数时延细化；not_detected不返回真值，search_boundary_peak无效。门限为离散搜索AWGN并集界，高斯等效预白化下使用假定PSD，interferenceFalseAlarmCalibrated=false。

geometry.estimateRanges针对**单个用户**，bsPosM[K×2]，rangesM/valid[K×1]，opts取sampledReceiver.position。返回positionM[1×2]、valid/status、usedBs、iterations、residualM、residualRmsM和conditionNumber。至少3站且非共线；缺站、锚点/解秩亏、求解失败、残差超限、域外分别标记，域外不裁剪。没有真实用户坐标输入。

experiment.sampledReceiver的opts取cfg.models.sampledReceiver，packet是可省略的第7参数；提供时直接复用，不重新生成。opts.userIds为空时取全部用户，Luser=numel(userIds)。输出rangeM/delayS/detected/validRanges[K×Luser]，receivers[K×Luser cell]、fixes[Luser×1 cell]、estimatedPositionM[Luser×2]、validPosition/positionErrorM[Luser×1]以及rangeErrorM[K×Luser]。真实坐标/时延只在传播和事后误差使用，无效误差NaN。packet元数据保留实际chips/通信符号/种子但去除大FFT缓存；txAudit去掉spectrum。sampledReceiverImplemented=true、steadyStateDllImplemented=false、fullyMatched=false。

可选whitening为none、known_average_C_psd或默认known_average_CP_psd。后者使用理想CSI/已分配功率的平均C及其他P功率谱，排除当前目标；不是含码交叉协方差的完整似然。预白化后noiseOnlyRangeCrbM2置NaN，另存gaussianEquivalentRangeCrbM2，不把干扰等效界误标为纯噪声界。

positionFeedback的context由runSingle构造，包含scene/ch/num/A/qos/rangeFactor；feasibility可复用，因为反馈只改目标几何，不改变约束。首次lambda[K×M]全部1；每轮求解并认证，再实际接收/测距/定位；全部用户和全部K×M测距有效才重算估计几何。支持central/dual/cp4a_subgradient_candidate；history每项保存iteration、geometrySource、lambda、solution、acquisition、positionChangeM、relativePowerChange。无oracle回填，无删用户重定义目标。失败返回power_*、acquisition_or_position_failure或estimated_geometry_rank_failure；收敛需最少fix数后位置/功率变化同时达标；迭代上限为fix_limit_reached。固定CSI假设和commonRandomNumbers单独记录。原六策略的oracle对象继续保留。

stepSensitivity由cfg.experiment.stepSensitivity控制，multipliers同时缩放内外初始步长（递减模式还缩放diminishingInitialStep）。每trial保存multiplier、innerInitialStep、outerInitialStep、完整solution；返回verifiedCount及objectiveRelativeSpread。authorStepsKnown=false；仅以真实证书认定收敛，失败试验也保留。CSV中的lastAcceptedInnerStep向后检索整个求解的实际正接受步长，避免最终warm start零步更新被错记为步长0；未发生接受则NaN。

新输出sampled_ranges/positions/packet_audit、feedback_iterations/geometry_powers/ranges/positions、step_sensitivity共8类CSV；禁用或没有结果时不创建对应输出。run.sampledReceiver、run.positionFeedback、run.stepSensitivity和manifest保留实际状态及边界。PRS原解析审计内的sampledReceiverImplemented=false仍描述该支路自身，不被MS有限采样支路覆盖。run_sampled_tests含17项检查；阶段7整体验证入口run_stage7_checks，结果见STAGE7_VERIFICATION.json。具体假设、失败与当前数值见SAMPLED_RECEIVER_AND_FEEDBACK.md和PROGRESS阶段7。

## 估计P干扰消除及仿真消融（阶段8基础接口）

```matlab
bank = waveform.positioningBank(packet,powerW)
fit = receiver.estimatePChannels(receivedSpectrum,bank,packet,opts,noisePsd)
prediction = receiver.reconstructPositioning(bank,fit,packet)
audit = experiment.interferenceAblation(ch,num,powerW,opts,noisePsd,userIds)
report.exportPChannels(study,folder,prefix)
figures = report.plotCancellation(run,cfg.io)
```

bank[K×fftLength]是每BS全部P已知码及sqrt(powerW)相加、通过既有前端的FFT；不读取通信数据。estimatePChannels输入receivedSpectrum[1×fftLength]、bank、公开采样/频率元数据、opts取sampledReceiver且另含speedMps、假定高斯PSD标量。**没有真实时延/增益/坐标接口**，但该接收模式额外假设全部P码及相对功率可知。

fit.delayS/complexGain[K×1]来自逐次捕获初值及联合最小二乘/投影Gauss-Newton；coarseOrder记录BS顺序；history列为[iteration,weightedResidualEnergy,proposedMaximumStepSamples,previousAcceptedStep]；iterations/rankRatio/relativeResidualEnergy保留诊断。valid表示局部优化达到容差且非边界，不是正确峰保证。status支持signal_absent、coarse_acquisition_failed、zero_template_energy、template_rank_deficient、delay_information_rank_deficient、line_search_failed、iteration_limit、search_boundary_peak和converged。无效fit禁止reconstructPositioning，失败时不消除。

拟合有效后另保存gaussianDelayStdSamples/gaussianRelativeGainStd[K×1]、cancellationEligible和eligibilityStatus。延迟代理用投影信息逆矩阵，复增益代理是**固定时延条件下**LS方差，未传播时延估计误差；门控默认全部站std≤0.25采样且相对增益std≤0.25。它们依赖假定PSD，不用真实误差判定，不能保证低预算/强C时改善。所有配置均在sampledReceiver.interferenceCancellation，默认enabled=false。

reconstructPositioning返回sum_k fittedGain_k×bank_k×exp(-j2πf fittedDelay_k)，只依据估计值。sampledReceiver在有效且通过精度门控后减去全部估计P、加回目标P，再以C+噪声白化重测目标；由此同时改变样本、权重及检测分母，实际whitening写入每链路对象。没有传播消除估计误差到虚警/CRB，estimatedCancellationUncertaintyPropagated=false；noiseOnlyRangeCrbM2置NaN，gaussianEquivalentRangeCrbM2只为忽略消除不确定性的代理。门控或拟合失败则用原CP谱/原样本接收，非oracle填补。

run.sampledReceiver新增estimatedPChannels[Luser×1 cell]、allPositioningCodesKnown、communicationDataUsedByEstimator=false和estimatedCancellationUncertaintyPropagated=false；receivers新增pCancellationApplied/pChannelFitStatus。反馈迭代复用同一逻辑，全部有效要求保持不变。run.schemaVersion当前8.1，旧cfg从base重建；消除默认关闭，原11/13图路径保留。

interferenceAblation为**受控发射仿真诊断**，userIds非空；内部把已分配功率按desired_only、same_carrier_P、same_BS_P、all_P、all_P_C、all_P_C_AWGN六种掩码生成接收。返回audit.rows，行数6×K×用户数，字段包含positioningUser/baseStation、transmitterMask、通信/噪声是否存在、valid、estimatedRangeM、rangeErrorForEvaluationM、normalizedPeak、threshold、sameCarrierNearFarRatio及status。真实距离只在传播和事后误差用；不能将消融的“只发目标”当业务可选接收器。PSD核随发送掩码变化，故非固定权重的单因素误差分解。runSingle/main通过cfg.experiment.interferenceAblation(enabled=false,userIds=18)调用，缺参考返回reference_solution_unavailable。

exportPChannels保存sampled_estimated_P_channels.csv或feedback_XX_estimated_P_channels.csv；拟合即使未通过精度门控也保留。sampled_ranges/feedback_ranges增加消除状态，另存interference_ablation.csv。plotCancellation由plotSampled调用，消除和消融各加1图；saveRun记录接收器附加知识、消除开关和未传播不确定性。run_cancellation_tests含15项检查，完整验收run_stage8_checks：112项、三个同随机数配对接收对照、原全区域无C无AWGN控制，共七次完整main，源码哈希校验。结果和限制见STAGE8_VERIFICATION.json、ESTIMATED_P_CANCELLATION.md与PROGRESS阶段8。

## 阶段8.1 残差主导门控

新增 `gate = msnoma.receiver.cancellationGate(fit,opts)`，opts取sampledReceiver.interferenceCancellation。gate包含eligible/status、precisionPassed、residualPassed、maximumRelativeResidualEnergy。接口没有真值、误差、场景或通信数据。fit新增precisionGatePassed/residualGatePassed/maximumRelativeResidualEnergy并保留原估计及局部valid。仅全部门控通过才消除；拒绝时原样执行原接收器。

maximumRelativeResidualEnergy默认0.5，是加权残差/接收能量上限；设1回放阶段8仅精度门控分支。它是保守工程工作区而非统计检验或精度保证。CSV和manifest新增门限/门控字段；图增加残差上限线。消除仍默认关闭，其他接口不变。验收入口run_stage8_1_checks（原112+新5项、六个main及同场景0.25/0.5/0.75门限），详见STAGE8_1_RESIDUAL_GATE.md；历史阶段8文件不覆盖。