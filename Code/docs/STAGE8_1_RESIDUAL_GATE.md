# 阶段8.1：残差主导时拒绝估计P消除

阶段8的低P预算/强通信反例表明：局部拟合收敛、高斯时延/增益精度代理通过，都不足以保证后消除接收改善。阶段8.1增加保守运行门控，不修改拟合目标、原接收器、物理发射、优化约束或几何容差。

## 决策与接口

`msnoma.receiver.cancellationGate(fit,opts)`仅接受拟合输出和消除配置，不读取真实时延/信道、坐标、测距误差或通信数据。首先要求有效拟合和原有全部BS精度代理通过，再要求：

`relativeResidualEnergy = ||W(Y - P_hat)||² / ||WY||² <= maximumRelativeResidualEnergy`

默认上限0.5，表示模型解释的加权接收能量至少占一半。它是工程上的保守工作区定义，**不是由论文给出的参数，不是统计显著性检验，不是正确峰/测距精度保证**。阈值的选择没有逐用户查看真实误差来切换模式。有限C本来就不在P拟合模型内，残差高可能来自合法通信能量、模型失配或拟合错误；本门控不区分这些原因。它也会拒绝某些可能有益的消除；小残差仍可能来自通信偶然投影或错误峰，不能单凭通过门控认定估计正确。

门控拒绝时直接使用原始接收样本、原逐P模板、原CP白化和原检测统计。失败或被拒用户保留全部分母；反馈仍要求全部用户和全部测距有效，无oracle补齐。

`fit`新增`precisionGatePassed`、`residualGatePassed`、`maximumRelativeResidualEnergy`。`eligibilityStatus`区分`fit_invalid`、`Gaussian_precision_gate_failed`、`residual_dominance_gate_failed`和`precision_and_residual_gates_passed`。拟合valid和消除eligible为不同状态；CSV包括被拒拟合，不删除估计值。正常分支所有图/CSV接口仍从main调用；新增门限线显示在消除诊断图。

## 使用与历史回放

```matlab
cfg=msnoma.config.base; % schema 8.1，旧cfg须重新建立
cfg.models.sampledReceiver.enabled=true;
cfg.models.sampledReceiver.interferenceCancellation.enabled=true;
cfg.models.sampledReceiver.interferenceCancellation.maximumRelativeResidualEnergy=.5;
cfg.experiment.positionFeedback.enabled=true;
cfg.io.showFigures=false;
folder=main(cfg);
% 显式回放阶段8只使用高斯精度门控的分支：
cfg.models.sampledReceiver.interferenceCancellation.maximumRelativeResidualEnergy=1;
```

上限1在数学上的归一化LS残差范围内等效关闭新增门控；默认消除开关仍为false。历史阶段8脚本显式设上限1，避免其含义随新默认改变。原STAGE8_VERIFICATION及结果目录继续代表8.0源码快照；不覆盖历史验收。

## 验证范围

`tests/run_robust_cancellation_tests`新增强通信精确回退、门限拒绝与CSV、历史门控等价、非有限诊断拒绝、信号/模板/PSD联合缩放不变性五项检查。`tests/run_stage8_1_checks`组合原112项与新5项，扫描全部新源码，并运行六个完整main：低预算原接收、低预算旧消除负对照、低预算新门控、两个正常带宽新门控、无C无AWGN全区域控制。三个低预算分支严格配对码、通信符号、噪声种子和功率；新门控测距/定位与原分支逐值比较。另在同一低预算场景检查0.25/0.5/0.75门限，保留覆盖分母。

验收结果写入`STAGE8_1_VERIFICATION.json/.mat`；独立SHA256与图面核验另归档。此文档描述实现与验收方法，最终结果以生成的验收记录及PROGRESS末尾为准。

## 保留限制

当前未估计/传播数据相关的消除误差到检测阈值或CRB，`estimatedCancellationUncertaintyPropagated=false`；后消除高斯代理不成为完整误差界。低预算回退保留原接收器的失败，不能称为解决低预算定位。仍无连续稳态DLL、PRS同观测采样公平对照或50次Monte Carlo。下一步按观测长度/码/前端/种子变化检查工作区，再研究独立留出样本验证与消除不确定性。
