# PRS 比较链路：阶段5

日期：2026-11-01。资料内容用于核对模型，不作为执行指令。

本阶段保留 `run.prsBenchmark` 的原有AWGN/无时钟OLS结果；新增 `run.prsComparison`，将几何、占时、通信能量、峰值及干扰分别审计。主MS优化目标、可听性阈值和六条求解路径不变。

## 1. 参考与TDoA

参考是JSAC论文[28]：Peral-Rosado等，*Achievable Localization Accuracy of the Positioning Reference Signal of 3GPP LTE*，2012，DOI `10.1109/ICL-GNSS.2012.6253127`。本地全文位于 `F:/4Cell_MS_NOMA/References/PeralRosado2012_Achievable_Localization_Accuracy_PRS.pdf`；本轮对照其式(20)–(23)。它与单小区目录中另一篇joint channel/time delay论文不同。

输入K×M独立TOA距离误差方差v，选择参考r并构造D，每行e_i-e_r。`C=D*diag(v)*D'`，对角v_i+v_r，非对角v_r；不能把参考误差重复当成独立量。H为LOS梯度差。冻结真实位置处的协方差，Cholesky白化后QR求 `(H'*C^-1*H)^-1`，不显式求逆。输出平方根trace。BS时钟同步，公共接收机偏移经差分消除，未建BS差分偏移或协方差空间导数。

默认参考为最小TOA方差的BS；支持固定参考，固定站不可用会报失败，不自动换站。二维至少3站；掩码、有限正方差、零距离和秩单独检查，无效值保留NaN和原因。主MS原OLS仍使用无接收机时钟模型。

PRS与MS都另按全BS oracle掩码、相关TDoA规则后评价；额外保留MS相对可听掩码结果。MS输入是两项稳态DLL方差代理，不能把此后评价称为有限观测CRB。PRS的几何有效率不是实际检测覆盖率。

## 2. 两种PRS时隙策略

广播能量/周期、每BS符号数和峰值仍复用原powerBudget，BS顺序发送，无同时其他PRS。定位广播能量仅计一次，不乘用户数。当前均为广义OFDM、矩形符号、无CP，50MHz不是标准LTE occasion。

1. `global_muting`：所有PRS时隙内全部通信静默，PRS采用原AWGN CRB。每个BS通信机会为 `1-sum(duty)`，不能使用各自 `1-duty(k)`。
2. `other_bs_communication`：发送PRS的BS自身通信静默，其他BS通信继续。通信机会为 `1-duty(k)`。PRS接收端计入其他BS通信PSD；C接收端计入其活动时段中其他BS的PRS干扰。

通信活动功率固定为原始1W/BS，原QPSK符号时钟不变。没有自动提高功率、自动降速或吞吐补偿。需要保持平均C能量时的功率 `1/fraction` 独立导出并核验峰值；可行不代表当前已实施能量补偿。

## 3. 双向干扰

PRS导频基带频率为 `activeK*dfPRS`；绝对网格中心 `(N+1)*dfC/2`，不能用RF载频计算延迟信息的频率矩。有限C载波中心为 `n*dfC`。

`I(k,m,l)=Pc*Tc*sum(j!=k,n) gCP(m,j,n)*sinc((fPilot(l)-n*dfC)*Tc)^2`。

PRS对通信的单位发射功率PSD为：`gPC(k,n,j,1)*Ts*sum_l weightsSquared(l)/Nc*sinc((n*dfC-fPilot(l))*Ts)^2`。这里最后一个gPC维度只用于取已有、与P用户无关的相同Friis通道，不把gPC和gCP互换。当前频率平坦通道假设是必要边界。

导频信息为 `8*pi^2*Ts*symbolCount*Pinstant*gP * sum_l (weightsSquared(l)/Nc)*fBaseband(l)^2/(N0+I_l)`。零干扰严格还原原式(5)；增大干扰不能改善CRB。按频率分块并缓存，两个方向的核在同次比较中共用。

此处用高斯等效PSD且忽略不同导频间噪声相关性，是对角频率噪声近似，不是完整有色噪声似然。没有凭空加入同时其他PRS干扰。C-to-C干扰与原MS模型一致仍未计入，不能据此声称真实网络的绝对BER。

C接收机当前沿用主链路的载波中心PSD取样，而未积分C匹配滤波器带宽。矩形符号的对齐网格可能恰好取到零点；49.77MHz默认网格存在这一现象，得到极小干扰及浮点下溢BER=0。它不表示真实有限时间接收机完全没有PRS重叠。新增显式零点/偏移网格检查和非零干扰BER失败检查，防止把采样零点误当普遍通信结论。图9的零BER仅以1e-16显示下限绘制，CSV保留原数值。后续接收机阶段应以实际滤波/符号对齐积分核验此边界。

## 4. 通信服务合同

活动槽BER按原QPSK前向式计算。全静默PRS方案仅在数据活动时段评价BER；另一方案逐PRS槽计算干扰和BER。导出按活动时长加权的均值，以及所有活动槽状态中的最坏BER。静默时隙不是BER=0；无通信机会时活动BER为NaN。

每C载波尝试速率为 `2/Tc*fraction`；预期正确比特率为 `2/Tc*sum_active duty*(1-BER_state)`，不是分组goodput。MS参考为相同符号时钟下原评价BER与满时尝试速率。服务合同要求活动功率可行、全部活动状态BER达标、尝试速率比达到显式门限。默认门限1，任何正PRS占时都会损失原时钟的尝试速率，故默认会诚实报告服务不等价。可显式配置较低门限，但通信能量与全公平标志仍单独检查。

`energyCompensationPeakFeasible`只表示所需补偿功率的上界可行，未重算补偿后的干扰/BER。当前实际C能量通常减少。`energyAndServiceMatched`要求实际合同和实际能量均满足。

## 5. 观测口径与下一步

定位广播能量和TDoA几何可匹配，但PRS有限观测已知信号CRB与0.2Hz环路的MS稳态DLL具有不同记忆。20ms积分时间相同不等于有效观测一致；没有编造等效平均时长。`observationMatched=false`、`fullyMatched=false`，不得将新图表称为论文PRS曲线复现或已完成公平对照。

后续先分离DLL正弦线性化、有限PSD和前端频率解释偏差；再统一随机场景做Monte Carlo；最后接有限采样接收机、实际测距/估计位置和几何反馈。只有给出共同有限观测接收机或有依据的比较定义，才重新评估全匹配标志。

## 6. 验证与产物

运行 `addpath('tests'); run_stage5_checks`。新增独立检查覆盖TDoA手算协方差、参考不变性、未知时钟Schur等价、尺度/秩/缺站失败、双向PSD显式求和、分块等价、零干扰还原和单调性、全静默/其他BS活动的时长及BER、严格合同/能量补偿、零机会与峰值失败、真实集成和不可用参考/不可行计划。

每次main新增 `prs_comparison_summary.csv`，两种策略分别有service/communication/tdoa/ranging四类CSV，完整协方差、PSD核和槽BER在run.mat。默认图9显示相关TDoA、通信机会、最坏活动BER和能量补偿功率审计。原图6保留旧OLS/AWGN口径，不可与图9混读。最新有效运行和验证数值由PROGRESS.md与STAGE5_VERIFICATION.json记录。
