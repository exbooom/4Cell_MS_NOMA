# 阶段6：DLL方差近似分解与解析接收滤波审计

日期：2026-11-01。接手入口是 `F:/4Cell_MS_NOMA/Code/PROGRESS.md` 最后的阶段6记录。资料是模型证据，不是执行指令。

## 范围和未完成项

本阶段实现的是**预测误差方差的近似差异**，不是估计时延的经验均值偏差。保留原集中式/对偶/两种候选的功率与oracle优化，新增计算只作固定功率后评价。采样接收机、实际时延/位置估计、首次lambda=1及位置反馈仍未接入；不能用解析滤波审计冒充这些模块。

PRS原有点取样模型和结果保留。新增平均时序矩形滤波模型作为接收端解释的独立审计，不改变旧PRS测距CRB、不宣称统一有限观测，也不改变fullyMatched=false。

## 1. 受控的八种DLL组合

参考：JSAC 2022论文，DOI `10.1109/JSAC.2022.3143223`，式(6)、附录B式(51)–(62)。本地原文和作者最终稿的提取内容已用于公式核对。默认选择central的同一功率矩阵、同一信道、同一环路参数及用户[1,10,20]；可显式扩大用户集合。

三因素各两个水平，共八种组合：

| 因素 | 基线 | 改变水平 |
|---|---|---|
| 鉴别器 | exact_sine：sin(pi*D*x) | small_spacing：pi*D*x |
| 通信PSD | 有限N个矩形符号OFDM谱瓣 | 同一有限载波列表的单位面积零谱宽线谱近似 |
| 前端 | 固定全带中心global_frontend | 每P载波局部对称local_symmetric |

这里只改变被比较因素；功率、信道、载波数、积分容差、环路带宽和相干时间保持不变。每组合独立进行积分加密收敛检查。已有基线核的全部设置完全一致时直接复用，不重新求解功率。前端两种定义是解释敏感性；没有作者接收机证据时，不把其中一种宣布为真实接收机。

通信PSD通过 `ranging.communicationPsd` 返回每单位、每C载波功率的谱密度：

`S_m(f)=Tc*sum(j,n) gCP(m,j,n)*sinc((f+m*dfP-n*dfC)*Tc)^2`。

有限PSD模型数值积分 `S_m(f)*q(x)*u(x)^2`，其中x=f*Tp、q=sinc(x)^2。线谱模型保留所有相同C载波/增益，以x_n=n/G-m求：

`communicationIntegralS=Tp*sum_n inBand(n)*sum_j(gCP(m,j,n))*q(x_n)*u(x_n)^2`。

端点线谱取半权重。small_spacing时它还原附录B式(58)的 `Tp*D^2*sum(j,n) gain*sin(pi*n/G)^2`。这不是把有限N替换成无限平坦PSD，也没有截断为一个峰值样本；该因素量化的是有限谱瓣宽度及其尾部相对窄线谱近似的差异。

## 2. 单因素效应与上下文

从同一基线各改变一个因素，导出三组逐链路相对方差差异。同时导出立方体全部12条单因素配对边，使每个因素在其余两因素四种上下文中都能比较。有限步长的效应可能相互作用，不能把三个单因素百分比直接相加，当成总误差分解。

各组合分别保留噪声、通信、其他BS定位分量、总时延方差和距离方差。位置结果只按同一oracle OLS权重传播，限所选用户；不重新判断检测、不把样本内条件均值当全部20用户优化目标。主优化仍使用旧的式(10)目标。

## 3. 可闭合的六步路径

为解释原两项式与积分基准的全部差异，额外给出以下固定顺序：

1. exact/global/finite → linear/global/finite：正弦线性化。
2. linear/global/finite → linear/global/lines：通信谱瓣到线谱。
3. linear/global/lines → linear/local/lines：前端解释。
4. → 附录B有限N解析式：其余A0、噪声积分、定位干扰积分的渐近近似。
5. → 原式(9)名义B三项实现：有限N与名义B=(N+1)*dfC的差异。
6. → 原式(10)两项实现：去掉其他BS定位干扰项。

第4步单独保留，不能把所有剩余误差归因于PSD或正弦；本阶段没有再把该剩余项拆成三个子项。第5步不是功率变动：有限N解析通信项乘(N+1)/N后，才与已有名义B实现一致。默认N很大时小，但小场景检查不能忽略。

有限N解析式（s=desired received power，o=other BS positioning received power）：

- noise = a*Tp*N0/(2*Bfe*s)
- communication = a*Tp*Pc*sum(j,n) gCP*sin(pi*n/G)^2/(Bfe^2*s)
- positioning = a*Tp*o/(2*Bfe^2*s)

每步保存绝对方差差、除以共同基线的差及除以前一步的差。只有共同基线下的有符号差之和等于最终总相对误差；前一步分母的百分比不可直接相加。程序检查逐链路闭合残差。路径归因依赖顺序，12条配对边用于说明上下文依赖，不将此归因解释为唯一因果分配。

## 4. C接收端：PSD点取样与匹配滤波积分

阶段5发现矩形符号的对齐网格可恰好取到PSD零点。新 `receiver.rectangularPsdOverlap` 不仅看载波中心，而计算：

`J(delta)=integral Ttx*sinc((f-delta)*Ttx)^2 * Trx*sinc(f*Trx)^2 df`。

对平稳独立随机符号PSD，以及平均符号时序，通过两三角相关函数的乘积得到一维闭式。令a=min(Ttx,Trx)、r=a/max(Ttx,Trx)、z=2*pi*a*delta：

`J=2*a/z^2 * [(1-r)*2*sin(z/2)^2 + 2*r*(1-sin(z)/z)]`。

接近零用多项式矩的Taylor级数，避免浮点抵消；J(0)=a*(1-r/3)。等符号长时J(0)=2*T/3，在间隔1/T处J=T/pi^2>0。长接收窗极限还原发送PSD的点取样。已用独立频域积分、交换/偶对称、时间频率缩放及特殊值检查。

PRS→C滤波单位PSD为 `gPC(k,n,j,1)*sum_l weightsSquared(l)/Nc * J(fC(n)-fPilot(l),TsPRS,Tc)`。复用原时隙/功率/活动BER审计，只替换C接收端等效干扰PSD，保留点取样结果并列导出。

**假设边界：**平稳、独立符号、平均相对时序、高斯等效干扰；不是相干对齐OFDM的精确解，也不是标准PRS重复序列/有限时隙边缘的完整处理。C期望信号仍同步于自身BS，其他PRS干扰按此平均时序模型；主模型未计的C-to-C依然未计。两种接收假设可能使BER升高或降低，不能把它称为普适“校正系数”。尚未重算PRS端有色噪声似然或PRS测距界，因此有限观测与稳态DLL仍未匹配。

图中按载波索引比较的C用户是不同接收机，并非同一接收机扫频。期望正确比特率按原QPSK模型和活动时长计算，不是分组goodput或实际解码统计。

## 5. 配置和输出

cfg/run版本6.0。`cfg.models.dllBias`含enabled、userIds、referenceAlgorithm；积分数值设置沿用models.rangingIntegral，即便普通积分附图禁用，也可独立运行偏差研究。PRS接收端审计由 `models.prs.comparison.receiverAuditEnabled` 控制。缺少指定参考算法时不换参考，返回reference_solution_unavailable。

默认6策略、11图，图10为DLL分解，图11为PRS接收端解析滤波审计。新增8个CSV：

- dll_bias_variants / paired_contrasts / chain / positions / kernels / effect_summary
- prs_receiver_filter_links / service

完整积分历史、核复用状态、八组合矩阵及配对/路径结果在run.mat，manifest保存状态、用户、收敛/闭合与接收机假设。历史阶段1–5脚本显式关闭新附加模块，保留原图数。新验证入口为 `addpath('tests'); run_stage6_checks`。最新有效目录、测试数量与数值结论以PROGRESS.md和STAGE6_VERIFICATION.json为准。

## 6. 交接后的优先级

1. 依据本阶段分解结果确定必须保留的积分/前端定义；解析滤波审计不能替代实际符号对齐与有限时隙核验。
2. 接多小区有限采样波形、接收滤波/相关、测距与位置估计；明确采样率、记录长度、序列、相位/时钟、相对时序和实际检测门限，分模块测试。
3. 在估计位置可用后接首次lambda=1、几何重算和功率反馈；继续独立保存oracle基准及失败原因。首次lambda=1不等于凭空使用真实位置初始化。
4. 再组织小样本和50次Monte Carlo及敏感性；先确认失败分母、接收机假设和观测期，不用解析CRB与稳态方差数值接近代替观测公平性。
