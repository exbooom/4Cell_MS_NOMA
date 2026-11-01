# 单小区 MS-NOMA 快报复现实施指南

**目标论文**：Lu Yin, Jiameng Cao, Kaiqin Lin, Zhongliang Deng, Qiang Ni, “A Novel Positioning-Communication Integrated Signal in Wireless Communication Systems,” *IEEE Wireless Communications Letters*, vol. 8, no. 5, pp. 1353-1356, October 2019. DOI: 10.1109/LWC.2019.2917670。

**用途**：给负责实现的 AI/研究者作为工程说明。本文定义复现范围、代码结构、公式推导、检查方法和不确定项；不是声称已经复现论文的结果。

**复现策略**：单小区、无凸优化。先复用当前 `MS_NOMA` 的信号生成；新增解析性能计算和图2；先完成 AWGN，再单独加入 ETU 衰落。代码输出、参数和图表都要能够追溯到这里的配置与公式。

## 0. 执行边界

1. 只复现2019年快报的单小区结果：图1的波形结构、图2(a)平均 BER、图2(b)每通信子载波 BER、图2(c)测距精度。此阶段不做多基站位置解算、CP4A 或功率优化。
2. 保留现有代码作为波形生成原型；不要重写掉原文件。先建立新入口和新计算函数，复现结果稳定后再判断是否需要调整波形生成器。
3. 将**解析性能计算**与**采样波形接收机仿真**分开。论文图2主要由解析 BER 和误差公式得到；先实现解析图，不要把接收机链路当作复现所有曲线的前置条件。
4. 公式原样、公式推导候选、工程假设必须分别标记。不能通过调未披露参数来强行贴近图形。
5. 论文给出的“range measurement accuracy”是测距标准差/误差精度，单位米；由时延/码相位方差换算时使用 $c_0^2$，标准差使用 $c_0$。

## 1. 工作区现状与复用判断

目标工作区：`F:\OneCell_MS_NOMA`。当前源文件集中在 `F:\OneCell_MS_NOMA\MS_NOMA`：

| 文件 | 当前用途 | 复用建议 |
|---|---|---|
| `generateSignal.m` | 生成 QPSK 通信、BPSK/±1 定位码片和叠加信号 | 第一阶段直接调用。验证功率、码片/符号时间尺度、实际带宽 |
| `calPSD.m` | Welch 风格分段谱估计 | 用于波形 sanity check；不替代理论 PSD |
| `calBER.m` | 当前长文版的单小区 sinc² 干扰和近似 BER | 可复用 sinc² 重叠计算思路；BER 系数、能量定义和 CPR 定义需为快报单独实现 |
| `main.m` | 波形/PSD/码相关性和后续长文 BER 展示 | 新建独立 `run_letter2019.m`，不要把快报逻辑继续堆入长文入口 |
| `L1CPcodes.mat` | 32条、每条10230 chip 的 ±1 码表 | 选码用于波形演示；论文图2解析结果不应依赖某一 PRN，除非评估码相关/接收机性能 |

已有生成器使用 `G=80, M=20, dfC=30 kHz`，因而 $N=(M+1)G-1=1679$，有效带宽为50.4 MHz。快报使用 `G=50, M=20, N=1049`。第一项实际配置改动应是将复现配置设为 `G=50` 并自动校验 $N$，不可沿用 `main.m` 里目前的 `G=80`。

快报定义 $\kappa=2G-1$，用它计算 CPR；当前长文代码使用 $2G$。建立 `letter2019` 公式模块，避免覆盖长文实现。

## 2. 论文仿真参数台账

### 2.1 原文明确值

| 参数 | 快报设置 | 单位/说明 |
|---|---:|---|
| $G$ | 50 | $T_c/T_p$，定位码片比通信符号快50倍 |
| $M$ | 20 | 定位用户/定位子载波数 |
| $N$ | 1049 | $N=(M+1)G-1$ |
| $\Delta f_c$ | 未在数值段明确重述 | 沿用目标系统配置时需显式写入 config；不得假称此快报给了值 |
| $\Delta f_p$ | $G\Delta f_c$ | Hz |
| 通信调制 | QPSK | 快报数值段明确 |
| 定位调制/码 | BPSK、长度 $L$ 的扩频码 | 码长度 $L$ 未在性能段给出 |
| $B_{fe}$ | $2B$ | 双边前端带宽 |
| $B_L$ | 0.2 | Hz，DLL loop noise bandwidth |
| $T_{coh}$ | 0.02 | s，predetection/coherent integration time |
| $D$ | 0.02 | chip，early-late spacing |
| 快报图2(c) $P_p/N_0$ | 45 dB-Hz | 定位载波与单边 PSD 之比 |
| 信道 | AWGN/flat 与 ETU fading | ETU 分支参数需声明 |
| 多径路径 | 快报未给出 | 使用标准 ETU 实现时记录工具箱版本、采样率、多普勒等 |
| 基站/位置 | 单小区 | 不产生二维位置几何 |

由环路参数：

$$
a=B_L(1-0.5B_LT_{coh})=0.2(1-0.002)=0.1996\ \mathrm{s}^{-1}.
$$

### 2.2 配置文件应明确的工程假设

以下值/定义没有被论文唯一指定，放进 `config_letter2019.m` 并在结果旁输出：

```matlab
cfg.seed = 2026;                 % 自选并记录
cfg.bandwidthMode = "integer_multiscale"; % 默认严格满足整数多尺度
cfg.deltaFcHz = ...;             % 需从系统配置/原波形设置确定
cfg.noiseConvention = ...;       % PSD 单/双边定义
cfg.berVariant = ...;            % paper_literal / standard_qpsk_candidate
cfg.Gamma = ...;
cfg.gamma = ...;
cfg.codeLength = ...;
cfg.ETUmaxDopplerHz = ...;
cfg.ETUnormalization = ...;
cfg.numMonteCarlo = ...;
cfg.prsModel = ...;              % external_reference / comparison_disabled
```

若从现有长文代码继承 `dfC=30 kHz`，则 $T_c=33.333\ \mu s$、$T_p=0.6667\ \mu s$、$B=(N+1)\Delta f_c=31.5$ MHz、$B_{fe}=63$ MHz。这是**沿用现有项目的工程配置**，不是快报数值段明确给出的载波间隔。先以此做内部对照可以，但图注须注明。若尚未确定物理带宽，默认采用这一具名分支 `legacy_30kHz_assumption` 开工，并额外扫描其他带宽；不必因此阻塞全部实现。

【可核验线索】图2(a)上轴标 $P_p/N_0+CPR$，图2(c)上轴标 $E_b/N_0$。线性量恒等式为 $\mathrm{CN0}\times\mathrm{CPR}=\kappa\rho/T_c$，因此上、下两种坐标轴的 dB 数字只差 $10\log_{10}(\kappa/T_c)$。若从图2(c)估读 CPR=10 dB、CN0=45 dB-Hz 时 $\rho$ 约−7 dB，则该差值约62 dB，提示 $T_c\approx99/10^{6.2}\approx62.47\ \mu s$、$\Delta f_c\approx16.01$ kHz，而非30 kHz。该数值来自读图，不能标为正文参数；可将其作为 `axis_inferred` 敏感性分支，并记录读图误差。绝对测距米数依赖这一时间尺度，所以选择30 kHz后不能要求图2(c)逐点吻合。

快报关系 $N=B/\Delta f_c-1$ 和 $M=B/\Delta f_p-1$ 配上整数 G 给出 $N+1=(M+1)G$。按本文建议用 $B_{eff}=(N+1)\Delta f_c$ 反算有效带宽，使离散仿真严格一致。

## 3. 公式符号与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| $B$ | MS-NOMA 总带宽 | Hz |
| $B_{fe}$ | 双边前端带宽 | Hz |
| $\Delta f_c,\Delta f_p$ | 通信、定位子载波间隔 | Hz |
| $G=\Delta f_p/\Delta f_c$ | 多尺度比 | 无量纲正整数 |
| $N,M$ | 通信、定位子载波/用户数 | 整数 |
| $T_c=1/\Delta f_c$ | 通信符号周期 | s |
| $T_p=1/\Delta f_p=T_c/G$ | 定位码片周期 | s |
| $P_c$ | 每条通信子载波功率 | W |
| $P_{p,m}$ | 第 m 条定位子载波功率 | W |
| $N_0$ | 论文定义的环境噪声单边 PSD | W/Hz |
| $I(n)$ | 定位能量泄漏到通信子载波 n 的等效 PSD | W/Hz |
| $E_b$ | **论文记号** $P_cT_c$，文本称 communication symbol energy | J；并非通常的 QPSK 每比特能量 |
| $\mathrm{CPR}_m$ | 通信/定位功率比例 | 线性比值；绘图用 dB |
| $\sigma^2_{LB,m}$ | AWGN 下码相位估计误差下界 | 快报式子实际具有 s² 量纲 |
| $\sigma^2_{t,m}$ | DLL tracking error 方差 | 快报式子实际具有 s² 量纲 |
| $\sigma_{range}$ | 测距标准差 | m |

定义 $\operatorname{sinc}(x)=\sin(\pi x)/(\pi x)$，且 $\operatorname{sinc}(0)=1$。所有功率计算用线性单位；仅出图转 dB。

## 4. 多尺度信号与现有生成器

快报第II节的结构是：通信符号按 $T_c$ 更新；每个定位用户的 BPSK 伪随机码按 $T_p$ 更新并调制到间隔为 $\Delta f_p$ 的独立定位子载波；定位和通信波形在时域叠加。定位信号持续发送。

便于实现的复基带形式：

$$
s(t)=\sum_{n=1}^{N}\sqrt{P_c}\,d_n(t)e^{j2\pi f_{c,n}t}
+\sum_{m=1}^{M}\sqrt{P_{p,m}}\,c_m(t)e^{j2\pi f_{p,m}t},
$$

其中 $d_n$ 为单位平均能量 QPSK 符号流，$c_m$ 为 ±1 定位码片流。定位子载波与通信子载波不要求正交；定位信号的 PSD 旁瓣造成跨尺度干扰。定位子载波彼此正交，理想同步模型下不同 P-User 间无干扰。

现有 `generateSignal.m` 已按通信符号/定位码片两个时间尺度构造复基带信号，适合作为结构演示原型。先加以下运行时断言，不要先改其数学结构：

```matlab
assert(cfg.G == 50);
assert(signal.N == (cfg.M + 1) * cfg.G - 1);
assert(abs(signal.Tc / signal.Tp - cfg.G) < 1e-12);
assert(size(signal.chips,1) >= cfg.M);
assert(all(ismember(unique(signal.chips(:)),[-1,1])));
```

之后对每类分量核对：平均通信功率、每路定位功率、总平均功率、有效子载波频率位置和 PSD 面积。**复用波形生成器不等于 BER/测距曲线已复现**：图2可先由理论模型直接产生。

### 4.1 复用生成器必须加的适配层

原生成器按 $P_p=2GP_c/10^{CPR_{legacy}/10}$ 设功率，快报按 $P_p=\kappa P_c/10^{CPR_{letter}/10}$。为了保持原文件不动且输出相同目标功率，应传入：

$$
CPR_{legacy,dB}=CPR_{letter,dB}+10\log_{10}(2G/\kappa).
$$

G=50时偏移约0.04365 dB。这个差很小，但应显式适配，而不是只改G就声称定义一致。`cfg.communicationPower` 在旧生成器中是**总通信功率**；如果新性能模块使用每子载波 $P_c$，适配层应传 `cfg.communicationPower=N*Pc`，否则相差1049倍。

旧生成器要求 CPR 为有限值，不能直接传 `Inf`；无定位基线复用其通信分量并令定位分量为零，或在适配层明确定义此分支。当前生成器只支持所有定位子载波等功率；图2(b) DPp 属于解析计算，不必立即改生成器。需要 DPp 采样波形时再增加独立扩展函数。

其每个定位码片采样数为 `NFFT/G`。应断言为整数；当前 `NFFT=(N+1)*oversampling` 且 `N+1=(M+1)*G` 时自动成立。L1CP码表不是论文已确认的码族，有限重复码的实际PSD也不必逐点等于独立随机码的sinc²统计PSD。

## 5. 归一化 PSD 与频率重叠

### 5.1 定位、通信单子载波 PSD

单位面积的矩形脉冲 PSD：

$$
G_{p,m}(f)=T_p\operatorname{sinc}^2[(f-m\Delta f_p)T_p],
$$
$$
G_{c,n}(f)=T_c\operatorname{sinc}^2[(f-n\Delta f_c)T_c].
$$

对全频率积分面积都为1。实际功率谱分别是 $P_{p,m}G_{p,m}(f)$、$P_cG_{c,n}(f)$。

### 5.2 定位到通信子载波的等效干扰【快报式（2）】

在通信子载波中心 $f=n\Delta f_c$ 采样定位第 m 个子载波 PSD：

$$
\bar P_{p,m}(n)=P_{p,m}G_{p,m}(n\Delta f_c).
$$

代入 $G_{p,m}$：

$$
\begin{aligned}
G_{p,m}(n\Delta f_c)
&=T_p\operatorname{sinc}^2[(n\Delta f_c-m\Delta f_p)T_p]\\
&=T_p\operatorname{sinc}^2\left(\frac nG-m\right)\\
&=T_p\operatorname{sinc}^2\left(m-\frac nG\right),
\end{aligned}
$$

因为 sinc² 是偶函数。因此：

$$
I(n)=\sum_{m=1}^{M}\bar P_{p,m}(n)
=T_p\sum_{m=1}^{M}P_{p,m}W_{nm},\quad
W_{nm}=\operatorname{sinc}^2\left(m-\frac nG\right).
$$

`W` 是预计算的 `N × M` 干扰重叠矩阵。注意这是论文的等效干扰采样模型；不等同于对某个具体接收滤波器输出干扰方差的普遍定义。

等功率场景：$P_{p,m}=P_p$，故 $I(n)=P_pT_p\sum_mW_{nm}$。不等功率场景使用相同公式，将定位功率向量代入即可。

## 6. 通信 BER 推导【快报式（1）】

### 6.1 论文印刷公式

$$
\mathrm{BER}(n)=K\operatorname{erfc}\left(\frac{\lambda P_cT_c}{I(n)+2N_0}\right).\tag{1 原文}
$$

论文说明 $K,\lambda$ 由调制和编码方式决定；本篇没有给出足以唯一恢复它们的数值表。公式的 `erfc` 自变量按印刷版看不到平方根。更重要的是，标准相干 Gray QPSK 无干扰 BER 是

$$
P_b=\frac12\operatorname{erfc}\sqrt{E_{bit}/N_{0,std}}.
$$

这里明确约定物理实通带 AWGN 双边 PSD 为 $N_{0,std}/2$（对应单边 $N_{0,std}$），复基带噪声自相关为 $N_{0,std}\delta(t-t')$。不使用含混的 `N0_2-sided` 变量名把 PSD 取值与其参数混为一谈。

快报沿用“$E_b=P_cT_c$”作为通信符号能量的文字，与 QPSK 每比特能量又相差2。因此实施者必须把下面两个分支分别保留，不能静默加根号或换能量定义。

### 6.2 分支 A：逐字论文模型

```text
BER_paper(n) = K * erfc(lambda * Pc * Tc / (I(n) + 2*N0))
```

按此分支画出的干净曲线是校验 $K,\lambda,N_0$ 的第一手依据。因为论文未列出具体系数，记录拟采用的 $K,\lambda$ 以及来源；若只根据图2曲线反推，明确标记为“从图像估读/拟合”，不称为论文明确参数。

### 6.3 分支 B：标准 QPSK 物理候选

若根据自洽噪声定义、调制推导而验证标准形式，可使用：

$$
P_b(n)=\frac12\operatorname{erfc}\sqrt{\frac{E_{bit,rx}}{N_{0,std}+I_{eff}(n)}}.
$$

这里 $E_{bit,rx}=gP_cT_c/2$（单位平均 QPSK 符号能量 $gP_cT_c$ 除以2 bit），但 $I_{eff}$ 如何从论文的 PSD 采样 $I(n)$ 换算到匹配滤波后的等效噪声，取决于滤波器、积分时间、单/双边 PSD 约定。默认具名候选可取 `N0_std=N0_paper, Ieff=I`，与现有代码口径一致，但必须标记为待验证的物理近似。

为了能立即实现并诊断图2，额外提供 `sqrt_axis_calibrated` 候选：令论文横轴 $\rho=P_cT_c/N_0$，并使无定位基线满足 $0.5\operatorname{erfc}\sqrt\rho$，则可定义

$$
P_b(n)=\tfrac12\operatorname{erfc}\sqrt{\frac{2P_cT_c}{I(n)+2N_0}}.
$$

这个形式的根号、系数2属于基线校准候选，并非论文明确给出的修正。三个候选可并行计算，无须等待参数确认才开始写程序。物理比特口径和论文标称横轴口径必须独立标注。

**不能**把上述标准公式只加进 `calBER`，却继续用原文式（12）推优化阈值。本文无优化，但 BER 分支仍需数学自洽。

### 6.4 Eb/N0 扫描与噪声反算

横轴记号按快报的定义为 $E_b=P_cT_c$。给定横轴 $\rho_b=10^{(E_b/N_0)_{dB}/10}$，可令

$$
N_0=\frac{P_cT_c}{\rho_b}
$$

来扫描 BER，并同时报告定位比值 $P_p/N_0$。对等定位功率：

$$
\mathrm{CPR}_m=\frac{\kappa P_c}{P_{p,m}},\qquad\kappa=2G-1=99,
$$
$$
P_p=\frac{\kappa P_c}{10^{\mathrm{CPR}_{dB}/10}}.
$$

无定位信号对应 $P_{p,m}=0$、$I(n)=0$、CPR=$+\infty$。

高信噪比下应观察定位干扰导致的 BER 平台；定位功率越小（CPR越大），平台越低。干净 BER 曲线可用于检查噪声与 $K,\lambda$ 约定。在论文定义不闭合前，不要只凭“这是 QPSK”就强设 $K=1/2$、$\lambda=\sqrt2$ 或其他组合。

## 7. 快报测距误差推导

快报考虑 BPSK 定位分量的码跟踪。图2(c)报告测距米数；分析从码相位/时延方差开始。

### 7.1 AWGN 码相位/时延下界【式（3）】

论文给出的形式：

$$
\sigma^2_{LB,m}=\frac{a}{(2\pi)^2}\left\{
\int_{B_0-B_{fe}/2}^{B_0+B_{fe}/2}
f^2\left[\frac{P_{p,m}G_{p,m}(f)}{N_0+NP_cG_c(f)}\right]df
\right\}^{-1},\tag{3}
$$

其中 $G_c(f)=\frac1N\sum_{n=1}^NG_{c,n}(f)$。论文数值分析用 OFDM 总通信 PSD 在整个带宽大致平坦的近似：

$$
G_c(f)\approx\frac1B\,\mathbf 1_{|f-B_0|<B/2}.
$$

【建模区别】从RF频率换到本地基带不能简单把式（3）里所有 $f^2$ 机械替换为 $x^2$：原式本身已经含有中心频率简写/时延信息定义问题。为匹配其式（5）—（8）的码包络时延模型，后续统一用本地基带 $x$ 计算频率二阶矩，不能把3.5GHz等RF载频加入该二阶矩。

采用论文推导中“通信平坦频带以目标P-Sub为中心”的对称近似，令 $Q(x)=G_{p,m}(x+m\Delta f_p)=T_p\operatorname{sinc}^2(xT_p)$。此谱没有有限支撑，**不能忽略定位信号在通信带外的旁瓣**；通信带外只含噪声的部分正是下界积分的一部分。

若严格保留通信实际公共频带，其在目标P-Sub本地频率中通常不是 $[-B/2,B/2]$，会随m平移。需另外提供 `physical_band` 分支，以有限载波和定义 $G_c$，不能将其结果不加说明地称为原式（5）对称近似。

PDF式（5）第一段印刷成 $\int_{B_{fe}/2}^{B/2}$，由于 $B_{fe}>B$，该段为负。将其作为 `paper_literal_reversed_limits` 仅供诊断；物理上带外能量贡献为正，应使用如下正向形式：

$$
\begin{aligned}
\sigma^2_{LB,m}=\frac{a}{(2\pi)^2}\Bigg\{&2\int_{B/2}^{B_{fe}/2}
\frac{x^2P_{p,m}G_{p,m}(x+m\Delta f_p)}{N_0}dx\\
&+\int_{-B/2}^{B/2}
\frac{x^2P_{p,m}G_{p,m}(x+m\Delta f_p)}{N_0+NP_c/B}dx
\Bigg\}^{-1}.
\end{aligned}\tag{5 正向积分候选}
$$

这等价于先定义总通信PSD平台 $U=NP_c/B$，再计算本地基带的正值信息积分

$$
J=\int_{-B_{fe}/2}^{B_{fe}/2}\frac{P_p x^2Q(x)}{N_0+U\mathbf1_{|x|<B/2}}dx,\qquad v_{LB}=\frac{a}{4\pi^2J}.
$$

`rangingLowerBoundExact` 的“exact”应说明是**给定平坦PSD近似模型后的精确数值积分**，不代表保留所有有限OFDM子载波的真实物理模型。

### 7.2 式（6）的近似解析式及 CPR 推导

原文近似用 $\operatorname{sinc}(BT_p)\approx0$、$\operatorname{sinc}(B_{fe}T_p)\approx0$，且 $B_{fe}>B$，将矩形 chip PSD 的有限带积分拆成带外（只受噪声）和带内（受噪声+通信）两部分，得到：

$$
\sigma^2_{LB,m}\approx0.25aT_p^2\frac{\mathrm{CPR}_m}{E_b/N_0}
\left(B_{fe}-\frac{E_b/N_0}{E_b/N_0+1}B\right)^{-1}.\tag{6}
$$

式（6）里的 $E_b/N_0$ 采用线性值，不是 dB 数字。下面给出从正值积分推到此式的中间步骤，并保留有限N/κ的误差。

定义可直接编码的偶函数积分：

$$
L(F)=\int_{-F/2}^{F/2}x^2Q(x)dx
=\frac{1}{\pi^2T_p}\int_{-F/2}^{F/2}\sin^2(\pi xT_p)dx
=\frac{1}{2\pi^2T_p}\left[F-\frac{\sin(\pi FT_p)}{\pi T_p}\right].
$$

因此不用数值积分也能得到平坦PSD模型的有限带精确结果：

$$
J=P_p\left[\frac{L(B_{fe})-L(B)}{N_0}+\frac{L(B)}{N_0+U}\right],\qquad v_{LB}=\frac{a}{4\pi^2J}.
$$

这可作为积分器独立验证。令 $Q_0=P_p/N_0$、$r=U/(N_0+U)$，忽略有限带的正弦边缘项后：

$$
J\approx\frac{Q_0}{2\pi^2T_p}(B_{fe}-rB),\qquad
v_{LB}\approx\frac{aT_p}{2Q_0(B_{fe}-rB)}.
$$

又由 $B=(N+1)/T_c$、$\rho=P_cT_c/N_0$：

$$
\frac{U}{N_0}=\frac{N}{N+1}\rho\equiv\zeta\rho,\quad
r=\frac{\zeta\rho}{1+\zeta\rho},\quad\zeta=\frac{N}{N+1}.
$$

取大N近似 $\zeta\approx1$，并采用 $\kappa\approx2G$，就能把上式写成论文式（6）。这不是对任意有限参数都严格相等。定义

$$
\kappa=2G-1,\qquad E_b=P_cT_c,\qquad
\mathrm{CPR}_m=\frac{\kappa P_c}{P_{p,m}}.
$$

令 $\rho=E_b/N_0=P_cT_c/N_0$，可核对：

$$
\frac{\mathrm{CPR}_m}{\rho}
=\frac{\kappa P_c/P_{p,m}}{P_cT_c/N_0}
=\frac{\kappa N_0}{T_cP_{p,m}}.
$$

这说明第一项由定位载噪比控制，且式（6）最终量纲为秒平方。变成距离标准差：

具体地，$\mathrm{CPR}/\rho=\kappa/(T_cQ_0)$。论文式（6）前面的系数为 $aT_p^2\kappa/(4T_cQ_0)$；正向积分近似的前系数为 $aT_p/(2Q_0)$，两者比值为 $\kappa/(2G)=99/100$。应同时输出 `paper_eq6` 和 `finite_sum_consistent`，避免把这约1%的方差差别当成积分失败。

$$
\sigma_{range,LB,m}=c_0\sqrt{\sigma^2_{LB,m}}.
$$

若画下界图，横轴 CPR 从 dB 转线性；固定 $P_p/N_0$ 时可用 $\mathrm{CPR}=\kappa(P_c/N_0)/(P_p/N_0)$ 求 $E_b/N_0$，但要维持论文横轴/固定量定义一致。

### 7.3 DLL tracking 误差精确积分【式（7）】

定义目标定位子载波本地基带 $x=f-m\Delta f_p$，则其 PSD 为 $T_p\operatorname{sinc}^2(xT_p)$。相干 early-late DLL 的误差方差写成：

$$
\sigma^2_{t,m}=\frac{a\displaystyle\int_{-B_{fe}/2}^{B_{fe}/2}
 [N_0+S_c^{local}(x)]G_{p,m}(x+m\Delta f_p)
 \sin^2(\pi xDT_p)dx}
{(2\pi)^2P_{p,m}\left[\displaystyle\int_{-B_{fe}/2}^{B_{fe}/2}
 xG_{p,m}(x+m\Delta f_p)\sin(\pi xDT_p)dx\right]^2}.\tag{7}
$$

其中 $S_c^{local}(x)$ 是接收定位用户处的本地基带通信PSD。快报式（7）印刷写作 $G_s(f)$ 而把目标定位PSD写成 $G_{p,m}(f+m\Delta f_p)$；实现统一用具名本地PSD接口，避免重复或遗漏频移。对称平坦模型为 $S_c^{local}(x)=U\mathbf1_{|x|<B/2}$。若使用绝对频率函数，则先做 $S_c^{local}(x)=G_s(x+f_{p,m})$。

AWGN、平坦信道、单位增益下，绝对频率通信PSD模型为：

$$
G_s(f)=NP_cG_c(f).
$$

论文为窄带早迟间隔 $D\to0$ 使用 $\sin(\pi xDT_p)\approx\pi xDT_p$。数值复现应同时保留 exact-sine 与 small-D 两个分支，检查给定 $D=0.02$ chip 时近似误差，不要默认 D 无关紧要。

### 7.4 DLL 近似式【式（8）】

原文在 OFDM 通信 PSD 近似平坦、$B_{fe}\gg B_p$、小早迟间隔后得到：

$$
\sigma^2_{t,m}\approx0.25aT_p^2\left[
\frac{2}{B_{fe}T_p(P_{p,m}/N_0)}
+\frac{B}{B_{fe}^2}\mathrm{CPR}_m
\left(\beta+\bar P_{c,multi-path}\right)\right].\tag{8}
$$

AWGN 与单一路径分支使用 $\beta=1$、$\bar P_{c,multi-path}=0$，所以括号最后因子为1。论文定义 $\beta=1$ 对应 Rician、$\beta=0$ 对应 Rayleigh，并由归一化传播幅度期望给出 $E(\alpha^2)=\beta+\bar P_{c,multi-path}$；衰落复现必须按该归一化定义配置，不能把工具箱默认 fading gain 不经归一化地代入。

同样按距离换算：$\sigma_{range,t}=c_0\sqrt{\sigma^2_{t,m}}$。

### 7.4.1 式（8）的中间推导与可计算核

定义 $u(x)=\sin(\pi xDT_p)$、$V(x)=N_0+S_c^{local}(x)$：

$$
A_0=\int xQ(x)u(x)dx,\quad A_N=N_0\int Q(x)u(x)^2dx,\quad A_C=\int S_c^{local}(x)Q(x)u(x)^2dx,
$$

积分域均为 $[-B_{fe}/2,B_{fe}/2]$，于是 $v_t=a(A_N+A_C)/(4\pi^2P_pA_0^2)$。无量纲变量 $y=xT_p$ 时，节点域是 $[-B_{fe}T_p/2,B_{fe}T_p/2]$，可复用 `sinc(y)` 与 `sin(pi*D*y)`，减少单位和采样错误。

小D近似后，利用上面的 $L(F)$：

$$
A_0\approx\pi DT_pL(B_{fe}),\quad
A_N\approx N_0\pi^2D^2T_p^2L(B_{fe}),\quad
A_C\approx U\pi^2D^2T_p^2L(B).
$$

故小D、但保留有限带边缘项的结果为：

$$
v_t\approx\frac{a\,[N_0L(B_{fe})+UL(B)]}{4\pi^2P_pL(B_{fe})^2}.
$$

再忽略正弦边缘项：

$$
v_t\approx\frac{aT_p}{2P_pB_{fe}^2}(N_0B_{fe}+UB)
=\frac{aT_p}{2B_{fe}Q_0}+\frac{aT_pNP_c}{2P_pB_{fe}^2}.
$$

论文式（8）的AWGN通信项是 $aT_p^2 B\kappa P_c/(4P_pB_{fe}^2)$。它与上式通信项的比值为 $\kappa(N+1)/(2GN)$，在本篇参数下约0.991。保留该有限N差异；衰落统计近似把通信项乘 $E(\alpha^2)$。

### 7.4.2 lower bound 与 DLL 并不普遍相等

相同噪声+干扰PSD下，下界积分权重为 $Q/V$，DLL的分子却含 $VQ$。Cauchy-Schwarz给出

$$
\left(\int xQ(x)u(x)dx\right)^2\le
\left(\int V(x)Q(x)u(x)^2dx\right)
\left(\int \frac{x^2Q(x)}{V(x)}dx\right),
$$

所以自洽积分模型应满足 $v_t\ge v_{LB}$。无通信时V是常数且D足够小时等号趋近成立；有带限通信干扰时V不再是常数，不能因为是AWGN传播环境就要求两条曲线完全相等。论文关于“converge to the lower bound”的文字应结合其图中条件理解。

特别是 $B_{fe}=2B$ 且通信功率很大时，带外仍只有环境噪声，lower bound保留带外信息，其误差可以趋于有限平台；固定 early-late DLL 的误差则随通信项增长。这正是图2(c)应检查的形状，而不是强行把两条曲线重合。

### 7.5 积分和近似公式实施细节

- 在 `f` 域直接积分时，载频 $B_0$ 极大，而有效基带带宽只有 MHz 量级。实现宜改用局部频率 $x=f-f_{p,m}$，避免浮点分辨率问题。
- 数值单位先统一 SI：Hz、秒、W、W/Hz。$D=0.02$ chip 以 chip 为单位输入公式时，要注意 $xDT_p$ 的定义中 $D$ 是 chip 数、$T_p$ 为每 chip 秒，乘积单位秒。
- 计算完整下界和 DLL 时使用自适应积分或按带宽/谱瓣分段加密网格，比较容差收紧后的结果。
- 公式以方差形式返回 `variance_s2`；绘图函数单独乘 $c_0$ 并开方。严禁把秒标准差直接标成米。
- 检查功率标度：其他参数固定时增大 $P_{p,m}/N_0$，误差应下降；增加通信干扰时 DLL 误差应增加。

## 8. 衰落信道实现边界

论文结果包含 flat/AWGN 与 ETU fading。先完整跑通 AWGN 解析结果，再做 ETU。

### 8.1 BER 的衰落分支

原文文字假设调制符号/定位扩频序列独立，并讨论反射定位分量的干扰。在 Rician 下称其较弱可忽略，在 Rayleigh 多径下可能需要计入。单小区图2(a)至少要区分：

- flat：通信链路增益固定为1；定位信号干扰按式（2）。
- ETU：对每个通信载波应用一致定义的信道功率增益；若定位用户接收的多径定位干扰也建模，按论文假设决定保留/忽略。

快报对 ETU 标准、移动速度/最大 Doppler、信道增益归一化、同步、均衡方式没有在四页内完整指定。实施 ETU 时在结果元数据中记录 MATLAB 版本、Communications Toolbox 版本、模型路径功率、最大 Doppler、采样率、过滤瞬态方式、平均方式。不要把“用了 ETU 对象”当成已复现论文 fading 图。

### 8.2 测距 fading 分支

在 tracking 式（7）中以接收通信 PSD $G_s(f)$ 表示多径通信功率谱；若采用论文简化慢衰落模型，可用

$$
G_s(f)\approx NP_cG_c(f)E(\alpha^2),\qquad
E(\alpha^2)=\beta+\bar P_{c,multi-path}.
$$

要复现 exact fading 曲线则需要具体多径 tap 和信号模拟。可以先复现快报式（8）的统计近似；不得称为逐路径精确仿真。

## 9. 图2逐图实现配方

### 9.1 图2(a)：平均 BER

输入：$G=50,M=20,N=1049$；横轴 $E_b/N_0=0:20$ dB；图中 flat 与 fading；CPR 曲线含15/20/30 dB以及无定位基线。按原图读数为主，最终用像素数字化/曲线核验时记录读图误差。

流程：

1. 由 $E_b/N_0$ dB 转成线性 $\rho$。
2. 指定 $P_c$ 的归一化（例如 $P_c=1$ W 只影响归一化噪声单位），$N_0=P_cT_c/\rho$。
3. CPR dB 转线性，按 $P_p=\kappa P_c/\mathrm{CPR}$ 计算等定位功率。
4. 每个横轴点、每个 CPR 计算向量 $I(n)=T_pW P_p\mathbf1$。
5. 用指定 BER 分支得到 $N$ 个载波 BER，做算术平均。
6. 输出 flat；另按 ETU 分支计算 fading 平均（独立 Monte Carlo 或统计期望，配置明确）。
7. 对数纵轴绘制，横轴标论文 $E_b/N_0$ dB，图例说明 CPR 与信道。

验收：干净 BER 随 $E_b/N_0$ 下降；定位越弱，曲线越接近无定位；低 CPR 高 SNR 有干扰平台；定位载波泄漏的平均对称/边缘效应可解释。

### 9.2 图2(b)：每通信子载波 BER

固定 $E_b/N_0=5$ dB、CPR=15 dB。画三条：

- `NP`：无定位，所有通信载波 BER 相同（若平坦、相同功率）。
- `IPp`：所有定位用户功率相同，沿频带的 BER 大致相近。
- `DPp`：定位功率不相同，定位子载波泄漏重叠使不同通信子载波 BER 起伏。

论文未给 `DPp` 的定位功率向量，不能唯一还原原蓝线。用注明种子的示例向量（例如在总定位功率相同条件下对用户功率作有界非均匀分配）画机制演示，并标为 `illustrative DPp, not original vector`。保存向量。BER 纵轴采用线性范围以便对应原图。

### 9.3 图2(c)：测距精度

固定 $P_p/N_0=45$ dB-Hz。横轴 CPR=10…30 dB。使用式（6）测距下界和式（7）DLL tracking 精确积分，分别叠加式（8）近似。图例包含：MS-NOMA lower bound exact/approx、DLL exact/approx、no communication baseline，以及 PRS 对照（若参考文献模型完成）。

固定定位载噪比下，改变 CPR 等价于改变通信功率/通信定位比例。按定义核算：

$$
\mathrm{CPR}=\kappa\frac{P_c}{P_p}
=\kappa\frac{P_c/N_0}{P_p/N_0},
$$

其中 $P_p/N_0=10^{45/10}$ Hz。由 CPR 反推 $P_c/N_0$。然后用式（6）、（8）算解析近似；精确式用自身积分所需的 $P_c,P_p,N_0$ 同比缩放，绝对共同尺度会消去（在模型线性成立时）。

PRS 下界引用快报参考文献[10]，即 *Joint channel and time delay estimation for LTE positioning reference signals*，不是期刊长文所引的参考文献[28]。如果参考文献[10]尚未核对，则将 PRS 曲线标为“待复现/未绘制”，不要用别的模型冒充原曲线。

## 10. 建议代码架构

原文件保持不动。新建以下结构：

```text
F:/OneCell_MS_NOMA/
  REPRODUCTION_GUIDE.md
  MS_NOMA/                         # 现有作者/整理代码
  reproduction/
    config_letter2019.m            # 快报参数和显式工程假设
    run_letter2019.m               # 编排、保存结果、图表
    +letter2019/
      numerology.m                 # Tc/Tp/N/B/Bfe/κ 与一致性检查
      makeSignalFromExisting.m     # 调用既有 generateSignal，不复制实现
      sincSquared.m                # 统一 sinc 定义及零点处理
      interferenceMatrix.m         # N×M重叠矩阵W
      positioningInterference.m    # 根据W和Pp向量算I(n)
      berPaperLiteral.m            # 快报式(1)
      berStandardQpskCandidate.m   # 对照候选，独立命名
      berSweep.m                   # 图2(a)、图2(b) 数据
      rangingLowerBoundExact.m     # 式(3)/(5) 数值积分
      rangingLowerBoundApprox.m    # 式(6)
      dllTrackingExact.m           # 式(7) 数值积分
      dllTrackingApprox.m          # 式(8)
      rangeUnits.m                 # variance s² -> std m
      etuChannel.m                 # 可选；记录所有模型设置
      prsLowerBound.m              # 完成参考文献[10]后才实现
      plotFigure2.m                # 只读结果数据并绘图
      saveRun.m                    # config、seed、数组、版本、日志
    tests/
      test_numerology.m
      test_interference.m
      test_ber.m
      test_ranging.m
```

### 10.1 函数接口建议

```matlab
[num, report] = letter2019.numerology(cfg)
W = letter2019.interferenceMatrix(num)
[I, terms] = letter2019.positioningInterference(num, W, Pp)
ber = letter2019.berPaperLiteral(Pc, Tc, I, N0, cfg)
varS2 = letter2019.rangingLowerBoundExact(num, cfg, power)
varS2 = letter2019.dllTrackingExact(num, cfg, power, channelPsd)
rangeStdM = letter2019.rangeUnits(varS2)
```

所有性能函数返回线性单位数组和辅助分项；绘图不负责重新计算物理量。

## 11. 验收测试

### 数学/数值单元测试

1. Numerology：`G=50,M=20` 得 `N=1049`；`Tc/Tp=50`；有效带宽等于 `(N+1)*dfC`。
2. sinc²：输入0得1；非零整数得0（数值容差）；W非负且尺寸 N×M。
3. 干扰：矩阵乘法 `Tp*W*Pp(:)` 与显式 m 循环逐元素一致。
4. CPR：`Pp=κ Pc/CPR` 反算 CPR 误差小于浮点容差；CPR=Inf 时定位功率为0且 BER 走无定位支路。
5. BER：扫 Eb/N0 后 clean BER 单调下降；对固定噪声，定位功率增加时 BER 不下降。检查 BER 处于[0,0.5]。
6. Measurement：`variance_s2 >= 0`；输出米数严格按 `c0*sqrt(variance_s2)` 转换；数值积分加密后应收敛到同一积分模型的解析核结果。small-D/大带宽近似的模型误差不会因为积分网格加密自动消失，必须单独报告。
7. Eq. (8) AWGN：令 $\beta=1,\bar P_{c,multi-path}=0$；关闭通信分量时 tracking 误差应趋近定位单独基线。
8. 所有配置和随机种子重复运行，可复现数据数组。

建议初始积分相对容差 $10^{-8}$；这是数值工程设定，不是论文参数。按收敛测试调整。

### 结果核对（不是拟合许可证）

- 图2(a)：比较曲线排序、平坦区和数量级；图像像素估值需标出估读误差。
- 图2(b)：clean 基线平直；等定位功率时曲线近似平直；不等功率体现频谱重叠结构。DPp 原始功率未知，因此只验机制，不声称逐点再现。
- 图2(c)：原文报告 CPR<30 dB 时 MS-NOMA 误差优于 PRS，MS-NOMA 结果低于1 m；固定 $P_p/N_0=45$ dB-Hz；近似曲线应与精确曲线相符。
- 输出精确/近似的最大相对差，不能只凭肉眼判断“很接近”。

## 12. 公式中必须保留的疑点

### 12.1 BER 的 erfc 自变量

2019快报式（1）和2022长文式（2）、（34）均按出版 PDF 可见形式在 erfc 内**没有平方根**；没有证据支持将其静默改成标准形式。另一方面，论文没有给 $K,\lambda$ 明确数值，且其能量记号为符号能量，因此干净 QPSK 曲线不能唯一从正文复原。

执行要求：同时运行 `paper_literal` 与 `standard_qpsk_candidate`，将 $K,\lambda$ 和噪声定义列为 config。只有依据参考[6]/[7]调制 BER 表、作者 errata 或干净曲线证据作出选择后，才指定主分支；保留另一分支结果供审计。

### 12.2 误差方差到距离单位

式（3）、（5）—（8）由码跟踪时延理论得到，量纲是 $s^2$；图2(c)的纵轴以 m 表示标准差。转换为：

$$
v_r=c_0^2v_\tau,\qquad\sigma_r=c_0\sqrt{v_\tau}.
$$

论文简写中没有把该光速换算明确写在这些公式旁。代码变量名必须写明 `_s2`、`_m2`、`_m`，杜绝重复或遗漏换算。

### 12.3 CPR 的近似与快报/长文差异

快报使用 $\kappa=2G-1$（跨一个 P-Sub 的通信用户数）；长文后续近似用 $2G$。本阶段以快报原值 $2G-1$ 为主，并可额外输出 $2G$ 敏感性，不能混用。

### 12.4 PRS 曲线

快报[10]才是本篇 PRS 模型来源。没有该文精度公式和相同测量参数前，图2(c)中的 PRS 曲线是未闭合依赖。可先完成 MS-NOMA 的 exact/approx 曲线。

### 12.5 ETU 配置

只知道使用 ETU，不足以唯一复现衰落数值。先完成 AWGN，再显式记录最大 Doppler、路径功率、归一化、接收处理和 Monte Carlo 次数；不确定时称“ETU sensitivity run”，不要称 exact reproduction。

### 12.6 式（5）积分方向、模型频移和有限项

式（5）印刷带外积分从 $B_{fe}/2$ 到 $B/2$，对正被积函数给负值，正向物理贡献应反向。优先实现7.1—7.2节的正值信息积分，并保留印刷式诊断结果。通信平坦频带以目标定位子载波中心对称是另一个近似，不是无条件的频率变量代换。式（6）、（8）还含大N、$\kappa\approx2G$、小D等近似，不能把其约1%的有限载波误差归咎于积分器。

## 13. 分阶段执行清单

### 阶段A：配置和波形 sanity check

- [ ] 新建快报独立配置：G=50、M=20、N自动算1049。
- [ ] 调用现有 `generateSignal.m`，检查 Tc/Tp、功率、频谱和有效带宽。
- [ ] 结果文件保存线性值、配置、随机种子和 MATLAB/工具箱版本。

### 阶段B：解析通信结果

- [ ] 建立 W 矩阵和单小区干扰 I(n)。
- [ ] BER paper-literal 和标准 QPSK candidate 分支都可运行。
- [ ] 先复现 flat 图2(a)、图2(b)，明确说明 DPp 示例向量并非原始作者数据。
- [ ] 依据干净基线和参考文献决定汇报分支。

### 阶段C：测距和 DLL

- [ ] 实现式（3）直接频域积分作为 lower-bound 基准。
- [ ] 实现式（5）分段形式并数值核对。
- [ ] 实现式（6）下界近似。
- [ ] 实现式（7）DLL 精确积分与式（8）近似。
- [ ] 以 `c0` 把方差换成米标准差，检查数值积分收敛。

### 阶段D：衰落与论文完整图2

- [ ] 明确 ETU 模型所有参数与接收处理后，复现 fading BER/tracking。
- [ ] 核对/实现快报参考[10]的 PRS 下界；否则如实留空该曲线。
- [ ] 生成统一绘图入口与逐图误差报告。

## 14. 给 AI 的执行任务文本

```text
请在 F:/OneCell_MS_NOMA 工作区内，严格依照 REPRODUCTION_GUIDE.md
复现 2019 单小区 MS-NOMA 快报。先读现有 MS_NOMA 源码，保留原文件，
新建 reproduction/ 模块和独立 run_letter2019.m。

按阶段 A、B、C、D 实施。每阶段完成后运行 MATLAB 测试，保存结果配置、
随机种子、MATLAB/工具箱版本、曲线底层数据和验证日志。优先实现解析
模型；将 ETU 和 PRS 作为可选后续模块。

硬性约束：
1. 快报使用 G=50、M=20、N=1049、κ=2G-1。不要沿用当前主程序的 G=80
   或长文 κ=2G。
2. BER 的论文原式与标准 QPSK 候选公式分支并存；不得无记录地加根号、
   选择 K/λ 或改变 Eb 定义。
3. 测距公式返回时延方差 s²；转成米标准差只能用 c0*sqrt(var_s2)。
4. 完整积分和近似式分开实现并报告相对误差、积分收敛。
5. 图2(b)不等定位功率原始向量没有披露；保存并标明所用示例向量。
6. ETU 参数和 PRS 参考文献[10]未闭合时，继续完成不依赖它们的 AWGN
   与 MS-NOMA 测距结果；不得伪称图2(c)的 PRS 已复现。
7. 不做 CP4A/多基站优化，本任务是无凸优化的单小区快报。

交付：可运行入口、每个公式的独立函数、必要单元测试、图2结果、
README/运行命令、未解决参数清单和逐图复现状态。仅当曲线定义、单位、
参数及数值核验均通过时才标记对应子图为“已复现”。
```

## 15. 来源和范围说明

- 主依据：用户提供的2019年快报 PDF，全文4页；公式（1）—（8）、数值设置和图2逐页核对。
- 项目依据：`F:\OneCell_MS_NOMA\MS_NOMA` 现有 MATLAB 代码与码表。
- 长文只用于确认 CPR 约定后续发生变化；当前复现以2019快报为准。
- 快报参考[8] Betz & Kolodziejski (2009) 支持 early-late DLL 理论；参考[10] Peral-Rosado et al. (2012) 是 PRS 下界外部依赖；参考[6]、[7]与 BER 调制/干扰背景有关。没有取得这些全文并核对前，不把具体系数或 PRS 公式标成已验证原文设置。
- 作者论文图是核对目标，不是原始数据。对曲线逐点比较需要图像数字化并估算读图误差；若无法确认唯一系统参数，应报告“模型已实现、数值尚未闭合”，不要用参数拟合伪装精确复现。

## 16. 公式覆盖与文档自检记录

| 快报公式 | 内容 | 本文实现位置 |
|---|---|---|
| (1) | 通信 BER | 第6节：印刷式、标准候选、横轴校准候选 |
| (2) | 定位到通信的频谱重叠 | 第5节：W矩阵与干扰I |
| (3) | AWGN码时延下界 | 7.1：基带含义与信息积分 |
| (4) | 通信PSD矩形近似 | 7.1：平坦平台及频带指示函数 |
| (5) | 分段下界积分 | 7.1：保留印刷反向积分、给出正向候选 |
| (6) | 下界解析近似 | 7.2：L(F)核、有限项与大N/κ近似 |
| (7) | DLL tracking积分 | 7.3、7.4.1：本地PSD、A0/AN/AC核 |
| (8) | DLL tracking近似 | 7.4、7.4.1：小D极限及有限载波差异 |

2026-10-31，针对指南推导做了独立数值检查：

- 用G=50、M=20、继承30kHz间隔的测试配置，平坦对称PSD模型的下界积分与L(F)闭式核一致，相对误差约 $9\times10^{-16}$。
- 以D=$10^{-5}$ chip验证DLL积分趋近小D有限带闭式式子，相对误差约 $1.3\times10^{-8}$。
- 对D=0.02、0.005、0.001 chip的同一测试配置，DLL方差不低于对应测距下界。
- 原生成器CPR适配的0.04365dB偏移能恢复快报定义下的目标定位功率。

这些检查只验证推导内部一致性，不表示图2的作者原始参数已经确定，也不表示已运行MATLAB复现程序。后续AI仍需实现文档所列模块、测试和逐图结果。
