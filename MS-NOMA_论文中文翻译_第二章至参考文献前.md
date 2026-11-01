# 面向未来通信-定位一体化系统的多尺度 NOMA 设计与性能分析

## 第二章 系统模型

考虑图 2 所示的典型空间定位场景，其中包含 `K` 个基站（BS）。通信用户（C-User）和定位用户（P-User）可以位于不同位置。定义 `h_{k_y\leftarrow k_x}` 为子载波 `k_x` 被用户 `k_y` 接收时的瞬时信道增益。典型定位场景下 MS-NOMA 波形的瞬时信道增益见表 I。

每个通信用户只连接一个基站，而每个定位用户连接多个基站，因此通信用户使用索引 `kn`，定位用户使用索引 `m`。

假设各基站完全同步，定位用户使用基于时间的算法估计位置。定位用户 `m` 的水平定位精度可以表示为：

$$
\Psi_m=\sqrt{\sum_{k\in K}(\lambda_{km}\sigma_{\rho}^{km})^2}
$$

其中，`lambda_{km}` 表示几何衰减因子，`sigma_{rho}^{km}` 表示测距误差方差（文中符号写作测距误差标准差相关量）。公式（1）同时也是定位误差的下界。其详细推导见附录 A。

## 第三章 MS-NOMA 波形的特性

### A. 定位信号对通信的干扰

为了更准确地评估定位信号对通信的干扰，假设通信分量之间的小区间干扰可以被理想消除。通信子载波 `kn` 的误码率为：

$$
\mathrm{BER}_{kn}=\Gamma\operatorname{erfc}\left(\sqrt{\frac{\gamma |h_{kn\leftarrow kn}^{c}|^2P_cT_c}{I_{kn}+2N_0}}\right)
$$

其中，`Gamma` 和 `gamma` 由调制编码方案决定；`P_c` 是通信分量功率；`T_c` 是通信符号周期；`N_0` 是环境噪声的单边功率谱密度；`I_{kn}` 是定位分量对通信用户 `kn` 造成的干扰。

由于不同基站发出的定位分量可能具有相近功率，邻近基站定位分量造成的干扰不能忽略。因此：

$$
I_{kn}=\sum_{k'\in K}\sum_{m\in M}\overline P_{kn\leftarrow k'm}^{p}
$$

其中，`\overline P_{kn\leftarrow k'm}^{p}` 是基站 `k'` 的定位子载波 `k'm` 到达通信用户 `kn` 后的功率：

$$
\overline P_{kn\leftarrow k'm}^{p}=|h_{kn\leftarrow k'm}^{p}|^2P_{k'm}^{p}G_m^p(n\Delta f_c)
$$

定位子载波 `m` 的归一化功率谱密度为：

$$
G_m^p(f)=T_p\operatorname{sinc}^2[(f-m\Delta f_p)T_p]
$$

将 `Delta f_p=G Delta f_c` 代入，可写成：

$$
\overline P_{kn\leftarrow k'm}^{p}=|h_{kn\leftarrow k'm}^{p}|^2P_{k'm}^{p}T_p\operatorname{sinc}^2\left(m-\frac nG\right)
$$

这说明定位信号对通信信号的影响由三部分决定：传播信道增益、定位发射功率，以及两个子载波频谱之间的重叠程度。

### B. MS-NOMA 波形的测距精度

接收机可以使用延迟锁定环（DLL）跟踪定位分量。以相干早迟鉴别器为例，定位子载波 `km` 的跟踪/测距精度可以由公式（6）表示。其物理结构是：噪声和干扰位于分子，有用定位信号强度和相关峰斜率位于分母。

其中，`a=B_L(1-0.5B_LT_{coh})`，由环路噪声带宽 `B_L` 和相干积分时间 `T_{coh}` 决定；`B_0` 是 MS-NOMA 波形中心频率；`B_{fe}` 是双边前端带宽；`D` 是 DLL 的早迟间隔。通信分量接收功率谱密度为：

$$
G_s^m(f)=\sum_{k'\in K}\sum_{n\in N}|h_{m\leftarrow k'n}^{c}|^2P_cG_n^c(f)
$$

其中：

$$
G_n^c(f)=T_c\operatorname{sinc}^2[(f-n\Delta f_c)T_c]
$$

其他基站定位分量的功率谱密度为：

$$
G_q^{km}(f)=\sum_{k'\in K_k}|h_{m\leftarrow k'm}^{p}|^2P_{k'm}^{p}G_m^p(f)
$$

将噪声、通信分量和其他定位分量代入 DLL 误差表达式，并进行近似，可得到：

$$
(\sigma_\rho^{km})^2\approx \frac{aT_p^2}{2}\left[\frac{1}{B_{fe}T_p(C/N_0)_{km}}+\frac{B(CPR)_{km\leftarrow k'}}{2B_{fe}^2}+\frac{(PPR)_{km\leftarrow k'm}}{2B_{fe}^2}\right]
$$

三项分别由环境噪声、所有基站通信分量、其他基站定位分量造成。由于 MS-NOMA 中定位分量被设计得远弱于通信分量，定位分量之间的干扰可以忽略，于是：

$$
(\sigma_\rho^{km})^2\approx \frac{aT_p^2}{2}\left[\frac{1}{B_{fe}T_p(C/N_0)_{km}}+\frac{B(CPR)_{km\leftarrow k'}}{2B_{fe}^2}\right]
$$

定义测距因子：

$$
(\widetilde\sigma_\rho^{km})^2=(\sigma_\rho^{km})^2P_{km}^{p}
$$

该量与定位发射功率无关或近似变化较小，便于后续功率优化。

## 第四章 MS-NOMA 波形的功率分配

### A. 约束条件

#### 1）QoS 约束下的 BER 阈值

为了保证通信用户的服务质量，所有通信用户的 BER 都必须低于阈值：

$$
\mathrm{BER}_{kn}\leq \Xi_{th},\quad \forall k\in K,\ \forall n\in N
$$

将 BER 表达式代入并整理，可以得到通信干扰阈值：

$$
I_{kn}\leq I_{kn}^{th}
=\frac{\gamma |h_{kn\leftarrow kn}^{c}|^2P_cT_c}{\operatorname{erfc}^{-1}(\Xi_{th}/\Gamma)}-2N_0
$$

#### 2）总功率限制

每个基站的总发射功率受到限制：

$$
\sum_{m\in M}P_{km}^{p}+NP_c\leq P_k^T
$$

定义基站 `k` 的定位功率预算：

$$
P_k^{th}=P_k^T-NP_c
$$

则：

$$
\sum_{m\in M}P_{km}^{p}\leq P_k^{th}
$$

#### 3）消除远近效应

为了让定位用户尽可能接收到多个基站的定位信号，不同基站的接收定位功率必须满足：

$$
\frac{|h_{m\leftarrow km}^{p}|^2P_{km}^{p}}
{|h_{m\leftarrow k'm}^{p}|^2P_{k'm}^{p}}
\geq \Gamma\Omega,
\quad \forall m,\ k,\ k'\neq k
$$

其中，`Omega` 是定位伪随机码自相关与互相关之比，`Gamma` 是由接收机性能决定的裕量，通常大于 1。对于特定的定位子载波 `km`，只要相对于最强的其他定位子载波满足约束，就能保证相对于所有更弱信号也满足约束。因此约束可改写为：

$$
|h_{m\leftarrow km}^{p}|^2P_{km}^{p}
\geq \Gamma\Omega|h_{m\leftarrow k'_{km}m}^{p}|^2P_{k'_{km}m}^{p}
$$

### B. 联合功率分配模型

目标是在 QoS、总功率和可听性约束下，使网络中所有定位用户的平均水平定位误差下界最小：

$$
\mathrm{OP1}:\quad \max_{\{P_{km}^{p}\}}-\frac1M\sum_{m\in M}(\Psi_m)^2
$$

满足：

$$
I_{kn}\leq I_{kn}^{th}
$$

$$
\sum_{m\in M}P_{km}^{p}\leq P_k^{th}
$$

$$
|h_{m\leftarrow km}^{p}|^2P_{km}^{p}
\geq \Gamma\Omega|h_{m\leftarrow k'_{km}m}^{p}|^2P_{k'_{km}m}^{p}
$$

最大化负误差与最小化误差完全等价。该问题可用拉格朗日对偶方法求解。

### C. 联合功率分配算法

对通信干扰约束、基站功率约束和远近效应约束分别引入非负对偶变量 `mu_{kn}`、`nu_k` 和 `beta_{km}`。拉格朗日对偶函数为：

$$
g(\mu,\nu,\beta)=\max_{\{P_{km}^{p}\}}L(\{P_{km}^{p}\},\mu,\nu,\beta)
$$

对应的对偶问题为：

$$
\min g(\mu,\nu,\beta)
$$

$$
\mu\geq0,\quad \nu\geq0,\quad \beta\geq0
$$

定义与定位子载波 `m` 相关的通信子载波集合：

$$
N_m=\{(2G-1)(m-1)+1,\ldots,(2G-1)m\}
$$

利用 `sum_{n in N}=sum_m sum_{n in N_m}`，对偶函数可以分解为：

$$
g(\mu,\nu,\beta)=\sum_{k\in K}g_k(\mu,\nu,\beta)
$$

$$
g_k(\mu,\nu,\beta)=\sum_{m\in M}g_{km}(\mu,\nu,\beta)+\nu_kP_k^{th}
$$

其中，每个局部函数只与一个 `P_{km}^{p}` 相关。单个局部问题为：

$$
\mathrm{OP2}:\quad \max_{P_{km}^{p}}
-\frac1M(\lambda_{km}\sigma_\rho^{km})^2-\nu_kP_{km}^{p}
$$

满足：

$$
I_{kn}\leq I_{kn}^{th},\quad n\in N_m
$$

$$
|h_{m\leftarrow km}^{p}|^2P_{km}^{p}
\geq \Gamma\Omega|h_{m\leftarrow k'_{km}m}^{p}|^2P_{k'_{km}m}^{p}
$$

局部对偶变量满足：

$$
\widetilde\mu_{kn}\geq0,\quad \widetilde\beta_{km}\geq0
$$

用 KKT 条件求解，可得闭式功率：

$$
\widetilde P_{km}^{p}=\lambda_{km}\widetilde\sigma_\rho^{km}
\left[\frac{M}{\widetilde\beta_{km}|h_{m\leftarrow km}^{p}|^2-\nu_k-\sum_{n\in N_m}\widetilde\mu_{kn}J_{kn\leftarrow m}}\right]^{-1/2}
$$

其中：

$$
J_{kn\leftarrow m}=\sum_{k'\in K}|h_{kn\leftarrow k'm}^{p}|^2T_p\operatorname{sinc}^2\left(m-\frac nG\right)
$$

定位用户的位置未知，而几何衰减因子 `lambda_{km}` 需要用户位置。因此第一次迭代时将 `lambda_{km}` 设置为 1，得到近似位置和初始功率；随后交替更新位置、功率和对偶变量，最终收敛到较优解。算法 CP4A 的过程是：初始化 `nu_k`；循环更新 `mu_{kn}`、`beta_{km}`；依据闭式解计算功率；用次梯度更新乘子；直到局部乘子和基站功率乘子变化小于 `epsilon`。

## 第五章 性能评估

首先在单小区网络中评估通信性能、资源消耗和测距精度，然后在四基站网络中综合考察定位性能。

仿真采用：通信分量 QPSK，定位分量 BPSK；载波频率 3.5 GHz；通信子载波间隔 30 kHz；总带宽分别为 20 MHz 和 50 MHz；定位用户数 `M=20`；前端带宽为总带宽的两倍；环路噪声带宽 0.2 Hz；相干积分时间 0.02 s；早迟间隔 0.02 chips。

### A. 通信性能

当所有定位用户的定位功率相等时，`Eb/N0` 和 `C/N0` 增大，平均 BER 下降。当 CPR 较小时，随着 `Eb/N0` 增大，BER 曲线会变平，因为此时定位干扰超过环境噪声，成为主要限制。定位分量越弱，BER 越接近只有噪声时的曲线。没有定位分量时，BER 最低。

### B. 定位资源消耗

MS-NOMA 连续发送定位波形，其单位时间通信能量和定位能量分别近似等于通信功率总和与定位功率总和。PRS 的定位功率可以接近通信功率，但由于其间断发送，能量由每个测量周期内占用的资源元素决定。

仿真显示，MS-NOMA 在大多数情况下消耗的等效资源元素少于 PRS。带宽越大，定位波形的符号周期越短，定位干扰越小，MS-NOMA 的等效资源消耗越少。PRS 的资源消耗随着测量周期缩短而显著增加；MS-NOMA 由于连续发送，与测量周期无关，因此更适合自动驾驶等高频定位场景。

在相同总能量和带宽下，MS-NOMA 的测距精度始终优于 PRS。MS-NOMA 中定位能量占比越低，通信干扰相对越大，测距精度会逐渐下降；PRS 中通信与定位采用时分方式，功率不变时测距精度基本不受能量比例影响。

### C. 定位性能

定位功率越强，测量精度通常越高，但其上限受通信 QoS 约束。因此必须通过 CP4A 合理分配功率。仿真采用四个基站，位置为 `(0,0)`、`(200,0)`、`(200,200)`、`(0,200)`，20 个定位用户随机分布，采用自由空间传播模型，每次仿真进行 50 次蒙特卡洛试验。通信总功率 `NP_c=1 W`，每个基站定位功率预算 `P_k^{th}=0.1 W`，BER 阈值为 `10^{-2}`，远近效应裕量 `Gamma=2`。

CP4A 能够显著降低远近效应。等功率传输时，约 50.8% 的定位用户无法同时接收到足够多的基站信号；CP4A 下所有用户都有定位结果。信道条件较差的用户通常获得更大的定位功率；位于区域边缘、几何条件较差的用户，即使信道增益相近，也可能获得更多功率。

与 PRS 相比，MS-NOMA 的平均定位精度提高约 93.0%（20 MHz）和 93.6%（50 MHz）。PRS 在覆盖边缘的误差更大，因为路径损耗和几何衰减同时恶化；CP4A 联合考虑测距误差和几何衰减，使 MS-NOMA 的定位误差分布更平坦。

## 第六章 结论

本文研究了面向下一代通信-定位一体化系统的多尺度非正交多址（MS-NOMA）技术。MS-NOMA 将功率可配置的定位分量叠加到通信分量上，从而在较少资源消耗下获得较高测距精度。

本文分析了通信分量和定位分量同时存在时的通信 BER 与定位测距精度。结果表明，两者之间的相互影响较小，MS-NOMA 相比传统 PRS 能显著提高测距精度。资源和能量分析表明，MS-NOMA 尤其适合需要高频定位的场景，因为它连续发送定位波形，不需要周期性占用大量 PRS 资源。

此外，本文将多用户定位功率分配建模为凸优化问题，在通信 QoS、功率预算和远近效应约束下最小化定位误差。提出的 CP4A 算法能够显著提升定位精度和定位覆盖率。

## 附录 A：水平定位精度的推导

定义定位用户 `m` 的测距误差向量：

$$
\boldsymbol\varepsilon_\rho^m=[\varepsilon_\rho^{1m},\varepsilon_\rho^{2m},\ldots,\varepsilon_\rho^{Km}]^T
$$

线性化后，位置误差为：

$$
\boldsymbol\varepsilon_X^m=((G^m)^TG^m)^{-1}(G^m)^T\boldsymbol\varepsilon_\rho^m=H^m\boldsymbol\varepsilon_\rho^m
$$

几何矩阵为：

$$
G^m=\begin{bmatrix}\iota_x^{1m}&\iota_y^{1m}&\iota_z^{1m}\\\iota_x^{2m}&\iota_y^{2m}&\iota_z^{2m}\\\vdots&\vdots&\vdots\\\iota_x^{Km}&\iota_y^{Km}&\iota_z^{Km}\end{bmatrix}
$$

方向余弦为：

$$
\iota_x^{km}=\frac{x_p^m-x_b^k}{\|X_b^k-X_p^m\|},\quad
\iota_y^{km}=\frac{y_p^m-y_b^k}{\|X_b^k-X_p^m\|},\quad
\iota_z^{km}=\frac{z_p^m-z_b^k}{\|X_b^k-X_p^m\|}
$$

假设测距误差相互独立且无偏，则协方差矩阵为：

$$
\Sigma_\rho^m=\operatorname{diag}((\sigma_\rho^{1m})^2,\ldots,(\sigma_\rho^{Km})^2)
$$

位置误差协方差为：

$$
\Sigma_X^m=H^m\Sigma_\rho^m(H^m)^T
$$

取水平面两个方向的误差，得到：

$$
\Psi_m=\sqrt{\sum_{k\in K}\left[\left((\eta_{1k}^m)^2+(\eta_{2k}^m)^2\right)(\sigma_\rho^{km})^2\right]}
$$

定义：

$$
\lambda_{km}=\sqrt{(\eta_{1k}^m)^2+(\eta_{2k}^m)^2}
$$

即可得到正文公式（1）。

## 附录 B：测距误差方差的推导

定义 DLL 相关积分：

$$
A_0^m=\int_{B_0-B_{fe}/2}^{B_0+B_{fe}/2}fG_m^p(f+m\Delta f_p)\sin(\pi fDT_p)df
$$

$$
A_1^m=\int_{B_0-B_{fe}/2}^{B_0+B_{fe}/2}N_0G_m^p(f+m\Delta f_p)\sin^2(\pi fDT_p)df
$$

$$
A_2^m=\int_{B_0-B_{fe}/2}^{B_0+B_{fe}/2}G_s^m(f+m\Delta f_p)G_m^p(f+m\Delta f_p)\sin^2(\pi fDT_p)df
$$

$$
A_3^{km}=\int_{B_0-B_{fe}/2}^{B_0+B_{fe}/2}G_q^{km}(f+m\Delta f_p)G_m^p(f+m\Delta f_p)\sin^2(\pi fDT_p)df
$$

DLL 测距误差方差为：

$$
(\sigma_\rho^{km})^2=\frac{a(A_1^m+A_2^m+A_3^{km})}{(2\pi)^2|h_{m\leftarrow km}^p|^2P_{km}^p(A_0^m)^2}
$$

由于定位波形的带宽、前端带宽和早迟间隔满足相应近似条件，使用泰勒展开和 sinc 函数积分，可得：

$$
A_0^m\approx\frac{1}{2\pi DB_{fe}}
$$

$$
A_1^m=\pi DT_pN_0A_0^m
$$

$$
A_2^m\approx D^2T_pP_c\sum_{k'\in K}\sum_{n\in N}|h_{m\leftarrow k'n}^c|^2\sin^2\left(\frac{n\pi}{G}\right)
$$

$$
A_3^{km}\approx \pi^2D^2T_p^4\sum_{k'\in K_k}|h_{m\leftarrow k'm}^p|^2P_{k'm}^p\overline A_3
$$

其中：

$$
\overline A_3\approx\frac{1}{2\pi^2T_p^3}
$$

代回即可得到正文公式（9）和（10）。定义等效通信信道增益：

$$
|h_{m\leftarrow k'}^c|^2=\frac{2}{N}\sum_{n\in N}|h_{m\leftarrow k'n}^c|^2\sin^2\left(\frac{n\pi}{G}\right)
$$

## 附录 C：闭式功率解的推导

局部问题的 KKT 互补松弛条件为：

$$
\sum_{n\in N_m}\widetilde\mu_{kn}(I_{kn}^{th}-I_{kn})=0
$$

$$
\widetilde\beta_{km}\left(|h_{m\leftarrow km}^p|^2P_{km}^p-\Gamma\Omega|h_{m\leftarrow k'_{km}m}^p|^2P_{k'_{km}m}^p\right)=0
$$

驻点条件为：

$$
\frac{\partial\widetilde L(\{P_{km}^p\},\widetilde\mu,\widetilde\beta)}{\partial P_{km}^p}=0
$$

利用：

$$
(\lambda_{km}\sigma_\rho^{km})^2=\frac{(\lambda_{km}\widetilde\sigma_\rho^{km})^2}{P_{km}^p}
$$

以及：

$$
J_{kn\leftarrow m}=\sum_{k'\in K}|h_{kn\leftarrow k'm}^p|^2T_p\operatorname{sinc}^2\left(m-\frac nG\right)
$$

将一阶条件整理并令其为零，就得到公式（33）的闭式解。这个闭式解随后被 CP4A 算法反复调用，乘子则用次梯度方法更新。
