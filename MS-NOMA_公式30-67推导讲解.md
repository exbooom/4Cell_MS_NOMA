# MS-NOMA 论文公式（30）-（67）推导与意义

本文档面向刚接触通导一体化的读者，解释论文从公式（30）开始到附录公式（67）的来源、推导逻辑和物理意义。论文中的 `C-Sub` 是通信子载波，`P-Sub` 是定位子载波；`P_{km}^p` 表示基站 `k` 给定位用户 `m` 的定位发射功率。

## 0. 先记住总路线

原问题 OP1 要在通信 QoS、总功率和远近效应约束下，最小化所有定位用户的平均定位误差。论文采用两层处理：

1. 对 OP1 做拉格朗日对偶分解，把所有功率变量拆成每个 `(k,m)` 的局部子问题 OP2。
2. 对 OP2 再引入局部乘子，用 KKT 条件得到 `P_{km}^p` 的闭式表达式（33）。

之后的公式（34）-（44）用于单小区性能和资源比较；公式（45）-（50）解释定位误差如何由测距误差和基站几何关系组成；公式（51）-（62）把 DLL 测距误差化简成可用于功率优化的形式；公式（63）-（67）说明闭式功率解来自哪些 KKT 条件。

---

## 1. 公式（30）-（33）：从局部问题到闭式功率

### 1.1 公式（30）-（32）：OP2 的对偶问题

对固定的基站 `k` 和定位用户 `m`，论文把局部问题写为 OP2：

$$
\mathrm{OP2}:\quad
\max_{P_{km}^{p}}
 -\frac{1}{M}(\lambda_{km}\sigma_{\rho}^{km})^2-\nu_kP_{km}^{p}
$$

约束为：

$$
I_{kn}\leq I_{kn}^{\mathrm{th}},\quad n\in N_m
$$

$$
|h_{m\leftarrow km}^{p}|^2P_{km}^{p}
\geq \Gamma\Omega |h_{m\leftarrow k'_{km}m}^{p}|^2P_{k'_{km}m}^{p}
$$

这里 `ν_k` 是基站 `k` 的总定位功率价格。原来的总功率约束已经在第一次分解时通过 `ν_k` 进入目标函数。OP2 仍保留两类局部约束：通信干扰约束和远近效应约束。

再对 OP2 做一次拉格朗日处理，给通信干扰约束引入 `\tilde\mu_{kn}`，给远近效应约束引入 `\tilde\beta_{km}`。局部对偶问题是：

$$
\min \tilde g_{km}(\tilde\mu_{kn},\tilde\beta_{km})
$$

满足：

$$
\tilde\mu_{kn}\geq 0,\quad n\in N_m
$$

$$
\tilde\beta_{km}\geq 0
$$

乘子非负是因为原约束都写成“可行方向大于等于 0”的形式。例如通信约束写成 `I_{kn}^{th}-I_{kn} >= 0`，远近约束写成“目标接收功率 - 要求的竞争功率 >= 0”。

### 1.2 公式（33）：KKT 一阶条件得到闭式解

论文给出的局部最优功率是：

$$
\widetilde P_{km}^{p}
=\lambda_{km}\widetilde\sigma_{\rho}^{km}
\left[
\frac{M}{
\widetilde\beta_{km}|h_{m\leftarrow km}^{p}|^2
-\nu_k
-\sum_{n\in N_m}\widetilde\mu_{kn}J_{kn\leftarrow m}}
\right]^{1/2}
$$

其中：

$$
J_{kn\leftarrow m}
=\sum_{k'\in K}|h_{kn\leftarrow k'm}^{p}|^2T_p
\operatorname{sinc}^2\left(m-\frac{n}{G}\right)
$$

推导的关键只有三步。先把测距误差写成“常数除以功率”的形式：

$$
(\sigma_{\rho}^{km})^2
=\frac{(\widetilde\sigma_{\rho}^{km})^2}{P_{km}^{p}}
$$

于是 OP2 中与 `P_{km}^p` 有关的目标部分变成：

$$
-\frac{1}{M}
\frac{(\lambda_{km}\widetilde\sigma_{\rho}^{km})^2}{P_{km}^{p}}
-\nu_kP_{km}^{p}
$$

通信干扰对功率的导数为：

$$
\frac{\partial I_{kn}}{\partial P_{km}^{p}}
=J_{kn\leftarrow m}
$$

远近效应约束中目标信号功率的导数为：

$$
\frac{\partial}{\partial P_{km}^{p}}
\left(|h_{m\leftarrow km}^{p}|^2P_{km}^{p}
-\Gamma\Omega |h_{m\leftarrow k'_{km}m}^{p}|^2P_{k'_{km}m}^{p}\right)
=|h_{m\leftarrow km}^{p}|^2
$$

因此，对局部拉格朗日函数求导并令其为零：

$$
\frac{1}{M}
\frac{(\lambda_{km}\widetilde\sigma_{\rho}^{km})^2}
{(P_{km}^{p})^2}
-\nu_k
-\sum_{n\in N_m}\widetilde\mu_{kn}J_{kn\leftarrow m}
+\widetilde\beta_{km}|h_{m\leftarrow km}^{p}|^2=0
$$

移项并开平方即可得到公式（33）。

### 1.3 公式（33）的直观意义

公式（33）是“定位收益”和“资源/约束代价”的平衡：

- `λ_{km}` 越大，几何条件越差，需要更大功率；
- `\widetildeσ_{ρ}^{km}` 越大，基础测距越差，需要更大功率；
- `\widetildeβ_{km}` 越大，远近效应约束越紧，倾向于增加该信号功率；
- `ν_k` 越大，基站功率越紧张，分给该定位信号的功率越小；
- `\widetildeμ_{kn}` 越大，通信 QoS 越紧，定位功率越受抑制；
- `J_{kn\leftarrow m}` 越大，该定位功率对通信用户造成的干扰越强，功率越小。

括号内必须为正；实际算法通过更新乘子来保证可行性。

---

## 2. 公式（34）-（38）：单小区性能特例

单小区时没有其他基站，因此不考虑小区间干扰，通信和定位信道增益可省略基站索引。通信 BER 变为：

$$
\mathrm{BER}_n
=\Gamma\operatorname{erfc}
\left(\sqrt{\frac{\gamma P_cT_c}{I_n+2N_0}}\right)
$$

其中 `I_n` 是所有定位子载波对通信子载波 `n` 的总干扰：

$$
I_n
=\sum_{m=1}^{M}P_m^pT_p
\operatorname{sinc}^2\left(m-\frac{n}{G}\right)
$$

单小区测距误差近似为：

$$
(\sigma_{\rho}^{m})^2
\approx \frac{aT_p^2}{2}
\left[
\frac{1}{B_{\mathrm{fe}}T_p(C/N_0)_m}
+\frac{B(CPR)_m}{2B_{\mathrm{fe}}^2}
\right]
$$

其中：

$$
(C/N_0)_m=\frac{P_m^p}{N_0}
$$

$$
(CPR)_m=\frac{2GP_c}{P_m^p}
$$

第一项是噪声造成的误差，第二项是通信信号造成的误差。增大 `P_m^p` 会同时提高定位信号相对噪声的强度、降低通信信号相对定位信号的比例，因此测距误差下降；但它会增加 `I_n`，使通信 BER 变差。

---

## 3. 公式（39）-（44）：资源消耗比较

定义通信能量、定位能量和总能量：

$$
E_{\mathrm{total}}=E_c+E_p
$$

对于传统 PRS，定位信号是间断发送的，单位测量周期内的能量比例约为：

$$
\left(\frac{E_c}{E_p}\right)_{\mathrm{PRS}}
=\frac{N_{\mathrm{RE,com}}}{N_{\mathrm{RE,PRS}}}
\approx \frac{N_{\mathrm{RE,total}}}{N_{\mathrm{RE,PRS}}}
$$

为获得与 MS-NOMA 相同的扩频增益，PRS 所需资源元素数为：

$$
N_{\mathrm{RE,PRS}}=\Delta f_pT_{\mathrm{coh}}
$$

测量周期内总资源元素数为：

$$
N_{\mathrm{RE,total}}=\frac{NT_{\mathrm{meas}}}{T_c}
$$

代入得到：

$$
\left(\frac{E_c}{E_p}\right)_{\mathrm{PRS}}
=\frac{N}{G}\frac{T_{\mathrm{meas}}}{T_{\mathrm{coh}}}
$$

MS-NOMA 不直接占用独立资源块，而是通过增加 BER 来“间接占用”通信资源，所以定义等效资源元素消耗：

$$
\Lambda_{\mathrm{MS-NOMA}}
=\frac{1}{N}\sum_n(\mathrm{BER}_n-\mathrm{BER}_0)
$$

其中 `BER_0` 是没有定位信号时的 BER。PRS 直接占用资源元素，因此：

$$
\Lambda_{\mathrm{PRS}}
=\frac{N_{\mathrm{RE,PRS}}}{N_{\mathrm{RE,total}}}
=\left(\frac{N}{G}\frac{T_{\mathrm{meas}}}{T_{\mathrm{coh}}}\right)^{-1}
$$

关键结论：PRS 的测量越频繁，`T_meas` 越小，直接占用的资源越多；MS-NOMA 连续发送，测量频率不直接增加额外资源元素。

---

## 4. 公式（45）-（50）：测距误差如何变成位置误差

令定位用户 `m` 的各条测距误差组成向量：

$$
\boldsymbol\varepsilon_\rho^m
=\left[\varepsilon_\rho^{1m},\varepsilon_\rho^{2m},\ldots,\varepsilon_\rho^{Km}\right]^T
$$

线性化定位方程后，位置误差为：

$$
\boldsymbol\varepsilon_X^m
=\left((G^m)^TG^m\right)^{-1}(G^m)^T
\boldsymbol\varepsilon_\rho^m
=H^m\boldsymbol\varepsilon_\rho^m
$$

矩阵 `G^m` 的每一行是从基站指向用户的单位方向向量：

$$
\iota_{x}^{km}=\frac{x_p^m-x_b^k}{\|X_b^k-X_p^m\|},\quad
\iota_{y}^{km}=\frac{y_p^m-y_b^k}{\|X_b^k-X_p^m\|},\quad
\iota_{z}^{km}=\frac{z_p^m-z_b^k}{\|X_b^k-X_p^m\|}
$$

如果各基站测距误差独立且无偏，测距误差协方差是对角矩阵：

$$
\Sigma_\rho^m=\operatorname{cov}(\boldsymbol\varepsilon_\rho^m,\boldsymbol\varepsilon_\rho^m)
=\operatorname{diag}\left((\sigma_\rho^{1m})^2,\ldots,(\sigma_\rho^{Km})^2\right)
$$

位置误差协方差为：

$$
\Sigma_X^m=H^m\Sigma_\rho^m(H^m)^T
$$

取水平面两个方向的误差方差并求和，得到：

$$
\Psi_m
=\sqrt{\sum_{k\in K}(\lambda_{km}\sigma_\rho^{km})^2}
$$

其中：

$$
\lambda_{km}=\sqrt{(\eta_{1k}^m)^2+(\eta_{2k}^m)^2}
$$

`λ` 只由基站与用户的几何位置决定，描述几何放大；`σ_ρ` 描述该条链路的测距质量。

---

## 5. 公式（51）-（55）：DLL 测距误差的积分形式

DLL 的 Early-Late 鉴别器通过比较提前和延迟相关输出估计码延迟。定义四个积分：

$$
A_0^m=\int_{B_0-B_{\mathrm{fe}}/2}^{B_0+B_{\mathrm{fe}}/2}
fG_m^p(f+m\Delta f_p)\sin(\pi fDT_p)\,df
$$

$$
A_1^m=\int_{B_0-B_{\mathrm{fe}}/2}^{B_0+B_{\mathrm{fe}}/2}
N_0G_m^p(f+m\Delta f_p)\sin^2(\pi fDT_p)\,df
$$

$$
A_2^m=\int_{B_0-B_{\mathrm{fe}}/2}^{B_0+B_{\mathrm{fe}}/2}
G_s^m(f+m\Delta f_p)G_m^p(f+m\Delta f_p)\sin^2(\pi fDT_p)\,df
$$

$$
A_3^{km}=\int_{B_0-B_{\mathrm{fe}}/2}^{B_0+B_{\mathrm{fe}}/2}
G_q^{km}(f+m\Delta f_p)G_m^p(f+m\Delta f_p)\sin^2(\pi fDT_p)\,df
$$

四项的意义分别是：有用相关斜率、噪声、通信干扰、其他定位信号干扰。综合得到：

$$
(\sigma_\rho^{km})^2
=\frac{a(A_1^m+A_2^m+A_3^{km})}
{(2\pi)^2|h_{m\leftarrow km}^p|^2P_{km}^p(A_0^m)^2}
$$

分母是“有用定位信号强度 × 相关峰斜率的平方”，越大越容易跟踪；分子是噪声和干扰，越大误差越大。

---

## 6. 公式（56）-（62）：从积分到可优化的近似式

当前端带宽足够大、Early-Late 间隔 `D` 较小时，使用 sinc 函数和泰勒展开近似，得到：

$$
A_0^m\approx \frac{1}{2\pi DB_{\mathrm{fe}}}
$$

$$
A_1^m=\pi DT_pN_0A_0^m
$$

通信干扰项近似为：

$$
A_2^m\approx D^2T_pP_c
\sum_{k'\in K}\sum_{n\in N}
|h_{m\leftarrow k'n}^c|^2
\sin^2\left(\frac{n\pi}{G}\right)
$$

其他定位信号干扰项近似为：

$$
A_3^{km}\approx \pi^2D^2T_p^4
\sum_{k'\in K_k}|h_{m\leftarrow k'm}^p|^2P_{k'm}^p\,\overline A_3
$$

其中：

$$
\overline A_3\approx \frac{1}{2\pi^2T_p^3}
$$

把这些结果代回公式（55），得到：

$$
(\sigma_\rho^{km})^2
\approx \frac{aT_p^2}{2}
\left(
\frac{N_0}{B_{\mathrm{fe}}|h_{m\leftarrow km}^p|^2P_{km}^p}
+\frac{BG P_c\sum_{k'\in K}|h_{m\leftarrow k'}^c|^2}
{B_{\mathrm{fe}}^2|h_{m\leftarrow km}^p|^2P_{km}^p}
+\frac{\sum_{k'\in K_k}|h_{m\leftarrow k'm}^p|^2P_{k'm}^p}
{B_{\mathrm{fe}}^2T_p|h_{m\leftarrow km}^p|^2P_{km}^p}
\right)
$$

定义等效通信信道增益：

$$
|h_{m\leftarrow k'}^c|^2
=\frac{2}{N}\sum_{n\in N}|h_{m\leftarrow k'n}^c|^2
\sin^2\left(\frac{n\pi}{G}\right)
$$

这样就把许多通信子载波的影响压缩成一个等效量。公式（61）的三项清晰对应：环境噪声、通信干扰、其他定位信号干扰。

如果 MS-NOMA 的定位信号足够弱，定位信号之间的第三项可以忽略，就得到正文公式（10）：

$$
(\sigma_\rho^{km})^2
\approx \frac{aT_p^2}{2}
\left[
\frac{1}{B_{\mathrm{fe}}T_p(C/N_0)_{km}}
+\frac{B(CPR)_{km\leftarrow k'}}{2B_{\mathrm{fe}}^2}
\right]
$$

其中：

$$
(C/N_0)_{km}=\frac{|h_{m\leftarrow km}^p|^2P_{km}^p}{N_0}
$$

$$
(CPR)_{km\leftarrow k'}
=\frac{2G|h_{m\leftarrow k'}^c|^2P_c}
{|h_{m\leftarrow km}^p|^2P_{km}^p}
$$

再定义：

$$
(PPR)_{km\leftarrow k'm}
=\frac{|h_{m\leftarrow k'm}^p|^2P_{k'm}^p}
{|h_{m\leftarrow km}^p|^2P_{km}^p}
$$

`C/N0` 是目标定位信号相对噪声的强度，`CPR` 是通信信号相对目标定位信号的强度，`PPR` 是其他定位信号相对目标定位信号的强度。

---

## 7. 公式（63）-（67）：KKT 条件如何回到公式（33）

局部问题的 KKT 条件包括互补松弛：

$$
\sum_{n\in N_m}\widetilde\mu_{kn}
\left(I_{kn}^{\mathrm{th}}-I_{kn}\right)=0
$$

$$
\widetilde\beta_{km}
\left(|h_{m\leftarrow km}^p|^2P_{km}^p
-\Gamma\Omega|h_{m\leftarrow k'_{km}m}^p|^2P_{k'_{km}m}^p\right)=0
$$

以及对局部拉格朗日函数的一阶驻点条件：

$$
\frac{\partial \widetilde L}{\partial P_{km}^p}=0
$$

展开后：

$$
-\frac{1}{M}\frac{\partial(\lambda_{km}\sigma_\rho^{km})^2}{\partial P_{km}^p}
-\nu_k
-\sum_{n\in N_m}\widetilde\mu_{kn}J_{kn\leftarrow m}
+\widetilde\beta_{km}|h_{m\leftarrow km}^p|^2=0
$$

根据测距误差与功率的反比关系：

$$
(\lambda_{km}\sigma_\rho^{km})^2
=\frac{(\lambda_{km}\widetilde\sigma_\rho^{km})^2}{P_{km}^p}
$$

所以：

$$
-\frac{1}{M}\frac{\partial(\lambda_{km}\sigma_\rho^{km})^2}{\partial P_{km}^p}
=\frac{(\lambda_{km}\widetilde\sigma_\rho^{km})^2}
{M(P_{km}^p)^2}
$$

干扰对功率的导数为：

$$
J_{kn\leftarrow m}
=\sum_{k'\in K}|h_{kn\leftarrow k'm}^p|^2T_p
\operatorname{sinc}^2\left(m-\frac{n}{G}\right)
$$

将这些代入驻点条件并令其为零，整理、开平方，就重新得到公式（33）。因此公式（33）并不是猜出来的，而是局部拉格朗日函数的 KKT 一阶条件的显式解。

---

## 8. 公式（30）以后整段的核心意义

这部分可以浓缩为四句话：

1. 公式（30）-（32）把每个局部功率分配问题写成对偶问题。
2. 公式（33）用 KKT 条件直接给出一个定位子载波应分多少功率。
3. 公式（34）-（44）说明 MS-NOMA 的 BER、测距误差和资源消耗如何随功率、带宽、测量频率变化。
4. 公式（45）-（67）证明“测距误差 + 几何关系”如何构成位置误差，并把 DLL 的复杂积分化成可用于优化的功率倒数形式。

算法实际运行时并不是一次就知道最优的 `μ`、`ν`、`β`。它先给乘子初值，根据公式（33）算功率，再看约束是否紧张，用次梯度更新乘子，循环直到功率和约束收敛。

