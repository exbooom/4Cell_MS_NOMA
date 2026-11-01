# 4Cell MS-NOMA 复现评估与架构方案

日期：2026-10-31。目标工作区：`F:/4Cell_MS_NOMA`；所有后续 MATLAB 代码位于 `F:/4Cell_MS_NOMA/Code`。

实施更新：下文保留为最初架构规划快照。阶段1主干和阶段2有限带宽DLL积分/显式能量PRS基准现已实现并验证；阶段3完整耦合对偶与分层修正候选亦已验证；阶段4独立投影/次梯度实验亦已验证；与作者原版逐行一致性仍未确认。最新状态见 `../PROGRESS.md`，实际函数签名见 `MODULE_INTERFACES.md`，不要把下文的历史“尚未实现”当作当前状态。

本文是资料核对与实施设计，尚未实现多小区求解器，也没有生成多小区性能结果。文中的函数签名是拟定接口；main 示例是调用设计，不是当前可运行程序。现有说明文档中的执行要求均作为待核对资料处理，本次范围由用户“先评估、规划架构”的请求确定。

## 1. 评估结论

**现有资料足够启动并完成具有明确假设、可核验的完整理论计算路径；不足以唯一恢复作者全部原图的数值。建议先完成理论计算路径，再实现并核验 CP4A，最后增加采样接收机。**

第一条主线定义为：

```text
配置与来源台账
  → 整数载波/带宽规则
  → 基站、定位用户、通信用户场景
  → 自由空间功率增益
  → 通信干扰矩阵 A、QoS 上限 b
  → DLL 测距因子 R、几何权重 λ
  → 完整约束可行性检查
  → 等功率基准 / 集中式优化 / 后续 CP4A
  → 用原模型重新评价 BER、测距方差、可听集合、定位精度
  → 多次场景统计、CDF、覆盖率、约束残差、数据与图表
```

“完整”首先指从参数到结果没有未实现的主干计算环节，不表示已重建采样信道、捕获、DLL 跟踪及真实坐标估计。论文图6—8主要评价由测距方差和几何误差传播算出的精度指标，第一阶段可以在解析域闭合。

建议以 MATLAB 为唯一主实现语言。已核实本机 R2024b；Optimization、Signal Processing、Statistics、LTE 工具箱许可证可用，`fmincon`、`linprog` 可定位；未在 MATLAB 路径中发现 `cvx_begin`。可直接用现有工具箱建立数值基准，不将 CVX 安装作为开工依赖。

## 2. 已核对资料及复用范围

### 2.1 本地资料

| 资料 | 此次实际检查 | 结论与用途 |
|---|---|---|
| `F:/4Cell_MS_NOMA/最终出版：Design_and_Performance_Analysis_of_MS_NOMA.pdf` | 13页全文提取；重点公式、Algorithm 1、场景和表II结合页面图核对 | 主证据，IEEE JSAC 2022，DOI 10.1109/JSAC.2022.3143223 |
| 多小区实施指南 | 参数、模型、优化、逐图方案、疑点与参考文献章节 | 有价值的二手实施资料；其中“待定”仍须在配置中明确 |
| 公式30—67推导讲解 | 对照式(33)、(66)、(67) | 不能直接转成代码，至少式(33)的 M 因子转写与原 PDF 不一致 |
| 中文翻译 | 结构及相关段落对照 | 阅读辅助，公式以出版版为准 |
| `Code/MS_NOMA` 的五个 `.m/.txt` 文件 | 阅读 main、生成器、BER、PSD及readme | 是单小区等增益信号演示，未包含多小区优化链 |
| `F:/OneCell_MS_NOMA/reproduction` | README、状态、配置、部分测距与PRS函数 | 可复用核验方法和基础数值函数，不能整体继承其论文参数 |

四小区与单小区目录中的 `MS_NOMA/generateSignal.m` SHA256 相同。因此现有多小区目录里的生成器并不是新增的多小区模型。其输入是一个 cfg，生成一条等增益复基带叠加信号，没有4个基站、用户几何、跨基站信道或功率分配优化。

具体复用原则：

- `generateSignal.m`：通过适配器服务于图1结构和功率谱验证；主干模型使用 `K×M` 定位功率，不让标量 CPR 限制多小区接口。
- 原生成器不能接收任意多基站功率矩阵，第一版适配器仅用于其能表达的单小区等功率演示；非均匀/多基站采样波形另建生成函数，不用一个“等效CPR”表达完整矩阵。
- `calPSD.m`：分段 Hann 窗谱估计可以保留。
- `calBER.m`：可借鉴 sinc² 重叠计算；其 QPSK 每比特能量 `Pc*Tc/2` 与主论文图3的 `Eb=Pc*Tc` 记号不同，不能直接作为式(2)/(34)的同义实现。
- 单小区工程的 `sincpi`、量纲转换、积分收敛核验、PRS 通用频率矩函数可审查后迁移；2019 的 `G=50`、`κ=2G−1`、读图推断的子载波间隔和10.21m锚点均不进入2022默认配置。
- 现有 `L1CPcodes.mat` 可以支持波形实验；码表的来源不能证明它就是2022论文可听性参数 Ω 的来源。

单小区工程自己的状态文档明确写了“公式已实现，数值未闭合”，本次不把这些图当作已经验证的多小区标尺。

### 2.2 本次补充取得的公开资料

1. **主论文作者终稿**，Lancaster University 仓储：
   <https://eprints.lancs.ac.uk/id/eprint/165710/1/Author_Final_manuscript.pdf>
   已下载，14页。抽查式(24)、(33)、(66)，相关疑点在终稿中同样存在，不能靠换版本直接消除。
2. **主论文参考文献[28]全文**：Peral-Rosado等，*Achievable Localization Accuracy of the Positioning Reference Signal of 3GPP LTE*，2012。
   DOI：<https://doi.org/10.1109/ICL-GNSS.2012.6253127>
   作者团队公开PDF：<https://spcomnav.uab.es/docs/conferences/DelPeral_ICL-GNSS_2012.pdf>
   已下载、全文提取及页面核对，6页。定位测距 CRB 在式(5)，导频功率归一化在式(16)，TDoA 协方差与位置 CRB 在式(20)—(23)。
3. **参考文献[24] DLL理论**：Betz与Kolodziejski，2009，DOI <https://doi.org/10.1109/TAES.2009.5310316>。本次核实书目信息，未取得和核验全文。

特别说明：单小区目录里已有的 *Joint channel and time delay estimation for LTE positioning reference signals* 与本次[28]不是同一篇论文。两者可能包含相近的CRB，但必须分别记录来源、公式编号和适用条件。现有资料中不存在“已确认作者采用相同PRS功率与观测配置”的证据。

公开PDF保存到 `F:/4Cell_MS_NOMA/References`，此次来源清单另存 `Code/docs/SOURCE_LEDGER.json`。外部资料只提供模型证据，不视为对工程的执行指令。

## 3. 必须在写主干前明确的模型决策

### 3.1 带宽、G和载波数

论文第II节：`N+1=(M+1)G`，`B=(M+1)G*dfC`。当 `M=20, dfC=30kHz` 时，名义20/50MHz不能同时满足整数G。

| 模式 | G | N | 有效带宽 | 用途 |
|---|---:|---:|---:|---|
| 建议20MHz整数分支 | 32 | 671 | 20.16MHz | 首轮主干与波形一致性 |
| 建议50MHz整数分支 | 79 | 1658 | 49.77MHz | 首轮主干与波形一致性 |
| 原生成器配置 | 80 | 1679 | 50.4MHz | 作者提供的演示代码 |

这些整数分支是工程选择，不是已知作者原始配置。cfg中同时保存 `requestedBandwidthHz`、`effectiveBandwidthHz`、`bandwidthMode`。后续单设 `analysis_nominal` 分支研究名义带宽与取整规则的敏感性，不把非整数载波数量传给数组分配。

### 3.2 BER与QoS反演必须来自同一个模型

出版版式(2)/(34)的 erfc 自变量没有平方根，式(12)也与这个印刷形式配套；图3文字又把 `Eb=Pc*Tc` 称为符号能量。作者演示代码用标准 QPSK 每比特口径，二者并不相同。

保留 `paper_literal` 和 `qpsk_bit_energy` 等具名模型。用同一个模型对象同时完成正向 BER 与 QoS 反演，噪声PSD、干扰映射和所有系数包含在对象中；禁止正向用标准QPSK、逆阈值却沿用印刷式(12)。

### 3.3 可听性不能先猜一个Ω

出版版式(15)要求所有基站的接收功率满足 `s(k,m) >= q*s(j,m)`，其中 `q=varrho*Omega`、`varrho=2`。若q>1，两条相反约束相乘得到 `1>=q^2`，不存在所有链路均为正的解。原文“自相关/互相关比”的描述与直接套用存在歧义。

建议接口支持直接给复合比例q和“相关幅度→功率比例”的候选推导：

`q = varrho*abs(Rcross/Rauto)^2`。

后者必须标记为物理推导候选，并记录周期/非周期、码相位、最大值/RMS、码分配及幅度比/功率比。先验证数学路径可用，后续对q做敏感性和码相关性核验。

更深一层的限制：仅相对强弱条件不能保证真实接收机检测成功。若 `0<q<=1`、所有信道增益为正、QoS余量和预算都严格为正，则可令各基站对同一用户接收功率相等，再统一缩小功率，满足所有齐次近远约束和上界约束。因此这种模型下，经过优化的“全部用户可听”在一定程度上是约束所保证的，不能当作真实捕获成功率证据。绝对检测门限、最小C/N0属于后续接收机扩展，不能混入原文主结果而不说明。

### 3.4 几何、时钟与失效点

附录A用三维G矩阵，但仿真只给出四个二维角点。全部点共面且无高度差时，三维矩阵秩不足。第一阶段建议明确采用 `2d_known_height` 和论文的独立测距误差传播；这是一项实现假设。同步基站不自动等于已知接收机时钟，故要单列 `clockModel=paper_no_clock_state`，后续TDoA/未知钟差另建模型。

“more than 3 BSs”在K=4时字面意味着4站；一些常规二维定位讨论则取至少3站。主分支建议用 `minHearableBs=4` 对应原文字面，另报告3站敏感性。可听集合确定后重算G和λ，不能始终使用全4站λ。

`positionStdM=NaN` 表示失败，`valid=false`单独保存；原图将失败填0只用于显示。平均精度同时报告覆盖率和样本数，绝不把0计入平均。保存优化目标 `mean(Psi.^2)` 与报表均值 `mean(Psi(valid))`，二者不是同一个量。

### 3.5 噪声与通信用户位置

主论文未完整公开噪声PSD、接收机噪声系数、天线增益、高度、通信用户空间分布及随机种子。尤其QoS约束作用于通信用户，不能只生成20个定位用户就用它们替代全部通信用户。

第一阶段可建立一个显式假设的 `smoke_assumed` 配置，例如热噪声参考值、单位天线增益、二维自由空间、固定种子和按最近基站分区抽样通信用户。数值必须带 `assumed` 来源，不能自动成为 `paper_reproduction` 配置。正式配置允许 unresolved，运行前按所选实验的依赖检查，不依赖某项的单小区实验仍可运行。

### 3.6 CP4A需单独核验

- 式(24)把每组通信载波数设成 `2G−1`。`M*(2G−1)`通常不等于N；G=80时3180不等于1679。不能越界截断后仍声称原分解成立。
- 对固定正c，`d(-c/p)/dp=+c/p^2`，式(66)第二行却写负号；式(33)也受影响。先有限差分检验，再推导实现。
- 独立变量 `p(k,m)` 对某通信干扰的导数应是A的一列；式(67)对所有发射基站求和不能直接视为该独立变量偏导。
- 最强干扰基站会切换，完整对偶还包括自身作为其他链路干扰的贡献。首次λ=1也不能替代实际位置估计流程。

建议实现顺序是“原文目标与完整约束的集中式求解 → 全耦合对偶核验 → 论文双层CP4A流程研究”。集中式结果不命名为CP4A；经推导修正的算法标为 `cp4a_consistent_candidate`，未经核验的字面算法仅作对照。

另：现有公式推导讲解把式(33)写为含 `sqrt(M/denominator)` 的形式；出版PDF实际是 `[M*(...)]^(-1/2)`，还差一个M因子关系。这说明二手推导不能直接转写实现。

## 4. 工程目录与依赖方向

```text
F:/4Cell_MS_NOMA/
  Code/
    main.m                         唯一主入口，参数选择与调度
    README.md
    docs/                          本设计、来源台账、后续决策记录
    MS_NOMA/                       现有作者演示，保留独立入口
    +msnoma/
      +config/                     基础配置、实验配置、依赖验证
      +model/                      载波、频谱、单位与数值约定
      +scenario/                   坐标、用户抽样、信道
      +comm/                       干扰、BER与QoS反演
      +ranging/                    DLL完整积分、近似与测距因子
      +geometry/                   可听性、几何与误差传播
      +optim/                      约束组装、可行性、求解与KKT
      +prs/                        导频配置、功率映射与CRB
      +experiment/                 单次、多次与各类扫参
      +report/                     统计、绘图、导出
      +waveform/                   现有生成器适配，后续接收机入口
    tests/                         模型不变量与求解器交叉核验
  References/                      公开补充论文与来源信息
  Results/
    run_YYYYMMDD_HHMMSS_<id>/       所有实验输出，避免相互覆盖
```

只对 `Code` 使用 `addpath`，MATLAB包目录由命名空间解析；不使用 `addpath(genpath(...))` 把所有旧main、旧函数和结果目录混入路径。作者演示函数由适配器临时加其精确目录，并通过 `onCleanup` 恢复路径。

依赖只向下：`experiment → optim/comm/ranging/geometry → scenario/model`；`report`消费结果，不参与计算；`config`不调用绘图和求解器。一个文件一个主要函数，不把模型公式、求解过程或绘图塞进main。

第一阶段不需要类继承、全局变量或大型对象框架。使用有版本号的struct、显式参数传入和返回值即可。

## 5. 统一数据接口

### 5.1 维度、单位和向量化

| 数据 | 字段及形状 | 约定 |
|---|---|---|
| `cfg` | `paper, assumed, models, numerics, experiment, io, provenance` | 原文参数、假设与数值容差分开 |
| `num` | `K,M,N,G,TcS,TpS,dfCHz,dfPHz,requestedBandwidthHz,effectiveBandwidthHz,BfeHz,PcW,c0Mps` | 所有量在线性SI单位计算 |
| `scene` | `bsPosM[K×D]`, `pUserPosM[M×D]`, `cUserPosM[K×N×D]`, `seed` | 每条通信载波对应一个通信用户；快照可重放 |
| `ch` | `gCServing[K×N]`, `gPC[K×N×K×M]`, `gP[K×M]`, `gCP[M×K×N]` | 全是功率增益abs(h)^2，不是复振幅 |
| `kernels` | `W[N×M]`、有限频谱核、模型版本 | 同一载波配置跨随机场景复用 |
| `geom` | `lambda[K×M]`, `rank[M]`, `valid[M]`, `mode` | 对未采用的BS保留索引映射 |
| `rangeFactor` | `RtauS2W[K×M]`, `RrangeM2W[K×M]` | 两项近似：v_tau=Rtau/p，v_r=Rrange/p |
| `problem` | `c[KM×1], Aineq, bineq, lb, scale, rowMap, variableMap, provenance` | 求解与诊断共享一个模型 |
| `solution` | `powerW[K×M], status, objectiveM2, exitflag, multipliers, diagnostics, history` | 不可行时没有伪造可用功率 |
| `evaluation` | `ber[K×N]`, `delayVarianceS2[K×M]`, `rangeVarianceM2[K×M]`, `positionStdM[M×1]`, `hearable[K×M]`, `valid[M×1]` | 独立重算，不信任求解器缓存 |
| `summary` | `coverage, conditionalMeanM, objectiveMeanM2, nTotal, nValid, solverFailureRate, cdf` | 失败原因及统计分母显式保留 |

向量化统一采用MATLAB原生列优先：

```matlab
p = powerW(:);                % powerW: K-by-M
variableId(k,m) = k+(m-1)*K;
communicationRow(k,n) = k+(n-1)*K;
powerW = reshape(p,K,M);
```

这与旧指南建议的按基站连续编号不同，是本设计为减少MATLAB reshape错误做的明确选择。所有约束生成、导数、导出用相同映射；CSV额外存k、m、n列，不依赖读者推断展平顺序。

`gPC`与`gCP`不能随意互换：前者接收点是通信用户，后者是定位用户。自由空间平坦场景里gPC可沿m重复、gCP可沿n重复，内部可以压缩存储，但逻辑接口保持清楚。完整4×1658×4×20张量约53万个double、4.25MB，不需要为第一版过早引入复杂的分布式存储。

### 5.2 配置来源约定

每项未公开参数至少记录 `value / unit / sourceLevel / source / rationale`。sourceLevel采用 `paper`、`public_reference`、`derived`、`assumed`、`unresolved`。`strict`模式按实验所需字段拒绝unresolved输入；`smoke_assumed`允许有来源记录的工程假设并在导出标签中显示。

所有随机性由局部 `RandStream`生成；用场景编号导出可复现子流。20/50MHz与各算法使用同一定位用户场景，通信用户按统一抽样规则派生并记录；不得用不同随机场景混淆算法差异。

## 6. 函数规划及签名

以下均为拟实现接口，包前缀 `msnoma.` 在表内省略。核心计算函数不绘图、不写文件、不修改cfg或全局随机状态。每个输出包含必要的模型标识与诊断，避免以裸数组隐藏单位。

### 6.1 配置、基础模型、场景

| 拟定函数 | 输入→输出 | 责任 |
|---|---|---|
| `config.base(profile)` | profile字符串→cfg | 主论文参数与显式假设初值 |
| `config.forExperiment(cfg,name)` | cfg、实验名→cfg | 配置小场景/Monte Carlo/扫参，保持原文参数来源 |
| `config.validate(cfg)` | cfg→validation | 按依赖校验字段、分支、单位和工具箱；阻断项附原因 |
| `model.numerology(cfg)` | cfg→num | 唯一计算N、G、带宽与时间尺度的地方 |
| `model.spectralKernels(num,opts)` | num、积分策略→kernels | W、目标频移和可缓存谱核 |
| `model.sincpi(x)` | 数值数组→同形数组 | sin(pi*x)/(pi*x)，零点稳定 |
| `scenario.make(cfg,num,stream)` | 配置、载波、局部流→scene | 基站、定位用户和通信用户坐标 |
| `scenario.fromSnapshot(snapshot,num)` | 已保存坐标→scene | 重放相同场景，检查维度与用户编号 |
| `scenario.freeSpaceGain(txPosM,rxPosM,prop)` | 坐标、传播参数→g及距离 | Friis功率增益，最小距离策略显式 |
| `scenario.channels(scene,num,prop)` | scene、num、传播→ch | 组装各种接收点的增益；默认无额外Rayleigh衰落 |

### 6.2 通信与测距

| 拟定函数 | 输入→输出 | 对应原文/约束 |
|---|---|---|
| `comm.makeBerModel(spec)` | 模型与噪声约定→berModel | 打包模型ID及正反算系数 |
| `comm.interferenceMatrix(ch,num,kernels)` | ch、num、W→A及rowMap | 式(3)—(5)，A为KN×KM |
| `comm.interference(powerW,A,num)` | 功率、A→interferenceWPerHz[K×N] | 矩阵乘法及reshape |
| `comm.ber(gCServing,PcW,TcS,N0WPerHz,I,berModel)` | 链路、功率、噪声与I→BER[K×N] | 式(2)/(34)与候选分支 |
| `comm.qosLimit(gCServing,PcW,TcS,N0WPerHz,xi,berModel)` | 相同BER模型→b[K×N]及cleanFeasible | 式(12)或同模型推导；不静默裁负阈值 |
| `ranging.equivalentCommunicationGain(ch,num)` | gCP、num→gCeq[M×K] | 式(62)有限载波求和 |
| `ranging.factorTwoTerm(ch,num,noise,dll)` | 信道、数值、噪声、环路→rangeFactor | 式(10)，与定位功率无关的R |
| `ranging.varianceTwoTerm(powerW,rangeFactor)` | p、R→时延/距离方差 | 优化所用模型 |
| `ranging.varianceThreeTerm(powerW,ch,num,noise,dll)` | 完整功率、信道等→分项方差 | 式(9)/(61)，保留其他基站定位干扰 |
| `ranging.integralKernels(ch,num,dll,quadOpts)` | ch、num等→A0及噪声/通信/定位线性核 | 式(51)—(54)，核与p分离，返回收敛诊断 |
| `ranging.varianceIntegral(powerW,ch,kernel,noise)` | p、积分核→方差及各分项 | 式(6)/(55)，精确正弦判别器 |

频域积分须指定前端频率范围是相对目标P-Sub还是相对整个RF/基带信道；载波平移必须同时作用于积分上下限和所有PSD。论文与旧指南的这一步仍值得单独验证。先实现具名 `local_symmetric` 解释，并为 `global_frontend` 解释留分支，不把3.5GHz代入本地基带正弦。加密频率节点和延迟带宽近似比较只用于选定样本，不在每个优化迭代中重复积分。

### 6.3 几何、可听性与优化

| 拟定函数 | 输入→输出 | 责任 |
|---|---|---|
| `geometry.weights(scene,bsMask,geometrySpec)` | 坐标、K×M掩码、模型→geom | QR/SVD求H、λ、秩与条件数；不把秩亏pinv结果当有效 |
| `geometry.hearability(powerW,ch,hearSpec)` | 接收功率、q→mask及功率比 | 式(15)/(16)，零功率不算听到 |
| `geometry.positionAccuracy(rangeVarianceM2,geom)` | 测距方差与几何→position结果 | 式(1)、(45)—(50) |
| `optim.buildProblem(A,b,rangeFactor,geom,ch,cfg,num)` | 完整模型→problem | 目标、QoS、预算、全对可听、正定义域、映射与缩放 |
| `optim.checkFeasibility(problem,opts)` | problem→feasibility | 解析预检+线性可行性/正裕量检查，返回证据 |
| `optim.objective(z,scaledProblem)` | 缩放功率→f、梯度、Hessian | sum(c./p)，独立于具体求解器 |
| `optim.solveCentral(problem,opts)` | problem→solution | linprog起点+fmincon，解析导数与尺度归一化 |
| `optim.solveDual(problem,opts)` | 同一problem→solution | 完整耦合对偶候选；与集中式交叉验证 |
| `optim.solveCp4a(problem,opts)` | 同一problem→solution | 后期实现原文双层/修正候选，mode必须显式 |
| `optim.equalPower(problem,policy)` | problem、策略→solution | budget_equal与qos_scaled_equal；诊断全部约束 |
| `optim.diagnostics(problem,solution)` | 原始模型及结果→diagnostics | SI残差、归一化残差、KKT、互补松弛、目标差 |

等功率策略不施加会改变其本意的近远平衡；它可以违反可听性，结果通过统一评价反映覆盖差异。但必须报告QoS是否满足。`qos_scaled_equal`仅在预算与通信QoS下选择共同等功率尺度，不用可听约束预先把等功率变成优化方案。

### 6.4 PRS、实验、输出与波形适配

| 拟定函数 | 输入→输出 | 责任 |
|---|---|---|
| `prs.allocation(prsSpec)` | 资源配置→activeK、symbolIds、relativePowerWeights | LTE标准配置与50MHz推广分开 |
| `prs.powerBudget(energySpec,allocation,ch)` | 总能量、占空比等→每测量的接收导频功率和噪声比 | 公平比较的唯一映射入口 |
| `prs.delayCrb(allocation,receivedBudget,noiseSpec)` | 导频、接收功率→时延方差 | [28]式(5)，式(16)归一化；多符号按Fisher信息相加 |
| `prs.positionBenchmark(scene,rangeVarianceM2,geometrySpec)` | PRS方差、共同几何→精度 | 主比较沿用主论文几何，独立TDoA分支另标 |
| `experiment.evaluate(powerW,context,evalSpec)` | 功率及固定场景模型→evaluation | 从原始功率重算所有性能，模型外误差也单列 |
| `experiment.runSingle(cfg,algorithmIds,snapshot)` | cfg、算法、可选快照→run | 一次场景的完整路径；各算法共享context |
| `experiment.runMonteCarlo(cfg,algorithmIds)` | cfg、算法→batch | 重复runSingle，记录成功/失败场景及原因 |
| `experiment.runSingleCell(cfg)` | cfg→BER、资源、测距、PRS实验 | 图3—5输入数据；不混入多小区功率求解 |
| `experiment.runSensitivity(cfg,paramSpec)` | 配置、参数轴→对照集 | q、噪声、通信分布、几何和带宽规则 |
| `experiment.runSpatialMap(cfg,mapSpec)` | 固定19用户、查询点协议→地图结果 | 明确每点重新优化；不把所有网格点当M个用户 |
| `report.aggregate(batch,statisticsSpec)` | batch→summary | 条件均值、覆盖率、无条件CDF及分母 |
| `report.plot(runOrSummary,plotSpec)` | 数据→figure句柄 | 只可视化；visible由配置控制 |
| `report.saveRun(result,cfg,outputRoot)` | 完整结果→输出目录与manifest | MAT、CSV、PNG、配置与源码版本/哈希 |
| `waveform.authorAdapter(cfg,num,powerSpec)` | 统一模型→原生成器输入与信号 | 演示路径的尺度适配与功率断言 |

这些接口构成后续函数清单；不要求第一版把每种扩展都实现。`solveCp4a`、空间曲面、采样接收机应在前置验收通过后增加，不能用空函数返回伪正常结果。

## 7. 优化模型与集中式基准

固定场景、信道、几何和两项测距近似后：

```text
p = vec(P),  c(k,m) = lambda(k,m)^2 * Rrange(k,m) / M

minimize     F(p) = sum_i c_i / p_i
subject to   A*p <= b                         通信QoS
             Bbudget*p <= budgetW             每基站定位预算
             Chear*p <= 0                     所有(k,j,m)近远约束
             p > 0                           有效链路的定义域
```

其中一条 `Chear` 的两个非零元是：目标变量(k,m)处 `-gP(k,m)`、干扰变量(j,m)处 `q*gP(j,m)`。

K=4、M=20时只有80个功率变量、240条全对可听约束、4条预算约束。50MHz整数分支N=1658，通信约束6632条。A通常是稠密矩阵，因为sinc²旁瓣不严格为0；规模足以先用集中式求解，无须先用未经验证的载波截断获得“可解性”。

F在正定义域内凸，若所有c_i>0则严格凸。`fmincon`是数值实现工具，不因为其返回success就自动认定得到可用结果；必须独立检查可行性和KKT。近似求解器变量下限 `pMinW` 是工程措施，需要减小下限检查目标稳定，并对c=0情况单独处理。

缩放建议：`p=PscaleW*z`，每行约束做正尺度归一化，目标另用固定正尺度缩放。保存转换因子以恢复原始SI残差、乘子和KKT。不要以原始极小的自由空间增益直接套统一1e-6绝对容差。

无数值下限的内点对偶关系为：

```text
d = A'*mu + Bbudget'*nu + Chear'*beta
p_i = sqrt(c_i/d_i),  要求d_i>0
```

对偶更新或回溯必须维持有效域。不能用sqrt(abs(d))掩盖负分母；若使用有界内层，则求解该有界标量最小化并保持其数学定义。停止条件包括原始可行性、对偶残差、互补松弛和目标变化，不能只看相邻乘子变化。

## 8. PRS对照的实现边界

[28]式(5)给出单OFDM符号平坦AWGN时延CRB：

`v_tau = Ts^2 / (8*pi^2*SNR*sum_k(w_k*k^2))`，其中 `w_k=p_k^2`、`sum(w_k)=Nc`，`SNR=(C/N0)/B`。式(16)在等导频分配时规定 `w_k=Nc/Na`。C是该定义下的信号功率，不能无说明地用单载波Pc替代。

因此可以复用旧工程通用频率矩的计算思想，但必须重新验证功率归一化。多个符号、PRS占空比、积分时间、总测量能量、数据/PRS之间的时频资源关系都由 `prs.powerBudget`统一记录。2019论文的10.21m锚定曲线不进入2022独立预测。

两个关键分支：

1. **几何一致比较**：从[28]只借用时延CRB，转成距离方差后使用主论文式(45)—(50)计算PRS位置指标。这样MS-NOMA与PRS使用同一几何与时钟模型。
2. **[28]完整TDoA模型**：按式(20)—(23)引入参考站造成的相关测量协方差和TDoA方向矩阵。这个结果单独标记，不能与主论文OLS误差传播无说明混比。

[28]的标准LTE配置最高20MHz、15kHz子载波间隔；其公开示例为6RB、1.08MHz占用带宽、12导频、单符号。这不是主论文50MHz配置。后续 `prs.mode` 至少区分 `lte_standard` 与 `generalized_ofdm_50mhz`，后者是理论推广，或另明确载波聚合方式；不能向LTE Toolbox传一个无效的50MHz单载波配置。

[28]自身使用的43/46dBm发射功率、15dBi天线、9dB噪声系数、750m站间距及宏小区路损，都不是主论文4角点自由空间实验的默认参数。

## 9. main调用骨架

拟定实现后，桌面使用方式为：

```matlab
cd('F:/4Cell_MS_NOMA/Code');
main
```

main主体建议保持为配置与实验分派：

```matlab
% main.m -- 设计示意，目前尚未落地
codeRoot = fileparts(mfilename('fullpath'));
addpath(codeRoot);

profile = "smoke_assumed";
experimentName = "single_scene";
algorithmIds = ["budget_equal", "qos_scaled_equal", "central"];

cfg = msnoma.config.base(profile);
cfg = msnoma.config.forExperiment(cfg,experimentName);
cfg.io.resultsRoot = fullfile(fileparts(codeRoot),'Results');
cfg.io.showFigures = true;
validation = msnoma.config.validate(cfg);
assert(validation.canRun,validation.message);

switch experimentName
    case "single_scene"
        result = msnoma.experiment.runSingle(cfg,algorithmIds,[]);
    case "monte_carlo"
        result = msnoma.experiment.runMonteCarlo(cfg,algorithmIds);
    case "single_cell"
        result = msnoma.experiment.runSingleCell(cfg);
    otherwise
        error('Unknown experiment: %s',experimentName);
end

figures = msnoma.report.plot(result,cfg.io);
runDir = msnoma.report.saveRun(result,cfg,cfg.io.resultsRoot);
fprintf('Results: %s\n',runDir);
```

runSingle内部才是完整主干：

```matlab
num = msnoma.model.numerology(cfg);
kernels = msnoma.model.spectralKernels(num,cfg.numerics);
scene = ... scenario.make(...) 或 scenario.fromSnapshot(...);
ch = msnoma.scenario.channels(scene,num,cfg.assumed.propagation);
berModel = msnoma.comm.makeBerModel(cfg.models.ber);
[A,rowMap] = msnoma.comm.interferenceMatrix(ch,num,kernels);
qos = msnoma.comm.qosLimit(...,berModel);
rangeFactor = msnoma.ranging.factorTwoTerm(...);
geom = msnoma.geometry.weights(scene,true(num.K,num.M),cfg.models.geometry);
problem = msnoma.optim.buildProblem(A,qos,...);
feasibility = msnoma.optim.checkFeasibility(problem,cfg.numerics);
% 对同一problem运行各算法；不可行场景按原因退出，不生成伪正常结果
% 每个解均调用diagnostics和experiment.evaluate，最后计算PRS与统计
```

绘图统一接收结果对象，只改变显示方式不重新随机生成数据。默认桌面图窗打开，同时保存；批处理设置showFigures=false。

## 10. 分阶段实施与可检查交付

### 阶段A：数据契约与完整单场景计算路径

实现配置、载波、场景、信道、A与BER正反算、两项测距因子、二维几何、可听性、完整约束、集中式求解及等功率对照；main可以一次跑通4×20小实例。

此阶段明确使用 `geometry.source=oracle`，由仿真真实坐标计算几何，验证固定几何下的理论指标；并不声称已经实现论文所述未知位置到估计位置的反馈过程。

输出至少包括：固定坐标、全部80个功率、通信BER、80条测距方差、20个位置精度、覆盖率、目标与约束残差、未解决参数标签。第一版使用一个带宽和一个随机场景；这是完整主干的最小交付，不先追求所有图都出现。

### 阶段B：数值模型核验与单小区2022分支

增加精确DLL积分、三项近似、噪声/通信/定位分项比较，完成BER无干扰基线和反演核验；对照2022式(34)—(44)实现图3—5的MS-NOMA数据及PRS资源关系。它与2019快报属于不同配置。

### 阶段C：PRS与双带宽统计闭合

实现[28]时延CRB和公平功率映射，明确50MHz扩展，补上PRS主比较。完成相同场景下20/50MHz、50次Monte Carlo、无条件CDF与覆盖率，并导出图7快照与表II口径的统计。

### 阶段D：CP4A及公开疑点研究

使用同一problem完成全耦合对偶候选，再研究原文双层流程、载波分组与导数分支。与集中式在相同输入上比较目标、全约束残差和KKT，而非仅比较曲线外观。将“原文算法复现”和“经推导修正的实现”分别标记。

### 阶段E：空间曲面与采样链

固定可重放场景，明确空间网格查询协议；最后才扩展多基站波形、传播时延、接收信号、解调、DLL与位置估计。首次λ=1、估计位置反馈的闭环要在这个阶段有实际观测与解算接口，不能把真实坐标伪装成估计坐标。

资料缺口采用显式假设+敏感性推进，不把等待作者回复设为任何阶段的前置条件。当前没有依据承诺所有绝对数值均能逐点恢复。

## 11. 验收与结果存储

核心核验只覆盖会影响数学结论的内容，不为目录和简单包装函数堆测试。

| 检查 | 通过条件/记录 |
|---|---|
| 功率和单位 | 每载波/总功率不混用；时延s²到距离m²只乘一次c0² |
| 展平顺序与A | 随机场景下A*p与直接求和一致，含跨基站非对称增益 |
| BER正反算 | 任意测试点BER≤xi与I≤b判别一致，包含负b情形 |
| 两项测距 | p翻倍、其他量固定时方差减半；R不随p变化 |
| DLL积分 | 加密积分后稳定，再记录近似模型误差；不把近似偏差算成积分失败 |
| 几何 | 对称场景权重对称；秩亏返回无效，不产生伪有限精度 |
| 优化手算例 | 单预算无其他约束时p_i=P*sqrt(c_i)/sum(sqrt(c))；有限差分核验梯度/Hessian |
| 可行性 | q>1且全链路正功率时明确拒绝；噪声已超QoS时明确拒绝 |
| 求解质量 | 原始/归一化残差、KKT、互补松弛与目标变化同时报告 |
| CP4A候选 | 同一问题相对集中式目标差和约束残差合格；未收敛有明确状态 |
| PRS | 固定导频和能量时CRB满足SNR的倒数标度，功率归一化可追溯 |
| 统计 | 无条件CDF终点等于覆盖率；无效点不计0误差；求解失败不被静默丢弃 |
| 可重放 | 相同配置、种子和快照得到相同输入及相符数值结果 |

初始容差可将归一化约束残差1e-6、缩放KKT残差1e-4、不同求解器相对目标差1e-3作为工程目标，实际须在尺度与条件数试验后定稿，不称论文容差。

每次运行建议保存：

```text
manifest.json           模型分支、来源、工具箱版本、源文件哈希、终止原因
config.mat / config.json
scene.mat               全部坐标、种子、用户编号
model.mat               num、A/b、R、lambda、约束映射与缩放
solutions.mat           各算法功率、乘子、残差、目标轨迹
evaluation.mat          各条链路和用户的重新评价
powers.csv / users.csv / constraints.csv / summary.csv
figures/*.png           明确带模型、单位与假设标签
```

覆盖不足、优化不可行、求解未收敛和几何秩亏是不同状态。配置级失败单列；算法比较使用共同有效场景并另外报告失败率；不能把不可行样本静默删除后当作原始50次平均。

## 12. 论文目标与当前状态

原文表II的定位精度是20MHz下MS-NOMA 0.57m、PRS 8.11m，50MHz下MS-NOMA 0.21m、PRS 3.28m；这是待核对的论文结果，不是此次计算输出。原图等功率不能定位比例50.8%，也不应拿来反向拟合q、噪声或通信用户分布。

| 当前能确定的事项 | 当前不能确定的事项 |
|---|---|
| 完整理论主干可以实现，优化规模适中 | 作者的全部随机实例、通信用户布置与未披露噪声设置 |
| 多小区需新建场景、信道、几何和优化模块 | Ω的确切相关性定义及数值 |
| 原作者演示代码可服务波形侧检查 | 印刷CP4A在原作者实现中如何处理公式疑点 |
| PRS[28]全文已取得，可独立实现CRB | 主论文使用的PRS资源、50MHz映射和时间积累方式 |
| 本机具备第一阶段所需MATLAB工具箱 | 全部绝对精度能否与论文逐点重合 |

下一步应先实现阶段A，以可复现单场景输出和残差报告作为第一轮检查对象，再扩展模型核验与论文图表。
