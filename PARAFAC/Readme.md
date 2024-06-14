@[TOC](目录)
# 一、引言
&emsp; 内容简介： 本文主要讲解了**均匀面阵中基于PARAFAC分解的二维DOA估计算法**。详细推导了均匀面阵的阵列流型的构成。然后介绍了PAFARAC模型，通过切片示意图说明分解的方式。接着又分析了如何利用PARAFAC分解获得的阵列流型进行二维的DOA估计，顺便提及与PAFARAC相关的张量分解的部分内容。最后逐段给出该模型的仿真代码，搭建出与前文公式推导相匹配的基于PARAFAC分解的大体框架，对代码估计的效果进行了验证，并指出其中可以优化的方向。
# 二、均匀面阵模型

&emsp; 如图1所示的均匀面阵，该面阵共有 $M \times N$ 个阵元，阵元均匀分布，相邻阵元的间距是 $d$，$d \leqslant \lambda/2$（$\lambda$ 是信号波长）。
![ 图1.1 二维均匀面阵模型 ](https://i-blog.csdnimg.cn/direct/5ce4cbf3be314269a4387c6b95016f08.png#pic_center =500x300)
<center><p>图1 二维均匀面阵模型</p></center>


&emsp; 假设空间中有 $K$ 个信号入射到此均匀面阵上，其二维DOA为 $(\theta_k, \phi_k)$，$k = 1,2,\cdots,K$，其中 $\theta_k$ 和 $\phi_k$ 分别代表第 $k$ 个信源的仰角和方位角。定义 $u_k = \sin\theta_k \sin\phi_k$，$v_k = \sin\theta_k \cos\phi_k$，示意图如图2所示。  

![图1.2](https://i-blog.csdnimg.cn/direct/13d21bca8a514978b02418f0ea24eca8.png#pic_center)
<center><p>图2 二维均匀面阵角度对应关系</p></center>

&emsp; $x$ 轴和 $y$ 轴上阵元的方向向量分别为：  

$$
\begin{aligned}
	\boldsymbol{a}_x\left( \theta _k,\phi _k \right) &=\left[ \begin{array}{c}
	1\\
	{e}^{\text{j}2\pi d\sin \theta _k\cos \phi _k/\lambda}\\
	\vdots\\
	{e}^{\text{j}2\pi \left( M-1 \right) d\sin \theta _k\cos \phi _k/\lambda}\\
\end{array} \right] ,\ \boldsymbol{a}_y\left( \theta _k,\phi _k \right) =\left[ \begin{array}{c}
	1\\
	{e}^{\text{j}2\pi d\sin \theta _k\sin \phi _k/\lambda}\\
	\vdots\\
	{e}^{\text{j}2\pi \left( N-1 \right) d\sin \theta _k\sin \phi _k/\lambda}\\
\end{array} \right]\\
\end{aligned}  \tag{1}
$$

&emsp; 子阵 1 的接收信号可表示为：  
$$
\boldsymbol{x}_1(t) = \boldsymbol{A}_x \boldsymbol{s}(t) + \boldsymbol{n}_1(t)
$$  
式中，$\boldsymbol{A}_x = [\boldsymbol{a}_x(\theta_1, \phi_1), \boldsymbol{a}_x(\theta_2, \phi_2), \cdots, \boldsymbol{a}_x(\theta_K, \phi_K)]$ 为子阵 1 的方向矩阵；$\boldsymbol{n}_1(t)$ 为子阵 1 的加性高斯白噪声；$\boldsymbol{s}(t) \in \mathbb{C}^{K \times 1}$ 为信源向量，这与一维均匀线阵的标准形式相同。  
&emsp; 接下来，子阵 2 的接收信号为： 
$$
\boldsymbol{x}_2(t) = \boldsymbol{A}_x \boldsymbol{\varPhi}_y^{1} \boldsymbol{s}(t) + \boldsymbol{n}_2(t)
$$  
式中，$\boldsymbol{A}_x \boldsymbol{\varPhi}_y^{1}$ 为子阵 2 的方向矩阵，$\boldsymbol{\varPhi}_y = \mathrm{diag}(\mathrm{e}^{\mathrm{j}2\pi d \sin\theta_1 \sin\phi_1 / \lambda}, \cdots, \mathrm{e}^{\mathrm{j}2\pi d \sin\theta_K \sin\phi_K / \lambda})$；$\boldsymbol{n}_2(t)$ 为子阵 2 的加性高斯白噪声。  


><font color=#FF0000 size=4>解释 : 为何$\boldsymbol{A}_x \boldsymbol{\varPhi}_y^{1}$ 是子阵2的方向矩阵？</font>
><font size=3>
>- 对于同一直线上的两个阵元，如图(a)，信源入射的波程差为： $\Delta x=dcos\alpha$ 
>- 将这两个阵元在平面上并排放置，如图(b)，波程差为：$\Delta x=d_y cos \alpha$
></font>
>![在这里插入图片描述](https://i-blog.csdnimg.cn/direct/d2f2959c5f7b49749a257bb8c9ee5b7c.png#pic_center =720x250)
>- 再将这两个阵元在平面中任意放置，如图\(c\)和(d)，波程差为：$\Delta x=\Delta x_1+\Delta x_2=d_x cos\beta +d_y cos \alpha$
>- 由于**子阵2**相对于**子阵1**，在 $y$ 轴方向上只相差1个阵元间距 $d$，因此 $\boldsymbol{A}_x$ 每个列向量将分别乘以指数因子  $\mathrm{e}^{\mathrm{j}2\pi \textcolor{red}{d} \sin\theta_k \sin\phi_k / \lambda}$ ,  $k = 1,2,\cdots,K$ ，同理，子阵3的每个列向量将分别乘以指数因子 $\mathrm{e}^{\mathrm{j}2\pi \textcolor{red}{2d} \sin\theta_k \sin\phi_k / \lambda}$。![在这里插入图片描述](https://i-blog.csdnimg.cn/direct/387095fc85d04390945cb481f9eb0a19.png#pic_center)



所以推广至子阵 $n$ 时的接收信号为：  
$$
\boldsymbol{x}_n(t) = \boldsymbol{A}_x \boldsymbol{\varPhi}_y^{n - 1} \boldsymbol{s}(t) + \boldsymbol{n}_n(t) \tag{2}
$$  
由此可得，整个均匀面阵的接收信号为：  
$$
\boldsymbol{x}(t) = \begin{bmatrix}
\boldsymbol{x}_1(t) \\
\boldsymbol{x}_2(t) \\
\vdots \\
\boldsymbol{x}_N(t)
\end{bmatrix} = \begin{bmatrix}
\boldsymbol{A}_x \\
\boldsymbol{A}_x \boldsymbol{\varPhi}_y \\
\vdots \\
\boldsymbol{A}_x \boldsymbol{\varPhi}_y^{N - 1}
\end{bmatrix} \boldsymbol{s}(t) + \begin{bmatrix}
\boldsymbol{n}_1(t) \\
\boldsymbol{n}_2(t) \\
\vdots \\
\boldsymbol{n}_N(t)
\end{bmatrix}   \tag{3}
$$  

由于$y$ 轴上 $N$ 个阵元对应的方向矩阵为 $\boldsymbol{A}_y = [\boldsymbol{a}_y(\theta_1, \phi_1), \boldsymbol{a}_y(\theta_2, \phi_2), \cdots, \boldsymbol{a}_y(\theta_K, \phi_K)]$，具体表示为：  
$$
\boldsymbol{A}_y = \begin{bmatrix}
1 & 1 & \cdots & 1 \\
\mathrm{e}^{\mathrm{j}2\pi d \sin\theta_1 \sin\phi_1 / \lambda} & \mathrm{e}^{\mathrm{j}2\pi d \sin\theta_2 \sin\phi_2 / \lambda} & \cdots & \mathrm{e}^{\mathrm{j}2\pi d \sin\theta_K \sin\phi_K / \lambda} \\
\vdots & \vdots & & \vdots \\
\mathrm{e}^{\mathrm{j}2\pi (N - 1)d \sin\theta_1 \sin\phi_1 / \lambda} & \mathrm{e}^{\mathrm{j}2\pi (N - 1)d \sin\theta_2 \sin\phi_2 / \lambda} & \cdots & \mathrm{e}^{\mathrm{j}2\pi (N - 1)d \sin\theta_K \sin\phi_K / \lambda}
\end{bmatrix} \tag{4}
$$  

对比可见 $\boldsymbol{\varPhi}_y^{N - 1}$ 实际上就是 $\boldsymbol{A}_y$ 的第$N$行的元素组成的对角矩阵，因此上式也可表示为：

$$
\boldsymbol{x}(t) = [\boldsymbol{A}_y \odot \boldsymbol{A}_x] \boldsymbol{s}(t) + \boldsymbol{n}(t)  \tag{5}
$$


式中，$\boldsymbol{x}$ 是$MN\times1$维度的矩阵，$\odot$ 表示的Khatri-Rao积[(解释的链接)](https://blog.csdn.net/xuehuitanwan123/article/details/104291475?ops_request_misc=%257B%2522request%255Fid%2522%253A%2522ccf64cf36119580514c86f8a26091d90%2522%252C%2522scm%2522%253A%252220140713.130102334..%2522%257D&request_id=ccf64cf36119580514c86f8a26091d90&biz_id=0&utm_medium=distribute.pc_search_result.none-task-blog-2~all~baidu_landing_v2~default-6-104291475-null-null.142^v102^pc_search_result_base3&utm_term=khatri-rao%E7%A7%AF&spm=1018.2226.3001.4187)。根据 Khatri-Rao 积的定义，接收信号可以写为：
$$
\boldsymbol{x}(t) = [\boldsymbol{a}_y(\theta_1, \phi_1) \otimes \boldsymbol{a}_x(\theta_1, \phi_1), \cdots, \boldsymbol{a}_y(\theta_K, \phi_K) \otimes \boldsymbol{a}_x(\theta_K, \phi_K)] \boldsymbol{s}(t) + \boldsymbol{n}(t) 
$$  
式中，$\otimes$ 代表 Kronecker 积。假设对于 $L$ 次采样，$\boldsymbol{a}_x(\theta_k, \phi_k)$ 和 $\boldsymbol{a}_y(\theta_k, \phi_k)$ 不变，定义 $\boldsymbol{X} = [\boldsymbol{x}(1), \boldsymbol{x}(2), \cdots, \boldsymbol{x}(L)] \in \mathbb{C}^{MN \times L}$，均匀面阵的接收信号可表示为：

$$
\boldsymbol{X} = [\boldsymbol{A}_y \odot \boldsymbol{A}_x] \boldsymbol{S}^{\mathrm{T}} + \boldsymbol{N} = \begin{bmatrix}
\boldsymbol{X}_1 \\
\boldsymbol{X}_2 \\
\vdots \\
\boldsymbol{X}_N
\end{bmatrix} = \begin{bmatrix}
\boldsymbol{A}_x \boldsymbol{D}_1(\boldsymbol{A}_y) \\
\boldsymbol{A}_x \boldsymbol{D}_2(\boldsymbol{A}_y) \\
\vdots \\
\boldsymbol{A}_x \boldsymbol{D}_N(\boldsymbol{A}_y)
\end{bmatrix} \boldsymbol{S}^{\mathrm{T}} + \begin{bmatrix}
\boldsymbol{N}_1 \\
\boldsymbol{N}_2 \\
\vdots \\
\boldsymbol{N}_N
\end{bmatrix} \tag{6} 
$$  
式中，$\boldsymbol{S} = [\boldsymbol{s}(1), \boldsymbol{s}(2), \cdots, \boldsymbol{s}(L)]^{\mathrm{T}} \in \mathbb{C}^{L \times K}$ 由 $L$ 次采样的信号向量组成；$\boldsymbol{D}_m(\cdot)$ 为由矩阵第 $m$ 行构造的对角矩阵；$\boldsymbol{N} = [\boldsymbol{n}(1), \boldsymbol{n}(2), \cdots, \boldsymbol{n}(L)]$ 为加性高斯白噪声矩阵。$\boldsymbol{X}$矩阵的示意图如图3所示，它可以视为是将一个**3维的立方体数据**铺展成一系列的**2维的平面数据**。
![在这里插入图片描述](https://i-blog.csdnimg.cn/direct/a10ee20ab6ef4b5d9b7c3adda72848e3.png#pic_center =761x317)
<center><p>图3 矩阵X的切片示意图</p></center>

---

# 三、PAFARAC模型及分解
## 3.1 数据模型

&emsp; 如果不将图3中的三维立方体展开，则其中的元素可按照如下方式描述： 
$$
x_{m,n,l} = \sum_{k=1}^{K} \boldsymbol{A}_x(m,k) \boldsymbol{A}_y(n,k) \boldsymbol{S}(l,k), \ m = 1,2,\cdots,M, \ n = 1,2,\cdots,N, \ l = 1,2,\cdots,L  \tag{7}
$$  
该形式称为阵列接收信号的PARAFAC模型， 式中，$\boldsymbol{A}_x(m,k)$、$\boldsymbol{A}_y(n,k)$ 和 $\boldsymbol{S}(l,k)$分别为 $x$ 轴方向矩阵 $\boldsymbol{A}_x$ 的第 $(m,k)$ 个元素、$y$ 轴方向矩阵 $\boldsymbol{A}_y$ 的第 $(n,k)$ 个元素和信源矩阵 $\boldsymbol{S}$ 的第 $(l,k)$ 个元素。
&emsp;上一节中得到的 $\boldsymbol{X}_n \ (n = 1,2,\cdots,N)$ 可以看作沿三个空间维度中Y轴的维度对三维矩阵切片得到的，每个矩阵可用如下形式表示：
$$
\boldsymbol{X}_n = \boldsymbol{A}_x \boldsymbol{D}_n(\boldsymbol{A}_y) \boldsymbol{S}^{\mathrm{T}} + \boldsymbol{N}_n, \ n = 1,2,\cdots,N \tag{8.1}
$$  

由PARAFAC模型的对称性可以得到三维矩阵沿另外两个维度的切片形式，如图4和图5，公式如下： 

$$
\boldsymbol{Y}_m = \boldsymbol{S} \boldsymbol{D}_m(\boldsymbol{A}_x) \boldsymbol{A}_y^{\mathrm{T}} + \boldsymbol{N}, \ m = 1,2,\cdots,M \tag{8.2}
$$  

$$
\boldsymbol{Z}_l = \boldsymbol{A}_y \boldsymbol{D}_l(\boldsymbol{S}) \boldsymbol{A}_x^{\mathrm{T}} + \boldsymbol{N}_t, \ l = 1,2,\cdots,L \tag{8.3}
$$ 
![在这里插入图片描述](https://i-blog.csdnimg.cn/direct/1c583d8c47fe498aa02fb96c7eece87b.png#pic_center =761x317)
<center><p>图4 矩阵Y的切片示意图</p></center>


![在这里插入图片描述](https://i-blog.csdnimg.cn/direct/ac628f39674941919825c3c5a93aa2f4.png#pic_center =761x317)

<center><p>图5 矩阵Z的切片示意图</p></center>

>- 注：上图中，矩阵$\boldsymbol{X}$沿$y$轴切片，矩阵$\boldsymbol{Y}$沿着$x$轴切片，矩阵$\boldsymbol{Z}$沿着时域切片。<font color=#FF0000 size=3>这里可能疑惑为什么矩阵$\boldsymbol{X}$不沿$x$轴切片</font>。实际上$\boldsymbol{X}$代表样本值，与$x$并无直接关系，$\boldsymbol{Y}$和$\boldsymbol{Z}$只是顺着$\boldsymbol{X}$向后选了两个字母，因此 $\boldsymbol{X,Y,Z}$三个矩阵与$x,y,z$轴无直接关系。如果非要大小写字母对应，可以选择$\boldsymbol{Y}$作为信号接收阵列，另外两者可逐一对应即可。
>- 原始的数学模型可见硕士论文`《基于平行因子分析的阵列参数估计》`的 3.3.1 小节，这其中是使用$\boldsymbol{A,B,C}$三个矩阵进行的解释。

按照公式(6)对 $\boldsymbol{X}$ 的处理，矩阵$\boldsymbol{X}$、$\boldsymbol{Y}$和$\boldsymbol{Z}$可以表示为：

$$
\boldsymbol{X} = [\boldsymbol{X}_1 \boldsymbol{X}_2 \dots \boldsymbol{X}_N]^\mathrm{T} 
= [\boldsymbol{A}_y \odot \boldsymbol{A}_x] \boldsymbol{S}^{\mathrm{T}} + \boldsymbol{N}_x \tag{9.1}
$$  

$$
\boldsymbol{Y} = [\boldsymbol{Y}_1 \boldsymbol{Y}_2 \dots \boldsymbol{Y}_M]
 = [\boldsymbol{A}_x \odot \boldsymbol{S}] \boldsymbol{A}_y^{\mathrm{T}} + \boldsymbol{N}_y \tag{9.2}
$$  

$$
\boldsymbol{Z} = [\boldsymbol{Z}_1 \boldsymbol{Z}_2 \dots \boldsymbol{Z}_L]
= [\boldsymbol{S} \odot \boldsymbol{A}_y] \boldsymbol{A}_x^{\mathrm{T}} + \boldsymbol{N}_z \tag{9.3}
$$  

## 3.2 PARAFAC 分解  
&emsp; 交替最小二乘法(Alternating Least Squares , ALS)是常用于 PARAFAC 模型的一种方法。ALS 通过交替固定其他因子矩阵，仅对某一个因子矩阵进行最小二乘优化，逐步逼近最优解。运行时首先进行初始化，为各因子矩阵（如信源矩阵、阵列流型）赋予初始值，然后根据式(9.1)计算 $\boldsymbol{S}$ 的LS解。式中 $[\cdot]^{\dagger}$ 表示广义逆运算。
$$
\boldsymbol{S}^{\mathrm{T}} = [\boldsymbol{A}_y \odot \boldsymbol{A}_x]^{\dagger} \boldsymbol{X} \tag{10.1}
$$
接着由式(9.2)得到 $y$ 轴的阵列流型
$$
\boldsymbol{A}_y^{\mathrm{T}} = [\boldsymbol{A}_x \odot \boldsymbol{S}]^{\dagger} \boldsymbol{Y} \tag{10.2}
$$
接着由式(9.3)得到 $x$ 轴的阵列流型
$$
\boldsymbol{A}_x^{\mathrm{T}} = [\boldsymbol{S} \odot \boldsymbol{A}_y]^{\dagger} \boldsymbol{Z} \tag{10.3}
$$

重复交替上述优化步骤，直至因子矩阵收敛（比如迭代前后变化量小于设定阈值）。


---

# 四、二维DOA估计
 
&emsp; 利用 PARAFAC 分解估计出 $\boldsymbol{A}_x$、$\boldsymbol{A}_y$ 后，设 $\boldsymbol{A}_x$ 的某一列为 $\boldsymbol{a}_x$，取 $\text{angle}(\boldsymbol{a}_x)$，因其范围为 $[-\pi, \pi]$，通过对某些项加上 $2k\pi$ 使其成为递增的序列。按以上方法对 $\boldsymbol{A}_x$ 的每一列进行调整后，利用阵列之间的相位差可以估计其 DOA。具体过程如下：

$$
\boldsymbol{a}_x(\theta_k, \phi_k) = [1 \ \mathrm{e}^{\mathrm{j}2\pi d \sin\theta_k \cos\phi_k / \lambda} \ 
 \dots \ \mathrm{e}^{\mathrm{j}2\pi (M - 1)d \sin\theta_k \cos\phi_k / \lambda} ]^{T}
$$  
使用 $\text{angle}$ 函数获得向量的角度值 $\boldsymbol{h}$ ，
$$
\begin{aligned}
\boldsymbol{h} &= \text{angle}(\boldsymbol{a}_x(\theta_k, \phi_k)) \\
  &= \left[0, 2\pi d \sin\theta_k \cos\phi_k / \lambda, \cdots, 2\pi (M - 1)d \sin\theta_k \cos\phi_k / \lambda\right]^{\mathrm{T}}
\end{aligned} \tag{11}
$$
由于$\boldsymbol{h}$的每个元素之间的差值相同，所以通过自相减可得 $\sin\theta_k \cos\phi_k$ 的值，然后取平均以减小误差。
$$
v_k = \sin\theta_k \cos\phi_k = mean\left[\frac{(\boldsymbol{h}[2:end] - \boldsymbol{h}[1:M-1]) }{2\pi d / \lambda}\right]  \tag{12}
$$  
同理，依据$\boldsymbol{A}_y$可计算出$u_k$。综合两者可得目标方向的二维DOA估计为：
$$
\begin{aligned}
\hat{\theta}_k &= \sin^{-1}(\sqrt{\hat{u}_k^2 + \hat{v}_k^2}) \\
\hat{\phi}_k &= \tan^{-1}(\hat{u}_k / \hat{v}_k)
\end{aligned}  \tag{13}
$$

&emsp; 至此，均匀面阵中基于PARAFAC模型的二维DOA估计算法描述完毕。

---
# 五、引申-张量分解

&emsp; 前文中的立方体数据，在数学上有专用名称：<font color=#FF0000 size=4>张量</font>。简单理解，一维数据构成向量，二维数据构成矩阵，三维及更高维数据则称为张量。维数的增高意味着处理的复杂度随之升高，为简化运算，出现许多张量分解算法。<font color=#FF0000 size=4>CP分解</font>是张量分解的算法之一，它将张量分解为一系列的秩一张量之和，其中[秩一张量](https://zhuanlan.zhihu.com/p/382219345)是由三个向量张成的。CP分解的示意图如图6，它类似于矩阵的SVD分解，将矩阵分解为许多因子矩阵之和。从图6这个角度来看，CP分解中最小的单元是不同的向量，向量的外积组成了不同的张量，因子张量之和就是原始张量。
![在这里插入图片描述](https://i-blog.csdnimg.cn/direct/70f07cb78afd4635a6fc5089df496e87.png#pic_center =780x230)
<center><p>图6 CP分解示意图</p></center>
&emsp; 数学公式可由公式(14)表示：

$$
x_{ijk} \approx \sum_{r=1}^{R} a_{ir} b_{jr} c_{kr} \quad \text{for } i = 1, \ldots, I, \ j = 1, \ldots, J, \ k = 1, \ldots, K. \tag{14}
$$

&emsp; CP分解不止一种表示形式，它还可以以矩阵为最小单元，按照如下方式分解：
$$
\begin{aligned}
\boldsymbol{X}_{(1)} &\approx \boldsymbol{A}(\boldsymbol{C} \odot \boldsymbol{B})^{\mathrm{T}}, \\
\boldsymbol{X}_{(2)} &\approx \boldsymbol{B}(\boldsymbol{C} \odot \boldsymbol{A})^{\mathrm{T}}, \\
\boldsymbol{X}_{(3)} &\approx \boldsymbol{C}(\boldsymbol{B} \odot \boldsymbol{A})^{\mathrm{T}}.
\end{aligned} \tag{15}
$$
&emsp; 可以发现，这种分解方式与PARAFAC的分解是相同的，所以PARAFAC的分解实际属于CP分解的类型之一。从更大的范围讲，CP分解又属于Tucker分解的特殊类型，BTD分解则兼具CP分解和Tucker分解两者的特点......太多了写不完 [附上链接](https://blog.csdn.net/FDS99999/article/details/131234000)。

---
# 六、代码解析
以下是代码，运行环境是MATLAB R2020b，由于涉及张量分解的`cp_als 函数`等，需提前[安装Tensor Toolbox](https://blog.csdn.net/qq_37637914/article/details/116016157)工具包。

## 6.1 参数配置
```matlab
%% 参数设置
M = 6;           % x轴阵元个数
N = 5;           % y轴阵元个数
K = 512;        % 快拍数
fc = 100e6;     % 载波频率
fs = 300e6;     % 采样频率
Pn = 1;          % 噪声功率

fines = [20, 40]; % 方位角（度）
thetas = [5, 60];    % 俯仰角（度）
signal_SNR = [30, 30]; % 信噪比（dB）
signal_f = [15e6, 30e6]; % 信号频率

m = (0:M-1)';    % x轴阵元索引
n = (0:N-1)';    % y轴阵元索引
c = 3e8;         % 光速
lambda = c / fc;  % 波长
dx = lambda / 2;  % x轴阵元间距
dy = lambda / 2;  % y轴阵元间距
num_signals = length(fines); % 信号数目
```

## 6.2 生成信号
两路信号分别进行载波调制；
```matlab
%% 生成信号
t = (0:K-1)/fs; % 时间向量
S = zeros(num_signals, K);
for k = 1:num_signals
    A_k = sqrt(Pn)*10^(signal_SNR(k)/20); % 信号幅度
    S(k,:) = A_k * exp(1j*2*pi*signal_f(k)*t); % 载波调制
end
```
## 6.3 构造阵列流型
按照第一节中公式(1)进行构造；
```matlab
%% 构造阵列流型
A = zeros(M, num_signals);
B = zeros(N, num_signals);
for k = 1:num_signals
    phi = deg2rad(fines(k));
    theta = deg2rad(thetas(k));
    
    % 空间频率计算
    u = (dx/lambda)*sin(theta)*cos(phi);
    v = (dy/lambda)*sin(theta)*sin(phi);
    
    A(:,k) = exp(-1j*2*pi*m*u); % x方向导向矢量
    B(:,k) = exp(-1j*2*pi*n*v); % y方向导向矢量
end
```

## 6.4 构造接收张量
使用`ktensor 函数`创建秩为 1 的 Kruskal对象 component，其内容为$\boldsymbol{A}$的列向量、$\boldsymbol{B}$的列向量、$\boldsymbol{S}$的行向量张成的张量，循环遍历每个信号，最终叠加在一起形成张量$\boldsymbol{X}$，参考图6理解。
```matlab
%% 构造接收张量
X = tensor(zeros(M,N,K));
noise = (randn(M,N,K) + 1j*randn(M,N,K)) * sqrt(Pn/2);

for k = 1:num_signals
    component = ktensor(1, A(:,k), B(:,k), S(k,:).');
    X = X + tensor(component);
end
X = X + tensor(noise); 
X_normalized = X / norm(X); % 张量归一化
```
## 6.5 ALS分解

结构体`struct`存储 PARAFAC 分解的相关参数。初始化方法为`nvecs`，这是种较为稳定的初始化方法（也选`random`）；`printitn` 置为1在分解过程中会输出每次迭代的相关信息；`tol` 设置收敛容差，当相邻两次迭代的误差变化小于这个容差时，就认为分解已经收敛，迭代过程会停止；`maxiters` 设置最大迭代次数为 150。若迭代次数达到这个上限，即便还未满足收敛条件，分解过程也会停止。
><font color=#FF0000 size=3>待优化：迭代次数与收敛容差的平衡</font>
>- 当收敛容差设置较大时，迭代会提前停止；设置较小时，迭代次数过多，出现过拟合；
>- 当迭代次数设置较小时，迭代会提前停止；设置较大时，迭代次数过多，出现过拟合；
>- 两者都会导致出现误差，感兴趣的朋友可自行优化。
```matlab
%% ALS分解
R = num_signals; % 分解秩
options = struct;
options.init = 'nvecs';  % 使用nvecs初始化
options.printitn = 1;    % 显示迭代过程
options.tol = 1e-4;      % 设置收敛容差
options.maxiters = 150;  % 增加最大迭代次数

[Factors, ~] = cp_als(X_normalized, R, options);

A_est = Factors{1}; 
B_est = Factors{2};
```
## 6.6 参数估计
此处按照公式(11-13)编写代码，可从估计得到的`A_est`和`B_est`中获得角度 $(\theta_k, \phi_k)$。`unwrap 函数`对应第四节提到的“**通过对某些项加上 $2k\pi$ 使其成为递增的序列**”。需要说明的是，这段代码有个<font color=#FF0000>小问题</font>：计算`phi_est`的`atan2d 函数`前添加了负号`-`。这与公式(13)是相悖的，但是不添加的话，估计的角度则刚好是真实值的相反数，我还没想明白 TT，欢迎评论区讨论。
```matlab
%% 参数估计
estimated_angles = zeros(num_signals, 2);
for d = 1:R
    % x方向估计
    a = A_est(:,d);
    phase_x = unwrap(angle(a));
    u_est = -(phase_x(2:end) - phase_x(1:end-1))/(2*pi*(dx/lambda));
    u_est_avg = mean(u_est);
    
    % y方向估计
    b = B_est(:,d);
    phase_y = unwrap(angle(b));
    v_est = -(phase_y(2:end) - phase_y(1:end-1))/(2*pi*(dy/lambda));
    v_est_avg = mean(v_est);
    
    % 角度解算
    phi_est = -atan2d(v_est_avg , u_est_avg);
    theta_est = asind(sqrt(u_est_avg.^2 + v_est_avg.^2));
    
    estimated_angles(d,:) = [phi_est, theta_est];
end
```
## 6.7 结果展示
```matlab
%% 二维可视化
figure(10);
% 真实位置三维散点图
scatter(fines, thetas, 'r', 'filled');hold on;
% 估计位置三维散点图
scatter(estimated_angles(:,1), estimated_angles(:,2), 'b', 'd');
xlabel('方位角 (度)'); ylabel('俯仰角 (度)'); zlabel('幅度');
title('信号分布');
axis([0,90,0,90]);
grid on; 
% 性能评估
RMSE_phi = sqrt(mean((fines - estimated_angles(:,1)').^2));
RMSE_theta = sqrt(mean((thetas - estimated_angles(:,2)').^2));
disp(['方位角RMSE: ', num2str(RMSE_phi), ' 度']);
disp(['俯仰角RMSE: ', num2str(RMSE_theta), ' 度']);
```
&emsp; 估计结果如图7所示，设定估计的连个方向分别为 $(20°,5°)$ 、 $(40°,60°)$，RMSE的误差也展示在图中。如果你更改了方位角和俯仰角，重新估计的话，大概率会出现图8所示情况，这说明出现了**欠拟合和过拟合**哈哈哈，如果想让代码稳定，可以去6.5小节优化分解过程，暂时就写到这里吧。
![avb](https://i-blog.csdnimg.cn/direct/e71be01a9ee147249d1cd8f1eb74274d.png#pic_center =800x)
<center><p>图7 正常估计结果</p></center>

![在这里插入图片描述](https://i-blog.csdnimg.cn/direct/1f1f754d6eb149ccb0f47c0b8d26f2a4.png#pic_center =800x)
<center><p>图8 异常估计结果</p></center>


==更新中......==

---


<font  size=5>参考文献</font >
- 学位论文 [基于平行因子分析的阵列参数估计](https://d.wanfangdata.com.cn/thesis/D053542)
- 知乎 [张量基础| 张量的秩与秩一张量](https://zhuanlan.zhihu.com/p/382219345)
- 知乎 [深入理解 | CP、Tucker分解](https://zhuanlan.zhihu.com/p/302453223)
- Github [二维DOA估计算法](https://github.com/yashcao/Graduate-Study-Report/blob/master/PhD-1%20study/%E4%BA%8C%E7%BB%B4DOA%E4%BC%B0%E8%AE%A1%E7%AE%97%E6%B3%95.pdf)
-  图书 [阵列信号处理及MATLAB实现 第三版]()
- 工具 [DeepSeek](https://chat.deepseek.com/) & [豆包](https://www.doubao.com/chat/)

