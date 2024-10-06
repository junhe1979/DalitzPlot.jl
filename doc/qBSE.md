# Quasipotential approximation

The general form of the BSE for the scattering amplitude is in a from

$$
\begin{align}
{\cal M}(k'_1k'_2,k_1k_2;P)&={\cal
V}(k'_1k'_2,k_1k_2;P)+\int\frac{d^4
k''_2}{(2\pi)^4}
{\cal
V}(k'_1k'_2,k''_1k''_2;P)G(k''_1k''_2;P){\cal
M}(k''_1k''_2,k_1k_2;P),\quad
\end{align}
$$

where ${\cal V}$ is the potential kernel and $G$ is the propagators for two constituent particles. Here the momentum of the system $P=k_1+k_2=k'_1+k'_2=k''_1+k''_2$. It can be abbreviated as

$$
\begin{align}
{\cal M}&={\cal V}+{\cal V}{G}{\cal M}.
\end{align}
$$

The Gross form of proposed quasipotential propagators for particles 1 and 2 with mass $m_1$ and $m_2$ written down in
the center of mass frame where $P=(W,{\bm 0})$ with particle 2 being on shell are

$$
\begin{align}
g=2\pi i\frac{\delta^+(k_2^2-m_2^2)}{k_1^2-m_1^2}=2\pi
i\frac{\delta(k^0_2-E_2)}{2E_2[(W-E_2)^2-E_1^2]},
\end{align}
$$

where $k_1=(k_1^0,\bm k)=(E_1,\bm k)$, $k_2=(k_2^0,-\bm k)=(W-E_1,-\bm k)$ with $E_1=\sqrt{m_1^2+|\bm k|^2}$.

With the define of $G_0=g/(2\pi i)$, the four-dimensional BSE can be reduced to a three-dimensional equation in center of mass frame

$$
\begin{align}
i{\cal M}({\bm k}',{\bm k})&=i{\cal
V}({\bm k}',{\bm k})+\int\frac{d
{\bm k}''}{(2\pi)^3}
i{\cal
V}({\bm k}',{\bm k}'')G_0({\bm k}'')i{\cal
M}({\bm k}'',{\bm k}),\quad
\end{align}
$$

Note: the $i{\cal M}$ and $i{\cal V}$ are usually real. In the center of mass frame. We choose ${\bm k}_2={\bm k}$ and ${\bm k}_1=-{\bm k}$.

## Partial-wave expansion

To reduce the equation to one-dimensional equation, we apply the partial wave expansion,

$$
\begin{align}
	&{\cal V}_{\lambda'\lambda}({\bm k}',{\bm k})=
	\sum
	 _{J\lambda_R}{\frac{2J+1}{4\pi}}D^{J*}_{\lambda_R,\lambda'}(\phi',\theta',-\phi'){\cal
	V}^J_{\lambda'\lambda,\lambda_R}({\rm k}',{\rm k})D^{J}_{\lambda_R,\lambda}(\phi,\theta,-\phi)
\nonumber\\
&\to{\cal V}_{\lambda'\lambda,\lambda_R}^J({\rm k}',{\rm k})=
{\frac{2J+1}{4\pi}}\int d\Omega' d\Omega
D^{J}_{\lambda_R,\lambda'}(\phi',\theta',-\phi'){\cal V}_{\lambda'\lambda}({\bm k}',{\bm k})D^{J*}_{\lambda_R,\lambda}(\phi,\theta,-\phi)
\nonumber\\&\to{\cal V}_{\lambda'\lambda}^J({\rm k}',{\rm k})=2\pi\int d\cos\theta_{k,k'} d^{J}_{\lambda,\lambda'}(\theta_{k',k})
{\cal V}_{\lambda'\lambda}({\bm k}',{\bm k})
\end{align}
$$

where the momenta are chosen as $k_2=(E_2,0,0,{\rm k})$, $k_1=(W-E_2,0,0,-{\rm k})$  and $k'_2=(E'_2,{\rm k}'\sin\theta_{k,k'},0,{\rm k}'\cos\theta_{k,k'})$, $k=(W-E_2,-{\rm k}'\sin\theta_{k,k'},0,-{\rm k}'\cos\theta_{k,k'})$ with ${\rm k}=|{\bm k}|$ and ${\rm k}'=|{\bm k}'|$. NOTE: Which particle is chosen to parallel to $z$ axis is related to the order of $\lambda$ and $\lambda'$ in $d^{J}_{\lambda'\lambda}(\theta_{k,k'})$, so it can not be chosen arbitrarily. And the definition of helicity is also dependent of the definition of ${\bm k}_{1,2}$. Here, $\lambda=\lambda_2-\lambda_1$  and $\lambda_1=-s_1$, $\lambda_2=s_2$. The scattering amplitudes ${\cal M}$ has analogous relations.

Now we have the partial wave BS equation,

$$
\begin{align}
i{\cal M}^J_{\lambda',\lambda}({\rm k}',{\rm k})&=
{\frac{2J+1}{4\pi}}\int d\Omega' d\Omega D^{J*}_{\lambda_R,\lambda'}(\phi',\theta',-\phi')D^{J}_{\lambda_R,\lambda}(\phi,\theta,-\phi)\nonumber\\
&\cdot\left[i{\cal V}_{\lambda'\lambda}({\bm k}',{\bm
k})+\int\frac{d{\bm k}''}{(2\pi)^3}i{\cal V}_{\lambda'\lambda''}({\bm
k}',{\bm k}'') G_0({\bm k}'')i{\cal M}_{\lambda''\lambda}({\bm k}'',{\bm k})\right]\nonumber\\
&=i{\cal V}^J_{\lambda'\lambda}({\rm k}',{\rm k})+\int\frac{{\rm k}''^2d{\rm k}''}{(2\pi)^3}i{\cal V}^J_{\lambda'\lambda''}({\rm k}',{\rm k}'')
G_0({\bm k}'')i{\cal M}^J_{\lambda''\lambda}({\rm k}'',{\rm k})
\end{align}
$$

where $\int d\Omega D^{J*}_{\lambda_1,\lambda_2}(\phi,\theta,-\phi)D^{J'*}_{\lambda'_1,\lambda'_2}(\phi,\theta,-\phi)={\frac{4\pi}{2J+1}}$ is used.

## Fixed parity

For a helicity state $|J,\lambda\rangle=|J,\lambda_1\lambda_2\rangle$ fulfill the party property,

$$
\begin{align}
	P|J,\lambda\rangle=P|J,\lambda_1\lambda_2\rangle
	=\eta_1\eta_2(-1)^{J-s_1-s_2}|J,-\lambda_1-\lambda_2\rangle
	\equiv\tilde{\eta}|J-\lambda\rangle
\end{align}
$$

The construction of normalized states with parity $\pm$ is now straightforward:

$$
\begin{align}
	 |J,\lambda;\pm\rangle&=\frac{1}{\sqrt{2}}(|J,+\lambda\rangle\pm\tilde{\eta}
	|J,-\lambda\rangle\nonumber\\
	\Rightarrow P|J,\lambda;\pm\rangle&=
	 \frac{1}{\sqrt{2}}(\tilde{\eta}|J,-\lambda\rangle\pm|J,+\lambda\rangle)
	 =\pm\frac{1}{\sqrt{2}}(\pm\tilde{\eta}|J,-\lambda\rangle+|J,\lambda\rangle)
        =\pm|J,\lambda;\pm\rangle
\end{align}
$$

$$
\begin{align}
	{\cal M}^{J\pm}_{\lambda'\lambda}=\langle J,\lambda';\pm|{\cal M}|J,\lambda;\pm\rangle
	={\cal M}^{J}_{\lambda'\lambda}\pm \tilde{\eta} {\cal M}^J_{\lambda'-\lambda},\quad
{\cal M}^{J\pm}_{\lambda'-\lambda}=\pm\tilde{\eta} {\cal M}^{J\pm}_{\lambda'\lambda}\equiv\eta {\cal M}^{J\pm}_{\lambda'\lambda},\ \ {\cal M}^{J\pm}_{-\lambda'\lambda}=\pm\tilde{\eta}' {\cal M}^{J\pm}_{\lambda'\lambda}\equiv\eta' {\cal M}^{J\pm}_{\lambda'\lambda}
\end{align}
$$

with $\eta=PP_1P_2(-1)^{J_1+J_2-J}=P(-1)^{1/2+J}$.
The potential ${\cal V}^{J^P}_{\lambda'\lambda}$ has analogous relations.

$$
\begin{align}
&	 i{\cal M}_{\lambda\lambda'}=i{\cal V}_{\lambda\lambda'}+\sum_{\lambda''}i{\cal V}_{\lambda\lambda''}Gi{\cal M}_{\lambda''\lambda'},\quad	 \eta'i{\cal M}_{\lambda-\lambda'}=\eta'i{\cal V}_{\lambda-\lambda'}+\sum_{\lambda''}i{\cal V}_{\lambda\lambda''}G\eta'i{\cal M}_{\lambda''-\lambda'}\nonumber\\
&\Rightarrow i{\cal M}^{J^P}_{\lambda\lambda'}=i{\cal V}^{J^P}_{\lambda\lambda'}+\sum_{\lambda''}i{\cal V}_{\lambda\lambda''}Gi{\cal M}^{J^P}_{\lambda''\lambda'},\quad
i{\cal M}^{J^P}_{\lambda\lambda'}=i{\cal V}^{J^P}_{\lambda\lambda'}+\sum_{\lambda''}i{\cal V}_{\lambda-\lambda''}Gi{\cal M}^{J^P}_{-\lambda''\lambda'}
=iV^{J^P}_{\lambda\lambda'}+\sum_{\lambda''}i{\cal V}_{\lambda-\lambda''}G\eta''i{\cal M}^{J^P}_{-\lambda''\lambda'}\nonumber\\
&\Rightarrow  i{\cal M}^{J^P}_{\lambda\lambda'}=i{\cal V}^{J^P}_{\lambda\lambda'}+\frac{1}{2}\sum_{\lambda''}i{\cal V}^{J^P}_{\lambda\lambda''}Gi{\cal M}^{J^P}_{\lambda''\lambda'},
\end{align}
$$

As shown in Eq.~(9), the amplitudes are not independent. If we only keep the independent amplitudes, the equation for definite parity can be written as

$$
\begin{align}
&i{\cal M}^{J^P}_{ij}=i{\cal V}^{J^P}_{ij}+\frac{1}{2}\sum_{\lambda''}i{\cal V}^{J^P}_{i\lambda''}Gi{\cal M}^{J^P}_{\lambda''j}
=i{\cal V}^{J^P}_{ij}+\frac{1}{2}i{\cal V}^{J^P}_{i0}Gi{\cal M}^{J^P}_{0j}+\sum_{k\neq0}i{\cal V}^{J^P}_{ik}Gi{\cal M}^{J^P}_{kj},\nonumber\\
\Rightarrow&i\hat{{\cal M}}^{J^P}_{ij}=i\hat{\cal V}^{J^P}_{ij}+\sum_{k}i\hat{\cal V}^{J^P}_{ik}Gi\hat{\cal M}^{J^P}_{kj}.
\end{align}
$$

where $i$, $j$, $k$ are the indix for the independent amplitudes and $\hat{\cal M}=f_if_j {\cal M}$ with $f_0=\frac{1}{\sqrt{2}}$ and $f_{i\neq0}=1$ with $0$ for the amplitudes with $\lambda_1=\lambda_2=0$.

The Bethe-Saltpeter equation for partial-wave amplitude with fixed spin-parity $J^P$ reads ,

$$
\begin{align}
i\hat{\cal M}^{J^P}_{\lambda'\lambda}({\rm p}',{\rm p})
&=i\hat{\cal V}^{J^P}_{\lambda',\lambda}({\rm p}',{\rm
p})+\sum_{\lambda''}\int\frac{{\rm
p}''^2d{\rm p}''}{(2\pi)^3}~
i\hat{\cal V}^{J^P}_{\lambda'\lambda''}({\rm p}',{\rm p}'')
G_0({\rm p}'')i\hat{\cal M}^{J^P}_{\lambda''\lambda}({\rm p}'',{\rm
p}).
\end{align}
$$

Note that the sum extends only over nonnegative helicity $\lambda''$. The partial wave potential is defined as

$$
\begin{align}
i\hat{\cal V}^{J^P}_{\lambda'\lambda''}({\rm p}',{\rm p}'')=f_if_j i{\cal V}_{\lambda'\lambda}^{J^P}({\rm p}',{\rm p})
&=f_if_j 2\pi\int d\cos\theta
~[d^{J}_{\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'\lambda}({\bm p}',{\bm p})
+\eta d^{J}_{-\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'-\lambda}({\bm p}',{\bm p})],
\end{align}
$$

## Transformation to a matrix equation

Now We have a integral equation with singularity in $G_0=\frac{1}{2 E((s-E)^2-\omega^2)}$  at $W=E_1+E_2$. This singularity can be isolated as,

$$
\begin{align}
i\hat{\cal M}^{J^P}({\rm p}',{\rm p})
&=i\hat{\cal V}^{J^P}({\rm p},{\rm p}')+\int\frac{{\rm p}''^2d {\rm p}''}{(2\pi)^3}i\hat{\cal V}^{J^P}({\rm p},{\rm p}'') G_0({\rm p}'')i\hat{\cal M}^{J^P}({\rm p}'', {\rm p}')\nonumber\\
&=i\hat{\cal V}^{J^P}+\mathcal{P}\int\frac{{\rm p}''^2d {\rm p}''}{(2\pi)^3}i\hat{\cal V}^{J^P}
G_0\hat{\cal M}^{J^P}-\pi i\int\frac{{\rm p}''^2d {\rm p}''}{(2\pi)^3}\hat{\cal V}^{J^P}_o\delta G_0\hat{\cal M}^{J^P}_o\nonumber\\
&=i\hat{\cal V}^{J^P}+\int\frac{d {\rm p}''}{(2\pi)^3}\left[{\rm p}''^2 G_0 i\hat{\cal V}^{J^P} \hat{\cal M}^{J^P}-\frac{M\hat{\cal V}^{J^P}_o\hat{\cal M}^{J^P}_o}{{\rm p}''^2-\bar{{\rm p}''}^2}\right]
-i\frac{\bar{{\rm p}''}^2\delta\bar{G}_0}{8\pi^2}\hat{\cal V}^{J^P}_o\hat{\cal M}^{J^P}_o\nonumber\\
&=i\hat{\cal V}^{J^P}+\int\frac{d {\rm p}''}{(2\pi)^3}{\rm p}''^2 G_0 i{\cal V}^{J^P} \hat{\cal M}^{J^P}
-[\int\frac{d {\rm p}''}{(2\pi)^3}\frac{M}{{\rm p}''^2-\bar{{\rm p}''}^2}+i\frac{\bar{{\rm p}''}^2\delta\bar{G}_0}{8\pi^2}]\hat{\cal V}^{J^P}_o\hat{\cal M}^{J^P}_o\nonumber\\
\end{align}
$$

with $\delta G_0=\frac{\delta(s-E-\omega))}{2E(s-E+\omega)}\theta(s-m_1-m_2)=\delta \bar{G}_0\delta({\rm p}''-\bar{{\rm p}''})=\frac{1}{4W\bar{{\rm p}''}}\delta({\rm p}''-\bar{{\rm p}''})\theta(s-m_1-m_2)$, $M=[{\rm p}''^2({\rm p}''^2-\bar{{\rm p}''}^2)G_0]_{{\rm p}''\to\bar{{\rm p}''}}\theta(s-m_1-m_2)=-\frac{\bar{{\rm p}''}^2}{2W}\theta(s-m_1-m_2)$.

We have

$$
\begin{align}
Im~G=-\rho/2=-\frac{\bar{{\rm p}''}}{32\pi^2 W}.
\end{align}
$$

It suggests the unitary is satisfied if the potential $i{\cal V}$ is real.

$$
\begin{align}
-T^\dag \rho T=2 T^\dag~ ImG~T=2T^\dag(-Im T^{-1})T=2T^\dag\frac{1}{2i}(T^{\dag-1}- T^{-1})T=i(T-T^\dag)
\end{align}
$$

where $T=i{\cal M}$.

With the Gauss discretization, the one-dimensional equation can be transformed as a matrix equation as

$$
\begin{align}
i\hat{\cal M}^{J^P}_{ik}
&=&i\hat{\cal V}^{J^P}_{ik}+\sum_{j=0}^N i\hat{\cal V}^{J^P}_{ij}G_ji\hat{\cal M}^{J^P}_{jk}\Rightarrow \hat{M}^{J^P}=\hat{V}^{J^P}+\hat{V}^{J^P}G\hat{M}^{J^P}
\end{align}
$$

$$
\begin{align}
	G_j=\left\{\begin{array}{cl}-\frac{i\bar{q}}{32\pi^2 W}+\sum_j
\left[\frac{w(q_j)}{(2\pi)^3}\frac{\bar{q}^2}
{2W{(q_j^2-\bar{q}^2)}}\right] & {\rm for}\ j=0,\ {\rm if}\ Re(W)>m_1+m_2,\nonumber\\
\frac{w(q_j)}{(2\pi)^3}\frac{q_j^2}
	{2E(q_j)[(W-E(q_j))^2-\omega^2(q_j)]}& {\rm for}\ j\neq0
	\end{array}\right.
\end{align}
$$

where $\bar{q}=\frac{1}{2W}\sqrt{[W^2-(m_1+m_2)^2][W^2-(m_1-m_2)^2]}$. The indices $i, j, k$ is for discrete momentum values, independent helicities, and coupled channels.

The default dimension is $[G] = 1$. Recalling that a factor of $2m$ should be included if a constituent particle is a fermion, we have $[G] = \text{GeV}^{n_f} \to [V] = [M] = \text{GeV}^{-n_f}$, with $n_f$ being the number of fermions. Therefore, under our convention where $\bar{u}u = 1$, the dimension of the potential must satisfy the above requirements. This criterion can be employed to verify the consistency of the Lagrangian and the derived potential.

Hence, for the channels above its thresholds, the matrix have an extra dimension.
We take two channel as example to explain the coupled-channel equation. The region of $W$ is divided as
$W<m_{1}$, $m_{1}<W<m_{2}$ and $W>m_{2}$.

$$
\begin{align}
%
V&=\left(\begin{array}{cc}
V^{NN}_{11}&V^{NN}_{12}\\
V^{NN}_{21}&V^{NN}_{22}
\end{array}\right),\quad
G=\left(\begin{array}{cc}
G^{N}_{1}&0\\
0&G^{N}_{2}
\end{array}\right),\quad
W<m_{1},\\
%
V&=\left(\begin{array}{cc}
V^{N+1N+1}_{11}&V^{N+1N}_{12}\\
V^{NN+1}_{21}&V^{NN}_{22}
\end{array}\right),\quad
G=\left(\begin{array}{cc}
G^{N+1}_{1}&0\\
0&G^{N}_{2}
\end{array}\right),\quad
 m_{1}<W<m_{2}
\\
V&=\left(\begin{array}{cc}
V^{N+1N+1}_{11}&V^{N+1N+1}_{12}\\
V^{N+1N+1}_{21}&V^{N+1N+1}_{22}
\end{array}\right),\quad
G=\left(\begin{array}{cc}
G^{N+1}_{1}&0\\
0&G^{N+1}_{2}
\end{array}\right),\quad W>m_{2}
\end{align}
$$

## Code

In code, we choose $V^{J^P}=\hat{V}^{J^P}/4\pi$, $G\to4\pi{G}$, and $M^{J^P}=\hat{M}^{J^P}/4\pi$. Hence, the qBSE becomes

$$
\begin{align}
M^{J^P}=V^{J^P}+V^{J^P}GM^{J^P}.
\end{align}
$$

Such convention is consistent with that in the chiral unitary approach.

$$
\begin{align}
V^{J^P}&=\hat{V}^{J^P}/4\pi=i\hat{\cal V}^{J^P}_{\lambda'\lambda''}({\rm p}',{\rm p}'')/4\pi=f_if_j i{\cal V}_{\lambda'\lambda}^{J^P}({\rm p}',{\rm p})/4\pi \nonumber\\
&=\frac{1}{2}f_if_j \int d\cos\theta
~[d^{J}_{\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'\lambda}({\bm p}',{\bm p})
+\eta d^{J}_{-\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'-\lambda}({\bm p}',{\bm p})],
\end{align}
$$


$$
\begin{align}
	G_j=\left\{\begin{array}{cl}-\frac{i\bar{q}}{8\pi W}+\sum_j
\left[\frac{w(q_j)}{2\pi^2}\frac{\bar{q}^2}
{2W{(q_j^2-\bar{q}^2)}}\right] & {\rm for}\ j=0,\ {\rm if}\ Re(W)>m_1+m_2,\nonumber\\
\frac{w(q_j)}{2\pi^2}\frac{q_j^2}
	{2E(q_j)[(W-E(q_j))^2-\omega^2(q_j)]}& {\rm for}\ j\neq0
	\end{array}\right.
\end{align}
$$
## Physical observable

After extend the energy in the center of mass frame $W$ into complex energy plane as $z$, the pole can be found by variation of $z$ to satisfy

$$
\begin{align}
|1-V(z)G(z)|=0
\end{align}
$$

with $z=E_R+i\Gamma_R/2$.

With the obtained amplitude $M^{J^P}$, we can also calculate the physical observable. Note that all physical observable are at real axis, we choose the onshell momentum as

$$
\begin{align}
M_{ij}(z)=\{[(1-VG)^{-1}]V\}_{ij}
\end{align}
$$

with $i$ and $k$ chosen as the onshell momentum, that is, $0$ dimension for $G$, and extra dimension for $V$.



### The cross section for the channel considered

For the open channel, the cross section can be obtained as
$$
\begin{align}
	\frac{d\sigma}{d\Omega}=\frac{1}{\tilde{j}_1\tilde{j}_2}\frac{1}{64\pi^2
	s}\frac{|{\bm k}'|}{|{\bm k}|}\sum_{\lambda\lambda'}|{\cal M}_{\lambda\lambda'}({\bm k})|^2.
\end{align}
$$
The total cross section can be written as

$$
\begin{align}
	\sigma&=\frac{1}{\tilde{j}_1\tilde{j}_2}\frac{1}{64\pi^2
	s}\frac{|{\bm k}'|}{|{\bm k}|}\int d\Omega\sum_{\lambda}|{\cal M}_{\lambda\lambda'}(|{\bm k}|)|^2
=\frac{1}{\tilde{j}_1\tilde{j}_2}\frac{1}{64\pi^2
	s}\frac{|{\bm k}'|}{|{\bm k}|}\sum_{J,\lambda}\frac{\tilde{J}}{4\pi}|{\cal M}^J_{\lambda\lambda'}(|{\bm k}|)|^2\nonumber\\
&=\frac{1}{16\pi
	s}\frac{|{\bm k}'|}{|{\bm k}|}\sum_{J^P,i\geq0 j\geq0}\frac{\tilde{J}}{\tilde{j}_1\tilde{j}_2}\left|\frac{\hat{{\cal M}}^{J^P}_{ij}(|{\bm k}|)}{4\pi}\right|^2
	=\frac{1}{16\pi
	s}\frac{|{\bm k}'|}{|{\bm k}|}\sum_{J^P,i\geq0 j\geq0}\frac{\tilde{J}}{\tilde{j}_1\tilde{j}_2}\left|M_{ij}^{J^P}\right|^2.
\end{align}
$$
NOTE: $4MM'$ should be multiplied if convention $\bar{u}u=1$ adopted.

$$
\begin{align}
	\sigma&\propto \sum_{J,\lambda'\lambda}
	|M^{J}_{\lambda'\lambda}|^2=\sum_{J,\lambda'j=0}|M^{J}_{\lambda'j}|^2+\sum_{J,\lambda'j>0}
	\left[|M^{J}_{\lambda'j}|^2+|M^{J}_{\lambda'-j}|^2\right]\nonumber\\
&=\sum_{J^P,\lambda'j=0}\delta_{\eta,1}|\frac{1}{2} M^{J^P}_{\lambda'j}|^2+\sum_{J,\lambda'j>0}
	 \left[\frac{1}{4}|M^{J^+}_{\lambda'j}+M^{J^-}_{\lambda'j}|^2+\frac{1}{4}|M^{J^+}_{\lambda'j}-M^{J^-}_{\lambda'j}|^2\right]\nonumber\\
&=\sum_{J^P,i=0j=0}\delta_{\eta,1}\delta_{\eta',1}\frac{1}{2}M^{J^P}_{ij}|^2+\sum_{J^P,i>0j=0}2\delta_{\eta,1}|\frac{1}{2} M^{J^P}_{ij}|^2+\sum_{J^P,\lambda'j>0}
	\frac{1}{2}|M^{J^P}_{\lambda'j}|^2\nonumber\\
&=\sum_{J^P,i=0j=0}\delta_{\eta,1}\delta_{\eta',1}|\frac{1}{2} M^{J^P}_{ij}|^2
+\sum_{J^P,i>0j=0}\delta_{\eta,1}|\frac{1}{\sqrt{2}} M^{J^P}_{ij}|^2
+\sum_{J^P,i=0j>0}\delta_{\eta',1}|\frac{1}{\sqrt{2}}M^{J^P}_{ij}|^2
+\sum_{J^P,i>0j>0}|M^{J^P}_{\lambda'j}|^2\nonumber\\
&=\sum_{J^P,i=0j=0}\delta_{\eta,1}\delta_{\eta',1}|\frac{1}{2} M^{J^P}_{ij}|^2
+\sum_{J^P,i>0j=0}\delta_{\eta,1}|\frac{1}{\sqrt{2}} M^{J^P}_{ij}|^2
+\sum_{J^P,i=0j>0}\delta_{\eta',1}|\frac{1}{\sqrt{2}}M^{J^P}_{ij}|^2
+\sum_{J^P,i>0j>0}|M^{J^P}_{ij}|^2\nonumber\\
&=\sum_{J^P,i\geq0j\geq0}|f_{i}f_{j}M^{J^P}_{ij}|^2=\sum_{J^P,i\geq0j\geq0}|\hat{M}^{J^P}_{ij}|^2
\end{align}
$$


### Argand plot

The amplitudes can be written as

$$
\begin{align}
	{\cal M}({\bm k})=-8\pi\sqrt{s}f({\bm
	k})=-\frac{8\pi\sqrt{s}}{|{\bm k}|}\sum
	 _{J\lambda_R}{\frac{2J+1}{4\pi}}D^{J*}_{\lambda_R,\lambda}(\phi,\theta,-\phi)
	a^J_{\lambda\lambda'}(|{\bm
	k}|)D^{J}_{\lambda_R,\lambda'}(\phi',\theta',-\phi').
\end{align}
$$
where $a^J=-\frac{|{\bm k}|}{8\pi\sqrt{s}}{\cal M}(|{\bm k}|)$, which can be displayed as a trajectory in an Argand plot.

# Three body decay

## kinematics

### Lorentz boost
Here, we consider an process $Y\to m_1X\to m_1m_2m_3$.  
To study a $1\to3$ decay with the qBSE, we need consider the center of mass frame (CMS) of $Y$ (which is also the laboratory frame in this issue) and the $m_1m_2$ where the qBSE is applied. The momenta of initial and final particles in the CMS of $Y$ are
$$
\begin{align}
P=(W,0,0,0);\  \ p_1=(E_1,{\bm p}_1);\ \ p_2=(E_2,{\bm p}_2);\ \ p_3=(E_3,{\bm p}_3)
\end{align}
$$

The momenta of particles 2 and 3 can be written with ${\bm p}^{cm}_3$ with Lorentz boost from $(m,{\bm 0})$ to $(E,{\bm k})$,
$$
\begin{align}
\Lambda^{\mu}_\nu=\frac{1}{m}\left(\begin{array}{cccc}
E({\bm k})&k_x&k_y&k_z\nonumber\\
k_x&m+\frac{k_x k_x}{E+m}&\frac{k_x k_y}{E+m}&\frac{k_x k_z}{E+m}\nonumber\\
k_y&\frac{k_y k_x}{E+m}&m+\frac{k_y k_y}{E+m}&\frac{k_y k_z}{E+m}\nonumber\\
k_z&\frac{k_z k_x}{E+m}&\frac{k_z k_y}{E+m}&m+\frac{k_z k_z}{E+m}\nonumber\\
\end{array}\right).
\end{align}
$$
 

With  the Lorentz boost  the momenta for particle 23 in the laboratory frame $(E_{23},-{\bm p}_1)$ can be written with the momenta in the c.m.s of particles 23 $(M_{23},{\bm 0})$ as
$$
\begin{align}
{\bm p}&={\bm p}^{cm}-\frac{{\bm p}_1}{M_{23}}\left[\frac{-{\bm p}_1\cdot {\bm p}^{cm}}{W-E_1({\rm p}_1)+M_{23}}+p^{0cm}\right],\nonumber\\
p^{0}&=\frac{1}{M_{23}}\left[(M-E_1({\rm p}_1))p^{0cm}-{\bm p}_1\cdot{\bm p}^{cm}\right],
\end{align}
$$
where the $p_{23}+p_1=P$ is applied, and $M_{23}=\sqrt{(p_2+p_3)^2}=\sqrt{(p^{cm}_2+p^{cm}_3)^2}=\sqrt{(P-p_1)^2}$.
Above boost is performed by the momenta of final particles, which are all onshell. If applied to 
final pariticles directly, $p_{1,2,3}$ can be written with $(\Omega_1,\Omega^{cm}_3,M_{23})$ beacuase $p^{0cm}$ can be expressed by $(\Omega_1,\Omega^{cm}_3,M_{23})$. However, if applied to intermediate particles, their momenta should be dependent on  $(\Omega_1,p^{cm}_3,\Omega^{cm}_3,M_{23})$ due to  offshellness.

The momenta in CMS of $23$ can also be written with the momentum in laboratory frame as
$$
\begin{align}
{\bm p}^{cm}&={\bm p}+\frac{{\bm p}_1}{M_{23}}\left[\frac{{\bm p}_1\cdot {\bm p}}{M-E_1({\rm p}_1)+M_{23}}+p^0\right], \nonumber\\
p^{0cm}&=\frac{1}{M_{23}}\left[(M-E_1({\rm p}_1))p^0+{\bm p}_1\cdot{\bm p}\right].
\end{align}
$$

### Phase space
The general expression of the decay width is given by
$$\begin{align}
d\Gamma=\frac{1}{2E}|{\cal M}|^2 d\Phi, \quad d\Phi=(2\pi)^4\delta^4(P-p_1-p_2-p_3) \frac{d^3p_1}{2E_1(2\pi)^3}\frac{d^3p_2}{2E_2(2\pi)^3}\frac{d^3p_3}{2E_3(2\pi)^3}
\end{align}$$

To study the invariant mass spectrum of the particles 2 and 3, it is convenient to rewrite the Lorentz-invariant phase space $d\Phi$ by taking as integration variables the direction of the momentum of particle 2 ${\bm p}_2^{cm}$ in the center-of-mass (cm) frame of particles 2 and 3. Thus, we first rewrite the phase factor as
$$\begin{align}
d\Phi&=(2\pi)^4\delta^4(P-p_1-p_2-p_3) \frac{d^3p_1}{2E_1(2\pi)^3}\frac{d^3p_2}{2E_2(2\pi)^3}\frac{d^3p_3}{2E_3(2\pi)^3}\nonumber\\
&=(2\pi)^4\delta(E^{cm}_2+E^{cm}_3-W_{23})\delta^3({\bm p}^{cm}_2+{\bm p}^{cm}_3) \frac{d^3p_1}{2E_1(2\pi)^3}\frac{d^3p^{cm}_2}{2E_2^{cm}(2\pi)^3}\frac{d^3p^{cm}_3}{2E_3^{cm}(2\pi)^3}
\end{align}$$
where $W_{23}^2=(M-E_1)^2-|{\rm p}_1|^2$. Here the Lorentz invariance of the  $\frac{d^3p}{2E(2\pi)^3}$ and $\delta^4(P-p_1-p_2-p_3)$ is used.

In this work, we use the momentum of the particle 3 in the center of mass system of two rescattering particles.
The momentum of the particle 3 has a relation ${\rm p}^{cm}_3=\frac{1}{2M_{23}}\sqrt{\lambda(M_{23}^2,m_3^2,m_2^2)}$ with $M_{23}=\sqrt{(p_2+p_3)^2}=\sqrt{(p^{cm}_2+p^{cm}_3)^2}$.

Owing to the three-momentum $\delta$ function, the integral over ${\bm p}^{cm}_2$ can be eliminated. Next, the quantity $d^3p^{cm}_3$ is converted to $dM_{23}$ by the relation,
$$\begin{align}
d^3p^{cm}_3=\frac{E^{cm}_2E^{cm}_3{\rm p}^{cm}_3}{M_{23}}dM_{23}d\Omega^{cm}_3,
\end{align}$$
where $M_{23}$($=E^{cm}_2+E^{cm}_3$) is the invariant mass of the $23$ system. We would like to integrate over the magnitude of the neutron momentum $p_1$, which is related to $W_{23}$. Hence, the energy-conserving $\delta$ function is substituted as,
$$\begin{align}
\delta(M_{23}-W_{23})=\frac{W_{23}}{|M{\rm p}_1/E_1|}\delta(\breve{\rm p}_1-{\rm p}_1)
\end{align}$$
where the $\breve{\rm p}_1$ satisfies $M_{23}^2=(M-\breve{E}_1)^2-\breve{\rm p}_1^2$.

Performing the integral over ${\rm p}_1$, we obtain the final expression of the decay width,
$$\begin{align}
d\Phi
&=\frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}^{cm}_3}{8M}d\Omega_1d\Omega^{cm}_3dM_{23}
=\frac{1}{(2\pi)^3}\frac{\breve{\rm p}_1{\rm p}^{cm}_3}{8M}d\cos\theta_1d\cos\theta^{cm}_3dM_{23},
\end{align}$$
Here, independence of the $\phi_1$ and $\phi^{cm}_2$ on the integrand is applied.

## Amplitude
Because the amplitude is covariant, the amplitude for the direct decay can be written as (we omit the helicity of particle 1, which is scalar meson)
$$\begin{align}
{\cal M}^{d}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3)&={\cal A}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3)={\cal A}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(\Omega_1,\Omega^{cm}_3,M_{23})\nonumber\\
%
&=\sum_{J\lambda_R}N_JD^{J*}_{\lambda_R,\lambda_{32}}( \Omega_3^{cm}){\cal A}^J_{\lambda_1;\lambda_2,\lambda_3;\lambda}(\Omega_1,M_{23}),\nonumber\\
%
%
{\cal M}^{Z}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3)&=\int \frac{d^4p'_3}{(2\pi)^4} {\cal T}_{\lambda_2,\lambda_3}(p_2,p_3;p'_2,p'_3)  G(p'_3){\cal A}_{\lambda_1;\lambda}(p_1,p'_2,p'_3)\nonumber\\
%
&=\sum_{\lambda'_2\lambda'_3}\int \frac{d^4p'_3}{(2\pi)^4} {\cal T}_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}(p_2,p_3;p'_2,p'_3)  G_0(p'_3){\cal A}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(p_1,p'_2,p'_3)\nonumber\\
%
&=\sum_{\lambda'_2\lambda'_3}\int \frac{d^3p'^{cm}_3}{(2\pi)^3} i{\cal T}_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}(\Omega_3^{cm},{\rm p}'^{cm}_3,M_{23})  G_0({\rm p}'^{cm}_3){\cal A}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(\Omega_1,{\rm p}'^{cm}_3,\Omega'^{cm}_3,M_{23})\nonumber\\
%
&=\sum_{J\lambda_R}D^{J*}_{\lambda_R,\lambda_{32}}( \Omega_3^{cm})\sum_{\lambda'_2\lambda'_3}\int \frac{{\rm p}'^{cm2}_3d{\rm p}'^{cm}_3d\Omega'^{cm}_3}{(2\pi)^3} N_J^2i{\cal T}^J_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}({\rm p}'^{cm}_3,M_{23})D^{J}_{\lambda_R,\lambda'_{32}}(\Omega'^{cm}_3)  \nonumber\\
&\ \cdot\ G_0({\rm p}'^{cm}_3)\sum_{J'\lambda'_R}N_{J'}D^{J'*}_{\lambda'_R,\lambda'_{32}}( \Omega'^{cm}_3){\cal A}^{J'}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(\Omega_1,{\rm p}'^{cm}_3,M_{23})\nonumber\\
%
&=\sum_{J\lambda_R}N_{J}D^{J*}_{\lambda_R,\lambda_{32}}( \Omega_3^{cm})\sum_{\lambda'_2\lambda'_3}\int \frac{{\rm p}'^{cm2}_3d{\rm p}'^{cm}_3}{(2\pi)^3} i{\cal T}^J_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}({\rm p}'^{cm}_3,M_{23}) G_0({\rm p}'^{cm}_3) {\cal A}^{J}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(\Omega_1,{\rm p}'^{cm}_3,M_{23})\nonumber\\
&\equiv{\cal M}^{Z}_{\lambda_2,\lambda_3;\lambda}({\Omega}_3^{cm},\Omega_1,M_{23}),\nonumber\\
%
%
{\cal M}^{Z}_{\lambda_2;\lambda_1,\lambda_3;\lambda}(p_1,p_2,p_3)%
&=\sum_{J\lambda_R}N_{J}D^{J*}_{\lambda_R,\lambda_{31}}( \tilde{\Omega}_3^{cm})\sum_{\lambda'_1\lambda'_3}\int \frac{\tilde{\rm p}'^{cm2}_3d\tilde{\rm p}'^{cm}_3}{(2\pi)^3} i{\cal T}^J_{\lambda_1,\lambda_3;\lambda'_1,\lambda'_3}(\tilde{\rm p}'^{cm}_3,M_{13}) G_0(\tilde{\rm p}'^{cm}_3) {\cal A}^{J}_{\lambda_2;\lambda'_1,\lambda'_3;\lambda}(\Omega_2,\tilde{\rm p}'^{cm}_3,M_{13})\nonumber\\
&\equiv{\cal M}^{Z}_{\lambda_2;\lambda_1,\lambda_3;\lambda}(\tilde{\Omega}_3^{cm},\Omega_2,M_{13}),\nonumber\\
%
%
{\cal M}^{Z}_{\lambda_3;\lambda_1,\lambda_2;\lambda}(p_1,p_2,p_3)%
&=\sum_{J\lambda_R}N_{J}D^{J*}_{\lambda_R,\lambda_{21}}( \dot{\Omega}_3^{cm})\sum_{\lambda'_1\lambda'_2}\int \frac{\dot{\rm p}'^{cm2}_2d\dot{\rm p}'^{cm}_2}{(2\pi)^3} i{\cal T}^J_{\lambda_1,\lambda_2;\lambda'_1,\lambda'_2}(\dot{\rm p}'^{cm}_2,M_{12}) G_0(\dot{\rm p}'^{cm}_2) {\cal A}^{J}_{\lambda_3;\lambda'_1,\lambda'_2;\lambda}(\Omega_3,\dot{\rm p}'^{cm}_2,M_{12})\nonumber\\
&\equiv{\cal M}^{Z}_{\lambda_3;\lambda_1,\lambda_2;\lambda}(\dot{\Omega}_2^{cm},\Omega_3,M_{12}).
\end{align}$$
where  ${\cal T}^{J}=({\cal T}^{J^-}+{\cal T}^{J^+})/2$ and $N_J=\sqrt{\frac{2J+1}{4\pi}}$. Note: when consider the rescattering of different particles, the different cm frames should be adopted.

In the amplitudes, 
$$\begin{align}
p'^\mu&\to \Lambda^\mu_\nu p^\nu,\nonumber\\
\epsilon'^\mu(p')&\to \Lambda^\mu_\nu \epsilon^\nu(p),\nonumber\\
\end{align}$$

## Decay width

### The case with one rescattering

Here, we consider the rescattering of particles 2 and 3. The distribution can be obtained as
$$\begin{align}
{d\Gamma\over dM_{23}}&=\int \frac{1}{6M}\sum_{\lambda_2,\lambda_3;\lambda}|{\cal M}_{\lambda_2,\lambda_3;\lambda}|^2 \frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}^{cm}_3}{8M}d\Omega_1d\Omega^{cm}_3
\end{align}$$

If the reflection effect is absent, for example, the $\pi DD^*$ final state. The results can be simplified further.
The amplitude can be written as 
$$\begin{align}
{\cal M}_{\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3)
&=\sum_{J\lambda_R}N_J^2D^{J*}_{\lambda_R,\lambda_{32}}( \Omega_3^{cm})D^J_{\lambda_R\lambda}(\Omega_1)\left[\hat{\cal A}^{d,J}_{\lambda_2,\lambda_3;\lambda}(M_{23})\right.\nonumber\\
%
&+\sum_{\lambda'_2\lambda'_3}\int \frac{{\rm p}'^{cm2}_3}{(2\pi)^3} {\cal T}^{{D^-D^{*0}},J}_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}({\rm p}'^{cm}_3,M_{23}) G^{D^-D^{*0}}_0({\rm p}'^{cm}_3) \hat{\cal A}^{{D^-D^{*0}},J}_{\lambda'_2,\lambda'_3;\lambda}({\rm p}'^{cm}_3,M_{23})\nonumber\\
%
&\left.+\sum_{\lambda'_2\lambda'_3}\int \frac{{\rm p}'^{cm2}_3}{(2\pi)^3} {\cal T}^{{D^{*-}D^{0}},J}_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}({\rm p}'^{cm}_3,M_{23}) G^{D^{*-}D^{0}}_0({\rm p}'^{cm}_3) \hat{\cal A}^{{D^{*-}D^{0}},J}_{\lambda'_2,\lambda'_3;\lambda}({\rm p}'^{cm}_3,M_{23})\right]\nonumber\\
&\equiv \sum_{J\lambda_R}N^2_JD^{J*}_{\lambda_R,\lambda_{23}}( \Omega_3^{cm})D^J_{\lambda_R\lambda}(\Omega_1)\hat{\cal M}^J_{\lambda_2,\lambda_3;\lambda}(M_{23})
\end{align}$$

Inserting the above amplitude to the definition of the invariant mass spectrum, we have
$$\begin{align}
{d\Gamma\over dM_{23}}
&=\frac{1}{6M}\frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}^{cm}_3}{8M}\sum_{\lambda_2,\lambda_3;\lambda}\int|\sum_{J\lambda_R}N^2_JD^{J*}_{\lambda_R,\lambda_{32}}( \Omega_3^{cm})D^J_{\lambda_R\lambda}(\Omega_1)\hat{\cal M}^J_{\lambda_2,\lambda_3;\lambda}(M_{23})|^2 d\Omega_1d\Omega^{cm}_3\nonumber\\
&=\frac{1}{6M}\frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}^{cm}_3}{8M}\sum_{\lambda_2,\lambda_3;\lambda;J}|\hat{\cal M}^J_{\lambda_2,\lambda_3;\lambda}(M_{23})|^2 \nonumber\\
\end{align}$$

Now we consider the amplitude with fixed parity,
$$\begin{align}
{\cal M}^{J}_{\lambda_{23};\lambda}&={\cal A}^{d,J}_{\lambda_{23};\lambda}+{\cal T}^{J}_{\lambda_{23},\lambda'_{23}}G_0{\cal A}^J_{\lambda'_{23};\lambda},\quad
\eta{\cal M}^{J}_{-\lambda_{23};\lambda}=\eta{\cal A}^{d,J}_{-\lambda_{23};\lambda}+\eta{\cal T}^{J}_{-\lambda_{23},\lambda'_{23}}G_0{\cal A}^J_{\lambda'_{23};\lambda}\nonumber\\
\Rightarrow {\cal M}^{J^P}_{\lambda_{23};\lambda}&={\cal A}^{d,J^P}_{\lambda_{23};\lambda}+{\cal T}^{J^P}_{\lambda_{23},\lambda'_{23}}G_0{\cal A}^{J}_{\lambda'_{23}
;\lambda}={\cal A}^{d,J^P}_{\lambda_{23};\lambda}+\eta'{\cal T}^{J^P}_{\lambda_{23},\lambda'_{23}}G_0{\cal A}^{J}_{-\lambda'_{23}
;\lambda}={\cal A}^{d,J^P}_{\lambda_{23};\lambda}+{1\over 2}{\cal T}^{J^P}_{\lambda_{23},\lambda'_{23}}G_0{\cal A}^{J^P}_{\lambda'_{23}
;\lambda}
\end{align}$$

Here ${\cal M}^{J^P}_{\lambda_{23};\lambda}={\cal M}^{J}_{\lambda_{23};\lambda}+\eta{\cal M}^{J}_{\lambda_{23};\lambda}$, ${\cal A}^{d,J^P}_{\lambda_{23};\lambda}={\cal A}^{d,J}_{\lambda_{23};\lambda}+\eta{\cal A}^{d,J}_{\lambda_{23};\lambda}$, ${\cal T}^{J^P}_{\lambda_{23},\lambda'_{23}}={\cal T}^{J}_{\lambda_{23},\lambda'_{23}}+\eta{\cal T}^{J}_{-\lambda_{23},\lambda'_{23}}={\cal T}^{J}_{\lambda_{23},\lambda'_{23}}+\eta'{\cal T}^{J}_{\lambda_{23},-\lambda'_{23}}$.
NOTE: For parity conserving interaction, we have ${\cal T}^J_{\lambda'\lambda}=\eta(\eta')^{-1}{\cal T}_{-\lambda'\lambda}$, which can be checked in the code by different definitions of ${\cal T}^{J^P}$.   


We summarize the results as following,
$$\begin{align}
{d\Gamma\over dM_{23}}&=\frac{1}{6M}\frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}^{cm}_3}{8M}\sum_{i\ge0;j\ge0;J^P}\frac{1}{N_J^2}|\hat{\cal M}^{J^P}_{i;j}(M_{23})|^2, \nonumber\\
\hat{\cal M}^{J^P}_{i;i}(M_{23})
&=\hat{\cal A}^{d,J^P}_{j;i}(M_{23})+\sum_{k}\int \frac{d{\rm p}'^{cm}_3{\rm p}'^{cm2}_3}{(2\pi)^3}\hat{\cal T}^{J^P}_{j;k}({\rm p}'^{cm}_3,M_{23}) G_0({\rm p}'^{cm}_3) \hat{\cal A}^{J^P}_{k;i}({\rm p}'^{cm}_3,M_{23}),
\end{align}$$
where $i$ and $j$ denote the independent $\lambda_{2,3}$ and $\lambda$, and the factors $f_{i=0}=1/\sqrt{2}$ and $f_{i\neq0}=1$ are inserted. The above equation can be abbreviated as $M={A}^{d}+{T}G { A}$, where  $T$ is solved by the Bethe-Salpeter equation $T=V+VGT$. NOTE: The amplitudes obtained by equation $T=V+TGV$  is the same, which can be checked by the code.

### The case with more than one rescattering

In the case with more than one resacttering, we should conisder Monte-Carlo method to generate the 
event. 

$$\begin{align}
d\Gamma=\frac{1}{2E}\sum|{\cal M}|^2 d\Phi=\frac{1}{2E}\sum|{\cal M}|^2 (2\pi)^{4-3n} dR
\end{align}$$

Here we consider a process with direct, 23 rescattering, 13 rescattering.
$$\begin{align}
{\cal M}^{d}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3)&={\cal A}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3),\nonumber\\
%
%
{\cal M}^{Z}_{\lambda_k;\lambda_i,\lambda_j;\lambda}(p_k,p_i,p_j)
&=\sum_{J\lambda_R}N_{J}D^{J*}_{\lambda_R,\lambda_{ji}}( \Omega_j^{cm})\sum_{\lambda'_i\lambda'_j}\int \frac{{\rm p}'^{cm2}_jd{\rm p}'^{cm}_j}{(2\pi)^3} i{\cal T}^J_{\lambda_i,\lambda_j;\lambda'_i,\lambda'_j}({\rm p}'^{cm}_j,M_{ij}) \nonumber\\
&\  \cdot \ G_0({\rm p}'^{cm}_j) {\cal A}^{J}_{\lambda_k;\lambda'_i,\lambda'_j;\lambda}(\Omega_k,{\rm p}'^{cm}_j,M_{ij})
\end{align}$$

$$
\begin{align}
&{\cal A}^J_{\lambda_k;\lambda'_i,\lambda'_j;\lambda}(\Omega_k,{\rm p}'^{cm}_j,M_{ij})=
{\frac{2J+1}{4\pi}}\int d\Omega'^{cm}_jD^{J}_{\lambda,\lambda_{ji}}( \Omega'^{cm}_j) {\cal A}_{\lambda_k;\lambda'_i,\lambda'_j;\lambda}(\Omega_k,\Omega'^{cm}_j,{\rm p}'^{cm}_j,M_{ij})
\end{align}
$$

- The GEN pacakge can generate the events with $(p_1, p_2, p_3)$ in the frame of $Y$.
- The $M_{ij}$ can be obtained as $(p_i+p_j)^2$, and used to solve the ${\cal T}^{J^P}_{\lambda_i,\lambda_j;\lambda'_i,\lambda'_j}({\rm p}'^{cm}_j,M_{ij})$ with qBSE. Here, the ${\rm p}'^{cm}_j$ is generated by the Gauss discretization.
- With Lorentz boost, the momentum of $p^{cm}$ can be calculated, which can further used to calculated $D^{J*}_{\lambda_R,\lambda_{ji}}( \Omega_j^{cm})$. 
- With the $p^{cm}$, ${\cal A}^{J}_{\lambda_k;\lambda'_i,\lambda'_j;\lambda}(\Omega_k,{\rm p}'^{cm}_j,M_{ij})$ can be calcualted from ${\cal A}_{\lambda_k;\lambda'_i,\lambda'_j;\lambda}(\Omega_k,\Omega'^{cm}_j,{\rm p}'^{cm}_j,M_{ij})$.

$$
\begin{align}
	\sum_{\lambda'}
	T^{J}_{\lambda\lambda'}A^{J}_{\lambda'}&=T^{J}_{\lambda0}A^{J}_{0}+\sum_{j>0}
	\left[T^{J}_{\lambda j}A^{J}_{j}+T^{J}_{\lambda -j}A^{J}_{-j}\right]\nonumber\\
&=\frac{1}{4}T^{J^P}_{\lambda0}A^{J^P}_{0}+\sum_{j>0}
	 \left[\frac{1}{4}(T^{J+}_{\lambda j}+T^{J-}_{\lambda j})(A^{J+}_{j}+A^{J-}_{j})
	 +\frac{1}{4}(T^{J+}_{\lambda j}-T^{J-}_{\lambda j})(A^{J+}_{j}-A^{J-}_{j})\right]\nonumber\\
&=\frac{1}{4}T^{J^P}_{\lambda0}A^{J^P}_{0}+\sum_{j>0}
	 \left[\frac{1}{2}T^{J+}_{\lambda j}A^{J+}_{j}	
	  +\frac{1}{2}T^{J-}_{\lambda j}A^{J-}_{j}\right]\nonumber\\
&=\sum_P T^{J^P}_{\lambda}A^{J^P}\nonumber\\
\end{align}
$$
$$\begin{align}
{\cal M}^{Z}_{\lambda}
&\propto  N_JD^{J*}_{\lambda_R,\lambda_{ji}}T^{J^P}_{\lambda}A^{J^P}\nonumber\\
&=N_JD^{J*}_{\lambda_R,0}T^{J^P}_{0}A^{J^P}+\sum_{l>0}
( N_JD^{J*}_{\lambda_R,l_{ji}}T^{J^P}_{l}A^{J^P}+N_JD^{J*}_{\lambda_R,-l_{ji}}\eta T^{J^P}_{l}A^{J^P})
\nonumber\\
&=N_JD^{J*}_{\lambda_R,0}T^{J^P}_{0}A^{J^P}+\sum_{l>0}
N_J( D^{J*}_{\lambda_R,l_{ji}}+\eta D^{J*}_{\lambda_R,-l_{ji}} )T^{J^P}_{l}A^{J^P}
\nonumber\\
&=\sum_l D^{J^P}_{\lambda_R l}T^{J^P}_{l}A^{J^P}
\end{align}$$
$$\begin{align}
{\cal M}^{Z}_{\lambda_k;\lambda_i,\lambda_j;\lambda}(p_k,p_i,p_j)
&=\sum_{J^P\lambda_R}N_{J}D^{J*}_{\lambda_R,\lambda_{ji}}( \Omega_j^{cm}){T}^{J^P}_{\lambda_i,\lambda_j;} G_0 A^{J^P}_{\lambda_k;\lambda}
\end{align}$$