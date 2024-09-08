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
D^{J*}_{\lambda_R,\lambda'}(\phi',\theta',-\phi'){\cal V}_{\lambda'\lambda}({\bm k}',{\bm k})D^{J}_{\lambda_R,\lambda}(\phi,\theta,-\phi)
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
