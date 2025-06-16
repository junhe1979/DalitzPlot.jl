# Outline

<!-- toc -->

- [Outline](#outline)
- [Quasipotential approximation](#quasipotential-approximation)
  - [Partial-wave expansion](#partial-wave-expansion)
  - [Fixed parity](#fixed-parity)
  - [Transformation to a matrix equation](#transformation-to-a-matrix-equation)
  - [Code](#code)
  - [Pole search](#pole-search)
  - [Physical observable](#physical-observable)
    - [The cross section for the channel considered](#the-cross-section-for-the-channel-considered)
    - [Argand plot](#argand-plot)
- [Three body decay](#three-body-decay)
  - [kinematics](#kinematics)
    - [Lorentz boost](#lorentz-boost)
  - [Amplitude](#amplitude)
  - [Decay width](#decay-width)
    - [The case with one rescattering](#the-case-with-one-rescattering)
      - [Phase space](#phase-space)
      - [Differential decay width](#differential-decay-width)
    - [The case with more than one rescattering](#the-case-with-more-than-one-rescattering)
- [qBSE package](#qbse-package)
  - [Data Structures for the Dimensions](#data-structures-for-the-dimensions)
    - [`structChannel`](#structchannel)
    - [`structIndependentHelicity`](#structindependenthelicity)
    - [`structDimension`](#structdimension)
  - [Data Structures for the Interactions](#data-structures-for-the-interactions)
    - [`structSys`](#structsys)
    - [`structInterAction`](#structinteraction)
  - [Additional Data Structures](#additional-data-structures)
    - [`structMomentum`](#structmomentum)
    - [`structHelicity`](#structhelicity)
    - [`structParticle`](#structparticle)
  - [Functions for the qBSE](#functions-for-the-qbse)
    - [`function preprocessing(Sys, channels, Ff, qn; Np=10, Nx=5, Nphi=5)`](#function-preprocessingsys-channels-ff-qn-np10-nx5-nphi5)
    - [`function res(Range, iER, qn, SYS, IA, CH, IH, fV)`](#function-resrange-ier-qn-sys-ia-ch-ih-fv)
    - [`function fV(k, l, SYS, IA0, CHf, CHi)`](#function-fvk-l-sys-ia0-chf-chi)
    - [`function propFF(k, ex, L, LLi, LLf; lregu=1, lFFex=0)`](#function-propffk-ex-l-lli-llf-lregu1-lffex0)
    - [`function simpleXsection(ER, resM2, CH, qn; Ep="cm")`](#function-simplexsectioner-resm2-ch-qn-epcm)
    - [`function lambda(m1, m2, m3)`](#function-lambdam1-m2-m3)
  - [Decay](#decay)
    - [`function setTGA(par, sij, k, tecm, i, j)`](#function-settgapar-sij-k-tecm-i-j)
    - [`TGA(cfinal, cinter, Vert, para)`](#tgacfinal-cinter-vert-para)
    - [`function Vertex14(p1, p2, l1, l2, Vert, k, P)`](#function-vertex14p1-p2-l1-l2-vert-k-p)
    - [auxiliary function](#auxiliary-function)
  - [Additional functions](#additional-functions)
    - [`function ch(pf, pin, amps)`](#function-chpf-pin-amps)
    - [`function LorentzBoostRotation(k, tecm, p1, p2)`](#function-lorentzboostrotationk-tecm-p1-p2)

<!-- tocstop -->

# Quasipotential approximation

The general form of the Bethe-Salpeter equation (BSE) for the scattering amplitude can be written as follows:

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

where ${\cal V}$ is the potential kernel and $G$ is the propagator for the two constituent particles. The total momentum of the system is denoted by $P = k_1 + k_2 = k'_1 + k'_2 = k''_1 + k''_2$.

The Bethe-Salpeter equation can be succinctly expressed as

$$
{\cal M} = {\cal V} + {\cal V} G {\cal M},
$$

The Bethe-Salpeter (BS) equation for the amputated scattering matrix with external legs, denoted as ${\cal M}_{[\mu_f][\mu_i]}$, is given by:

$$
\begin{align}
	{\cal M}_{[\mu_f][\mu_i]}={\cal V}_{[\mu_f][\mu_i]}
	+{\cal V}_{[\mu_f][\mu']}G^{[\mu'][\mu]}{\cal
	M}_{[\mu][\mu_i]},
\end{align}
$$

where the propagator for the two constituent particles is given by

$$
\begin{align}
G^{[\mu'][\mu]} = G_1^{[\mu'_1][\mu_1]} \, G_2^{[\mu'_2][\mu_2]}
= \frac{-P^{[\mu'][\mu]}}{(k_1^2 - m_1^2)(k_2^2 - m_2^2)}
= P^{[\mu'][\mu]} \, \tilde{G}_0,
\end{align}
$$

where $P^{[\mu'][\mu]}$ is a general tensor structure. For example, for two vector mesons, $P^{[\mu'][\mu]} = (-g^{\mu'_1\mu_1} + k^{\mu'_1}_1 k^{\mu_1}_1 / m_1^2)(-g^{\mu'_2\mu_2} + k^{\mu'_2}_2 k^{\mu_2}_2 / m_2^2)$; for two spin-$1/2$ baryons, $P^{[\mu'][\mu]} = (\gamma \cdot k_1 + m_1)(\gamma \cdot k_2 + m_2)$.
In general, such propagators are difficult to handle because the potential ${\cal V}$ and amplitude ${\cal M}$ cannot be factorized. However, since the dominant contribution comes from the region where both constituents are near their mass shells, a form factor or cutoff is usually introduced to restrict the propagator to the near on-shell region. Therefore, it is reasonable to approximate $P^{[\mu'][\mu]}$ by its on-shell value, which can be expressed as a sum over polarization vectors or spinors: $P^{[\mu'][\mu]} \approx \sum_{\lambda_1, \lambda_2} A_{\lambda_1}^{[\mu'_1]} A_{\lambda_1}^{[\mu_1]} A_{\lambda_2}^{[\mu'_2]} A_{\lambda_2}^{[\mu_2]}$, where $A$ denotes the polarization vector or spinor.

The Gross form of proposed quasipotential propagators for particles 1 and 2 with mass $m_1$ and $m_2$ written down in
the center of mass frame where $P=(W,{\boldsymbol 0})$ with particle 2 being on shell are

$$
\begin{align}
\tilde{G}_0=\frac{1}{(k_1^2-m_1^2)(k_2^2-m_2^2)}\to g=2\pi i\frac{\delta^+(k_2^2-m_2^2)}{k_1^2-m_1^2}=2\pi
i\frac{\delta(k^0_2-E_2)}{2E_2[(W-E_2)^2-E_1^2]},
\end{align}
$$

where $k_1=(k_1^0,\boldsymbol k)=(E_1,\boldsymbol k)$, $k_2=(k_2^0,-\boldsymbol k)=(W-E_1,-\boldsymbol k)$ with $E_1=\sqrt{m_1^2+|\boldsymbol k|^2}$.

With the define of $G_0=g/(2\pi i)$, the four-dimensional BSE can be reduced to a three-dimensional equation in center of mass frame

$$
\begin{align}
i{\cal M}({\boldsymbol k}',{\boldsymbol k})&=i{\cal
V}({\boldsymbol k}',{\boldsymbol k})+\int\frac{d
{\boldsymbol k}''}{(2\pi)^3}
i{\cal
V}({\boldsymbol k}',{\boldsymbol k}'')G_0({\boldsymbol k}'')i{\cal
M}({\boldsymbol k}'',{\boldsymbol k}),\quad
\end{align}
$$

Note: the $i{\cal M}$ and $i{\cal V}$ are usually real. In the center of mass frame. We choose ${\boldsymbol k}_2={\boldsymbol k}$ and ${\boldsymbol k}_1=-{\boldsymbol k}$.

## Partial-wave expansion

To reduce the equation to one-dimensional equation, we apply the partial wave expansion,

$$
\begin{align}
{\cal V}_{\lambda'\lambda}({\boldsymbol k}',{\boldsymbol k})&\equiv\langle\theta'\phi',\lambda'|{\cal V}|\theta\phi,\lambda\rangle
=\sum_{J'M'JM}\langle\theta'\phi'|J'M',\lambda'\rangle\langle J'M'|{\cal V}|JM,\lambda\rangle\langle JM,\lambda|\theta\phi,\lambda\rangle
\nonumber\\
&=\sum_{JM }N_J^2D^{J*}_{M\lambda'}(\phi',\theta',0){\cal
	V}^{JM}_{\lambda'\lambda}({\rm k}',{\rm k})D^{J}_{M\lambda}(\phi,\theta,0)
=\sum_{J }N_J^2d^{J}_{\lambda\lambda'}(\theta_{k,k'}){\cal
	V}^{J}_{\lambda'\lambda}({\rm k}',{\rm k})
	.
\nonumber\\
{\cal V}_{\lambda'\lambda}^{JM}({\rm k}',{\rm k})&=
N^2_J\int d\Omega' d\Omega
D^{J}_{M,\lambda'}(\phi',\theta',0){\cal V}_{\lambda'\lambda}({\boldsymbol k}',{\boldsymbol k})D^{J*}_{M,\lambda}(\phi,\theta,0).
\end{align}
$$

where $N_J=\sqrt{\frac{2J+1}{4\pi}}$, $\int d\Omega D^{J*}_{\lambda_1,\lambda_2}(\phi,\theta,0)D^{J'*}_{\lambda'_1,\lambda'_2}(\phi,\theta,0)=N_J^{-2}$.

To calculate ${\cal V}_{\lambda'\lambda}^{JM}({\rm k}',{\rm k})$, we adopt a special CMS frame. The momenta are chosen as $k_2=(E_2,0,0,{\rm k})$, $k_1=(W-E_2,0,0,-{\rm k})$  and $k'_2=(E'_2,{\rm k}'\sin\theta_{k,k'},0,{\rm k}'\cos\theta_{k,k'})$, $k=(W-E_2, -{\rm k}'\sin\theta_{k,k'},0,-{\rm k}'\cos\theta_{k,k'})$ with ${\rm k}=|{\boldsymbol k}|$ and ${\rm k}'=|{\boldsymbol k}'|$.

$$
\begin{align}
\to{\cal V}_{\lambda'\lambda}({\boldsymbol k}',{\boldsymbol k})&=
	\sum_{JM}N^2_JD^{J*}_{M\lambda'}(0,\theta_{k',k},0){\cal
	V}^{JM}_{\lambda'\lambda}({\rm k}',{\rm k})D^{J}_{M\lambda}(0,0,0).
\nonumber\\
&=\sum_{JM}N^2_Jd^{J}_{M\lambda'}(\theta_{k',k}){\cal
	V}^{JM}_{\lambda'\lambda}({\rm k}',{\rm k})\delta_{M\lambda}
	=\sum_{J}N^2_Jd^{J}_{\lambda,\lambda'}(\theta_{k',k}){\cal
	V}^J_{\lambda'\lambda}({\rm k}',{\rm k})
\nonumber\\
{\cal V}_{\lambda'\lambda}^J({\rm k}',{\rm k})&=2\pi\int d\cos\theta_{k,k'} d^{J}_{\lambda,\lambda'}(\theta_{k',k})
{\cal V}_{\lambda'\lambda}({\boldsymbol k}',{\boldsymbol k}).
\end{align}
$$

where $\int^1_{-1}d\cos\theta' d^{J'}_{\lambda,\lambda'}(\theta')d^{J}_{\lambda,\lambda'}(\theta')=N_J^{-2}/2\pi\delta_{JJ'}$ are used.

NOTE: Which particle is chosen to parallel to $z$ axis is related to the order of $\lambda$ and $\lambda'$ in $d^{J}_{\lambda'\lambda}(\theta_{k,k'})$, so it can not be chosen arbitrarily. And the definition of helicity is also dependent of the definition of ${\boldsymbol k}_{1,2}$. Here, $\lambda=\lambda_2-\lambda_1$  and $\lambda_1=-s_1$, $\lambda_2=s_2$. The scattering amplitudes ${\cal M}$ has analogous relations.

Now we have the partial wave BS equation,

$$
\begin{align}
i{\cal M}^J_{\lambda',\lambda}({\rm k}',{\rm k})&=
N^2_J\int d\Omega' d\Omega D^{J*}_{\lambda_R,\lambda'}(\phi',\theta',0)D^{J}_{\lambda_R,\lambda}(\phi,\theta,0)\nonumber\\
&\cdot\left[i{\cal V}_{\lambda'\lambda}({\boldsymbol k}',{\boldsymbol
k})+\int\frac{d{\boldsymbol k}''}{(2\pi)^3}i{\cal V}_{\lambda'\lambda''}({\boldsymbol
k}',{\boldsymbol k}'') G_0({\boldsymbol k}'')i{\cal M}_{\lambda''\lambda}({\boldsymbol k}'',{\boldsymbol k})\right]\nonumber\\
&=i{\cal V}^J_{\lambda'\lambda}({\rm k}',{\rm k})+\int\frac{{\rm k}''^2d{\rm k}''}{(2\pi)^3}i{\cal V}^J_{\lambda'\lambda''}({\rm k}',{\rm k}'')
G_0({\boldsymbol k}'')i{\cal M}^J_{\lambda''\lambda}({\rm k}'',{\rm k})
\end{align}
$$

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
{\cal M}^{J}_{-\lambda'-\lambda}=\tilde{\eta}\tilde{\eta}'{\cal M}^{J}_{\lambda'\lambda}
\end{align}
$$

$$
\begin{align}
	&{\cal M}^{J\pm}_{\lambda'\lambda}=\langle J,\lambda';\pm|{\cal M}|J,\lambda;\pm\rangle
	={\cal M}^{J}_{\lambda'\lambda}\pm \tilde{\eta} {\cal M}^J_{\lambda'-\lambda}
	={\cal M}^{J}_{\lambda'\lambda}\pm \tilde{\eta}'{\cal M}^J_{-\lambda'\lambda},\\
&{\cal M}^{J\pm}_{\lambda'-\lambda}=\pm\tilde{\eta} {\cal M}^{J\pm}_{\lambda'\lambda}\equiv\eta {\cal M}^{J\pm}_{\lambda'\lambda},\ \ {\cal M}^{J\pm}_{-\lambda'\lambda}=\pm\tilde{\eta}' {\cal M}^{J\pm}_{\lambda'\lambda}\equiv\eta' {\cal M}^{J\pm}_{\lambda'\lambda}
\end{align}
$$

with $\eta=PP_1P_2(-1)^{J_1+J_2-J}=P(-1)^{1/2+J}$.
The potential ${\cal V}^{J^P}_{\lambda'\lambda}$ has analogous relations.

$$
\begin{align}
	 i{\cal M}^J_{\lambda\lambda'}&=i{\cal V}^J_{\lambda\lambda'}+\sum_{\lambda''}i{\cal V}^J_{\lambda\lambda''}G_0i{\cal M}^J_{\lambda''\lambda'},\quad	 \eta'i{\cal M}^J_{\lambda-\lambda'}=\eta'i{\cal V}^J_{\lambda-\lambda'}+\sum_{\lambda''}i{\cal V}^J_{\lambda\lambda''}G_0\eta'i{\cal M}^J_{\lambda''-\lambda'}\nonumber\\
\Rightarrow i{\cal M}^{J^P}_{\lambda\lambda'}&=i{\cal V}^{J^P}_{\lambda\lambda'}+\sum_{\lambda''}i{\cal V}^J_{\lambda\lambda''}G_0i{\cal M}^{J^P}_{\lambda''\lambda'},\nonumber\\
i{\cal M}^{J^P}_{\lambda\lambda'}&=i{\cal V}^{J^P}_{\lambda\lambda'}+\sum_{\lambda''}i{\cal V}^J_{\lambda-\lambda''}G_0i{\cal M}^{J^P}_{-\lambda''\lambda'}
=iV^{J^P}_{\lambda\lambda'}+\sum_{\lambda''}i{\cal V}^J_{\lambda-\lambda''}G_0\eta''i{\cal M}^{J^P}_{-\lambda''\lambda'}\nonumber\\
\Rightarrow  i{\cal M}^{J^P}_{\lambda\lambda'}&=i{\cal V}^{J^P}_{\lambda\lambda'}+\frac{1}{2}\sum_{\lambda''}i{\cal V}^{J^P}_{\lambda\lambda''}G_0i{\cal M}^{J^P}_{\lambda''\lambda'},
\end{align}
$$

As shown in Eq. (14), the amplitudes are not independent. If we only keep the independent amplitudes, the equation for definite parity can be written as

$$
\begin{align}
&i{\cal M}^{J^P}_{ij}=i{\cal V}^{J^P}_{ij}+\frac{1}{2}\sum_{\lambda''}i{\cal V}^{J^P}_{i\lambda''}G_0i{\cal M}^{J^P}_{\lambda''j}
=i{\cal V}^{J^P}_{ij}+\frac{1}{2}i{\cal V}^{J^P}_{i0}Gi{\cal M}^{J^P}_{0j}+\sum_{k\neq0}i{\cal V}^{J^P}_{ik}G_0i{\cal M}^{J^P}_{kj},\nonumber\\
\Rightarrow&i{{\cal M}}^{J^P}_{ij}=i{\cal V}^{J^P}_{ij}+\sum_{k}i{\cal V}^{J^P}_{ik}G_0i{\cal M}^{J^P}_{kj}.
\end{align}
$$

where $i$, $j$, $k$ are the indices for the independent amplitudes. And we redefine

$$
f_if_j {\cal M}^{J^P}_{\lambda'\lambda}\to{\cal M}^{J^P}_{ij}
$$

with $f_0=\frac{1}{\sqrt{2}}$ and $f_{i\neq0}=1$ with $0$ for the amplitudes with $\lambda_1=\lambda_2=0$.

The Bethe-Saltpeter equation for partial-wave amplitude with fixed spin-parity $J^P$ reads ,

$$
\begin{align}
i{\cal M}^{J^P}_{\lambda'\lambda}({\rm k}',{\rm k})
&=i{\cal V}^{J^P}_{\lambda',\lambda}({\rm k}',{\rm
k})+\sum_{\lambda''}\int\frac{{\rm
k}''^2d{\rm k}''}{(2\pi)^3}~
i{\cal V}^{J^P}_{\lambda'\lambda''}({\rm k}',{\rm k}'')
G_0({\rm k}'')i{\cal M}^{J^P}_{\lambda''\lambda}({\rm k}'',{\rm
k}).
\end{align}
$$

Note that here the sum extends only over indices for the independent amplitudes. The partial wave potential is defined as

$$
\begin{align}
i{\cal V}_{\lambda'\lambda}^{J^P}({\rm k}',{\rm k})
&=2\pi\int d\cos\theta
~[d^{J}_{\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'\lambda}({\boldsymbol k}',{\boldsymbol k})
+\eta d^{J}_{-\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'-\lambda}({\boldsymbol k}',{\boldsymbol k})],
\end{align}
$$

Or, with the independent amplitudes as

$$
\begin{align}
i{\cal M}^{J^P}_{ij}({\rm k}',{\rm k})
&=i{\cal V}^{J^P}_{ij}({\rm k}',{\rm
k})+\sum_{k}\int\frac{{\rm
k}''^2d{\rm k}''}{(2\pi)^3}~
i{\cal V}^{J^P}_{ik}({\rm k}',{\rm k}'')
G_0({\rm k}'')i{\cal M}^{J^P}_{kj}({\rm k}'',{\rm
k}).
\end{align}
$$

$$
\begin{align}
i{\cal V}^{J^P}_{ij}({\rm k}',{\rm k}'')
&=f_if_j2\pi\int d\cos\theta
~[d^{J}_{\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'\lambda}({\boldsymbol k}',{\boldsymbol k})
+\eta d^{J}_{-\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'-\lambda}({\boldsymbol k}',{\boldsymbol k})].
\end{align}
$$

Note here the $f_if_j$ is also incorporated. Additionally, the form factors for the interacting particles are also included in the potential, modifying it as ${\cal V}\to f(k'){\cal V}f(k)$. Consequently, the resulting amplitude ${\cal M}$ also includes these form factors.

The potential ${\cal V}_{\lambda'-\lambda}({\boldsymbol k}',{\boldsymbol k})$ is introduced by `fV` function in main file as

```julia
fV(k, l, SYS, IA0, CHf, CHi) 
```

where `k` and `l` are for the momenta and helicities of final and initial particles. `SYS` is
for the system information. `IA0` is for the interaction information., `CHf` and `CHi` are for the information of final and initial channels.

Transition of ${\cal V}_{\lambda'-\lambda}({\boldsymbol k}',{\boldsymbol k})$ to ${\cal V}^{J^P}_{ij}({\rm k}',{\rm k}'')$ performed in `qBSE.fKernel` which is an internal function.

## Transformation to a matrix equation

Now We have a integral equation with singularity in $G_0=\frac{1}{2 E_2[(W-E_2)^2-E_1^2]}=\frac{1}{2 E_2[(W-E_2-E_1+i\epsilon)(W-E_2+E_1)]}$  at $W=E_1+E_2$. This singularity can be isolated as,

$$
\begin{align}
i{\cal M}^{J^P}({\rm k}',{\rm k})
&=i{\cal V}^{J^P}({\rm k},{\rm k}')+\int\frac{{\rm k}''^2d {\rm k}''}{(2\pi)^3}i{\cal V}^{J^P}({\rm k},{\rm k}'') G_0({\rm k}'')i{\cal M}^{J^P}({\rm k}'', {\rm k}')\nonumber\\
&=i{\cal V}^{J^P}+\mathcal{P}\int\frac{{\rm k}''^2d {\rm k}''}{(2\pi)^3}i{\cal V}^{J^P}
G_0{\cal M}^{J^P}-\pi i\int\frac{{\rm k}''^2d {\rm k}''}{(2\pi)^3}{\cal V}^{J^P}_o\delta G_0{\cal M}^{J^P}_o\nonumber\\
&=i{\cal V}^{J^P}+\int\frac{d {\rm k}''}{(2\pi)^3}\left[{\rm k}''^2 G_0 i{\cal V}^{J^P} {\cal M}^{J^P}-\frac{M{\cal V}^{J^P}_o{\cal M}^{J^P}_o}{{\rm k}''^2-\bar{{\rm k}''}^2}\right]
-i\frac{\bar{{\rm k}''}^2\delta\bar{G}_0}{8\pi^2}{\cal V}^{J^P}_o{\cal M}^{J^P}_o\nonumber\\
&=i{\cal V}^{J^P}+\int\frac{d {\rm k}''}{(2\pi)^3}{\rm k}''^2 G_0 i{\cal V}^{J^P} {\cal M}^{J^P}
-[\int\frac{d {\rm k}''}{(2\pi)^3}\frac{M}{{\rm k}''^2-\bar{{\rm k}''}^2}+i\frac{\bar{{\rm k}''}^2\delta\bar{G}_0}{8\pi^2}]{\cal V}^{J^P}_o{\cal M}^{J^P}_o\nonumber\\
\end{align}
$$

with $\delta G_0=\frac{\delta(W-E_2-E_1)}{2E(W-E_2+E_1)}\theta(s-m_1-m_2)=\delta \bar{G}_0\delta({\rm k}''-\bar{{\rm k}''})=\frac{1}{4W\bar{{\rm k}''}}\delta({\rm k}''-\bar{{\rm k}''})\theta(s-m_1-m_2)$, $M=[{\rm k}''^2({\rm k}''^2-\bar{{\rm k}''}^2)G_0]_{{\rm k}''\to\bar{{\rm k}''}}\theta(s-m_1-m_2)=-\frac{\bar{{\rm k}''}^2}{2W}\theta(s-m_1-m_2)$.

We have

$$
\begin{align}
Im~G=-\rho/2=-\frac{\bar{{\rm k}''}}{32\pi^2 W}.
\end{align}
$$

It should be noted that in the region $W < m_1 + m_2$, the potential remains real, since no imaginary component emerges in this energy range.

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
i{\cal M}^{J^P}_{ik}
&=&i{\cal V}^{J^P}_{ik}+\sum_{j=0}^N i{\cal V}^{J^P}_{ij}G_ji{\cal M}^{J^P}_{jk}\Rightarrow {M}^{J^P}={V}^{J^P}+{V}^{J^P}G{M}^{J^P}
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

The propagator is calculated in `qBSE.fProp` which is an internal function.

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

The informations about the dimensions are calculated in `qBSE.WORKSPACE`. The matrix $V$ and $G$ are calculated in `qBSE.srAB`.
Note that such function is internal function and not used by the user.

## Code

**Attention**: The following details are specific to the old version of the code, which includes Fortran code and Julia code versions prior to v0.2.4. In the new version, the treatment described below is obsolete.

In old code, we choose $\hat{V}^{J^P}=f(k'){V}^{J^P}f(k)/4\pi$, $\hat{G}=4\pi{G}$, and $\hat{M}^{J^P}=f(k'){M}^{J^P}f(k)/4\pi$. The form factors are also included in to the ampltudes and the potential kernel. Hence, the qBSE becomes

$$
\begin{align}
\hat{M}^{J^P}=\hat{V}^{J^P}+\hat{V}^{J^P}G\hat{M}^{J^P}.
\end{align}
$$

Such convention is consistent with that in the chiral unitary approach.

$$
\begin{align}
\hat{V}^{J^P}&={V}^{J^P}/4\pi=i{\cal V}^{J^P}_{\lambda'\lambda''}({\rm p}',{\rm p}'')/4\pi=f_if_j i{\cal V}_{\lambda'\lambda}^{J^P}({\rm p}',{\rm p})/4\pi \nonumber\\
&=\frac{1}{2}f_if_j \int d\cos\theta
~[d^{J}_{\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'\lambda}({\boldsymbol p}',{\boldsymbol p})
+\eta d^{J}_{-\lambda\lambda'}(\theta)
i{\cal V}_{\lambda'-\lambda}({\boldsymbol p}',{\boldsymbol p})],
\end{align}
$$

$$
\begin{align}
	\hat{G}_j=\left\{\begin{array}{cl}-\frac{i\bar{q}}{8\pi W}+\sum_j
\left[\frac{w(q_j)}{2\pi^2}\frac{\bar{q}^2}
{2W{(q_j^2-\bar{q}^2)}}\right] & {\rm for}\ j=0,\ {\rm if}\ Re(W)>m_1+m_2,\nonumber\\
\frac{w(q_j)}{2\pi^2}\frac{q_j^2}
	{2E(q_j)[(W-E(q_j))^2-\omega^2(q_j)]}& {\rm for}\ j\neq0
	\end{array}\right.
\end{align}
$$

## Pole search

After extend the energy in the center of mass frame $W$ into complex energy plane as $z$, the pole can be found by variation of $z$ to satisfy

$$
\begin{align}
|1-V(z)G(z)|=0
\end{align}
$$

with $z=E_R+i\Gamma_R/2$.

The $|1-V(z)G(z)|$ in complex enery plane is calculated in `qBSE.res` and the pole can be found with function `qBSE.showPoleInfo`.

As shown above, the propagator acquires an imaginary part in the region $W > m_1 + m_2$, specifically $-i\rho/2$. According to the Schwarz reflection principle, if a function $f(z)$ is analytic in a region of the complex plane that includes a portion of the real axis where $f$ is real, then it satisfies $[f(z^*)]^* = f(z)$. Since the propagator $G$ satisfies these conditions, we have, for ${\rm Re}(z) > m_1 + m_2$,

$$
G(z - i\epsilon) = G^*(z + i\epsilon) = G(z + i\epsilon) + i\rho.
$$

The value of the propagator at the lower edge of the branch cut on the first Riemann sheet, $G^I(z - i\epsilon)$, coincides with the upper edge of the second Riemann sheet, $G^{II}(z + i\epsilon)$. Therefore,

$$
G^{II}(z + i\epsilon) = G^I(z - i\epsilon) = G^I(z + i\epsilon) + i\rho.
$$

Before analytically continuing into the complex plane, we adopt the convention that on the first Riemann sheet, the propagator is real for $W < m_1 + m_2$, and has an imaginary part $-i\rho/2$ for $W > m_1 + m_2$. The second Riemann sheet is defined as the region where the propagator has an imaginary part of $i\rho/2$, i.e., shifted by $i\rho$ relative to the first sheet.

When searching for the poles, we use the first Riemann sheet for $W < m_1 + m_2$, and the second Riemann sheet for $W > m_1 + m_2$. In the code, the `lRm` in `qn` is chosen as `0` in this case. If only first or  second Riemann sheet is considered, `lRm` is set to `1` or `2`.

In the code, the `lRm` can be set to a `Int64` number `0`, `1`, or `2` for all channels, or `(n1,n2,n2,...)` for each channel, where `n1`, `n2`, etc. are the Riemann sheet numbers for each channel.

**Attention**: Physical observables are computed on the first Riemann sheet—more precisely, along the real axis (`lRm=1`). Therefore, for $W < m_1 + m_2$, the treatment is consistent with that used in pole searches. However, for $W > m_1 + m_2$, the pole search is carried out on the second Riemann sheet, where the propagator has an imaginary part of $i\rho/2$, in contrast to $-i\rho/2$ on the first sheet. This leads to a relative phase difference between the two Riemann sheets along the real axis.

## Physical observable

With the obtained amplitude $M^{J^P}$, we can also calculate the physical observable. Note that all physical observable are at real axis, we choose the onshell momentum as

$$
\begin{align}
M_{ij}(z)=\{[(1-VG)^{-1}]V\}_{ij}
\end{align}
$$

with $i$ and $k$ chosen as the onshell momentum, that is, $0$ dimension for $G$, and extra dimension for $V$.

The $|M|^2$ for each channel is calculated in `qBSE.res`.

### The cross section for the channel considered

For the open channel, the cross section can be obtained as

$$
\begin{align}
	\frac{d\sigma}{d\Omega}=\frac{1}{\tilde{j}_1\tilde{j}_2}\frac{1}{64\pi^2
	s}\frac{{\rm k}'}{{\rm k}}\sum_{\lambda\lambda'}|i{\cal M}_{\lambda\lambda'}({\boldsymbol k}',{\boldsymbol k})|^2.
\end{align}
$$

where $j_1$ and $j_2$ is the spin of the intitial particles, and we define $\tilde{j}=2j+1$.

The total cross section can be written as

$$
\begin{align}
	\sigma
=\frac{1}{\tilde{j}_1\tilde{j}_2}\frac{1}{64\pi^2
	s}\frac{{\rm k}'}{{\rm k}}\sum_{J,\lambda}N_J^2|i{\cal M}^J_{\lambda\lambda'}({\rm k}',{\rm k})|^2
=\frac{1}{\tilde{j}_1\tilde{j}_2}\frac{1}{64\pi^2
	s}\frac{{\rm k}'}{{\rm k}}\sum_{J^P,ij}N_J^2\left|{{ M}}^{J^P}_{ij}\right|^2.
\end{align}
$$

Here, ${\rm k}$' and ${\rm k}$ are onshell momenta, so we only choose $ij$ for the onshell momenta. Since we adopt ${\cal V}\to f(k'){\cal V}f(k)$, the amplitudes ${M}^{J^P}\to f(k'){M}^{J^P}f(k)$. The form factors vanish due to onshelness for initial and final states of a scattering.
NOTE: $4MM'$ should be multiplied due to convention $\bar{u}u=1$ adopted in code.

$$
\begin{align}
	\sigma&\propto \sum_{J,\lambda'\lambda}
	|M^{J}_{\lambda'\lambda}|^2=\sum_{J,\lambda'j=0}|M^{J}_{\lambda'j}|^2+\sum_{J,\lambda'j>0}
	\left[|M^{J}_{\lambda'j}|^2+|M^{J}_{\lambda'-j}|^2\right]\nonumber\\
&=\sum_{J^P,\lambda'j=0}\delta_{\eta,1}|\frac{1}{2} M^{J^P}_{\lambda'j}|^2+\sum_{J,\lambda'j>0}
	 \left[\frac{1}{4}|M^{J^+}_{\lambda'j}+M^{J^-}_{\lambda'j}|^2+\frac{1}{4}|M^{J^+}_{\lambda'j}-M^{J^-}_{\lambda'j}|^2\right]\nonumber\\
&=\sum_{J^P,i=0j=0}\delta_{\eta,1}\delta_{\eta',1}|\frac{1}{2}M^{J^P}_{ij}|^2+\sum_{J^P,i>0j=0}2\delta_{\eta,1}|\frac{1}{2} M^{J^P}_{ij}|^2+\sum_{J^P,\lambda'j>0}
	\frac{1}{2}|M^{J^P}_{\lambda'j}|^2\nonumber\\
&=\sum_{J^P,i=0j=0}\delta_{\eta,1}\delta_{\eta',1}|\frac{1}{2} M^{J^P}_{ij}|^2
+\sum_{J^P,i>0j=0}\delta_{\eta,1}|\frac{1}{\sqrt{2}} M^{J^P}_{ij}|^2
+\sum_{J^P,i=0j>0}\delta_{\eta',1}|\frac{1}{\sqrt{2}}M^{J^P}_{ij}|^2
+\sum_{J^P,i>0j>0}|M^{J^P}_{\lambda'j}|^2\nonumber\\
&=\sum_{J^P,i=0j=0}\delta_{\eta,1}\delta_{\eta',1}|\frac{1}{2} M^{J^P}_{ij}|^2
+\sum_{J^P,i>0j=0}\delta_{\eta,1}|\frac{1}{\sqrt{2}} M^{J^P}_{ij}|^2
+\sum_{J^P,i=0j>0}\delta_{\eta',1}|\frac{1}{\sqrt{2}}M^{J^P}_{ij}|^2
+\sum_{J^P,i>0j>0}|M^{J^P}_{ij}|^2\nonumber\\
&=\sum_{J^P,i\geq0j\geq0}|f_{i}f_{j}M^{J^P}_{ij}|^2=\sum_{J^P,i\geq0j\geq0}|{M}^{J^P}_{ij}|^2
\end{align}
$$

The cross section for certain channel is calculated with `qBSE.simpleXsection`

### Argand plot

The amplitudes can be written as

$$
\begin{align}
	i{\cal M}({\boldsymbol k}', {\boldsymbol k})=-8\pi\sqrt{s}f({\boldsymbol k}',{\boldsymbol
	k})=-\frac{8\pi\sqrt{s}}{|{\boldsymbol k}|}\sum
	 _{JM}N^2_J D^{J*}_{\lambda_R,\lambda}(\phi',\theta',0)
	a^J_{\lambda\lambda'}(|{\boldsymbol	k}'|,|{\boldsymbol k}|)D^{J}_{M,\lambda'}(\phi,\theta,0).
\end{align}
$$

where $a^J=-\frac{|{\boldsymbol k}|}{8\pi\sqrt{s}}{\cal M}^J(|{\boldsymbol k}|)$, which can be displayed as a trajectory in an Argand plot.

# Three body decay

## kinematics

### Lorentz boost

Here, we consider an process $Y\to m_1X\to m_1m_2m_3$.
To study a $1\to3$ decay with the qBSE, we need consider the center of mass frame (CMS) of $Y$ (which is also the laboratory frame in this issue) and the $m_1m_2$ where the qBSE is applied. The momenta of initial and final particles in the CMS of $Y$, remarked as $lab$,  are

$$
\begin{align}
P^{lab}=(W,0,0,0);\  \ p^{lab}_1=(E^{lab}_1,{\boldsymbol p}^{lab}_1);\ \ p^{lab}_2=(E^{lab}_2,{\boldsymbol p}^{lab}_2);\ \ p^{lab}_3=(E^{lab}_3,{\boldsymbol p}^{lab}_3)
\end{align}
$$

The  Lorentz boost from $(m,{\boldsymbol 0})$ to $(E,{\boldsymbol k})$,

$$
\begin{align}
\Lambda^{\mu}_\nu=\frac{1}{m}\left(\begin{array}{cccc}
E({\boldsymbol k})&k_x&k_y&k_z\nonumber\\
k_x&m+\frac{k_x k_x}{E+m}&\frac{k_x k_y}{E+m}&\frac{k_x k_z}{E+m}\nonumber\\
k_y&\frac{k_y k_x}{E+m}&m+\frac{k_y k_y}{E+m}&\frac{k_y k_z}{E+m}\nonumber\\
k_z&\frac{k_z k_x}{E+m}&\frac{k_z k_y}{E+m}&m+\frac{k_z k_z}{E+m}\nonumber\\
\end{array}\right).
\end{align}
$$

With  the Lorentz boost  the momenta for particle 23 in the laboratory frame $(E^{lab}_{23},-{\boldsymbol p}^{lab}_1)$ can be written with the momenta in the CMS of particles 23 $(M_{23},{\boldsymbol 0})$ as $p^{lab}=\Lambda(E^{lab}_{23},-{\boldsymbol p}^{lab}_1) p^{cm}$,

$$
\begin{align}
{\boldsymbol p}^{lab}&={\boldsymbol p}^{cm}-\frac{{\boldsymbol p}^{lab}_1}{M_{23}}\left[\frac{-{\boldsymbol p}^{lab}_1\cdot {\boldsymbol p}^{cm}}{W-E^{lab}_1({\rm p}^{lab}_1)+M_{23}}+p^{0cm}\right],\nonumber\\
p^{0lab}&=\frac{1}{M_{23}}\left[(W-E^{lab}_1({\rm p}^{lab}_1))p^{0cm}-{\boldsymbol p}^{lab}_1\cdot{\boldsymbol p}^{cm}\right],
\end{align}
$$

where the $p_{23}+p_1=P$ is applied, and $M_{23}=\sqrt{(p^{lab}_2+p^{lab}_3)^2}=\sqrt{(p^{cm}_2+p^{cm}_3)^2}=\sqrt{(P-p_1)^2}$.

The momenta in CMS of $23$ can also be written with the momentum in laboratory frame as $p=\Lambda(E^{lab}_{23},{\boldsymbol p}^{lab}_1) p^{lab}$

$$
\begin{align}
{\boldsymbol p}^{cm}&={\boldsymbol p}^{lab}+\frac{{\boldsymbol p}^{lab}_1}{M_{23}}\left[\frac{{\boldsymbol p}^{lab}_1\cdot {\boldsymbol p}^{lab}}{M-E^{lab}_1({\rm p}^{lab}_1)+M_{23}}+p^{0lab}\right], \nonumber\\
p^{0cm}&=\frac{1}{M_{23}}\left[(M-E^{lab}_1({\rm p}^{lab}_1))p^{0lab}+{\boldsymbol p}^{lab}_1\cdot{\boldsymbol p}^{lab}\right].
\end{align}
$$

For the outgoing particles 2 and 3, the ${\boldsymbol p}^{cm}$ should be set along the $z$ axis, hence, an additional rotation is adopted as $p=R(-\theta^{cm},-\phi^{cm})p^{cm}$

$$
\begin{align}
R(-\theta,-\phi)=\frac{1}{m}\left(\begin{array}{cccc}
1&0&0&0\\
0&\cos \theta \cos \phi&\cos\theta\sin\phi&-\sin\theta\nonumber\\
0&-\sin\phi&\cos\phi&0\nonumber\\
0&\sin\theta\cos\phi&\sin\theta\sin\phi&\cos\theta\nonumber\\
\end{array}\right).
\end{align}
$$

As described above, we perform both a Lorentz boost and a rotation. Through the Lorentz boost, the momentum $p^{lab}$ in the lab frame transforms into $p^{cm}$ in the center-of-mass system (CMS) of particles 2 and 3. The value of $p_1^{cm}$ depends on $p_1^{lab}$ and $M_{23}$, and the lab-frame energy $E_1^{lab}$ can also be derived from $M_{23}$. Thus, $p_1^{cm}$ is ultimately dependent only on $\Omega_1^{lab}$ and $M_{23}$. The momenta $p_{2,3}^{cm}$ depend on $p^{lab}_{2,3}$ and $\Omega_1^{lab}$ and can be expressed in terms of $M_{23}$ and the spherical angle $\Omega_3^{cm}$.

After rotation, $p^{cm}$ becomes $p$ in the new CMS frame of particles 2 and 3. The momenta $p_{2,3}$ of the final particles are aligned along the z-axis and depend solely on $M_{23}$. However, for intermediate particles, they are functions of ${\rm p}'_3$ and $\Omega'_3$. Finally, $p_1$ depends on $\Omega_1^{lab}$, $M_{23}$, and $\Omega_3^{cm}$.

The Lorentz boost is performed by `Xs.LorentzBoost`, and rotation by `Xs.Rotation`.
A function `qBSE.LorentzBoostRotation` is also provided to make both a Lorentz boost and a rotation.

## Amplitude

Because the $|{\cal M}|^2$ is invariant in different reference frame, the amplitude for the direct decay can be written with the momenta in cm frame of partilces 1 and 2 obtained with Lorentz boost and rotation, as

$$
\begin{align}
i{\cal M}^{d}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3)&=i{\cal A}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(\Lambda p_1,p_2,p_3)=i{\cal A}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(\Omega^{lab}_1,M_{23},\Omega_3^{cm}), \ \ \ {\rm for\ onshell}\nonumber\\
i{\cal M}^{d}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p'_2,p'_3)&=i{\cal A}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(\Omega^{lab}_1,\Omega'_3,{\rm p}'_3,M_{23},\Omega_3^{cm}) \nonumber\\
%
&=\sum_{JM}N_JD^{J*}_{ M\lambda_{32}}( \Omega'_3)i{\cal A}^{JM}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(\Omega^{lab}_1,{\rm p}'_3,M_{23},\Omega_3^{cm}),\ \ \ {\rm for\ offshell}\nonumber\\
%
i{\cal A}^{JM}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(\Omega^{lab}_1,{\rm p}'_2,M_{23},\Omega_3^{cm})&=
N_J\int d\Omega'_3D^{J}_{M,\lambda'_{32}}( \Omega'_3) i{\cal A}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(\Omega^{lab}_1,\Omega'_3,{\rm p}'_3,M_{23},\Omega_3^{cm})\nonumber\\
%
i{\cal M}^{Z}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3)&=i\int \frac{d^4p'_3}{(2\pi)^4} {\cal T}_{\lambda_2,\lambda_3}(p_2,p_3;p'_2,p'_3)  G(p'_3){\cal A}_{\lambda_1;\lambda}(p_1,p'_2,p'_3)\nonumber\\
%
&=\sum_{\lambda'_2\lambda'_3}\int \frac{d^3p'_3}{(2\pi)^3} i{\cal T}_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}(p_2,p_3;p'_2,p'_3)  G_0(p'_3)i{\cal A}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(p_1,p'_2,p'_3)\nonumber\\
%
&=\sum_{\lambda'_2\lambda'_3}\int \frac{d^3p'_3}{(2\pi)^3} i{\cal T}_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}(\Omega'_3,{\rm p}'_3,M_{23})  G_0({\rm p}'_3)i{\cal A}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(\Omega^{lab}_1,\Omega'_3,{\rm p}'_3,M_{23},\Omega_3^{cm})\nonumber\\
%
&=\sum_{\lambda'_2\lambda'_3}\int \frac{{\rm p}'^{2}_3d{\rm p}'^{}_3d\Omega'_3}{(2\pi)^3} \sum_{JM}N_J^2\delta_{M\lambda_{32}}i{\cal T}^{JM}_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}({\rm p}'_3,M_{23})D^{J}_{M\lambda'_{32}}(\Omega'_3)  \nonumber\\
&\ \cdot\ G_0({\rm p}'_3)\sum_{J'M'}N_{J'}D^{J'*}_{M'\lambda'_{32}}( \Omega'_3)i{\cal A}^{J'M'}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(\Omega^{lab}_1,{\rm p}'_3,M_{23},\Omega_3^{cm})\nonumber\\
%
&=\sum_{\lambda'_2\lambda'_3}\int \frac{{\rm p}'^{2}_3d{\rm p}'^{}_3d\Omega'_3}{(2\pi)^3} \sum_{J}N_Ji{\cal T}^J_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}({\rm p}'_3,M_{23})   G_0({\rm p}'_3)i{\cal A}^{J\lambda_{32}}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(\Omega^{lab}_1,{\rm p}'_3,M_{23},\Omega_3^{cm})\nonumber\\
%
&=\sum_{J}N_J\sum_{\lambda'_2\lambda'_3}\int \frac{{\rm p}'^{2}_3d{\rm p}'^{}_3}{(2\pi)^3} i{\cal T}^J_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}({\rm p}'_3,M_{23}) G_0({\rm p}'_3) i{\cal A}^{J\lambda_{32}}_{\lambda_1;\lambda'_2,\lambda'_3;\lambda}(\Omega^{lab}_1,{\rm p}'_3,M_{23},\Omega_3^{cm})\nonumber\\
&\equiv i{\cal M}^{Z}_{\lambda_2,\lambda_3;\lambda}(\Omega^{lab}_1,M_{23},\Omega_3^{cm}),\nonumber\\
%
%
i{\cal M}^{Z}_{\lambda_2;\lambda_1,\lambda_3;\lambda}(p_1,p_2,p_3)%
&=\sum_{J}N_J\sum_{\lambda'_1\lambda'_3}\int \frac{\tilde{\rm p}'^{2}_3d\tilde{\rm p}'_3}{(2\pi)^3} i{\cal T}^J_{\lambda_1,\lambda_3;\lambda'_1,\lambda'_3}(\tilde{\rm p}'_3,M_{13}) G_0(\tilde{\rm p}'_3) i{\cal A}^{J\lambda_{31}}_{\lambda_2;\lambda'_1,\lambda'_3;\lambda}(\tilde{\Omega}_3^{lab}, \tilde{\rm p}'_3,M_{13},\tilde{\Omega}_3^{cm})\nonumber\\
&\equiv i{\cal M}^{Z}_{\lambda_2;\lambda_1,\lambda_3;\lambda}(\tilde{\Omega}_3^{lab},M_{13},\tilde{\Omega}_3^{cm}).
\end{align}
$$

Note: when consider the rescattering of different particles, the different cm frames should be adopted.

## Decay width

### The case with one rescattering

#### Phase space

The phase space is given by

$$
\begin{align}
d\Phi=(2\pi)^4\delta^4(P-p_1-p_2-p_3) \frac{d^3p_1}{2E_1(2\pi)^3}\frac{d^3p_2}{2E_2(2\pi)^3}\frac{d^3p_3}{2E_3(2\pi)^3}
\end{align}
$$

To study the invariant mass spectrum of the particles 2 and 3, it is convenient to rewrite the Lorentz-invariant phase space $d\Phi$ by taking as integration variables the direction of the momentum of particle 2 ${\boldsymbol p}_2^{cm}$ in the center-of-mass (cm) frame of particles 2 and 3. Thus, we first rewrite the phase factor as

$$
\begin{align}
d\Phi
&=(2\pi)^4\delta(E^{cm}_2+E^{cm}_3-W_{23})\delta^3({\boldsymbol p}^{cm}_2+{\boldsymbol p}^{cm}_3) \frac{d^3p^{lab}_1}{2E^{lab}_1(2\pi)^3}\frac{d^3p^{cm}_2}{2E_2^{cm}(2\pi)^3}\frac{d^3p^{cm}_3}{2E_3^{cm}(2\pi)^3}
\end{align}
$$

where $W_{23}^2=(M-E^{lab}_1)^2-|{\rm p}^{lab}_1|^2$. Here the Lorentz invariance of the  $\frac{d^3p}{2E(2\pi)^3}$ and $\delta^4(P-p_1-p_2-p_3)$ is used.

Here, we use the momentum of the particle 3 in the center of mass system of two rescattering particles.
The momentum of the particle 3 has a relation ${\rm p}^{cm}_3=\frac{1}{2M_{23}}\sqrt{\lambda(M_{23}^2,m_3^2,m_2^2)}$.

Owing to the three-momentum $\delta$ function, the integral over ${\boldsymbol p}^{cm}_2$ can be eliminated. Next, the quantity $d^3p^{cm}_3$ is converted to $dM_{23}$ by the relation,

$$
\begin{align}
d^3p^{cm}_3=\frac{E^{cm}_2E^{cm}_3{\rm p}^{cm}_3}{M_{23}}dM_{23}d\Omega^{cm}_3,
\end{align}
$$

where $M_{23}$($=E^{cm}_2+E^{cm}_3$) is the invariant mass of the $23$ system. We would like to integrate over the magnitude of the neutron momentum $p_1$, which is related to $W_{23}$. Hence, the energy-conserving $\delta$ function is substituted as,

$$
\begin{align}
\delta(M_{23}-W_{23})=\frac{W_{23}}{|M{\rm p}^{lab}_1/E^{lab}_1|}\delta(\breve{\rm p}_1-{\rm p}^{lab}_1)
\end{align}
$$

where the $\breve{\rm p}_1$ satisfies $M_{23}^2=(M-\breve{E}_1)^2-\breve{\rm p}_1^2$.

Performing the integral over ${\rm p}^{lab}_1$, we obtain the final expression of the decay width,

$$
\begin{align}
d\Phi
&=\frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}^{cm}_3}{8M}d\Omega^{lab}_1d\Omega^{cm}_3dM_{23}
=\frac{1}{(2\pi)^3}\frac{\breve{\rm p}_1{\rm p}^{cm}_3}{8M}d\cos\theta^{lab}_1d\cos\theta^{cm}_3dM_{23},
\end{align}
$$

Here, independence of the $\phi^{lab}_1$ and $\phi^{cm}_2$ on the integrand is applied.

#### Differential decay width

Here, we consider the rescattering of particles 2 and 3. The distribution can be obtained as

$$
\begin{align}
{d\Gamma\over dM_{23}}&=\int \frac{1}{2M}\sum_{\lambda_2,\lambda_3;\lambda}|i{\cal M}_{\lambda_2,\lambda_3;\lambda}|^2 \frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}_3}{8M}d\Omega^{lab}_1d\Omega^{cm}_3
\end{align}
$$

If the $Y$ is not scalar, there shoold be an additional factor $1/\tilde{j}_Y$. Here we
use the distribution above in the frame with Lorentz boost but without rotation. However,
we calculate the $|i{\cal M}_{\lambda_2,\lambda_3;\lambda}|^2$ in the frame with both Lorentz boost and rotation due to it is independent on the frames.

If the reflection effect is absent, for example, the $\pi DD^*$ final state. The results can be simplified further.
The amplitude can be written as

$$
\begin{align}
i{\cal M}_{\lambda_1\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3)
&=\sum_{JM}N^3_JD^J_{\lambda_R\lambda_1}D^J_{\lambda_R\lambda_{32}}(\Omega^{cm}_3)\left[i{\cal A}^{d,J}_{\lambda_1,\lambda_2,\lambda_3;\lambda}(M_{23})\right.\nonumber\\
%
&\left.+\sum_{\lambda'_2\lambda'_3}\int \frac{{\rm p}'^{2}_3}{(2\pi)^3} i{\cal T}^{J}_{\lambda_2,\lambda_3;\lambda'_2,\lambda'_3}({\rm p}'^{}_3,M_{23}) G^{}_0({\rm p}'^{}_3) i{\cal A}^{J}_{\lambda'_2,\lambda'_3;\lambda}({\rm p}'_3,M_{23})\right]\nonumber\\
&\equiv \sum_{JM}N^3_JD^J_{\lambda_R\lambda_1}(\Omega^{lab}_1)D^J_{\lambda_R\lambda_{32}}(\Omega^{cm}_3)i{\cal M}^J_{\lambda_2,\lambda_3;\lambda}(M_{23})
\end{align}
$$

Inserting the above amplitude to the definition of the invariant mass spectrum, we have

$$
\begin{align}
{d\Gamma\over dM_{23}}
&=\frac{1}{2M}\frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}_3}{8M}\sum_{\lambda_2,\lambda_3;\lambda}\int|\sum_{JM}N^3_JD^J_{\lambda_R\lambda_1}(\Omega^{lab}_1)D^J_{\lambda_R\lambda_{32}}(\Omega^{cm}_3)i{\cal M}^J_{\lambda_2,\lambda_3;\lambda}(M_{23})|^2 d\Omega^{lab}_1d\Omega_3\nonumber\\
&=\frac{1}{2M}\frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}^{cm}_3}{8M}\sum_{\lambda_2,\lambda_3;\lambda;J}|N_Ji{\cal M}^J_{\lambda_2,\lambda_3;\lambda}(M_{23})|^2 \nonumber\\
\end{align}
$$

Now we consider the amplitude with fixed parity,

$$
\begin{align}
i{\cal M}^{J}_{\lambda_{23};\lambda}&=i{\cal A}^{d,J}_{\lambda_{23};\lambda}+i{\cal T}^{J}_{\lambda_{23},\lambda'_{23}}G_0i{\cal A}^J_{\lambda'_{23};\lambda},\quad
\eta i{\cal M}^{J}_{-\lambda_{23};\lambda}=\eta i{\cal A}^{d,J}_{-\lambda_{23};\lambda}+\eta i{\cal T}^{J}_{-\lambda_{23},\lambda'_{23}}G_0i{\cal A}^J_{\lambda'_{23};\lambda}\nonumber\\
\Rightarrow i{\cal M}^{J^P}_{\lambda_{23};\lambda}&=i{\cal A}^{d,J^P}_{\lambda_{23};\lambda}+i{\cal T}^{J^P}_{\lambda_{23},\lambda'_{23}}G_0i{\cal A}^{J}_{\lambda'_{23}
;\lambda}=i{\cal A}^{d,J^P}_{\lambda_{23};\lambda}+\eta'i{\cal T}^{J^P}_{\lambda_{23},\lambda'_{23}}G_0i{\cal A}^{J}_{-\lambda'_{23}
;\lambda}=i{\cal A}^{d,J^P}_{\lambda_{23};\lambda}+{1\over 2}i{\cal T}^{J^P}_{\lambda_{23},\lambda'_{23}}G_0i{\cal A}^{J^P}_{\lambda'_{23}
;\lambda}
\end{align}
$$

Here ${\cal M}^{J^P}_{\lambda_{23};\lambda}={\cal M}^{J}_{\lambda_{23};\lambda}+\eta{\cal M}^{J}_{\lambda_{23};\lambda}$, ${\cal A}^{d,J^P}_{\lambda_{23};\lambda}={\cal A}^{d,J}_{\lambda_{23};\lambda}+\eta{\cal A}^{d,J}_{\lambda_{23};\lambda}$, ${\cal T}^{J^P}_{\lambda_{23},\lambda'_{23}}={\cal T}^{J}_{\lambda_{23},\lambda'_{23}}+\eta{\cal T}^{J}_{-\lambda_{23},\lambda'_{23}}={\cal T}^{J}_{\lambda_{23},\lambda'_{23}}+\eta'{\cal T}^{J}_{\lambda_{23},-\lambda'_{23}}$.
NOTE: For parity conserving interaction, we have ${\cal T}^J_{\lambda'\lambda}=\eta(\eta')^{-1}{\cal T}_{-\lambda'\lambda}$, which can be checked in the code by different definitions of ${\cal T}^{J^P}$.

We summarize the results as following,

$$
\begin{align}
{d\Gamma\over dM_{23}}&=\frac{1}{2M}\frac{1}{(2\pi)^5}\frac{\breve{\rm p}_1{\rm p}_3}{8M}\sum_{ij;J^P}|N_J{\cal M}^{J^P}_{ij}(M_{23})|^2, \nonumber\\
\hat{\cal M}^{J^P}_{i;i}(M_{23})
&={\cal A}^{d,J^P}_{j;i}(M_{23})+\sum_{k}\int \frac{d{\rm p}'_3{\rm p}'^{2}_3}{(2\pi)^3}{\cal T}^{J^P}_{j;k}({\rm p}'_3,M_{23}) G_0({\rm p}'_3) {\cal A}^{J^P}_{k;i}({\rm p}'_3,M_{23}),
\end{align}
$$

where $i$ and $j$ denote the independent $\lambda_{2,3}$ and $\lambda$, and the factors $f_{i=0}=1/\sqrt{2}$ and $f_{i\neq 0}=1$ are inserted. The above equation can be abbreviated as $M={A}^{d}+{T}G { A}$, where  $T$ is solved by the Bethe-Salpeter equation $T=V+VGT$. NOTE: The $\hat{T}$ should be multiplied by $4\pi$ to be used as ${\cal T}$, and $\hat{G}_0$ should be divided by $4\pi$ to be used as ${G}_0$, which is cancelled by each other.

### The case with more than one rescattering

In the case with more than one resacttering, we should conisder Monte-Carlo method to generate the
event.

$$
\begin{align}
d\Gamma=\frac{1}{2E}\sum|{\cal M}|^2 d\Phi=\frac{1}{2E}\sum|{\cal M}|^2 (2\pi)^{4-3n} dR
\end{align}
$$

The distribution can be calculated with `qBSE.Xsection`.

Here we consider a process with direct, 23 rescattering, 13 rescattering.

$$
\begin{align}
i{\cal M}^{d}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3)&=i{\cal A}_{\lambda_1;\lambda_2,\lambda_3;\lambda}(p_1,p_2,p_3),\nonumber\\
%
%
i{\cal M}^{Z}_{\lambda_k;\lambda_i,\lambda_j;\lambda}(p_k,p_i,p_j)
&=\sum_{J}N_{J}\sum_{\lambda'_i\lambda'_j}\int \frac{{\rm p}'^{2}_jd{\rm p}'_j}{(2\pi)^3} i{\cal T}^J_{\lambda_i,\lambda_j;\lambda'_i,\lambda'_j}({\rm p}'_j,M_{ij})  \ G_0({\rm p}'_j) i{\cal A}^{J\lambda_{ji}}_{\lambda_k;\lambda'_i,\lambda'_j;\lambda}(\Omega^{lab}_k,{\rm p}'_j,M_{ij},\Omega^{cm}_j).
\end{align}
$$

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
&=\frac{1}{2}\sum_P T^{J^P}_{\lambda}A^{J^P}.\nonumber\\
\end{align}
$$

$$
\begin{align}
&{\cal A}^{J^P\lambda_{ji}}_{\lambda_k;\lambda'_i,\lambda'_j;\lambda}(...)=
\int d\Omega'_j \left[D^{J}_{\lambda_{ji},\lambda'_{ji}}(\phi'_j, \theta'_j,0) {\cal A}_{\lambda_k;\lambda'_i,\lambda'_j;\lambda}(...,\Omega'_j,...)+\eta'D^{J}_{\lambda_{ji},-\lambda'_{ji}}(\phi'_j, \theta'_j,0) {\cal A}_{\lambda_k;-\lambda'_i,-\lambda'_j;\lambda}(...,\Omega'_j,...)\right].
\end{align}
$$

$$
\begin{align}
i{\cal M}^{Z}_{\lambda_k;\lambda_i,\lambda_j;\lambda}(p_k,p_i,p_j)
&=\sum_{J^P}\frac{N_J^2}{2}\int \frac{{\rm p}'^{2}_jd{\rm p}'_j}{(2\pi)^3} \sum_{i'j'}i{\cal T}^{J^P}_{\lambda_i,\lambda_j;i'j'}({\rm p}'_j,M_{ij})  \ G_0({\rm p}'_j) i{\cal A}^{J^P\lambda_{ji}}_{\lambda_k;i'j';\lambda}(\Omega^{lab}_k,{\rm p}'_j,M_{ij},\Omega^{cm}_j).
\end{align}
$$

Here we move the $N_J$ in ${\cal A}^{J^P}$ to the definition of ${\cal M}$. The amplitude
are calculated in `qBSE.TGA`.

$$
\begin{align}
	\sum_{\lambda}|\frac{1}{2}T^{J^P}_{\lambda}A^{J^P}|^2&=|\frac{1}{2}T^{J^P}_{0}A^{J^P}|^2+2\sum_{j>0}|\frac{1}{2}T^{J^P}_{j}A^{J^P}|^2=\frac{1}{2}\sum_j|T^{J^P}A^{J^P}|^2.
\end{align}
$$

$$
\begin{align}
\sum_{\lambda_k;\lambda_i,\lambda_j;\lambda}|i{\cal M}^{Z}_{\lambda_k;\lambda_i,\lambda_j;\lambda}(p_k,p_i,p_j)|^2
&=\frac{1}{2}\sum_{J^P,\lambda_k,\lambda}N^2_{J}\sum_{j}|T^{J^P}A^{J^P}|^2.
\end{align}
$$

# qBSE package

The qBSE package is used to solve the Bethe-Salpeter equation with some auxiliary functions.

## Data Structures for the Dimensions

In the qBSE approach, matrix dimensions for a coupled-channel system are organized hierarchically: the outermost level corresponds to different channels, the next to independent helicities, and the innermost to discretized momentum points. The data structures are designed to clearly separate and manage these levels.

### `structChannel`

The `struct structChannel` structure encapsulates the properties of a physical channel in the qBSE package, often used as `CH` in the code. Its fields include:

- `p::Vector{Float64}`, `m::Vector{Float64}`, `J::Vector{Int64}`, `Jh::Vector{Int64}`, `P::Vector{Int64}`: Vectors specifying the momentum $p$, mass $m$, total angular momentum $J$=`J/Jh`, and parity (`P`) for each particle in the channel. Here, `Jh = 1` for integer spin and `Jh = 2` for half-integer spin.
- `cutoff::Float64` specifying the cutoff parameter for the channel.
- `name::String`, `name0::String`  for the channel name, with and without charge. `name0` is used to define channels with definite isospin.
- `IHb::Int64`, `IHe::Int64`, `IHn::Int64` indicating the starting and ending indices of independent helicities, and the total number of independent helicity states for the channel.

These definitions ensure that each channel's quantum numbers and relevant parameters are explicitly tracked, supporting efficient matrix construction and manipulation in qBSE calculations.

### `structIndependentHelicity`

The `mutable struct structIndependentHelicity` structure defines the independent helicity states in the qBSE package, commonly referenced as `IH` in the code. Its fields are:

- `iCH::Int64`: Index of the channel to which this independent helicity belongs.
- `hel::Vector{Int64}`, `helh::Vector{Int64}`: Vector of helicity values for the independent helicity states, with `helh` indicating whether each helicity state corresponds to a fermion or boson.
- `Dimb::Int64`, `Dime::Int64`, `Dimn::Int64` for the starting index, ending index, and total count of independent helicities (i.e., the dimension range).
- `Dimo::Int64`: Indicates whether an extra on-shell dimension is present (`+1` if $W > \sum m$, otherwise `0`). The total number of momentum discretization points for an independent helicity is `Np + Dimo`.

These fields help manage the indexing and dimensionality of helicities, especially when handling on-shell and off-shell contributions in the qBSE framework.

### `structDimension`

The `mutable struct structDimension` structure stores information about each dimension in total matrix. Its fields are:

- `iIH::Int64`: Index indicating which independent helicity this dimension belongs to.
- `k::ComplexF64`, `w::Float64`: The discretized momentum value and its corresponding integration weight for this dimension (used in numerical integration).
- `Dimo::Int64`: Indicates whether an extra on-shell dimension is present (`+1` if $W > \sum m$, otherwise `0`). The total number of discretization points for an independent helicity is `Np + Dimo`.

This structure is auxiliary one to above structures and ensures that each discretized momentum point is properly associated with its helicity and integration weight, supporting accurate matrix construction and numerical calculations in the qBSE framework.

## Data Structures for the Interactions

In the qBSE package, interaction and system information are encapsulated in dedicated data structures to facilitate efficient computation and clear organization.

### `structSys`

The `struct structSys` (often referenced as `SYS` in the code) stores information about the system and the generally used discretization and angular integration data. Its fields include:

- `Sys::String`: A label identifying the system.
- `kv::Vector{Float64}`,`wv::Vector{Float64}`: Discretized momentum points and weights  for the momentum discretization (used in numerical integration).
- `xv::Vector{Float64}`, `wxv::Vector{Float64}`: Discretized values and weight of $\cos\theta$.
- `d::Vector{Matrix{Float64}}`: Precomputed Wigner $d$-matrices of $\theta$.
- `pv::Vector{Float64}`, `wpv::Vector{Float64}`: Discretized azimuthal angles and weight of $\phi$ discretization.
- `sp::Vector{Float64}`, `cp::Vector{Float64}`: Sine values and Cosine values of the discretized $\phi$ angles.

These fields provide all necessary information for numerical integration over momentum and angular variables in the qBSE framework.

### `structInterAction`

The `struct structInterAction` structure (typically used as `IA` in the code) stores the properties of each interaction or exchange process. Its fields are:

- `Nex::Int64`: Total number of exchange particles or processes.
- `name_ex::Vector{String}`,`key_ex::Vector{Int64}`: Names and integer labels of the exchanged particles or interactions.
- `J_ex::Vector{Int64}`, `Jh_ex::Vector{Int64}`, `P_ex::Vector{Int64}`, `m_ex::Vector{Float64}`: Spin, parity, and mass of  exchanged particles.
- `dc::Vector{Int64}`: Indicates whether the exchange is direct or crossed.
- `Ff::Vector{Float64}`: Flavor factors associated with each exchange.

These definitions ensure that all relevant quantum numbers and parameters for each interaction are explicitly tracked, supporting flexible and accurate construction of the interaction kernel in qBSE calculations.

## Additional Data Structures

### `structMomentum`

The `mutable struct structMomentum` structure stores the four-momenta and related kinematic quantities for a $2 \to 2$ scattering process. Its fields are:

- `i1::SVector{5,ComplexF64}`, `f1::SVector{5,ComplexF64}`: Initial and final momentum of particle 1.
- `i2::SVector{5,ComplexF64}`, `f2::SVector{5,ComplexF64}`: Initial and final momentum of particle 2.
- `q::SVector{5,ComplexF64}`: Momentum of the exchanged particle.
- `q2::Complex{Float64}`: Invariant mass squared ($q^2$) of the exchange.
- `qt::Complex{Float64}`: Additional kinematic variable (e.g., $t$-channel momentum transfer).

### `structHelicity`

The `mutable struct structHelicity` structure stores the helicity information for the particles in a $2 \to 2$ interaction. Its fields are:

- `i1::Int64`, `i1h::Int64`: Helicity and denominator for initial particle 1.
- `f1::Int64`, `f1h::Int64`: Helicity and denominator for final particle 1.
- `i2::Int64`, `i2h::Int64`: Helicity and denominator for initial particle 2.
- `f2::Int64`, `f2h::Int64`: Helicity and denominator for final particle 2.

### `structParticle`

The `struct structParticle` structure defines the properties of a single particle. Its fields are:

- `name::String`, `name0::String`, `nameL::String`: Particle name (with charge and without charge), LaTeX representation of the particle name.
- `m::Float64`: Mass of the particle.
- `J::Int64`, `Jh::Int64`, `P::Int64`: Spin and Parity.

A function is also provided to read particle information from a file and populate the particle list and key dictionary:

`function particles!(particles::Vector{structParticle}, pkey::Dict{String,Int64}, filename::String)`

- `particles`: A global vector (`const p = structParticle[]`) storing all particle information.
- `pkey`: A global dictionary (`const pkey = Dict{String,Int64}()`) mapping particle names to their indices.

This setup allows efficient lookup and management of particle properties for use in qBSE calculations.

## Functions for the qBSE

### `function preprocessing(Sys, channels, Ff, qn; Np=10, Nx=5, Nphi=5)`

This function prepares the system and channel data structures for qBSE calculations by setting up the necessary discretization and quantum number information.

**Arguments:**

- `Sys::String`: Label identifying the system, stored in `SYS.Sys`.
- `channels`: List of channels to be included in the calculation, stored in `IA[]`.
- `Ff`: Flavor factors, stored in `IA[]`.
- `qn`: Quantum numbers for the process, and labels for Riemann sheets.
- `Np`, `Nx`, `Nphi`: Number of momentum discretization points (default: 10), $\cos\theta$ discretization points (default: 5), azimuthal angle discretization points (default: 5).

Example for arguments:

```julia
Sys == "KNcp"
qn = (I=1, Ih=1, J=1, Jh=2, P=-1, C=-1, lRm=1)

channels = (
    ("K_b0", "N_p", qBSE.p[qBSE.pkey["K"]].m + 0.22 * 1.63),
    ("pi_0", "Sigma_p", qBSE.p[qBSE.pkey["pi"]].m + 0.22 * 1.63),
    ("pi_p", "Sigma_0", qBSE.p[qBSE.pkey["pi"]].m + 0.22 * 1.63),
    ("pi_p", "Lambda", qBSE.p[qBSE.pkey["pi"]].m + 0.22 * 1.63),
    ("eta", "Sigma_p", qBSE.p[qBSE.pkey["eta"]].m + 0.22 * 1.63)
	)

Ff = Dict(
    "K_b0:N_p-->K_b0:N_p" => ([["V", 1], 1.0],),
    "K_b0:N_p-->pi_0:Sigma_p" => ([["V", 1], -sqrt(0.5)],),
    "K_b0:N_p-->pi_p:Sigma_0" => ([["V", 1], sqrt(0.5)],),
    "K_b0:N_p-->pi_p:Lambda" => ([["V", 1], -sqrt(1.5)],),
    "K_b0:N_p-->eta:Sigma_p" => ([["V", 1], sqrt(1.5)],),
    "pi_0:Sigma_p-->pi_p:Sigma_0" => ([["V", 1], -2.0],)
	)
```

**Returns:**

- `SYS::structSys`: The processed system object containing discretization and system information.
- `IA::Matrix{structInterAction}`: The processed interaction objects, organized as a matrix.
- `CH::Vector{structChannel}`: The processed channel objects, each describing a physical channel.
- `IH::Vector{structIndependentHelicity}`: The processed independent helicity objects, ready for use in qBSE calculations.

### `function res(Range, iER, qn, SYS, IA, CH, IH, fV)`

This function calculates the rescatering process by qBSE for a given range of energies.

**Arguments:**

- `Range`: The range of energies to be considered. For example `Range = (ERmin=1.2, ERmax=2.0, NER=200, EIt=0.200, NEI=20, Ep="cm")`
- `iER`: The index of the energy range, energy can be obtained as `ER = Range.ERmax - iER * (Range.ERmax - Range.ERmin) / (Range.NER - 1)`.
- `qn`: The quantum numbers for the process.
- `SYS`, `IA`, `CH`, `IH`:  Obtained by `preprocessing`.
- `fV`: The potential kernel for qBSE.

**Returns:**

- `Ect::Vector{ComplexF64}`: The complex energy values for the system.
- `reslogt::Vector{Float64}`: The values of $\log|1 - VG|$.
- `resM2::Matrix{Float64}`: The matrix $|M|^2$ for the channels.
- `IH::Vector{structIndependentHelicity}`: The independent helicity objects, recalculated from the input `IH`.
- `Dim::Vector{structDimension}`: The dimension objects, recalculated from the input `Dim`.
- `TG::Matrix{ComplexF64}`: The $TG$ matrix used for decays.

### `function fV(k, l, SYS, IA0, CHf, CHi)`

This function defines the potential kernel for the qBSE calculation. It is used to compute the interaction potential between particles in a scattering process.

**Arguments:**

- `k::structMomentum`,`l::structHelicity`: The momenta and helicities of initial and final particles.
- `SYS::structSys`: The same as above.
- `IA0::structInterAction`: The interaction object containing information about the interactions and their properties, one of element of `IA` above.
- `CHf::structChannel`,`CHi::structChannel`: The final and initial channel object containing information about the final state particles, one of element of `CH` above.

**Returns:**

A value representing the potential kernel, suitable for use in the qBSE package's calculations.

**Note:**
This is an internal function intended for use within the qBSE package. Users should not call this function directly.

### `function propFF(k, ex, L, LLi, LLf; lregu=1, lFFex=0)`

In the `fV` function, form factors can be included via the auxiliary function `propFF`, which is used to regulate the interaction kernel.

**Arguments:**

- `k::structMomentum`: The kinematic information for the process.
- `ex::Int64`: Index of the exchanged particle, usually set as `ex = IA0.key_ex[le]`, where `le` runs over the exchanged mesons in `1:IA0.Nex`.
- `L::Float64`: The cutoff parameter for the exchanged meson, often chosen as the average of the initial and final particle cutoffs: `L = (LLi + LLf) / 2`.
- `LLi::Float64`, `LLf::Float64`: Cutoff parameters for the initial and final particles, typically set as `LLi = CHi.cutoff`, `LLf = CHf.cutoff`.
- `lregu::Int64`: Regularization flag for constituent particles (Only `1`).
- `lFFex::Int64`: Type of form factor for the exchanged meson:
  - `0`: No form factor.
  - `1`: $\frac{L^2 - m^2}{L^2 - q^2}$
  - `2`: $\frac{L^4}{(m^2 - q^2)^2 + L^4}$
  - `3`: $\exp\left(-\frac{(m^2 - q^2)^2}{L^4}\right)$
  - `4`: $\frac{L^4 + (q_t - m^2)^2 / 4}{(q^2 - (q_t + m^2) / 2)^2 + L^4}$
  - `5`: $\frac{L^2}{L^2 - q^2}$

**Description:**

These options allow flexible control over the inclusion and type of form factors in the potential kernel, supporting different regularization schemes as needed for the physical system under study.

### `function simpleXsection(ER, resM2, CH, qn; Ep="cm")`

This function computes the total cross section for a $2 \to 2$ scattering process using the squared amplitude matrix and channel information.

**Arguments:**

- `ER::Vector{Float64}`: Array of energy values (typically from `Ect`).
- `resM2::Matrix{Float64}`: Matrix of squared amplitudes $|M|^2$ for the process, as obtained from the qBSE calculation.
- `CH`: Vector of channel structures, as described above.
- `qn`: Quantum numbers for the process.
- `Ep`: Energy frame to use (`"cm"` for center-of-mass, `"L"` for laboratory; default: `"cm"`).

**Returns:**
A `Vector{Matrix{Float64}}` value representing the total cross section at each energy point.

**Description:**
The function calculates the cross section by summing or integrating over the relevant momentum discretization points and applying the appropriate kinematic and symmetry factors, including spin averaging and flux normalization. It uses the channel and quantum number information to determine the initial and final state properties.

### `function lambda(m1, m2, m3)`

This function calculates the Kallen function $\lambda(m_1, m_2, m_3)$, which is used in the phase space integration.

## Decay

### `function setTGA(par, sij, k, tecm, i, j)`

set the frame and other things for calculating TGA, which should be usde before `TGA` function.

**Arguments:**

- `par`: The parameters for the TGA calculation.
- `sij`: The invariant mass squared of the two particles $i$ and $j$.
- `k`: The momentum of the particle in the center-of-mass (CM) frame.
- `tecm`: The total energy in the center-of-mass frame.
- `i`: The index of the first particle.
- `j`: The index of the second particle.

**Returns:**

`para=(E=sqrt(sij), par=par, IH=IH, Dim=Dim, k=kn, P=Pn)`

### `TGA(cfinal, cinter, Vert, para)`

This function calculates the transition amplitude for a given final state `cfinal` and interaction `cinter`.

**Arguments:**

- `cfinal`: The index of final state configuration (e.g., channel and particle information for the outgoing particles).
- `cinter`: The index of intermediate interaction.
- `Vert`: The vertex function and associated data, including the vertex structure, precomputed spinors or polarization vectors, and helicity assignments (e.g., `Vert=Vertex14, ULc=GA1 * FR.U(para14.P, lLc), l2=l2, l3=l3`).
- `para`: The calculation parameters, as returned by `setTGA`.

### `function Vertex14(p1, p2, l1, l2, Vert, k, P)`

This function computes the initial decay vertex for a process involving two rescattering particles.

**Arguments:**

- `p1`, `p2`: Four-momenta of the two rescattering particles.
- `l1`, `l2`: Helicity indices of the two rescattering particles.
- `Vert`: Vertex structure or function, as defined in the context of the calculation.
- `k`: Momenta of all final state particles.
- `P`: Four-momentum of the parent (initial) particle.

**Returns:**

A value representing the decay vertex amplitude, suitable for use in the qBSE package's transition amplitude calculations.

**Note:**
This is an internal function intended for use within the qBSE package. Users should not call this function directly.

### auxiliary function

`qBSE.showSYSInfo(Range, qn, IA, CH, IH)`

This function is used to display the system information, including the range of energies, quantum numbers, interaction information, and channel information.

`qBSE.showPoleInfo(qn, Ec, reslog, "data/output.txt")`

This function is used to display the pole information, including the energy, width, and other relevant parameters. The results are saved in a file named "output.txt" in the "data" directory.

## Additional functions

### `function ch(pf, pin, amps)`

This function constructs a channel tuple for use in decay width or cross section calculations, such as with `Xs.Xsection`.

**Arguments:**

- `pf`: List of final state particle names (e.g., `["pi_m", "pi_p", "pi_p", "Lambda"]`).
- `pin`: List of initial state particle names (can be empty if not needed).
- `amps`: Amplitude information or function.

**Returns:**

A tuple of the form `(pf=pf, namef=namef, mf=mf, pin=pin, namei=namei, mi=mi, amps=amps)`, where:

- `pf`, `pin`: Final and initial state particle names.
- `namef`, `namei`: Channel names for final and initial states.
- `mf`, `mi`: Masses of final and initial state particles.
- `amps`: Amplitude information.

**Example:**

To compute the decay width for $\Lambda_c \to \pi^- \pi^- \pi^+ \Lambda$:

```julia
ch = qBSE.ch(["pi_m", "pi_p", "pi_p", "Lambda"], [], amps)
```

The resulting `ch` can be passed directly to `Xs.Xsection`. See the `Xs.Xsection` documentation for further usage details.

### `function LorentzBoostRotation(k, tecm, p1, p2)`

This function performs a Lorentz boost and rotation to transform a particle's momentum from the lab (static) frame of the total system to the center-of-mass system (CMS) of two particles.

**Arguments:**

- `k`: The momentum of the particle in the CMS frame (after transformation).
- `tecm`: The total energy in the CMS frame.
- `p1`: The momentum of the first particle in the lab frame.
- `p2`: The momentum of the second particle in the lab frame.

**Returns:**

- `knew`: Array of momentum vectors after Lorentz boost and rotation.
- `Pnew`: Total four-momentum vector after the same transformations.

This function is useful for kinematic calculations where momenta need to be expressed in different reference frames, ensuring proper treatment of boosts and rotations.
