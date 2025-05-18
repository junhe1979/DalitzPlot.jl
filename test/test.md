The phase sapce is given as follows:

$d\Phi = (2\pi)^4\delta^4(\sum_{i}k_i-P)\prod_{i}\frac{d^3k_i}{(2\pi)^32E_i}$

Such definition is different from that in PDG by a factor of $(2\pi)^4$.

With such definition, the decay width in CMS is given as

$\Gamma = \frac{1}{2M}\int d\Phi |M|^2=\frac{ (2\pi)^{4-3n}}{2M}\int dR |M|^2$.

The cross section is given as

$d\sigma=\frac{1}{4\sqrt{(p_1\cdot p_2)^2-m_1^2m_2^2}}\frac{|p_1\cdot p_2|}{p_1^0p_2^0}d\Phi |M|^2=\frac{(2\pi)^{4-3n}}{4\sqrt{(p_1\cdot p_2)^2-m_1^2m_2^2}}\frac{|p_1\cdot p_2|}{p_1^0p_2^0}dR |M|^2$.

The phase space in the code is given as

$dR = (2\pi)^{3n-4} d\Phi = \prod_{i}\frac{d^3k_i}{2E_i}\delta^4(\sum_{i}k_i-P)$.

For the two-body final state, the phase space is given as

$dR_2 =\frac{d^3k_1}{2E_1}\frac{d^3k_2}{2E_2}\delta^4(k_1+k_2-P)=\frac{|{\bm p}|}{4\sqrt{s}}d\Omega$

$\int dR_2=\frac{\pi|{\bm p}|}{\sqrt{s}},$

with

$\frac{2|{\bm p}|}{\sqrt{s}}=\lambda^{1/2}(1,m_1^2/s,m^2_2/s)$,

and

$\lambda(x,y,z)=x^2+y^2+z^2-2xy-2yz-2xz$.

if $m_1=m_2=0$ then $\int dR_2=\pi/2$.

$dR_3 =\frac{d^3k_1}{2E_1}\frac{d^3k_2}{2E_2}\frac{d^3k_3}{2E_3}\delta^4(k_1+k_2+k_3-P)=\frac{(2\pi)^2}{16M^2}ds_{12}ds_{23}$

$\int dR_3=\int\frac{(2\pi)^2}{16M^2}ds_{12}d_{23}=\frac{\pi^2}{4M^2}\int^{(M-m_3)^2}_{(m_1+m_2)^2}ds_{12}\frac{\lambda^{1/2}(M^2,s_{12},m^2_3)\lambda^{1/2}(s_{12},m_1^2,m^2_2)}{s_{12}}$.

if $m_1=m_2=m_3=0$ then

$\int dR_3=\frac{\pi^2}{4M^2}\int^{M}_{0}ds_{12}\lambda^{1/2}(M^2,s_{12},m^2_2)=\frac{\pi^2}{4M^2}\int^{M^2}_{0}ds_{12}\sqrt{M^4+s_{12}^2-2M^2s_{12}}=\frac{\pi^2}{4M^2}\int^{M^2}_{0}ds_{12}(M^2-s_{12})=\frac{\pi^2M^2}{8}$.
