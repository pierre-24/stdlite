title: Theory

!!! note
    
    In the following, letters $i,j,\ldots$ ($a,b,\ldots$) are the label of occupied (virtual) molecular orbitals (MO), while $\mu, \nu, \ldots$ refers to atomic orbitals (AO).

## Introduction

Molecular properties are defined as a response of a (molecular) system to an external perturbation (change of geometry, application of an external electric field, etc).
The goal of `stdlite` is to compute response function/excited states.
It requires a ground state HF or DFT wavefunction, from which the orbitals (i.e., the MO energies, $\varepsilon$, and LCAO coefficients, $\mathbf C$) are extracted and used in the sTDA/sTD-DFT procedures.

## Response function theory

### Time-independent and dependent frameworks

In the time-independent (or static) case, a general expression for the electronic energy under perturbation is:

$$E(\{\kappa_i\}, \lambda) = \Braket{\Psi(\{\kappa_i\}) | \hat H(\lambda) | \Psi(\{\kappa_i\})},$$

where $\{\kappa_i\}$ are the set of parameters that minimize the energy (*e.g.*, the LCAO coefficients), while $\hat H(\lambda) = \hat H_0 + \lambda \hat V$, where $\hat H_0$ is the electronic Hamiltonian (for which solutions are known) and $\hat V$ is a perturbation.
In the framework of the [Rayleigh-Schrödinger perturbation theory](https://en.wikipedia.org/wiki/Perturbation_theory_(quantum_mechanics)), $E$ can be expanded in a power series, 

$$E= E^{(0)} + \lambda\,E^{(1)} + \lambda^2\,E^{(2)}+\ldots,$$

with:

$$E^{(k)} = \frac{1}{k!}\,\left.\frac{d^k E}{d \lambda^k}\right|_{\lambda=0}.$$

In particular, the first order derivative is:

$$\begin{aligned}
\frac{d E}{d\lambda} &= \frac{\partial E(\{\kappa_i\},\lambda)}{\partial\lambda} + \sum_i \left.\frac{ \partial E(\{\kappa_i \}, \lambda)}{\partial\kappa_i}\right|_{\lambda=0}\,\frac{\partial\kappa_i}{\partial\lambda} \\
&= \Braket{\Psi(\{\kappa_i \})|\frac{\partial\hat{H}(\lambda)}{\partial \lambda}|\Psi(\{\kappa_i\})},
\end{aligned}$$

which is the [Hellman-Feynman theorem](https://en.wikipedia.org/wiki/Hellmann%E2%80%93Feynman_theorem), obtained under the assumption that the $\{\kappa_i\}$ were obtained [variationaly](https://en.wikipedia.org/wiki/Variational_method_(quantum_mechanics)), given the stationary principle:

$$\forall\kappa_i: \left.\frac{\partial E(\{\kappa_i \}, \lambda)}{\partial\kappa_i}\right|_{\lambda=0} = 0.$$

The second-order derivative is obtained from differentiating the Hellman-Feynman theorem expression, yielding (under the assumption that the perturbation is linear):

$$\tag{1}\frac{d^2 E}{d\lambda^2} =
\sum_i \left.\frac{\partial^2 E(\{\kappa_i\},\lambda)}{\partial\lambda\partial\kappa_i}\right|_{\lambda=0}\,\frac{\partial\kappa_i}{\partial\lambda} 
\Leftrightarrow \braket{\braket{\hat V; \hat V}} = (\mathbf\kappa^\lambda)^T\,\eta^\lambda$$

where $\braket{\braket{\hat V; \hat V}}$ is the linear response function, $\mathbf\eta^\lambda$ is the perturbed electronic gradient, and $\mathbf\kappa^\lambda$, the first-order response vector to $\lambda$, containing all $\partial\kappa_i/\partial\lambda$.
$\kappa^\lambda$ is computed thanks to the following set of equations (derivative of the stationary condition):

$$\tag{2}\sum_j \left.\frac{\partial^2 E(\{\kappa_i\},\lambda)}{\partial\kappa_i\partial\kappa_j}\right|_{\lambda=0}\,\frac{\partial\kappa_j}{\partial\lambda} = -\left.\frac{\partial^2 E(\{\kappa_i\},\lambda)}{\partial\lambda\partial\kappa_i}\right|_{\lambda=0}  \Leftrightarrow \mathbf{E}^\lambda\,\mathbf\kappa^\lambda = -\mathbf\eta^\lambda$$
 
where $\mathbf E^\lambda$ is the electronic Hessian. Evaluating the elements of $\mathbf E^\lambda$ and $\mathbf\eta^\lambda$ actually result in new kinds of (bielectronic) integrals to evaluate (*vide supra*), while $\kappa^\lambda$ has to be computed by solving this equation.
Thanks to the [2n+1 rule](https://doi.org/10.1007/s10910-008-9497-x), the response vector can be used to access the quadratic response function (*i.e.*, the third derivative of the energy w.r.t. $\lambda$) as well.

The dynamic case is similar. 
Assuming a time-dependent perturbed Hamiltonian ([Heisenberg picture](https://en.wikipedia.org/wiki/Heisenberg_picture)), $\hat H(\lambda,t) = \hat H_0 + \hat V(\lambda, t)$, its development is governed by the time-dependent (TD) Shrödinger equation.
Following the [Floquet theory](https://en.wikipedia.org/wiki/Floquet_theory), it is proposed to extract a (position independent) phase factor, $F(t)$ from the wavefunction:

$$\Psi(\{\kappa_i\},t) = e^{i\,F(t)}\,\tilde\Psi(\{\kappa_i\},t),$$

where $\tilde\Psi$ is the phase isolated wavefunction. Inserting this expression in the TD Shrödinger equation results in:

$$\begin{aligned}
&\left[\hat H(\lambda,t)-i \frac{\partial}{\partial t}\right]\,\tilde\Psi(\{\kappa_i\},t) = Q(t)\,\tilde\Psi(\{\kappa_i\},t),\\
&\text{ with } Q(t) = \frac{\partial F(t)}{\partial t} = \Braket{\tilde\Psi(\{\kappa_i\},t)|\hat H(\lambda,t)-i \frac{\partial}{\partial t}|\tilde\Psi(\{\kappa_i\},t)}.
\end{aligned}$$

As a result, $Q(t)$, a real quantity, is called the TD quasienergy since it reduces to the energy $E_0$ in the time-independent case.
In Floquet theory, the perturbation is assumed to be periodic in time, of frequency $\omega$ (and period $T$), which means that $\hat V(\lambda,t)$ oscillates at a multiple of the fundamental frequency $\omega$.
Therefore, the Hamiltonian becomes periodic, which implies that $\tilde{\Psi}$ also oscillates with the same period.
As a consequence one can define the (time-averaged) quasienergy $\mathcal{Q}$ is given by:

$$\mathcal{Q} = \frac{1}{T}\,\int_0^T dt\,\Braket{\tilde\Psi(\{\kappa_i\},t)|\hat H(\lambda,t)-i \frac{\partial}{\partial t}|\tilde\Psi(\{\kappa_i\},t)}.$$

The quasienergy $\mathcal Q$ acts similarly to $E$ in the time-independent case, wich allows to obtain a TD version of the Hellman-Feynman theorem, and of the linear response equations.
Using a Fourier series of $\hat V$, one can obtain a TD version of the linear response equations [Eqs. (1) and (2)].

### Application to DFT: TD-DFT

Under the conditions of the [Runge and Gross theorem](https://en.wikipedia.org/wiki/Runge%E2%80%93Gross_theorem), this theory can be applied to DFT.
In the TD case and assuming the Hermicity of the different matrices and real orbitals, Eq. (2) can be written:

$$\tag{3}\left[\begin{pmatrix}
\mathbf A & \mathbf B \\
\mathbf B & \mathbf A
\end{pmatrix}-\omega\begin{pmatrix}
\mathbf 1 & \mathbf 0\\
\mathbf 0 & -\mathbf 1
\end{pmatrix}\right]\,\begin{pmatrix}
\mathbf x_\zeta^\omega\\ \mathbf y_\zeta^\omega
\end{pmatrix}=-\begin{pmatrix}
\mathbf \eta_\zeta\\ \mathbf \eta_\zeta
\end{pmatrix},$$

where $\mathbf x^\omega$ and $\mathbf y^\omega$ are the frequency-dependent linear response vector (to be determined) in direction $\zeta$.
The $\mathbf A$ and $\mathbf B$ are electronic Hessian (super-)matrices (related to orbital rotations). For a **global hybrid** density functional (independent of $\omega$),

$$E_{XC}= (1-a_x)\,E_X^{GGA}+a_x\,E_X^{HF}+E_C^{GGA},$$

their elements are:

$$\begin{aligned}
&A_{ia, jb} = \delta_{ij}\delta_{ab} (\epsilon_a - \epsilon_i) + 2(ia|jb) - a_x\,(ij|ab) + (1-a_x)\,(ia|f_{XC}|jb),\\
&B_{ia,jb} = 2(ia|bj) - a_x(ib|aj) + (1-a_x)\,(ia|f_{XC}|bj),
\end{aligned}$$

where, $\epsilon_i$ and $\epsilon_a$ are orbital energies, $a_x$ is the amount of non-local Fock exchange, $(ia|jb)$, $(ia|bj)$, and $(ib|aj)$ are exchange-type and $(ij|ab)$ Coulomb-type two-electron integrals, $(ia|f_{XC}|jb)$ and $(ia|f_{XC}|bj)$ are responses of the exchange-correlation functional.

Eq. (3) can be turned into a pseudo-eigenvalue problem:

$$[(\mathbf{A} + \mathbf{B}) - \omega^2(\mathbf{A}-\mathbf{B})^{-1}]\,[\mathbf x^\omega_\zeta + \mathbf y^\omega_\zeta] = -2\eta_\zeta,$$

which is solved to get the linear response vectors for a given time-dependent perturbation.

It is also customary to consider the case when $\eta = 0$, which lead to the following pseudo-hermitian problem:

$$\tag{4}\begin{pmatrix}
\mathbf A & \mathbf B \\
\mathbf B & \mathbf A
\end{pmatrix}\,\begin{pmatrix}
\mathbf x_\zeta^\omega\\ \mathbf y_\zeta^\omega
\end{pmatrix}=\omega\begin{pmatrix}
\mathbf 1 & \mathbf 0\\
\mathbf 0 & -\mathbf 1
\end{pmatrix}\,\begin{pmatrix}
\mathbf x_\zeta^\omega\\ \mathbf y_\zeta^\omega
\end{pmatrix}$$

which is generally referred to as the Casida equation [thought Eq. (3) may also be called that]. In this case, the $\omega$'s are the excitation energies while $\mathbf x^\omega$ and $\mathbf y^\omega$ might be seen as excitation and de-excitation vectors. 
Indeed, a common simplification to Eq. (4) is the **Tamm-Dancoff approximation**, where $\mathbf B = \mathbf 0$, and therefore:

$$\tag{5}\mathbf A\,\mathbf x^\omega = \omega\,\mathbf x^\omega,$$

where the $\mathbf x_{ia}^\omega$ is the coefficient associated to the $i \to a$ transition.

Solving Eq. (4) is done through another eigenvalue problem, namely:

$$(\mathbf{A}-\mathbf{B})^\frac{1}{2}\,(\mathbf{A}-\mathbf{B})\,(\mathbf{A}-\mathbf{B})^\frac{1}{2}\,\mathbf{Z} = \omega^2\,\mathbf{Z}, \text{ with } \mathbf{Z} = (\mathbf{A}-\mathbf{B})^\frac{1}{2} (\mathbf x + \mathbf y).$$

## The simplified approaches to TD-DFT

The simplified TD-DFT methods root in 3 approximations:

1. all integrals involving the XC-functionals are neglected,
2. the singly excited configuration space is truncated (see below), and
3. the [zero-differential overlap](https://en.wikipedia.org/wiki/Zero_differential_overlap) (ZDO) approximation is used for two-electron integrals which built $\mathbf A$ and $\mathbf B$. Different approximations defines different flavors of simplified TD-DFT (see below).

The truncation of the CI space is done in two steps:

1. an active MO space is defined by $\varepsilon_p \in [\varepsilon_{LUMO}-E_{w}, \varepsilon_{HOMO}+E_{w}]$, with $E_w = 2\,(1+0.8a_x)\,E_{thr}$, and then
2. configuration state functions (CSF) are selected within this active space, either if $A_{ia,ia} < E_{thr}$ (primary CSF, P-CSFs) or if $E^{(2)}_{jb} > E_{pt}$ (secondary CSF, S-CFSs), where
   $E^{(2)}_{jb} = \sum_{ia}^{P-CSFs} \frac{|A_{ia,jb}|}{A_{ia,ia}-A_{jb,jb}}$.


### sTD-DFT

To evaluate the integrals, the following formula is used:

$$(ia|jb) \approx \sum_{AB}^N Q_A^{ia}\,Q_B^{jb}(AA|BB), \text{ with } Q_A^{ia} = \sum_{\mu\in A} (C^{\perp}_{i\mu})^\star\,C^{\perp}_{a\mu},$$

where $Q_{ia}^A$ are the transition charges on atom A, computed from the Löwdin orthogonalized LCAO coefficients, $C^\perp = C\,S^{1/2}$.

??? note "Simplification of the integrals using the ZDO approximation"
    
    Starting from the definition of a 4-center integral:

    $$(ia|jb) = \sum_{\mu\nu\alpha\beta} (C_{i\mu})^\star\,C_{a\nu}\,(C_{j\alpha})^\star\,C_{b\beta}\,(\mu\nu|\alpha\beta),$$

    one can instead use Löwdin orthogonalized molecular orbitals, $C^\perp = C\,S^{1/2}$:

    $$(ia|jb) = \sum_{\mu\nu\alpha\beta} (C^\perp_{i\mu})^\star\,C^\perp_{a\nu}\,(C^\perp_{j\alpha})^\star\,C^\perp_{b\beta}\,(\lambda_\mu\lambda_\nu|\lambda_\alpha\lambda_\beta).$$

    Now, the [ZDO approximation](https://en.wikipedia.org/wiki/Zero_differential_overlap) impose that $\lambda_\mu\,\lambda_\nu = \delta_{\mu\nu}\lambda_\mu\lambda_\mu$, so:
    
    $$(ia|jb) \approx \sum_{\mu\nu} (C^\perp_{i\mu})^\star\,C^\perp_{a\mu}\,(C^\perp_{j\nu})^\star\,C^\perp_{b\nu}\,(\lambda_\mu\lambda_\mu|\lambda_\nu\lambda_\nu).$$

    Finally, it is assumed that $\lambda_\mu\lambda_\nu = \nu\nu$, and therefore:

    $$(ia|jb) \approx \sum_{\mu\nu} (C^\perp_{i\mu})^\star\,C^\perp_{a\mu}\,(C^\perp_{j\nu})^\star\,C^\perp_{b\nu}\,(\mu\mu|\nu\nu).$$

    Further simplification arise from population analysis, which defines a *transition charge*:

    $$Q_A^{ia} = \sum_{\mu\in A} (C^{\perp}_{i\mu})^\star\,C^{\perp}_{a\mu},$$

    which allows to arrive at the final expresssion.

The remaining $(AA|BB)$ integrals are Mataga–Nishimoto–Ohno–Klopman (MNOK) damped Coulomb operators, which are evaluated according to the type of bielectronic integral that they approximates:

+ For Coulomb-type integrals, $(ij|ab)$,
  
$$(AA|BB)_J = \left[\frac{1}{R_{AB}^{\gamma_J}+\left(a_x\,\eta_{AB}\right)^{-\gamma_J}}\right]^{1/\gamma_J}.$$

+ For exchange-type integrals, $(ia|jb)$,

$$(AA|BB)_K = \left[\frac{1}{R_{AB}^{\gamma_J}+\eta_{AB}^{-\gamma_J}}\right]^{1/\gamma_J}.$$

In both cases, $\eta_{AB} = \frac{1}{2}\,(\eta_A + \eta_B)$ where $\eta_A$ depends on the chemical hardness of A, while $\gamma_J$ and $\gamma_K$ are globally fitted parameters.  

### XsTD-DFT

!!! warning
    
    The publication describing the XsTD-DFT implementation is not yet available. Thus, this approach is not yet implemented.

## Sources and references

+ J. Toulouse, [Introduction to the calculation of molecular properties by response theory](https://www.lct.jussieu.fr/pagesperso/toulouse/enseignement/molecular_properties.pdf) (last consultation: January 2023). 
+ E. Fromager, [Linear response time-dependent density functional theory](https://quantique.u-strasbg.fr/lib/exe/fetch.php?media=fr:pageperso:ef:lecture_rctf_tddft_e_fromager.pdf)  (last consultation: January 2023).
+ M. E. Casida, Time-Dependent Density Functional Response Theory for Molecules. In D. E. Chong (ed.), *Recent Advances in Density Functional Methods*. World Scientific, **1995** ([10.1142/9789812830586_0005](https://doi.org/10.1142/9789812830586_0005)).
+ S. Hirata, M. Head-Gordon, Time-dependent density functional theory within the Tamm–Dancoff approximation. *Chem. Phys. Lett.*, **1999**, 314, 291 ([10.1016/S0009-2614(99)01149-5](https://doi.org/10.1016/S0009-2614(99)01149-5))
+ S. Löffelsender, P. Beaujean, M. de Wergifosse, Simplified quantum chemistry methods to evaluate non-linear optical properties of large systems. *WIREs Comput. Mol. Sci.* **2023**, 2023, e1695 ([10.1002/wcms.1695](https://dx.doi.org/10.1002/wcms.1695)). 
+ S. Grimme, A simplified Tamm-Dancoff density functional approach for the electronic excitation spectra of very large molecules. *J Chem Phys.* **2013**, 138, 244104 ([10.1063/1.4811331](https://doi.org/10.1063/1.4811331)).
+ C. Bannwarth, S. Grimme, A simplified time-dependent density functional theory approach for electronic ultraviolet and circular dichroism spectra of very large molecules. *Comput Theor Chem.* **2014**, 1040, 45 ([10.1016/j.comptc.2014.02.023](https://doi.org/10.1016/j.comptc.2014.02.023)).