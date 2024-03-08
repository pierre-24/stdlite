title: Theory

## Introduction

Molecular properties are defined as a response of a (molecular) system to an external perturbation (change of geometry, application of an external electric field, etc).
The goal of `stdlite` is to compute properties within the simplified TDA/TD-DFT framework.

!!! note "TL;DR:"

    The following steps are required to perform a sTD-DFT calculation:
    
    1. Extract data (*i.e.*, the atomic orbitals, as well as the MO energies, $\varepsilon$, and LCAO coefficients, $\mathbf C$) from a QC calculation;
    2. Select configuration state functions (CSFs, $\ket{\Psi_i^a}$, or singly excited determinants) using [a set of rules](#the-simplified-approaches-to-td-dft) and build the corresponding electronic hessian (super-)matrices $\mathbf A'$ and $\mathbf B'$ (might be zero);
    3. Using said matrices, compute the [linear response](#linear-response) or [amplitude](#excitations) vectors. This generally requires to compute extra [expectation values](#expectation-values-and-others) matrices in MO basis, such as the [dipole moment](#dipole-moment);
    4. Use said response/amplitude vectors to compute actual [properties](#properties).


## Response function theory

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

## Application to DFT: TD-DFT

!!! note

    In the following, $p$, $q$, $r$, $s$,... refer to molecular orbitals (MOs), $i$, $j$, $k$, $l$,... to occupied, $a$, $b$, $c$, $d$,... to unoccupied ones, $\ket{m}$, $\ket{n}$,... to electronic excited states, $\mu$, $\nu$, $\xi$... to atomic orbitals (AOs) and $\zeta$, $\sigma$, $\tau$, $\upsilon$, ... to cartesian direction. 
    All quantities are in [atomic units](https://en.wikipedia.org/wiki/Atomic_units), unless otherwise mentioned.


### Linear response

Under the conditions of the [Runge and Gross theorem](https://en.wikipedia.org/wiki/Runge%E2%80%93Gross_theorem), this theory can be applied to DFT.
For a time-dependent perturbation at frequency $\omega$, and assuming the Hermicity of the different matrices and real orbitals, Eq. (2) can be written (the superscripts $\lambda$ have been dropped for clarity):

$$\tag{3}\left[\begin{pmatrix}
\mathbf A & \mathbf B \\
\mathbf B & \mathbf A
\end{pmatrix}-\omega\begin{pmatrix}
\mathbf 1 & \mathbf 0\\
\mathbf 0 & -\mathbf 1
\end{pmatrix}\right]\,\begin{pmatrix}
\mathbf x_\zeta(\omega)\\ \mathbf y_\zeta(\omega)
\end{pmatrix}=-\begin{pmatrix}
\mathbf \eta_\zeta\\ \mathbf \eta_\zeta
\end{pmatrix},$$

where $\mathbf x_\zeta(\omega)$ and $\mathbf y_\zeta(\omega)$ are the frequency-dependent linear response vectors (to be determined) in direction $\zeta$.
The $\mathbf A$ and $\mathbf B$ are electronic Hessian (super-)matrices (related to orbital rotations).
The perturbed electronic gradient vector elements $\eta_{ia,\zeta} = \braket{i|\hat\eta_\zeta|a}$ are elements of the expectation value matrix corresponding to the perturbation (e.g., when the perturbation is an electric field, $\hat\eta$ corresponds to the dipole moment operator).

To solve this problem, Eq. (3) can be turned into a linear equation of the form:
    
$$[(\mathbf{A} + \mathbf{B}) - \omega^2(\mathbf{A}-\mathbf{B})^{-1}]\,[\mathbf x_\zeta(\omega) + \mathbf y_\zeta(\omega)] = -2\mathbf\eta_\zeta.$$

??? note "Detailed solution"
    
    The previous equation is easier seen as a linear system written in the following form:

    $$\mathbf L(\omega)\,\mathbf u_{\zeta}(\omega) = -2\mathbf\eta_\zeta,$$
    
    where $\mathbf L(\omega) = (\mathbf{A} + \mathbf{B}) - \omega^2(\mathbf{A}-\mathbf{B})^{-1}$ and $\mathbf u_{\zeta}(\omega)  = \mathbf x_\zeta(\omega) + \mathbf y_\zeta(\omega)$, and which is solved using any of the usual methods for linear systems (worst case scenario: $\mathbf u_\zeta(\omega) = -2\mathbf L^{-1}(\omega)\,\eta_\zeta$).
    Then, since, from Eq. (3),
    
    $$\mathbf x_\zeta(\omega) - \mathbf y_\zeta(\omega) = \omega\,(\mathbf A-\mathbf B)^{-1}\,[\mathbf x_\zeta(\omega) + \mathbf y_\zeta(\omega)],$$
    
    one can define $\mathbf v_{\zeta}(\omega)$ as:

    $$\mathbf v_{\zeta}(\omega)  =\mathbf x_\zeta(\omega) - \mathbf y_\zeta(\omega) =  \omega\,(\mathbf A-\mathbf B)^{-1}\,\mathbf u_{\zeta}(\omega),$$
    
    the response vector are obtained: $\mathbf x_{\zeta}(\omega) = \frac{1}{2}[\mathbf u_{\zeta}(\omega)  + \mathbf v_{\zeta}(\omega)]$ and $\mathbf y_{\zeta}(\omega) = \frac{1}{2}[\mathbf u_{\zeta}(\omega)  - \mathbf v_{\zeta}(\omega)]$.

The [Tamm-Dancoff approximation](https://doi.org/10.1016/S0009-2614(99)01149-5) (setting $\mathbf B = \mathbf 0$ in all previous equations) can also be used.

### Excitations

Since the pole of the response function corresponds to the electronic excited states, it is also customary to consider the case when $\eta = 0$, which lead to the following pseudo-hermitian problem:

$$\tag{4}\begin{pmatrix}
\mathbf A & \mathbf B \\
\mathbf B & \mathbf A
\end{pmatrix}\,\begin{pmatrix}
\mathbf x^m\\ \mathbf y^m
\end{pmatrix}=\omega_m\begin{pmatrix}
\mathbf 1 & \mathbf 0\\
\mathbf 0 & -\mathbf 1
\end{pmatrix}\,\begin{pmatrix}
\mathbf x^m\\ \mathbf y^m
\end{pmatrix}$$

which is generally referred to as the Casida equation. 
In this case, each $\omega_m$ (eigenvalue, corresponding to the transition energy bewteen $\ket{0}$ and $\ket{m}$) is associated to one $\mathbf x^m$ and one $\mathbf y^m$ ("eigenfunction"), might be seen as amplitude vectors associated to excitation and de-excitation, respectively. 
Solving this problem is done using two approaches.

On the one hand, Eq. (4) can be rewritten in a true eigenvalue problem, namely:

$$(\mathbf{A}-\mathbf{B})^\frac{1}{2}\,(\mathbf{A}+\mathbf{B})\,(\mathbf{A}-\mathbf{B})^\frac{1}{2}\,\mathbf{z}^m = \omega^2\,\mathbf{z}^m, \text{ with } \mathbf{z}^m = (\mathbf{A}-\mathbf{B})^{-\frac{1}{2}} (\mathbf x^m + \mathbf y^m).$$

??? note "Detailed solution"

    In this case, after $\mathbf z^m$ have been obtained, one extract using the following procedure. First, from previous expression, one can obtain:
    
    $$ \mathbf u^m = \mathbf x^m + \mathbf y^m = \frac{1}{\sqrt\omega}\,(\mathbf A-\mathbf B)^\frac{1}{2}\,\mathbf z^m.$$
    
    Now, since $(\mathbf A + \mathbf B)\,(\mathbf x^m+\mathbf y^m) = \omega\,(\mathbf x^m-\mathbf y^m)$ [obtained from Eq. (4)], one has:
    
    $$\mathbf v^m = \mathbf x^m-\mathbf y^m = \frac{1}{\omega}\,(\mathbf A + \mathbf B)\,(\mathbf x^m+\mathbf y^m)  = \frac{1}{\omega}\,(\mathbf A + \mathbf B)\,\mathbf u^m,$$
    
    and therefore the response vector are obtained: $\mathbf x^m = \frac{1}{2}(\mathbf u^m + \mathbf v^m)$ and $\mathbf y^m = \frac{1}{2}(\mathbf u^m - \mathbf v^m)$.

On the other hand, the [Tamm-Dancoff approximation](https://doi.org/10.1016/S0009-2614(99)01149-5) ($\mathbf B = \mathbf 0$ and thus $\mathbf y^m = 0$) simply leads to:

$$\mathbf A\,\mathbf x^m = \omega_m\,\mathbf x^m.$$

In this case, $\ket{m} = \sum_{ia}\,x^m_{ia}\,\ket{\Psi_i^a}$, where $\ket{\Psi_i^a}$ is a **configuration state function**, *i.e.*, a singly-excited determinant where an electron has been moved from $i$ to $a$.

In both cases, these amplitude vectors are linked to their linear response counterparts through the following spectral representation:

$$\begin{aligned}
x_{ia,\zeta}(\omega) &= \sum_{\ket{m}} \eta_{ia,\zeta}\,(x^{m}_{ia} + y^{m}_{ia})\,\left[\frac{x_{ia}^{m}}{\omega-\omega_m}-\frac{y_{ia}^{m}}{\omega+\omega_m}\right],\\
y_{ia,\zeta}(\omega) &= \sum_{\ket{m}} \eta_{ia,\zeta}\,(x^{m}_{ia} + y^{m}_{ia})\,\left[\frac{y_{ia}^{m}}{\omega-\omega_m}-\frac{x_{ia}^{m}}{\omega+\omega_m}\right],
\end{aligned}$$

where these expression involves a summation over the manifold $\{\ket{m}\}$ of excited states (and one can set $\mathbf y^m = 0$ to get the TDA version).
These representations lead to simplification when taking residue of response functions.

## The simplified approaches to TD-DFT

In the rest of this development a **global hybrid** density functional is assumed,

$$E_{XC}= (1-a_x)\,E_X^{GGA}+a_x\,E_X^{HF}+E_C^{GGA}.$$

Thanks to the [Slater-Condon rules](https://en.wikipedia.org/wiki/Slater%E2%80%93Condon_rules), one can evaluate the elements of $\mathbf A$ and $\mathbf B$, which are:

$$\begin{aligned}
&A_{ia, jb} = \delta_{ij}\delta_{ab} (\varepsilon_a - \varepsilon_i) + 2\,(ia|jb) - a_x\,(ij|ab) + (1-a_x)\,(ia|f_{XC}|jb),\\
&B_{ia,jb} = 2\,(ia|bj) - a_x\,(ib|aj) + (1-a_x)\,(ia|f_{XC}|bj),
\end{aligned}$$

where, $\varepsilon_i$ and $\varepsilon_a$ are orbital energies, $a_x$ is the amount of non-local Fock exchange, $(ia|jb)$, $(ia|bj)$, and $(ib|aj)$ are exchange-type and $(ij|ab)$ Coulomb-type two-electron integrals, $(ia|f_{XC}|jb)$ and $(ia|f_{XC}|bj)$ are responses of the exchange-correlation functional.

The simplified TD-DFT methods root in 3 approximations which simplify the content of $\mathbf A$ and $\mathbf B$:

1. all integrals involving the XC-functionals are neglected (referred to as the [random phase approximation (RPA)](https://en.wikipedia.org/wiki/Random_phase_approximation) approach),
2. the singly excited configuration space is truncated (see below), and
3. the [zero-differential overlap](https://en.wikipedia.org/wiki/Zero_differential_overlap) (ZDO) approximation is used for two-electron integrals which built $\mathbf A$ and $\mathbf B$. Different approximations for the remaining integrals define different flavors of simplified TD-DFT (see below).

Then, using those approximated $\mathbf A'$ and $\mathbf B'$ matrices, the linear response or Casida equations presented above are solved.

### Truncation of the active space

The truncation of the CI space is done in three steps:

1. An active MO space is defined by $\varepsilon_p \in [\varepsilon_{LUMO}-E_{w}, \varepsilon_{HOMO}+E_{w}]$, with $E_w = 2\,(1+0.8a_x)\,E_{thr}$.
2. From this active space, primary $\ket{\Psi_i^a}$ configuration state functions (P-CSFs) are selected, for which $A_{ia,ia} < E_{thr}$.
3. Then, from CSFs $\ket{\Psi_j^b}$ for which $A_{jb,jb} > E_{thr}$, a set of secondary CSFs (S-CSFs), for which $E^{(2)}_{jb} > E^{(2)}_{thr}$ is build (typically, $E^{(2)}_{thr} = 10^{-4}$). Other CSFs are discarded.

The selection of S-CSFs is based on a perturbative approach, where $E^{(2)}_{jb}$ measure the cumulative perturbative contributions of a given S-CSF $\ket{\Psi_j^b}$ to all the P-CSFs $\ket{\Psi_i^a}$:

$$E^{(2)}_{jb} = \sum_{ia}^{\text{P-CSFs}} \frac{|A_{ia,jb}|^2}{A_{jb,jb}-A_{ia,ia}}.$$

### sTD-DFT, or the monopole approximation

To evaluate the integrals, a formula of the following kind is used:

$$(ia|jb) \approx (ia|jb)' = \sum_{AB}^N q_A^{ia}\,q_B^{jb}(AA|BB), \text{ with } q_A^{ia} = \sum_{\mu\in A} (C^{\perp}_{i\mu})^\star\,C^{\perp}_{a\mu},$$

where the $q_{ia}^A$'s are the transition charges on atom A, computed using Löwdin orthogonalized LCAO coefficients, $C^\perp_{i\mu} = \sum_\nu C_{i\nu}\,S^{1/2}_{\nu\mu}$.

??? note "Details of the application of the ZDO approximation"
    
    Starting from the definition of a 4-center integral:

    $$(ia|jb) = \sum_{\mu\nu\alpha\beta} (C_{i\mu})^\star\,C_{a\nu}\,(C_{j\alpha})^\star\,C_{b\beta}\,(\mu\nu|\alpha\beta),$$

    one can instead use Löwdin orthogonalized molecular orbitals, $C^\perp = C\,S^{1/2}$:

    $$(ia|jb) = \sum_{\mu\nu\alpha\beta} (C^\perp_{i\mu})^\star\,C^\perp_{a\nu}\,(C^\perp_{j\alpha})^\star\,C^\perp_{b\beta}\,(\lambda_\mu\lambda_\nu|\lambda_\alpha\lambda_\beta),$$

    where $\lambda_\mu = \sum_\nu S^{-1/2}_{\mu\nu}\,\phi_\nu$.
    Now, the [ZDO approximation](https://en.wikipedia.org/wiki/Zero_differential_overlap) impose that $\lambda_\mu\,\lambda_\nu = \delta_{\mu\nu}\lambda_\mu\lambda_\mu$, so:
    
    $$(ia|jb) \approx \sum_{\mu\nu} (C^\perp_{i\mu})^\star\,C^\perp_{a\mu}\,(C^\perp_{j\nu})^\star\,C^\perp_{b\nu}\,(\lambda_\mu\lambda_\mu|\lambda_\nu\lambda_\nu).$$

    Finally, it is assumed that $\lambda_\mu\lambda_\mu \approx \mu\mu$, and therefore:

    $$(ia|jb) \approx \sum_{\mu\nu} (C^\perp_{i\mu})^\star\,C^\perp_{a\mu}\,(C^\perp_{j\nu})^\star\,C^\perp_{b\nu}\,(\mu\mu|\nu\nu).$$

    A simplification of the notation arise from population analysis, which allows to defines a *transition charge*:

    $$q_A^{ia} = \sum_{\mu\in A} (C^{\perp}_{i\mu})^\star\,C^{\perp}_{a\mu},$$

    leading to the final expresssion.

The remaining $(AA|BB)$ integrals are evaluated with Mataga–Nishimoto–Ohno–Klopman (MNOK) damped Coulomb operators, according to the type of bielectronic integral that they approximate:

+ For Coulomb-type integrals, $(ij|ab)$,
  
$$(AA|BB)_J = \left[\frac{1}{R_{AB}^{\gamma_J}+\left(a_x\,\eta_{AB}\right)^{-\gamma_J}}\right]^{1/\gamma_J}.$$

+ For exchange-type integrals, $(ia|jb)$,

$$(AA|BB)_K = \left[\frac{1}{R_{AB}^{\gamma_K}+\eta_{AB}^{-\gamma_K}}\right]^{1/\gamma_K}.$$

In both cases, $\eta_{AB} = \frac{1}{2}\,(\eta_A + \eta_B)$ where $\eta_A$ is the chemical hardness of A (obtained from [here](https://dx.doi.org/10.1002/qua.22202)), while $\gamma_J$ and $\gamma_K$ are globally fitted parameters.  

In practice, the remaining elements of the (now approximated) electronic Hessian matrices $\mathbf A'$ and $\mathbf B'$ are:

$$\begin{aligned}
A'_{ia,jb} =& \delta_{ij}\delta_{ab} (\varepsilon_a - \varepsilon_i) + 2\,(ia|jb)' - (ij|ab)',\\ 
B'_{ia,jb} =& 2\,(ia|bj)' -a_x\,(ib|aj)'.
\end{aligned}$$

??? note "Implementation detail"
    In order to be computationally efficient, these element can be evaluated by precomputing three kind of transition charges: $q_A^{ij}$, $q_A^{ia}$, and $q_A^{ab}$, and two intermediates:
    
    $$(ij|BB)_J = \sum_A^N q_A^{ij}\,(AA|BB)_J \text{ and } (ia|BB)_K = \sum_A^N q_A^{ia}\,(AA|BB)_K,$$
    
    so that a scalar product leads to the value of the different integrals. For example,
    
    $$(ia|jb)' \approx \sum_{B}^N (ia|BB)_K\,q_B^{jb}.$$
    
    However, these intermediates have a huge memory cost [$\mathcal O (N_{atm}\times N_{MO}^2)$ altogether], so a direct version (integrals are evaluated when requested) is also available.


### XsTD-DFT

!!! warning
    
    The publication describing the XsTD-DFT implementation is not yet available. Thus, this approach is not yet implemented.

## Expectation values and others

Things that can be obtained without solving the linear response equation.

!!! note
    
    When required, conversion of $X$ from the AO to the MO basis is done using:

    $$X^{MO}_{pq} = \sum^{AO}_{\mu\nu} C_{p\mu} X_{\mu\nu} C_{q\nu}.$$

### Density matrix

The density matrix elements are defined as:

$$P_{\mu\nu} = \sum^{MO}_p n_p\,C^\star_{p\mu}\,C_{p\nu},$$

where $n_p$ is the occupancy of $p$. If this is a closed shell calculation, occupied MO have an occupancy of 2, 0 otherwise.

### Dipole moment

The total dipole moment is obtained as:

$$\vec\mu = \vec\mu_e + \sum^{N_{nuc}}_A Z_A\,\vec r_A,$$

where $Z_A$ and $\vec r_A$ are the charge and position of nuclei $A$, respectively.
The electronic dipole moment, $\vec\mu_e$, is an expectation value of the wavefunction.
The electric dipole moment matrix elements (in AO basis) are computed as:

$$D_{\mu\nu} = \braket{\mu|\hat\mu|\nu} = -\,\braket{\mu|\vec{r}-\vec R_0|\nu},$$

where the dipole moment integral have been multiplied by the value of the electronic charge (-1 in atomic units) and $\vec R_0$ is the origin.
The electronic dipole moment is computed via:

$$\vec\mu_e = \sum^{AO}_{\mu\nu} P_{\mu\nu}\,D_{\nu\mu} = tr(\mathbf P\mathbf D),$$

where $\mathbf P$ is the density matrix. Alternatively, in MO basis:

$$\vec\mu_e = \sum^{MO}_p n_p\,\vec\mu_{pp},$$

where $\vec\mu_{pp}$ is a shorthand for $D^{MO}_{pp}$.

## Responses (properties)

Things that can be obtained thanks to the amplitude/linear response vectors.

!!! note
    
    Linear and quadratic response functions are generally noted $\braket{\braket{\hat A; \hat B}}_{\omega_B}$ and $\braket{\braket{\hat A; \hat B, \hat C}}_{\omega_B,\omega_C}$, which describes how the expectation value of $\hat A$ (at frequency $\omega_A = -\omega_B - \omega_C$) responds to a set of perturbation to first and second order in perturbation.
    Residue of the response functions provide information on the (excited states of the) unperturbed system.

    For example, the linear response function might be extended in the following spectral representation:

    $$\braket{\braket{\hat A; \hat B}}_{\omega_B} = \sum_{\ket{m}} \frac{\braket{0|\hat A|m}\braket{m|\hat B|0}}{\omega_B-\omega_m} + \frac{\braket{0|\hat B|m}\braket{m|\hat A|0}}{\omega_A+\omega_m},$$

    and a corresponding single residue might be:

    $$\lim_{\omega_B\to\omega_m} (\omega_B-\omega_m)\,\braket{\braket{\hat A; \hat B}}_{\omega_B} = \braket{0|\hat A|m}\braket{m|\hat B|0},$$

    which provides access to trasition matrix elements $\braket{0|\hat A|m}$ between the ground state $\ket{0}$ and an excited state $\ket{m}$. 
    In practice, such residues are thus evaluated thanks to amplitude vectors through the spectral representation of linear responses.

### Polarizability

The dynamic (electric) polarizability tensor elements are defined:

$$\alpha_{\zeta\sigma}(-\omega;\omega) = -\braket{\braket{\hat\mu_\zeta;\hat\mu_\sigma}}_\omega = -2\,\sum_{ia}^{CSF} \mu_{ia,\zeta}\,[x_{ia,\sigma}(\omega)+y_{ia,\sigma}(\omega)].$$

### Ground to excited transition dipole moment (and oscillator strength)

From the single residue of the electric polarizability, 

$$\lim_{\omega\to\omega_m} (\omega-\omega_m)\,\braket{\braket{\mu_\zeta;\mu_\sigma}}_\omega =  \braket{0|\hat\mu_\zeta|m}\,\braket{m|\hat\mu_\sigma|0},$$

the transition dipole moment elements (in the dipole length formalism) between ground $\ket{0}$ and excited state $\ket{m}$, associated with energy $\omega_{m}$ are given by:

$$\mu_{0m,\zeta} = \braket{0|\hat\mu_\zeta|m} =  \sqrt{2}\,\sum_{ia}^{CSF} \vec\mu_{ia,\zeta}\,(x^m_{ia}+y^m_{ia}),$$

where $\zeta$ is a cartesian direction and the $\sqrt 2$ factor is required for singlet excitations (it is 0 for triplet).
The associated [oscillator strength](https://en.wikipedia.org/wiki/Oscillator_strength) is defined by:

$$f_{0m} = \frac{2}{3}\,\omega_m\,|\vec\mu_{0m}|^2.$$

### First hyperpolarizability

By neglecting the response of the XC kernel and the Hartree XC kernel, the element of the first hyperpolarizability tensor in the sTD-DFT framework are:

$$\beta_{\zeta\sigma\tau}(\omega_\varsigma;\omega_1,\omega_2) = -\braket{\braket{\hat\mu_\zeta;\hat\mu_\sigma,\hat\mu_\tau}}_{\omega_1,\omega_2} = \mathcal A_{\zeta\sigma\tau}(\omega_\varsigma;\omega_1,\omega_2) - \mathcal B_{\zeta\sigma\tau}(\omega_\varsigma;\omega_1,\omega_2),$$

with $\omega_\varsigma = -\omega_1 - \omega_2$, and:

$$\begin{aligned}
\mathcal A_{\zeta\sigma\tau}(\omega_\varsigma;\omega_1,\omega_2) &= \sum_{\mathcal P} \sum_{ia,ja} x_{ia,\zeta}(\omega_\varsigma)\,[-\mu_{ij,\sigma}]\,y_{ja,\tau}(\omega_2), \\
\mathcal B_{\zeta\sigma\tau}(\omega_\varsigma;\omega_1,\omega_2) &= \sum_{\mathcal P} \sum_{ia,ib} x_{ia,\zeta}(\omega_\varsigma)\,[-\mu_{ab,\sigma}]\,y_{ib,\tau}(\omega_2),
\end{aligned}$$

where $\sum_{\mathcal P}$ is the sum over the sequence of permutations of the pairs of components and energies, $\{(\zeta,\omega_\varsigma), (\sigma, \omega_1), (\tau,\omega_2)\}$.

### Excited to excited transition dipole moment (and oscillator strength)

From the double residue of the electric first hyperpolarizability:

$$\begin{aligned}
\lim_{\omega_1\to-\omega_m,\omega_2\to\omega_n} &(\omega_1+\omega_m)\,(\omega_2-\omega_n)\,\braket{\braket{\hat\mu_\zeta;\hat\mu_\sigma,\hat\mu_\tau}}_{\omega_1,\omega_2}\\
&=  \braket{0|\hat\mu_\zeta|m}\,\braket{m|\hat\mu_\sigma - \delta_{mn}\,\braket{0|\hat\mu_\sigma|0}|n}\,\braket{n|\hat\mu_\tau|0},
\end{aligned}$$

one extract the (unrelaxed) singlet-to-singlet transition dipole moment from $\ket{m}$ to $\ket{n}$ (see [there](https://github.com/pierre-24/stdlite/blob/master/docs/assets/esa.pdf)):

$$\begin{aligned}
&\braket{m|\hat\mu_\zeta - \delta_{mn}\,\braket{0|\hat\mu_\zeta|0}|n}\\
&\hspace{2em}= \frac{1}{2}\left\{ \sum_{ia,ib} \mu_{ab,\zeta}\,[x^n_{ia}\,x^m_{ib} + y^m_{ia}\,y^n_{ib}]  - \sum_{ia,ja} \mu_{ij,\zeta}\,[x^n_{ia}\,x^m_{ja} + y^m_{ia}\,y^n_{ja}] \right\}.
\end{aligned}$$

which is equal to the fluctuation operator if $m = n$.
The corresponding oscillator strength is given by:

$$f_{mn} = \frac{2}{3}\,(\omega_n-\omega_m)\,(\vec\mu_{mn}\cdot\vec\mu_{nm}).$$

## Sources and references

+ J. Toulouse, [Introduction to the calculation of molecular properties by response theory](https://www.lct.jussieu.fr/pagesperso/toulouse/enseignement/molecular_properties.pdf) (last consultation: January 2023). 
+ E. Fromager, [Linear response time-dependent density functional theory](https://quantique.u-strasbg.fr/lib/exe/fetch.php?media=fr:pageperso:ef:lecture_rctf_tddft_e_fromager.pdf)  (last consultation: January 2023).
+ G. P. Chen, V. K. Voora, M. M. Agee, S. G. Balasubramani, and F. Furche, Random-Phase Approximation Methods. *Annu. Rev. Phys. Chem.* **2017**, 68, 421 ([10.1146/annurev-physchem-040215-112308](https://doi.org/10.1146/annurev-physchem-040215-112308))
+ M. E. Casida, Time-Dependent Density Functional Response Theory for Molecules. In D. E. Chong (ed.), *Recent Advances in Density Functional Methods*. World Scientific, **1995** ([10.1142/9789812830586_0005](https://doi.org/10.1142/9789812830586_0005)).
+ S. Hirata, M. Head-Gordon, Time-dependent density functional theory within the Tamm–Dancoff approximation. *Chem. Phys. Lett.*, **1999**, 314, 291 ([10.1016/S0009-2614(99)01149-5](https://doi.org/10.1016/S0009-2614(99)01149-5))
+ S. Löffelsender, P. Beaujean, M. de Wergifosse, Simplified quantum chemistry methods to evaluate non-linear optical properties of large systems. *WIREs Comput. Mol. Sci.* **2023**, 2023, e1695 ([10.1002/wcms.1695](https://dx.doi.org/10.1002/wcms.1695)) [and references therein]. 
+ S. Grimme, A simplified Tamm-Dancoff density functional approach for the electronic excitation spectra of very large molecules. *J Chem Phys.* **2013**, 138, 244104 ([10.1063/1.4811331](https://doi.org/10.1063/1.4811331)).
+ C. Bannwarth, S. Grimme, A simplified time-dependent density functional theory approach for electronic ultraviolet and circular dichroism spectra of very large molecules. *Comput Theor Chem.* **2014**, 1040, 45 ([10.1016/j.comptc.2014.02.023](https://doi.org/10.1016/j.comptc.2014.02.023)).