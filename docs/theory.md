title: Theory

## Introduction

Molecular properties are defined as a response of a (molecular) system to an external perturbation (change of geometry, application of an external electric field, etc).
The goal of `stdlite` is to compute properties within the simplified TDA/TD-DFT framework.

!!! note "TL;DR:"

    The following steps are required to perform a sTD-DFT calculation:
    
    1. Extract data (*i.e.*, the atomic orbitals, as well as the MO energies, $\varepsilon$, and LCAO coefficients, $\mathbf C$) from a QC calculation;
    2. Select configuration state functions (CSFs, $\ket{\Psi_i^a}$, or singly excited determinants) using [a set of rules](#the-simplified-approaches-to-td-dft) and build the corresponding electronic hessian (super-)matrices $\mathbf A'$ and $\mathbf B'$ (might be zero);
    3. Using said matrices, compute the [linear response](#linear-response) or [amplitude](#excitations) vectors. This generally requires to compute extra [AO integrals](#ao-integrals);
    4. Use said response/amplitude vectors to compute actual [properties](##linear-and-quadratic-response-functions).


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
\Leftrightarrow \braket{\braket{\hat V; \hat V}} = ({}^\lambda\kappa)^T\,{}^\lambda\eta$$

where $\braket{\braket{\hat V; \hat V}}$ is the linear response function, ${}^\lambda\eta$ is the perturbed electronic gradient, and ${}^\lambda\kappa$, the first-order response vector to $\lambda$, containing all $\partial\kappa_i/\partial\lambda$.
${}^\lambda\kappa$ is computed thanks to the following set of equations (derivative of the stationary condition):

$$\tag{2}\sum_j \left.\frac{\partial^2 E(\{\kappa_i\},\lambda)}{\partial\kappa_i\partial\kappa_j}\right|_{\lambda=0}\,\frac{\partial\kappa_j}{\partial\lambda} = -\left.\frac{\partial^2 E(\{\kappa_i\},\lambda)}{\partial\lambda\partial\kappa_i}\right|_{\lambda=0}  \Leftrightarrow {}^\lambda\mathbf E\,{}^\lambda\kappa = -{}^\lambda\eta$$
 
where ${}^\lambda\mathbf E$ is the electronic Hessian. Evaluating the elements of ${}^\lambda\mathbf E$ and ${}^\lambda\eta$ actually result in new kinds of (bielectronic) integrals to evaluate (*vide supra*), while ${}^\lambda\kappa$ has to be computed by solving this equation.
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
For a time-dependent perturbation at frequency $\omega$, and assuming  real orbitals, Eq. (2) can be written (the superscripts $\lambda$ have been dropped for clarity):

$$\tag{3}\left[\begin{pmatrix}
\mathbf A & \mathbf B \\
\mathbf B & \mathbf A
\end{pmatrix}-\omega\begin{pmatrix}
\mathbf 1 & \mathbf 0\\
\mathbf 0 & -\mathbf 1
\end{pmatrix}\right]\,\begin{pmatrix}
\mathbf x_\zeta(\omega)\\ \mathbf y_\zeta(\omega)
\end{pmatrix}=-\begin{pmatrix}
\mathbf \eta_\zeta\\ \mathbf \eta_\zeta^\star
\end{pmatrix},$$

where $\mathbf x_\zeta(\omega)$ and $\mathbf y_\zeta(\omega)$ are the frequency-dependent linear response vectors (to be determined) in direction $\zeta$.
The $\mathbf A$ and $\mathbf B$ are electronic Hessian (super-)matrices (related to orbital rotations), which are the occ-unocc parts of ${}^\lambda\mathbf E$.
The perturbed electronic gradient vector elements $\eta_{ia,\zeta} = \braket{i|\hat\eta_\zeta|a}$ are MO integral elements of the operator corresponding to the perturbation (e.g., when the perturbation is an electric field, $\hat\eta$ corresponds to the dipole moment operator, $\hat\mu$).

To solve this problem, Eq. (3) can be turned into a linear equation. 
Adding and subtracting the two lines of the system given in Eq. (3) results in:

$$
\begin{pmatrix}
\mathbf A+ \mathbf B & -\omega \\
-\omega & \mathbf A - \mathbf B
\end{pmatrix}\,\begin{pmatrix}
\mathbf t_\zeta(\omega)\\ \mathbf u_\zeta(\omega)
\end{pmatrix}=-2\,\begin{pmatrix}
\mathbf \eta_\zeta + \mathbf \eta_\zeta^\star \\ \mathbf \eta_\zeta - \mathbf \eta_\zeta^\star
\end{pmatrix},
$$

where $\mathbf t_\zeta(\omega) = \mathbf x_\zeta(\omega)+\mathbf y_\zeta(\omega)$ and $\mathbf u_\zeta(\omega) = \mathbf x_\zeta(\omega)-\mathbf y_\zeta(\omega)$.

In the case where $\eta^\star = \eta \Leftrightarrow \mathbf \eta_\zeta - \mathbf \eta_\zeta^\star = 0$ (hermitian operator), then:
    
$$[(\mathbf{A} + \mathbf{B}) - \omega^2(\mathbf{A}-\mathbf{B})^{-1}]\,\mathbf t_\zeta(\omega) = -2\mathbf\eta_\zeta,$$

which is to be solved, and then:

$$\mathbf u_{\zeta}(\omega) =  \omega\,(\mathbf A-\mathbf B)^{-1}\,\mathbf t_{\zeta}(\omega).$$


??? note "If $\hat\eta$ is anti-hermitian"

    On the other hand, in the case where $\eta = -\eta^\star \Leftrightarrow \mathbf \eta_\zeta + \mathbf \eta_\zeta^\star = 0$ (anti-hermitian operator), then:
    
    $$[(\mathbf{A} - \mathbf{B}) - \omega^2(\mathbf{A}+\mathbf{B})^{-1}]\,[\mathbf x_\zeta(\omega) - \mathbf y_\zeta(\omega)] = -2\Im(\mathbf\eta_\zeta).$$

The [Tamm-Dancoff approximation](https://doi.org/10.1016/S0009-2614(99)01149-5) (setting $\mathbf B = \mathbf 0$ in all previous equations) can also be used.


??? note "Implementation details"

    For efficiency reasons, $\mathbf t_\zeta(\omega)$ and $\mathbf u_\zeta(\omega)$ are actually used by `stdlite`, rather than $\mathbf x_\zeta(\omega)$ and $\mathbf y_\zeta(\omega)$.

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
They satisfy $(\mathbf x^m + \mathbf y^m)\cdot(\mathbf x^m - \mathbf y^m) = |\mathbf x^m|^2-|\mathbf y^m|^2=1$.
Solving this problem is done using two approaches.

On the one hand, Eq. (4) can be rewritten in a true eigenvalue problem, namely:

$$
\begin{pmatrix}
\mathbf A+ \mathbf B & -\omega \\
-\omega & \mathbf A - \mathbf B
\end{pmatrix}\,\begin{pmatrix}
\mathbf t^m\\ \mathbf u^m
\end{pmatrix}= \mathbf 0,
$$

with $\mathbf t^m = \mathbf x^m + \mathbf y^m$ and $\mathbf u^m = \mathbf x^m - \mathbf y^m$.
This turns, according to [10.1063/1.4867271](https://dx.doi.org/10.1063/1.4867271), in:

$$(\mathbf{A}-\mathbf{B})^\frac{1}{2}\,(\mathbf{A}+\mathbf{B})\,(\mathbf{A}-\mathbf{B})^\frac{1}{2}\,\mathbf{z}^m = \omega^2\,\mathbf{z}^m, \text{ with } \mathbf{z}^m = (\mathbf{A}-\mathbf{B})^{-\frac{1}{2}} \mathbf t^m,$$

to be solved.
After $\mathbf z^m$ have been obtained, one can then extract:

$$\begin{aligned}
\mathbf t^m &= \frac{1}{\sqrt\omega}\,(\mathbf A-\mathbf B)^\frac{1}{2}\,\mathbf z^m,\\
\mathbf u^m &= \frac{1}{\omega}\,(\mathbf A + \mathbf B)\,\mathbf t^m.
\end{aligned}
$$

On the other hand, the [Tamm-Dancoff approximation](https://doi.org/10.1016/S0009-2614(99)01149-5) ($\mathbf B = \mathbf 0$ and thus $\mathbf y^m = 0$) simply leads to:

$$\mathbf A\,\mathbf x^m = \omega_m\,\mathbf x^m.$$

In this case, $\ket{m} = \sum_{ia}\,x^m_{ia}\,\ket{\Psi_i^a}$, where $\ket{\Psi_i^a}$ is a **configuration state function**, *i.e.*, a singly-excited determinant where an electron has been moved from $i$ to $a$.

??? note "Implementation details"

    For efficiency reasons, $\mathbf t^m$ and $\mathbf u^m$ are actually used by `stdlite`, rather than $\mathbf x^m$ and $\mathbf y^m$.


### Spectral representations and their implications

These amplitude vectors are linked to their linear response counterparts through what is called their spectral representation. In particular, if $\hat\eta$ is hermitian:

$$\begin{aligned}
x_{ia,\zeta}(\omega) &= \sum_{\ket{m}} \eta_{ia,\zeta}\,(x^{m}_{ia} + y^{m}_{ia})\,\left[\frac{x_{ia}^{m}}{\omega-\omega_m}-\frac{y_{ia}^{m}}{\omega+\omega_m}\right],\\ 
y_{ia,\zeta}(\omega) &= \sum_{\ket{m}} \eta_{ia,\zeta}\,(x^{m}_{ia} + y^{m}_{ia})\,\left[\frac{y_{ia}^{m}}{\omega-\omega_m}-\frac{x_{ia}^{m}}{\omega+\omega_m}\right], 
\end{aligned}$$

which implies:

$$\mathbf x_\zeta(-\omega) = \mathbf y_\zeta(\omega) \land \mathbf y_\zeta(-\omega) = \mathbf x_\zeta(\omega) \Rightarrow \mathbf x_\zeta(0) = \mathbf y_\zeta(0).$$

These expression involves a summation over the manifold $\{\ket{m}\}$ of excited states (and one can set $\mathbf y^m = 0$ to get the TDA version).
Furthermore,

$$\begin{aligned}
t_{ia,\zeta}(\omega) &= x_{ia,\zeta}(\omega) + y_{ia,\zeta}(\omega) =  \sum_{\ket{m}} \eta_{ia,\zeta}\,[t^{m}_{ia}]^2\,\left[\frac{1}{\omega-\omega_m}-\frac{1}{\omega+\omega_m}\right],\\
u_{ia,\zeta}(\omega) &= x_{ia,\zeta}(\omega) - y_{ia,\zeta}(\omega) = \sum_{\ket{m}} \eta_{ia,\zeta}\,t^{m}_{ia}\,u_{ia}^{m}\,\left[\frac{1}{\omega-\omega_m}+\frac{1}{\omega+\omega_m}\right],
\end{aligned}$$

which implies that $\mathbf u_\zeta(\omega)$ is symmetric with respect to a change of the sign of the energy, while $\mathbf t_\zeta(\omega)$ is antisymmetric, since:

$$\mathbf t_\zeta(-\omega) = \mathbf t_\zeta(\omega) \land \mathbf u_\zeta(-\omega) = -\mathbf u_\zeta(\omega) \Rightarrow \mathbf u_\zeta(0) = \mathbf 0.$$

??? note "If $\hat\eta$ is anti-hermitian"

    $$\begin{aligned}
    x_{ia,\zeta}(\omega) &= \sum_{\ket{m}} \eta_{ia,\zeta}\,(x^{m}_{ia} - y^{m}_{ia})\,\left[\frac{x_{ia}^{m}}{\omega-\omega_m}+\frac{y_{ia}^{m}}{\omega+\omega_m}\right],\\
    y_{ia,\zeta}(\omega) &= \sum_{\ket{m}} \eta_{ia,\zeta}\,(x^{m}_{ia} - y^{m}_{ia})\,\left[\frac{y_{ia}^{m}}{\omega-\omega_m}+\frac{x_{ia}^{m}}{\omega+\omega_m}\right],
    \end{aligned}$$

    which implies:

    $$\mathbf x(-\omega) = -\mathbf y(\omega) \land \mathbf y(-\omega) = -\mathbf x(\omega).$$

These representations help in obtaining expressions when taking residues of response functions.

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

??? note "Implementation details"

    For efficiency reasons, $\mathbf A'+\mathbf B'$ and $\mathbf A'-\mathbf B'$ are actually used by `stdlite` rather than $\mathbf A'$ and $\mathbf B'$.

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

??? note "Implementation details"

    In order to be computationally efficient, these element can be evaluated by precomputing three kind of transition charges: $q_A^{ij}$, $q_A^{ia}$, and $q_A^{ab}$, and two intermediates:
    
    $$(ij|BB)_J = \sum_A^N q_A^{ij}\,(AA|BB)_J \text{ and } (ia|BB)_K = \sum_A^N q_A^{ia}\,(AA|BB)_K,$$
    
    so that a scalar product leads to the value of the different integrals. For example,
    
    $$(ia|jb)' \approx \sum_{B}^N (ia|BB)_K\,q_B^{jb}.$$


### XsTD-DFT

!!! warning
    
    The publication describing the XsTD-DFT implementation is not yet available. Thus, this approach is not yet implemented.

## Properties

### AO integrals

The AO integral elements of an operator, $\hat O$, are denoted $O_{\mu\nu} = \braket{\mu|\hat O|\nu}$.
Operator are either hermitian (which results in a [hermitian](https://en.wikipedia.org/wiki/Hermitian_matrix) matrix, $O_{\mu\nu} = O^\star_{\nu\mu}$ with real values on the diagonal) or anti-hermitian (which results in a [skew-symmetric](https://en.wikipedia.org/wiki/Skew-Hermitian_matrix) matrix, $O_{\mu\nu} = -O_{\nu\mu}^\star$).

A conversion of $A$ from the AO to the MO basis is done using:

$$O_{pq} = \sum^{AO}_{\mu\nu} C_{p\mu} O_{\mu\nu} C_{q\nu}.$$

For the moment, the following operators are handled by `stdlite`:

| Operator                  | Expression                                 | Dimensionality | Symmetry      | Hermitian | [Time-Reversal symmetry](https://en.wikipedia.org/wiki/T-symmetry) |
|---------------------------|--------------------------------------------|----------------|---------------|-----------|--------------------------------------------------------------------|
| Dipole length (`dipl`)    | $\hat\mu^L = e\,(\vec r - R_0)$            | 3 (x, y, z)    | Symmetric     | Yes       | Even                                                               |
| Dipole velocity (`dipv`)  | $\hat\mu^V = i\vec\nabla$                  | 3 (x, y, z)    | Antisymmetric | Yes       | Odd                                                                |
| Angular momentum (`angm`) | $\hat m = i(\hat r - R_0)\times\vec\nabla$ | 3 (x, y, z)    | Antisymmetric | Yes       | Odd                                                                |

An [even time-reversal symmetry](https://en.wikipedia.org/wiki/T-symmetry#Even) operator, $\hat A$, do not change upon time reversal.
Indeed, given $\hat\theta$ the so-called [time-reversal operator](https://bohr.physics.berkeley.edu/classes/221/9697/timerev.pdf) so that $\hat\theta\psi(t) = \psi(-t)^\star$,
one has:

$${}^A\theta = \hat\theta\hat A \hat\theta^\dagger = 1.$$

If, on the other hand, an operator $\hat B$ has an odd time-reversal symmetry, then:

$${}^B\theta = \hat\theta\hat B \hat\theta^\dagger = -1.$$

$\hat p = -i\vec\nabla$, being a purely imaginary operator, is an example of the latter. 

### Linear and quadratic response functions

Linear and quadratic response functions are generally noted $\braket{\braket{\hat A; \hat B}}_{\omega_B}$ and $\braket{\braket{\hat A; \hat B, \hat C}}_{\omega_B,\omega_C}$, which describes how the expectation value of $\hat A$ (at frequency $\omega_A = -\omega_B - \omega_C$) responds to a set of perturbation to first and second order in perturbation.

Linear responses result in a rank-2 tensor, ${}^{AB}\mathcal{T}$.
As discussed in [10.1016/S1380-7323(02)80033-4](https://dx.doi.org/10.1016/S1380-7323(02)80033-4), the time-reversal symmetry [and thus the symmetric or antisymmetric nature of $\mathbf t(\omega)$ and $\mathbf u(\omega)$ with respect to a change to the change of sign of $\omega$] can be exploited to compute the linear response:

$$\tag{5}
{}^{AB}\mathcal{T}_{\zeta\sigma} = -\braket{\braket{\hat A_\zeta;\hat B_\sigma}}_\omega
= -2 \sum_{ia}^{CSFs} A_{ia,\zeta}\,{}^B\kappa_{ia,\sigma}(\omega),$$

with:

$${}^B\kappa_{ia,\sigma}(\omega) = \begin{cases}
{}^Bt_{ia,\sigma}(\omega) & \text{ if } {}^A\theta{}^B\theta = 1, \\
{}^Bu_{ia,\sigma}(\omega) & \text{ otherwise},
\end{cases}$$

where ${}^B\mathbf \kappa(\omega)$ is computed from the linear response vectors ${}^B\mathbf x(\omega)$ and ${}^B\mathbf y(\omega)$, obtained using Eq. (3) with $\hat\eta = \hat B$.

??? note "Consequence of the spectral representation"

    If $\hat B$ is hermitian, then, from the spectral representations of the linear response vectors, it can be seen that ${}^B\kappa(\omega)$ is either:

    + hermitian (and real) if both $\hat A$ and $\hat B$ have the same symmetry with respect to time reversal, or 
    + anti-hermitian (and thus imaginary) if not. This further implies that $\lim_{\omega\to 0} {}^B\kappa(\omega)=0$.

A rank-3 tensor, ${}^{ABC}\mathcal{T}$ is the result of a quadratic response function.
Following [10.1021/acs.jctc.7b01008](https://dx.doi.org/10.1021/acs.jctc.7b01008) but neglecting the response of the XC kernel (*i.e.*, $f_{XC}$ and $g_{XC}$, which results in an "unrelaxed" expression), one gets:

$$\tag{6}{}^{ABC}\mathcal{T}_{\zeta\sigma\tau} = \braket{\braket{\hat A_\zeta;\hat B_\sigma,\hat C_\tau}}_{\omega_B,\omega_C} = {}^{ABC}\mathcal{M}_{\zeta\sigma\tau} + {}^{ABC}\mathcal{N}_{\zeta\sigma\tau},$$

where:

$${}^{ABC}\mathcal{M}_{\zeta\sigma\tau} = \sum_{[B,C]}\sum_{ia,ja} \frac{1}{4}\,A_{ij,\zeta}\,{}^{BC}\mathcal{K}_{ia,ja} 
- \frac{1}{2}\,B_{ij,\sigma}\,\left\{\begin{array}{}^At_{ia,\zeta}(\omega_B+\omega_C)\,{}^C\kappa^+_{ja,\tau}(\omega_C)\\ + {}^Au_{ia,\zeta}(\omega_B+\omega_C)\,{}^C\kappa^-_{ja,\tau}(\omega_C)\end{array}\right\},$$

with:

$${}^{BC}\mathcal{K}_{ia,ja} = {}^Bu_{ia,\sigma}(\omega_B)\,{}^Cu_{ja,\tau}(\omega_C) - {}^Bt_{ia,\sigma}(\omega_B)\,{}^Ct_{ja,\tau}(\omega_C) ,$$

and:

$${}^{ABC}\mathcal{N}_{\zeta\sigma\tau} =\sum_{[B,C]}\sum_{ia,ib} \frac{1}{4}\,A_{ab,\zeta}\,{}^{BC}\mathcal{K}_{ia,ib} 
+ \frac{1}{2}\,B_{ab,\sigma}\,\left\{\begin{array}{}^At_{ia,\zeta}(\omega_B+\omega_C)\,{}^C\kappa^+_{ib,\tau}(\omega_C)\\ + {}^Au_{ia,\zeta}(\omega_B+\omega_C)\,{}^C\kappa^-_{ib,\tau}(\omega_C)\end{array}\right\}$$

with:

$${}^{BC}\mathcal{K}_{ia,ib} = {}^Bt_{ia,\sigma}(\omega_B)\,{}^Ct_{ib,\tau}(\omega_C) - {}^Bu_{ia,\sigma}(\omega_B)\,{}^Cu_{ib,\tau}(\omega_C).$$

In both expressions, $\sum_{[B,C]}$ is a sum over the permutation of the pairs of a given operator and its corresponding frequency, $\{(\hat B_\sigma,\omega_B), (\hat C_\tau,\omega_C)\}$, and, using the time-reversal symmetry,

$${}^C\kappa^+_{jb,\tau}(\omega) = \begin{cases} {}^Ct_{jb,\tau}(\omega) & \text{if } \theta_B = 1,\\-{}^Cu_{jb,\tau}(\omega) & \text{otherwise}.\end{cases}$$ 

and

$${}^C\kappa^-_{jb,\tau}(\omega) = \begin{cases} {}^Cu_{jb,\tau}(\omega) & \text{if } \theta_B = 1,\\-{}^Ct_{jb,\tau}(\omega) & \text{otherwise}.\end{cases}$$

In the case of the SHG first hyperpolarizability (so, $\braket{\braket{\hat\mu;\hat\mu,\hat\mu}}_{\omega,\omega}$), these (rather long) expressions are equivalent to the one reported in [10.1002/wcms.1695](https://dx.doi.org/10.1002/wcms.1695).

### Residues

Residue of the response functions provide information on the (excited states of the) unperturbed system.

Assuming an exact wavefunction, the linear response function might be extended in:

$$\braket{\braket{\hat A; \hat B}}_{\omega_B} = \sum_{\ket{m}} \frac{\braket{0|\hat A|m}\braket{m|\hat B|0}}{\omega_B-\omega_m} - \frac{\braket{0|\hat B|m}\braket{m|\hat A|0}}{\omega_B+\omega_m},$$

and a corresponding single residue might be:

$$\lim_{\omega_B\to\omega_m} (\omega_B-\omega_m)\,\braket{\braket{\hat A; \hat B}}_{\omega_B} = \braket{0|\hat A|m}\braket{m|\hat B|0},$$

which provides access to transition matrix elements $\braket{0|\hat A|m}$ between the ground state $\ket{0}$ and an excited state $\ket{m}$. 

In practice, such residues are thus evaluated thanks to amplitude vectors through the spectral representation of linear responses, which gives (assuming that we have a singlet wavefunction):

$$A_{0m,\zeta} = \braket{0|\hat A_\zeta|m} = \sqrt{2}\,\sum_{ia}^{CSF} A_{ia,\zeta}\,\kappa_{ia}^m, \text{ with } \kappa^m_{ia} = \begin{cases}
t^m_{ia} & \text{if } {}^A\theta = 1, \\
u^m_{ia} & \text{otherwise}.
\end{cases}$$

??? note "Details"
    
    This expression was obtained using: 

    $$\begin{aligned}
    \lim_{\omega_B\to\omega_n}\,(\omega_B-\omega_m)\,{}^Bx_{ia,\sigma}(\omega_B) &= B_{ia,\sigma}\,(x^{m}_{ia} + y^{m}_{ia})\,(x_{ia}^m), \\
	\lim_{\omega_B\to\omega_n}\,(\omega_B-\omega_m)\,{}^By_{ia,\sigma}(\omega_B) &= B_{ia,\sigma}\,(x^{m}_{ia} + y^{m}_{ia})\,(y_{ia}^m),
    \end{aligned}$$

    or equivalently,

    $$\begin{aligned}
    \lim_{\omega_B\to\omega_m}\,(\omega_B-\omega_m)\,{}^Bt_{ia,\sigma}(\omega_B) &= B_{ia,\sigma}\,[t^{m}_{ia}]^2, \\
	\lim_{\omega_B\to\omega_m}\,(\omega_B-\omega_m)\,{}^Bu_{ia,\sigma}(\omega_B) &= B_{ia,\sigma}\,t^{m}_{ia}\,u^{m}_{ia},
    \end{aligned}$$
    
    on Eq. (5).

Different ground to excited moments, related to experimentally measurable properties, can be extracted (all given in atomic units):

| Property            | Dipole length                                                   | Dipole velocity                                                          |
|---------------------|-----------------------------------------------------------------|--------------------------------------------------------------------------|
| Oscillator strength | $f^L_{0m} = \frac{2}{3}\,\omega_m\, \|\vec\mu^L_{0m} \|^2$      | $f^L_{0m} = \frac{2}{3\,\omega_m}\, \|\vec\mu^V_{0m} \|^2$               |
| Rotatory strength   | $R^L_{0m} = -\frac{1}{2}\, \Im(\vec\mu^L_{0m}\cdot\vec m_{0m})$ | $R^V_{0m} = -\frac{1}{2\,\omega_m}\,\Re(\vec\mu^V_{0m}\cdot\vec m_{0m})$ | 

See, *e.g.*, [10.1016/j.comptc.2014.02.023](https://doi.org/10.1016/j.comptc.2014.02.023), for a discussion on the differences between the two formalisms.

For an exact wavefunction, quadratic responses are given by:

$$\braket{\braket{\hat A; \hat B, \hat C}} = \sum_{[A,B,C]} \sum_{\ket{m}, \ket{n}} \frac{\braket{0|\hat A|m}\braket{m|\hat B - \delta_{mn}\,\braket{0|\hat B|0}|n}\braket{n|\hat C|0}}{(\omega_A+\omega_m)\,(\omega_C-\omega_n)},$$

where $\omega_A = -\omega_B-\omega_C$ and $\sum_{[A,B,C]}$ is a sum over all permutation of the pairs of a given operator and its corresponding frequency, $\{(\hat A, \omega_A), (\hat B, \omega_B), (\hat C, \omega_C)\}$.

A possible double residue is:

$$\begin{aligned} \lim_{\omega_B\to-\omega_m,\omega_C\to\omega_n} &(\omega_B+\omega_m)\,(\omega_C-\omega_n)\,\braket{\braket{\hat A_\zeta;\hat B_\sigma,\hat C_\tau}}{\omega_B,\omega_C} \\
&= \braket{0|\hat B_\sigma|m}\,\braket{m|\hat A_\zeta - \delta_{mn}\,\braket{0|\hat A_\zeta|0}|n}\,\braket{n|\hat C_\tau|0}.
\end{aligned}$$

Therefore, one gets:

$$\begin{aligned}
&\braket{m|\hat\mu_\zeta - \delta_{mn}\,\braket{0|\hat\mu_\zeta|0}|n}\\
&\hspace{2em}= \frac{1}{4}\sum_{[m,n]}\left\{ \sum_{ia,ib} A_{ab,\zeta}\,[t^m_{ia}\,t^n_{ib} + u^m_{ia}\,u^n_{ib}]  - \sum_{ia,ja} A_{ij,\zeta}\,[t^m_{ia}\,t^n_{ja} + u^m_{ia}\,u^n_{ja}] \right\}.
\end{aligned}$$

which is equal to the fluctuation operator if $m = n$.

??? note "Details"

    This expression was obtained using

    $$\begin{aligned}
    \lim_{\omega_B\to-\omega_m}\,(\omega_B+\omega_m)\,{}^Bt_{ia,\sigma}(\omega_B) &= -B_{ia,\sigma}\,[t^{m}_{ia}]^2, \\
	\lim_{\omega_B\to-\omega_m}\,(\omega_B+\omega_m)\,{}^Bu_{ia,\sigma}(\omega_B) &= B_{ia,\sigma}\,t^{m}_{ia}\,u^{m}_{ia},\\
    \lim_{\omega_C\to\omega_n}\,(\omega_C-\omega_n)\,{}^Ct_{ia,\sigma}(\omega_C) &= C_{ia,\tau}\,[t^{m}_{ia}]^2, \\
	\lim_{\omega_C\to\omega_n}\,(\omega_C-\omega_n)\,{}^Cu_{ia,\sigma}(\omega_C) &= C_{ia,\tau}\,t^{m}_{ia}\,u^{m}_{ia},
    \end{aligned}$$
    
    on Eq. (6), and recognising the expression of $\braket{0|\hat B|m}$ and $\braket{n|\hat C|0}$ given above.

## Sources and references

+ J. Toulouse, [Introduction to the calculation of molecular properties by response theory](https://www.lct.jussieu.fr/pagesperso/toulouse/enseignement/molecular_properties.pdf) (last consultation: January 2023). 
+ E. Fromager, [Linear response time-dependent density functional theory](https://quantique.u-strasbg.fr/lib/exe/fetch.php?media=fr:pageperso:ef:lecture_rctf_tddft_e_fromager.pdf)  (last consultation: January 2023).
+ G. P. Chen, V. K. Voora, M. M. Agee, S. G. Balasubramani, and F. Furche, Random-Phase Approximation Methods. *Annu. Rev. Phys. Chem.* **2017**, 68, 421 ([10.1146/annurev-physchem-040215-112308](https://doi.org/10.1146/annurev-physchem-040215-112308))
+ S. M. Parker, D. Rappoport, and F. Furche, Quadratic Response Properties from TDDFT: Trials and Tribulations. *J. Chem. Theor. Comput.* **2018**, 14, 807 ([10.1021/acs.jctc.7b01008](https://dx.doi.org/10.1021/acs.jctc.7b01008))
+ M. E. Casida, Time-Dependent Density Functional Response Theory for Molecules. In D. E. Chong (ed.), *Recent Advances in Density Functional Methods*. World Scientific, **1995** ([10.1142/9789812830586_0005](https://doi.org/10.1142/9789812830586_0005)).
+ S. Hirata, M. Head-Gordon, Time-dependent density functional theory within the Tamm–Dancoff approximation. *Chem. Phys. Lett.*, **1999**, 314, 291 ([10.1016/S0009-2614(99)01149-5](https://doi.org/10.1016/S0009-2614(99)01149-5))
+ S. Löffelsender, P. Beaujean, M. de Wergifosse, Simplified quantum chemistry methods to evaluate non-linear optical properties of large systems. *WIREs Comput. Mol. Sci.* **2023**, 2023, e1695 ([10.1002/wcms.1695](https://dx.doi.org/10.1002/wcms.1695)) [and references therein]. 
+ S. Grimme, A simplified Tamm-Dancoff density functional approach for the electronic excitation spectra of very large molecules. *J Chem Phys.* **2013**, 138, 244104 ([10.1063/1.4811331](https://doi.org/10.1063/1.4811331)).
+ C. Bannwarth, S. Grimme, A simplified time-dependent density functional theory approach for electronic ultraviolet and circular dichroism spectra of very large molecules. *Comput Theor Chem.* **2014**, 1040, 45 ([10.1016/j.comptc.2014.02.023](https://doi.org/10.1016/j.comptc.2014.02.023)).