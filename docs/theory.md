title: Theory

## (Static) response function

Molecular properties are defined as a response of a (molecular) system to an external perturbation (change of geometry, application of an external electric field, etc).

A general expression for the electronic energy under perturbation is:

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

$$	\forall\kappa_i: \left.\frac{\partial E(\{\kappa_i \}, \lambda)}{\partial\kappa_i}\right|_{\lambda=0} = 0.$$

The second-order derivative is obtained from differentiating the Hellman-Feynman theorem, yielding (under the assumption that the perturbation is linear):

$$\frac{d^2 E}{d\lambda^2} =
\sum_i \left.\frac{\partial^2 E(\{\kappa_i\},\lambda)}{\partial\lambda\partial\kappa_i}\right|_{\lambda=0}\,\frac{\partial\kappa_i}{\partial\lambda} 
\Leftrightarrow \braket{\braket{\hat V; \hat V}} = (\mathbf\kappa^\lambda)^T\,\eta^\lambda$$

where $\braket{\braket{\hat V; \hat V}}$ is the linear response function, $\mathbf\eta^\lambda$ is the perturbed electronic gradient, and $\mathbf\kappa^\lambda$, the first-order response vector to $\lambda$, containing all $\partial\kappa_i/\partial\lambda$. 
The response vector can be optained thanks to the following set of equations (derivarive of the stationary condition):

$$	\sum_j \left.\frac{\partial^2 E(\{\kappa_i\},\lambda)}{\partial\kappa_i\partial\kappa_j}\right|_{\lambda=0}\,\frac{\partial\kappa_j}{\partial\lambda} = -\left.\frac{\partial^2 E(\{\kappa_i\},\lambda)}{\partial\lambda\partial\kappa_i}\right|_{\lambda=0}  \Leftrightarrow \mathbf{E}^\lambda\,\mathbf\kappa^\lambda = -\mathbf\eta^\lambda$$
 
where $\mathbf E^\lambda$ is the electronic Hessian. Both $\mathbf E^\lambda$ and $\mathbf\eta^\lambda$ actually result in new kinds of (bielectronic) integrals to evaluate (*vide supra*), while $\kappa^\lambda$ has to be computed by solving this equation.
Thanks to the [2n+1 rule](https://doi.org/10.1007/s10910-008-9497-x), the response vector can be used to access the quadratic response function (*i.e.*, the third derivative of the energy w.r.t. $\lambda$) as well.

## Sources

+ J. Toulouse, [*Introduction to the calculation of molecular properties by response theory*](https://www.lct.jussieu.fr/pagesperso/toulouse/enseignement/molecular_properties.pdf). 