# Theory Core

## 2. The Local Volatility Pricing Equation

Under the risk-neutral measure, the asset price $S_t$ satisfies

$$
dS_t = r S_t \, dt + \sigma(t,S_t) S_t \, dW_t.
$$

Let $x = \ln S$ and $\tau = T - t$. Define the local variance

$$
D(x,\tau) := \frac12 \sigma^2(t,S).
$$

Suppressing explicit $\tau$-dependence, the value $V(x,\tau)$ of a European claim satisfies

$$
V_\tau = D(x) V_{xx} - (D(x) - r) V_x - r V, \qquad V(x,0) = g(x),
$$

where $g(x)$ is the payoff, for example $(e^x-K)^+$ for a call.

Standing assumptions:

- $D(x)\in C^2(\mathbb{R})$ and $D(x)>0$ on the computational domain.
- Sufficient growth at infinity, or finite-domain boundary conditions, to guarantee well-posedness.
- Markovian local volatility.
- Constant interest rate $r\geq 0$.
- No jumps.

## 3. Gauge Transformation

**Theorem 3.1 (Gauge Transformation).** Let

$$
\Phi(x)=\frac{x}{2}-\int^x \frac{r}{2D(s)}\,ds.
$$

Set $V(x,\tau)=e^{\Phi(x)}\psi(x,\tau)$. Then $\psi$ satisfies

$$
\psi_\tau = D(x)\psi_{xx}+V_0(x)\psi,
$$

where

$$
V_0(x)=-\frac{D(x)}{4}-\frac{r}{2}-\frac{r^2}{4D(x)}+\frac{rD'(x)}{2D(x)}.
$$

**Proof.** Substitute $V=e^\Phi\psi$ into the pricing equation. Using

$$
V_x=e^\Phi(\psi_x+\Phi'\psi),
$$

and

$$
V_{xx}=e^\Phi\left[\psi_{xx}+2\Phi'\psi_x+(\Phi''+(\Phi')^2)\psi\right],
$$

division by $e^\Phi$ gives

$$
\psi_\tau
=D\psi_{xx}+\left[2D\Phi'-(D-r)\right]\psi_x
+\left[D(\Phi''+(\Phi')^2)-(D-r)\Phi'-r\right]\psi.
$$

The coefficient of $\psi_x$ vanishes if and only if

$$
2D\Phi'=D-r,
$$

or

$$
\Phi'=\frac12-\frac{r}{2D(x)}.
$$

Differentiating gives

$$
\Phi''=\frac{rD'}{2D^2}.
$$

Substitution and simplification yield the stated expression for $V_0$.

## 4. Rescaling to Self-Adjoint Form

**Theorem 4.1 (Main Result).** Let $\psi=\sqrt{D}\,\varphi$. Then $\varphi$ satisfies

$$
\varphi_\tau=-H\varphi,
$$

where

$$
H\varphi=-\partial_x\left(D(x)\partial_x\varphi\right)+U_{\rm eff}(x)\varphi
$$

and

$$
U_{\rm eff}(x)
=\frac{D}{4}+\frac{r}{2}+\frac{r^2}{4D}
-\frac{rD'}{2D}+\frac{(D')^2}{4D}-\frac{D''}{2}.
$$

**Proof.** Set $g(x)=\sqrt{D(x)}$, so $\psi=g\varphi$. Then

$$
\psi_x=g'\varphi+g\varphi_x,
$$

and

$$
\psi_{xx}=g''\varphi+2g'\varphi_x+g\varphi_{xx}.
$$

Substituting into the gauged equation and dividing by $g>0$ gives

$$
\varphi_\tau
=D\varphi_{xx}
+2D\frac{g'}{g}\varphi_x
+\left(D\frac{g''}{g}+V_0\right)\varphi.
$$

For $g=\sqrt{D}$,

$$
\frac{g'}{g}=\frac{D'}{2D},
$$

so the first-order coefficient is $D'$, matching

$$
\partial_x(D\partial_x\varphi)=D\varphi_{xx}+D'\varphi_x.
$$

The zero-order contribution is

$$
U_{\rm eff}=-V_0-D\frac{g''}{g}.
$$

Since

$$
g'=\frac{D'}{2\sqrt D},
$$

and

$$
g''=\frac{D''}{2\sqrt D}-\frac{(D')^2}{4D^{3/2}},
$$

we have

$$
\frac{g''}{g}=\frac{D''}{2D}-\frac{(D')^2}{4D^2}.
$$

Therefore

$$
D\frac{g''}{g}=\frac{D''}{2}-\frac{(D')^2}{4D}.
$$

Combining this with $V_0$ gives the stated effective potential.

The decomposition is

$$
U_{\rm eff}(x)
=
\underbrace{\frac{(D')^2}{4D}-\frac{D''}{2}}_{U_{\rm geom}(x)}
+
\underbrace{\frac{D}{4}+\frac{r}{2}+\frac{r^2}{4D}-\frac{rD'}{2D}}_{U_{\rm fin}(x)}.
$$

Here $U_{\rm geom}$ encodes volatility-gradient and geometric effects, while $U_{\rm fin}$ encodes discounting and no-arbitrage drift effects.

## 5. Riemannian Geometric Interpretation

The diffusion coefficient induces the metric

$$
ds^2=\frac{dx^2}{D(x)}.
$$

Equivalently,

$$
g_{xx}=D^{-1}, \qquad g^{xx}=D.
$$

High-volatility regions, where $D$ is large, shorten the metric distance and permit rapid probability transport. Low-volatility regions stretch the metric and behave like effective barriers.

## 6. Quadratic Normal Volatility and Regimes

For the quadratic normal volatility specification,

$$
\sigma(x)=\alpha x^2+\beta x+\gamma,
$$

and

$$
D(x)=\frac12\sigma(x)^2.
$$

The discriminant

$$
\Delta=\beta^2-4\alpha\gamma
$$

classifies the quadratic profile:

- $\Delta<0$: elliptic regime; no real roots.
- $\Delta=0$: parabolic transition; a double root.
- $\Delta>0$: hyperbolic regime; two real roots.

For numerical pricing on a finite domain, the hyperbolic regime is acceptable only when the roots remain outside the pricing window and $\sigma(x)>0$ on the grid.

The canonical hyperbolic parametrization is

$$
\sigma(x)=\sigma_0\cosh(a x),
$$

so

$$
D(x)=D_0\cosh^2(a x).
$$

For this family,

$$
D'=2aD_0\cosh(ax)\sinh(ax),
$$

and

$$
D''=2a^2D_0\left(\cosh^2(ax)+\sinh^2(ax)\right).
$$

The geometric potential simplifies to

$$
U_{\rm geom}
=a^2D_0\sinh^2(ax)
-a^2D_0\left(\cosh^2(ax)+\sinh^2(ax)\right)
=-a^2D_0\cosh^2(ax).
$$

After a coordinate change such as $z=\tanh(ax)$, the operator is related to a hyperbolic Pöschl--Teller form. A canonical shape is

$$
U_{\rm PT}(z)\propto-\frac{\lambda(\lambda-1)a^2}{1-z^2}.
$$

The idealized discrete spectrum is

$$
E_n=a^2(\lambda-1-n)^2,
\qquad
n=0,1,\dots,\lfloor \lambda-1\rfloor.
$$

## 7. Spectral Pricing Formula

By the spectral theorem, the transformed solution admits

$$
\varphi(x,\tau)
=
\sum_n c_n e^{-E_n\tau}\phi_n(x)
+\int_0^\infty c(k)e^{-E(k)\tau}\phi_k(x)\,dk.
$$

The coefficients are

$$
c_n=\langle \phi_n,\varphi(\cdot,0)\rangle_{L^2},
$$

with

$$
\varphi(x,0)=\frac{e^{-\Phi(x)}V(x,0)}{\sqrt{D(x)}}.
$$

The option price is recovered by

$$
V(x,\tau)=e^{\Phi(x)}\sqrt{D(x)}\,\varphi(x,\tau).
$$

For a truncated discrete expansion, one expects an error controlled by the neglected high-energy modes:

$$
\|\mathrm{remainder}\|\leq C e^{-E_N\tau}
$$

for fixed $\tau>0$, when the spectral assumptions and payoff regularity conditions are satisfied.
