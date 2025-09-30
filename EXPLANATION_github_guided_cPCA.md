# Complex PCA and Guided cPCA: Outputs and Interpretation 

## 1) Data and normalization
When `--normalize zscore` is used:
$$
z_p(t) = \frac{x_p(t) - \overline{x}_p}{\mathrm{sd}(x_p)}.
$$

## 2) cPCA via SVD
$$
X = U \Sigma V^{H}.
$$

Component time courses: $s_k(t) = (U\Sigma)_{t,k}$.  
Spatial loadings: $\ell_{k,p}$ from columns of $V$.  
Explained variance: $\lambda_k = \sigma_k^2/(T-1)$.

## 3) Reconstruction and phase binning
$$
\widehat{x}_k(t,p) = \Re\{ s_k(t)\,\ell_{k,p} \}.
$$

With $\phi_k(t) = \arg s_k(t)$ and $N$ bins over $[-\pi,\pi)$:
$$
M[b,p] = \mathrm{mean}_{\;t:\,\phi_k(t)\in \text{bin } b}\; \widehat{x}_k(t,p).
$$

## 4) Mapping bins to seconds
If dominant frequency is $f$ Hz and $N$ bins are used:
$$
\text{seconds per bin} = \frac{1}{fN}.
$$

## 5) Guided cPCA
$$
C = \frac{X_c^{H} X_c}{T-1}, \qquad
C_{\mathrm{guided}} = C + T\,\mathrm{diag}(\lambda)\,T^{H}.
$$

Eigendecomposition of $C_{\mathrm{guided}}$ yields guided components. Weights $\lambda$ control the guidance strength.
