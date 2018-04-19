# Kuramoto-Sivashinsky Module

### Basic Equations

This module solves the Kuramoto-Sivashinsky equation for $u(x,t)$:

$$\partial_t u + \partial_x^4 u + \partial_x^2 u + u\partial_x u = 0\ .$$

### Implementation

The equation is time-stepped forward in Fourier space:

$$\partial_t \widehat{u} + k_x^4 \widehat{u} - k_x^2 \widehat{u} + \widehat{ u\partial_x u } = 0\ .$$

Thus:

$$\mathcal{L} = k_x^2 - k_x^4\ ,$$
$$\mathcal{N}(\widehat{u}) = - \mathrm{FFT}(u \partial_x u)\ .$$

The function `calcN!` implements dealiasing to avoid energy piling up at the grid-scale.
