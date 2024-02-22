# IntFFT
Implementation of 6 FFTs - Radix 2, Radix 4, Split Radix - each one
decimated in time and frequency and their according integer
approximations implemented according to this paper: 
https://www.researchgate.net/publication/3908291_Integer_fast_Fourier_transform_IntFFT

Essentialy, the complex product $x_k \cdot e^{i \theta}$ in the FFT
is replaced with a <i>lifting scheme</i> (function `lift` in
`intfft.h`). It's using a <i>quantizator</i> (the `Q` template)
which converts floats to integers. Functions like `round`, `floor`
or `ceil` could be used as quantizators.
