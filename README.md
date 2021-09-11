# Time and Frequency Analysis Methods on GW Signals

###### Background Notes by X. Li

### `Part I: Basic Transform Methods`

#### 1.1 Revisit Continuous Fourier Transform and Discrete Fourier Transform

##### Continuous Fourier Transform:

​			
$$
X(f)=\int_{-\infty}^{\infty}x(t) e^{-j2\pi ft}dt
$$

##### Discrete Fourier Transform:

$$
X(k)=\sum_{n=1}^{N-1} x(n)e^{-j2\pi kn/N}
$$

where we let 
$$
\frac{k}{N}\equiv f
$$
and
$$
n \equiv t
$$

#### 1.2 Chirp Z Transform

The **chirp Z-transform** (**CZT**) is a generalization of the [discrete Fourier transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform) (DFT). While the DFT samples the [Z plane](https://en.wikipedia.org/wiki/Z-transform) at uniformly-spaced points along the unit circle, the chirp Z-transform samples along spiral arcs in the Z-plane, corresponding to straight lines in the [S plane](https://en.wikipedia.org/wiki/S_plane).[[1\]](https://en.wikipedia.org/wiki/Chirp_Z-transform#cite_note-Shilling-1)[[2\]](https://en.wikipedia.org/wiki/Chirp_Z-transform#cite_note-2) The DFT, real DFT, and zoom DFT can be calculated as special cases of the CZT. ( - Wikipedia)

To obtain DCZT from DFT, we introduce
$$
W^{kn}=e^{-j2\pi kn/N}
$$
where		
$$
W^{kn}=W^{k^{2}/2}+ W^{n^{2}/2} + W^{-(k-n)^{2}/2}
$$
thus, the DFT becomes
$$
X(k)=W^{k^{2}/2}\sum_{n=1}^{N-1} x(n)W^{n^{2}/2}W^{-(k-n)^{2}/2}
$$
where a complex chirp is introduced via 
$$
W^{n^{2}/2}
$$


##### Discrete Chirp Z Transform (*to be verified):

$$
X(k)=e^{-j\pi k^{2}/N}\sum_{n=1}^{N-1} x(n)e^{-j\pi n^{2}/N} e^{j\pi(k-n)^{2}/N}
$$

###### Compute with FFT & IFFT:

$$
X(k)=W^{k^{2}/2}IFFT[FFT(x(n)W^{n^{2}/2})FFT(W^{-(k-n)^{2}/2})]
$$

Some notes:

When computing using FFT, the signal in the block is treated as periodic**(to be verified), however, the input signal itself does not need to be.

Fresnel integral approximation?

#### 1.3 Gabor Transform

The **Gabor transform** is a special case of the [short-time Fourier transform](https://en.wikipedia.org/wiki/Short-time_Fourier_transform).  It is used to determine the [sinusoidal](https://en.wikipedia.org/wiki/Sine_wave) [frequency](https://en.wikipedia.org/wiki/Frequency) and [phase](https://en.wikipedia.org/wiki/Phase_(waves)) content of local sections of a signal as it changes over time. The function to be transformed is first multiplied by a [Gaussian function](https://en.wikipedia.org/wiki/Gaussian_function), which can be regarded as a [window function](https://en.wikipedia.org/wiki/Window_function), and the resulting function is then transformed with a Fourier transform to derive the [time-frequency analysis](https://en.wikipedia.org/wiki/Time-frequency_analysis). ( - Wikipedia)

##### Continuous Gabor Transform:

$$
G_x(\tau, \omega)=\int_{-\infty}^{\infty}x(t)e^{-\pi (t-\tau)^{2}} e^{-j\omega t}dt
$$

where the window function is a Gaussian in the form of
$$
e^{-\pi \alpha(t-\tau)^{2}}
$$
the width of which can be controlled via a scaling parameter α.

Thus the scaled Gabor is 
$$
G_x(\tau, \omega)=\sqrt[4]{\alpha}\int_{-\infty}^{\infty}x(t)e^{-\pi \alpha(t-\tau)^{2}} e^{-j\omega t}dt
$$


##### Inverse Gabor Transform:

$$
x(t)=\int_{-\infty}^{\infty}G_x(\tau, \omega)e^{j\omega t}e^{\pi(t-\tau)^{2}}d\omega/(2\pi)
$$

##### Gabor Expansion:

$$
y(t)=\sum_{m=-\infty}^{\infty}\sum_{n=-\infty}^{\infty}C_{nm}g_{nm}(t)
$$



##### Discrete Gabor Transform (*to be verified):

$$
y(k)=\sum_{m=1}^{M-1}\sum_{n=1}^{N-1}....
$$

#### 1.4 S Transform

The *S* transform is a generalization of the [short-time Fourier transform](https://en.wikipedia.org/wiki/Short-time_Fourier_transform) (STFT), extending the [continuous wavelet transform](https://en.wikipedia.org/wiki/Continuous_wavelet_transform) and overcoming some of its disadvantages. Modulation sinusoids are fixed with respect to the time axis; this localizes the scalable Gaussian window dilations and translations in *S* transform. The *S* transform doesn't have a cross-term problem and yields a better signal clarity than [Gabor transform](https://en.wikipedia.org/wiki/Gabor_transform). ( - Wikipedia)

##### Continuous S Transform:

$$
S_x(\tau, \omega)=\int_{-\infty}^{\infty}x(t)|\omega|e^{-\pi (t-\tau)^{2}\omega^{2}}e^{-j2\pi \omega \tau}d\tau
$$

where
$$
\omega = \frac{|f|}{\sqrt{2\pi}}
$$


##### Inverse S Transform:

$$
x(t)=\int_{-\infty}^{\infty}[\int_{-\infty}^{\infty}S_x(\tau, \omega)d\tau]e^{j2\pi \omega t}d\omega
$$

##### Spectral Form:

$$
S_x(\tau, \omega)=\int_{-\infty}^{\infty}X(\omega+\alpha)e^{-\pi \alpha^{2}/\omega^{2}}e^{j2\pi \alpha \tau}d\alpha
$$

where X(t) represents the Fourier Transform. (Convolution is used to get this spectral form)

Let
$$
\tau = n\triangle_T, \omega = p\triangle_F,\alpha=p\triangle_T
$$
We obtain the 

##### Discrete S Transform:

$$
S_x(n\triangle_T, p\triangle_F)=\sum_{p=0}^{N-1}X[(p+m)\triangle_F]e^{-\pi p^{2}/m^{2}}e^{j2\pi pm/N}
$$

##### Implementation:

```pseudocode
  Step1.Compute 
  
    
    {\displaystyle X[p\Delta _{F}]\,}
  
X[p\Delta _{{F}}]\, 

  loop{
     Step2.Compute
  
    
    {\displaystyle e^{-\pi {\frac {p^{2}}{m^{2}}}}}
  
e^{{-\pi {\frac  {p^{2}}{m^{2}}}}}for 
  
    
    {\displaystyle f=m\Delta _{F}}
  
f=m\Delta _{{F}}

     Step3.Move 
  
    
    {\displaystyle X[p\Delta _{F}]}
  
X[p\Delta _{{F}}]to
  
    
    {\displaystyle X[(p+m)\Delta _{F}]}
  
X[(p+m)\Delta _{{F}}]

     Step4.Multiply Step2 and Step3 
  
    
    {\displaystyle B[m,p]=X[(p+m)\Delta _{F}]\cdot e^{-\pi {\frac {p^{2}}{m^{2}}}}}
  
B[m,p]=X[(p+m)\Delta _{{F}}]\cdot e^{{-\pi {\frac  {p^{2}}{m^{2}}}}}

     Step5.IDFT(
  
    
    {\displaystyle B[m,p]}
  
B[m,p]).
  Repeat.}
```

#### 1.5 Constant Q Transform

**CQT** transforms a data series to the frequency domain. It is related to the [Fourier transform](https://en.wikipedia.org/wiki/Fourier_transform)[[1\]](https://en.wikipedia.org/wiki/Constant-Q_transform#cite_note-b91-1) and very closely related to the complex [Morlet wavelet](https://en.wikipedia.org/wiki/Morlet_wavelet) transform.[[2\]](https://en.wikipedia.org/wiki/Constant-Q_transform#cite_note-2)

The transform can be thought of as a series of filters *f**k*, logarithmically spaced in frequency, with the *k*-th filter having a [spectral width](https://en.wikipedia.org/wiki/Spectral_width) *δf**k* equal to a multiple of the previous filter's width:

![{\displaystyle \delta f_{k}=2^{1/n}\cdot \delta f_{k-1}=\left(2^{1/n}\right)^{k}\cdot \delta f_{\text{min}},}](https://wikimedia.org/api/rest_v1/media/math/render/svg/5d4620127b4c5a85de8a9849fe3a8e23274db0d9)

where *δf**k* is the bandwidth of the *k*-th filter, *f*min is the central frequency of the lowest filter, and *n* is the number of filters per [octave](https://en.wikipedia.org/wiki/Octave_(electronics)). ( - Wikipedia)

**Shifted (m) STFT:**
$$
X(k,m)=\sum_{p=0}^{N-1}W[n-m]x[n]e^{j2\pi pm/N}
$$
given data with sampling frequency fs = 1/T, we have filter width *δf_k* and quality factor Q such that
$$
Q = \frac{f_k}{\delta f_k}
$$
Our window length of the k_th bin is therefore
$$
N[k]=\frac{f_s}{\delta f_k}=\frac{f_s}{f_k}Q
$$
making the relative power for each bin decrease at higher frequencies. To compensate this effect, we normalize the bin length by itself N[k].

Using the Hamming Window as an example, the normalization is
$$
W[k,n]=\alpha-(1-\alpha)cos\frac{2\pi n}{N[k]-1}
$$
with
$$
\alpha=25/46, 1\le n \le N[k]-1
$$
Finally, the **Discrete Constant Q Transform:**
$$
X[k]=\frac{1}{N[k]}\sum_{n=0}^{N[k]-1}W[k,n]x[n]e^{\frac{-j2\pi Qn}{N[k]}}
$$
Notes:

CQT is well suited to music data as the out put of the transform is effectively amplitude/phase against log frequency which is useful when the frequency spans several octaves. However, in the case of gravitational waves

#### 1.6 Implementation on Chirp Signal (Windowed FFT - *Tukey window*)

In [1]:

```python
import numpy as np
import matplotlib.pyplot as plt
```

In [2]:

```python
dt = 0.001
t = np.arange(0,3,dt)
f0 = 50
f1 = 250
t1 = 2
x = np.cos(2*np.pi*t*(f0 + (f1 - f0)*np.power(t, 2)/(3*t1**2)))

fs = 1/dt

#Chirp signal
plt.plot(t, x)
```

![](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/chirp_signal.png)

In [3]:

```python
from scipy import signal
```

In [4]:

```python
freqs, times, spectrogram = signal.spectrogram(x)

#Scipy PSD Power Spectral Density
freqs, psd = signal.welch(x)

plt.figure(figsize=(5, 4))
plt.semilogx(freqs, psd)
plt.title('PSD: power spectral density')
plt.xlabel('Frequency')
plt.ylabel('Power')
plt.tight_layout()

#Scipy Spectrogram
freqs, times, spectrogram = signal.spectrogram(x)

plt.figure(figsize=(5, 4))
plt.imshow(spectrogram, aspect='auto', cmap='hot_r', origin='lower')
plt.title('Spectrogram')
plt.ylabel('Frequency band')
plt.xlabel('Time window')
plt.tight_layout()
```



![img](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/chirp_psd.png)



![img](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/chirp_spectrogram_scipy.png)



![chirp_spectrogram_matplot](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/chirp_spectrogram_matplot.png)

#### 1.7 On the implementation of Time-Frequency Analysis for CNN

![](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/compare_Gabor_CQT.png)

Constant Q Transform's high resolution at low frequencies is desirable in many cases for accurate analysis of the signal in both the time and frequency domain. However, when applying CQT for the purpose of identifying specific features on a spectrogram, high resolution in certain directions also means less obvious feature representation in certain cases (**this assumption is highly dependent on the type of data you have**). The Q range also affects the result significantly.

Some examples on GW chirp signal:

![GW_CQT_verylargeQ](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/GW_CQT_verylargeQ.png)

![GW_CQT_largeQ](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/GW_CQT_largeQ.png)

![GW_CQT_smallQ](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/GW_CQT_smallQ.png)

Thus, the time and frequency sampling width/length need to be adjusted accordingly for performance and resolution.

A **spectrogram** can also be **inverse transformed** to time-domain. Taking the generated chirp signal as an example:

![original_chirp](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/original_chirp.png)

is the original chirp signal.

![chirp_spectrogram_cqnsgt](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/chirp_spectrogram_cqnsgt.png)

is the spectrogram using invertible Constant Q - Nonstationary Gabor Transform (Matlab CQT).

![inverse_transform_chirp](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/inverse_transform_chirp.png)

is the inverse transformed chirp signal.

![original_inverse_difference](https://raw.githubusercontent.com/xli2522/GW-SignalGen/main/img/original_inverse_difference.png)

is the magnitude difference between the original chirp signal and the inverse transformed chirp signal.

#### 1.8 Resources on CQT

Invertible CQT: https://www.univie.ac.at/nonstatgab/cqt/index.php (Useful for the Spectrogram Denoising and Time-series parameter estimation)

Gabor Wavelets (bob package): https://pythonhosted.org/bob.ip.gabor/guide.html#gabor-wavelets

Window Functions (Wikipedia): https://en.wikipedia.org/wiki/Window_function

On the **discrete Gabor transform** and the **discrete Zak transform**: http://www.martinbastiaans.org/pdfs/sigpro.pdf

### `Part II: Transform Methods Targeting Chirp Signals`

#### 2.1 Fractional Fourier Transform (FrFT)

The Fractional Fourier Transform is a family of linear transforms generalizing the Fourier Transform. It can be thought of as the Fourier Transform to the n-th power, where n need not be an integer. Thus it can transform a function to any intermediate domain between time and frequency. Its applications range from *filter design* and signal analysis to phase retrieval and pattern recognition.

A completely different meaning for fractional Fourier Transform was introduced by Bailey and Swartztrauber as essentially another name for a *z transform*, and in particular for the case that corresponds to a discrete Fourier transform shifted by a fractional amount in frequency space (multiplying the input by a linear chirp) and evaluating at a fractional set of frequency points. ( - Wikipedia)

The FrFT provides a continuous representation of a signal from the time to the frequency domain at intermediate domains by means of the fractional order of the transform that changes from -pi/2 to pi/2. (- The DLCT and its Applications 2013)

**Continuous Fractional Fourier Transform**
$$
X_\alpha(u)=\int_{-\infty}^{\infty}x(t)K_\alpha(t,u)dt
$$
where alpha is the fractional order and K_alpha(t,u) is the kernel of the transformation which is defined as 
$$
X_\alpha(u)=\sqrt(\frac{1-j*cot(\alpha)}{2*\pi})exp(j*\frac{t^2+u^2}{2}cot(\alpha)-j*ut*csc(\alpha))
$$
When alpha = 0, the FrFT of the signal x(t) is the signal itself, while if alpha = -+pi/2, the FrFT becomes the Fourier Transform. 

**Inverse Continuous Fractional Fourier Transform**
$$
x(t)=\int_{\infty}^{\infty}X_\alpha(u)K^*_\alpha(t,u)du
$$
where * stands for the complex conjugate.

**Discrete Fractional Fourier Transform**
$$
X_\alpha(\rho)=\sum_{n=0}^{N-1}K^ˆ_\alpha(n, \rho)x(n)
$$
where the Kernel of the transformation has the following spectral expression
$$
X^ˆ_\alpha(n, \rho)=\sum_{k\in M}\nu^ˆ_\alpha(\rho)exp(-j*\alpha k)\nu_k(n)
$$
where nu_k(n) is the kth discrete Hermite - Gaussian function and M = {0, ..., N-2, N-N mod2}.

**Inverse Discrete Fractional Fourier Transform** 
$$
x(n)=\sum_{\rho=0}^{N - 1}K{^ˆ}_\alpha^*(n, \rho)X_\alpha(\rho)
$$
For a discrete signal x(n), we can define the connection between the fractional order (\alpha) and the chirp rate (\gamma) as 
$$
\alpha = - tan ^{-1}(\frac{1}{2\gamma})
$$
​																									More info see The DLCT and its Applications 2013

#### 2.2 Discrete Chirp-Fourier Transform 

**Discrete Chirp-Fourier Transform** 
$$
X_c(k,r) = \frac{1}{\sqrt{N}}\sum_{n=0}^{N-1}exp(-j*\frac{2\pi}{N}(rn^2+kn)), 0\le r,k\le N-1
$$
where k represents the frequencies and r is an arbitrarily fixed integer that represents the chirp rates. The DCFT is the same as the DFT when r = 0. (page 5)

**Inverse Discrete Chirp-Fourier Transform** 
$$
x(n) = exp(j*\frac{2\pi}{N}rn^2)\frac{1}{\sqrt{N}}\sum_{k=0}^{N-1}X_c(k,r)exp(j\frac{2\pi}{N}kn), 0\le n\le N-1
$$
The DCFT approximates the chirp rate by integer numbers r. Therefore, when using the DCFT to detect a chirp signal, the discrete chirp rate r0 of the signal should be an integer to guarantee that the parameter can be matched and that the peak will not be lost. This restriction affects the practical applications of the DCFT. (page 5)

​																									More info see The DCFT and its Applications 2013

#### **2.3 Linear Chirp Transform**

The DLCT uses discrete complex linear chirp bases. It is not a time-frequency but rather a *frequency chirp-rate* transformation, implementable using Fast Fourier Transform. The Discrete Fourier Transform is a special case of the DLCT which has the properties of modulation and duality in time and frequency. (page 6)

**Continuous Linear Chirp Transform**

(Other conditions see DCFT and its Applications 2013)
$$
X(\Omega,\gamma)=\int_{-\infty}^{\infty}x(t)exp(-j(\Omega t+\gamma^2))dt
$$
**Inverse Continuous Linear Chirp Transform**
$$
x(t)=\int_{-\infty}^{\infty}X(\Omega, \gamma) exp(j(\Omega t+\gamma^2))d\Omega
$$
The CLCT can provide us with the same bandwidth of sinusoid signals if we intent to filter linear chirps in the frequency chirp-rate space. Thus, the CLCT can eliminate the effect of the chirp rate on the channel bandwidth of chirp communication systems (page 11, more in chapter 6) if we filter the signal at the corresponding chirp rate. 

**Discrete Linear Chirp Transform**
$$
X(k,m)=\sum_{n=0}^{N-1}x(n)exp(-j \frac{2\pi}{N}(Cmn^2+kn)), 0\le k\le N-1, -L/2\le m\le(L/2)-1
$$
**Inverse Discrete Linear Chirp Transform**
$$
x(n)=\sum_{m=-L/2}^{L/2-1}\sum_{k=0}^{N-1}\frac{X(k,m)}{LN}exp(j\frac{2\pi}{N}(Cmn^2+kn)), 0\le k\le N-1, C=\frac{2\Lambda}{L}
$$
Note: One could think of the DLCT as a generalization of the DFT
$$
X(k,m)=\frac{1}{N}X(k)\circledast DFT[exp(-j\frac{2\pi}{N}Cm)]
$$
If m = 0, then X(k,0) is the DFT of x(n) or the representation using chirp base with zero chirp rates. (page 15)



**Properties of the Discrete Linear Chirp Transform (page 17):**

Properties of the DLCT are similar to those of the DFT. 

- Modulation property...
- Duality property...

***Implementation with the FFT*** (page 19)

The *DLCT* can be implemented using the FFT. Rewriting X(k,m) as
$$
X(k,m)=\sum_{n=0}^{N-1}[x(n)exp(j2\pi Cmn^2/N)]exp(j2\pi kn/N)
$$
and let
$$
h(n,m)=x(n)exp(j2\pi Cmn^2/N)
$$
then for each -L/2 ⩽ m0 ⩽ L/2 - 1 the X(k, m0) is the DFT of h(n,m0) which can be obtained by the FFT.

The IDLCT can be implemented using the inverse FFT. Rewriting the expression for x(n) as 
$$
x(n)=\frac{1}{L}\sum_{n=-L/2}^{L/2-1}[\sum_{n=0}^{N-1}\frac{X(k,m)}{N}exp(j2\pi kn/N)]exp(j2\pi Cmn^2/N)
$$
and let 
$$
g(n,m)=\sum_{n=0}^{N-1}\frac{X(k,m)}{N}exp(j2\pi kn/N)
$$
where g(n,m) is the inverse DFT for each  -L/2 ⩽ m0 ⩽ L/2 - 1. Then
$$
x(n)=\frac{1}{L}\sum_{n=-L/2}^{L/2-1}[g(n,m)]exp(j2\pi Cmn^2/N)
$$
If a vector X = [x(0)...x(N-1)]^T then
$$
x = \frac{1}{L}diag[G E]
$$
or the diagonal of the product of an N×L matrix.

G = {g(n, m)} with an L × N matrix E = {e(m, n) = exp(j2πCmn2/N )}

#### **2.4 Discrete Cosine Chirp Transform**

Page 22
$$

$$



### `Part III: Types of Chirp Signals and GW Signal Compression`

- Linear Chirp

  - A linear chirp is a function whose frequency changes linearly with time. For example, a linear chirp has an instantaneous frequency 
    $$
    \Omega_0+2\gamma_0t
    $$
    at time t. This is one of the reasons we need to use linear chirp bases instead of the classical Fourier bases because they are more suitable for representing the frequency changes of non-stationary signals. (page 10)

    

- Quadratic Chirp

- Logarithmic Chirp



