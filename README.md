[> Intro
--------
This provides a simple fixed-point FFT module.

The implementation is based on VFFT - [Ryuji Naitou's FPGA向きのFFTアルゴリズム(A FFT algorithm suitable for FPGA)](http://nahitafu.cocolog-nifty.com/nahitafu/2016/01/fpgafft-78b8.html) and written by amaranth HDL. The test portion uses the [amlib](https://github.com/amaranth-community-unofficial/amlib) library.

Since this module is for the display my radio system, its interface is not generic. It could be easily modified to adapt it to other purposes, I guess. The input/output is given with complex number which is a pair of two fixed-point numbers representing -1 to 1. Inputs should be synced with a strobe signal and outputs will put with the output strobe signal. The inputs are cooked by 4-term blackman-harris window and the equation

$$\frac{1}{N} \sum_{\rm k=0}^{N-1} x_{\rm k} e^{-i \frac{2\pi {\rm k}}{N}}$$

is used for DFT.

[> Updates
----------

[> Features
-----------
**TODO**

[> Getting started
------------------
**TODO**

[> Tests
--------
**TODO**

run-tests.sh generates a simple test results test_FixedPointFFTTest_test_fft.vcd which can be seen by gtkwave
```
gtkwave --rcvar 'enable_vcd_autosave yes' --rcvar 'do_initial_zoom_fit yes' test_FixedPointFFTTest_test_fft.vcd
```
for example.

[> Links
-------------

[1] Ryuji Naitou's article
* [FPGA向きのFFTアルゴリズム](http://nahitafu.cocolog-nifty.com/nahitafu/2016/01/fpgafft-78b8.html)
