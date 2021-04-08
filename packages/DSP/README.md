*DSP* is a [Datagrok package](https://datagrok.ai/help/develop/develop#packages) for digital signal processing. The package provides a set of core signal processing functions wrapped around fast [WebAssembly](https://webassembly.org/)-based implementations. These functions run directly in the browser on the data provided.

The following functions are currently implemented:

* SMA Filter
* Exponential Filter
* Kalman filter
* Min-Max Transform
* Z-score Transform
* Box-Cox Transform
* Trend Detection
* Detrending
* Fourier Filter
* Spectral Density
* Subsampling
* Averaging Downsampling

We would appreciate your contribution to this package with more filters!

The package is also used in the [BioSignals](https://github.com/datagrok-ai/public/tree/master/packages/BioSignals) package.