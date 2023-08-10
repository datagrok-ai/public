#name: Spectrogram
#description: Spectrogram plot
#reference: https://en.wikipedia.org/wiki/Spectrogram
#language: r
#tags: demo, viewers
#sample: sensors/ecg.csv
#input: dataframe data [Input data table]
#input: column signal {type:numerical; allowNulls:false} [Column with signal]
#input: double sampleRate = 256.0 [Signal sample rate, in Hz]
#input: int windowLength = 1024 [FFT window length, in samples]
#input: double timeStep = 0.1 [Analysis time step, in seconds]
#input: bool removeDc = TRUE [Flag to force DC component removal]
#output: graphics {name:spectrogram} [Spectrogram plot]

require(signal)
require(graphics)
require(colorRamps)

# Select signal column
signal = data[[signal]]

# Remove DC
if (removeDc == TRUE) {
  signal = signal - mean(signal)
}

# Generate spectrogram
spg <- specgram(signal, n=windowLength, Fs=sampleRate, window=hamming(windowLength),
                overlap=windowLength - ceiling(timeStep * sampleRate))
spect <- t(20*log10(abs(spg$S)))

# Plots
time <- spg$t
freq <- spg$f
image(time, freq, spect, col=matlab.like(256), xlab="Time, s", ylab="Frequency, Hz")
