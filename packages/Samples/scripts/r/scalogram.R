#name: Scalogram R
#description: CWT (Morlet wavelet) scalogram plot
#reference: https://en.wikipedia.org/wiki/Scaleogram
#language: r
#tags: demo, viewers
#sample: sensors/ecg.csv
#input: dataframe data [Input data table]
#input: column signal {type:numerical; allowNulls:false} [Column with signal]
#input: double sampleRate = 256.0 [Signal sample rate, in Hz]
#input: int octaves = 8 [Octaves in scale values]
#input: int voices = 16 [Voices in scale values]
#input: bool removeDc = TRUE [Flag to force DC component removal]
#output: graphics [Scalogram plot]

require(Rwave)
require(colorRamps)

# Select signal column
signal = data[[signal]]

# Remove DC
if (removeDc == TRUE) {
  signal = signal - mean(signal)
}

# Computes the scale values for a given number of octaves and voices
getScales <- function(numOctave, numVoice) {
  scale = vector(length=(numOctave * numVoice))
  scale.index = 1
  for (octave in 0:(numOctave - 1)) {
    for (voice in 0:(numVoice - 1)) {
      scale[scale.index] = 2^(octave + (voice / numVoice))
      scale.index = scale.index + 1
    }
  }
  return (scale)
}

# Computes the corresponding frequency values for a given vector of scales,
# center frequency, and sample rate
scaleToFreq <- function(scale, center, sampleRate) {
  return (center / (scale * sampleRate * 2.0))
}

# Perform CWT
coef = cwt(signal, noctave=octaves, nvoice=voices, w0=(2.0 * pi * 0.8125), twoD=TRUE, plot=FALSE)

# Plots
time = seq(0, (length(signal) - 1) / sampleRate, by=(1.0 / sampleRate))
scales = floor(scaleToFreq(getScales(octaves, voices), 0.8125, 1.0 / sampleRate) * 10) / 10
image(time, 1:dim(coef)[2], abs(coef)^2.0, yaxt="n", ylab="Frequency, Hz", xlab="Time, s",
      col=matlab.like(256))
axis(2, at=seq(1, (octaves * voices)), labels=scales)
