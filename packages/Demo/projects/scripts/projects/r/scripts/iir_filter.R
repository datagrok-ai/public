#name: IIR Filter
#description: IIR Butterworth filter (high or low pass)
#reference: https://en.wikipedia.org/wiki/Butterworth_filter
#language: r
#tags: demo
#sample: sensors/ecg.csv
#input: dataframe data [Input data table]
#input: column signal {type:numerical; allowNulls:false} [Column with signal]
#input: int order = 6 [Filter order]
#input: double frequency = 0.35 [Critical frequency of the filter, in Hz]
#input: double sampleRate = 256.0 [Signal sample rate, in Hz]
#input: string type = high {choices: ["low", "high"]} [Type of filter: low, high]
#output: dataframe filtered {action:join(data)} [Filtered signal vector]
#output: graphics [Raw and filtered signals]

require(signal)

# Select signal column
signal = data[[signal]]

# Filtering
filt <- butter(order, 2.0 * pi * frequency / sampleRate, type=type)
filtered = filtfilt(filt, signal)

# Plots
time = seq(0, length(signal)-1, by=1) / sampleRate
plot(time, signal, type="l", col="blue", ylab="Amplitude", xlab="Time, s",
     ylim=c(min(signal, filtered), max(signal, filtered)))
lines(time, filtered, type="l", col="red")
grid(10, 10)

filtered = data.frame(filtered)
