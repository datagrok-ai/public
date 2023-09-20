#name: Octave Graphics
#description: graphics output column input
#language: octave
#input: dataframe df
#input: column xName {type:numerical; allowNulls:false}
#input: column yName {type:numerical; allowNulls:false}
#output: graphics barchart

pkg load dataframe;

d = dataframe(df);

bar(d(1:3,[xName]));
