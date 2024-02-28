#name: Octave Exponential Operator Test
#description: Exponential Operator Test
#language: octave
#output: dataframe result
#test: OctaveExponentialOperatorTest().columns.length == 1

A = [-1.0000,0,0,0,0,0,0,0,1.0000,1.0000;1.0000,-1.0000,0,0,0,0,0,0,0,0;0,1.0000,-1.0000,0,0,0,0,0,0,0;0,0,1.0000,-1.0000,0,0,0,0,0,0;0,0,0,1.0000,-1.0000,0,0,0,0,0;0,0,0,0,1.0000,-1.0000,0,0,0,0;0,0,0,0,0,1.0000,-1.0000,0,0,0;0,0,0,0,0,0,1.0000,-1.0000,0,0;0.5000,0,0,0,0,0,0,0,0,-3.0000;0,0,0,0,0,0,0,1.0000,0,0];
b = [0;0;0;0;0;0;0;0;0;2000];
exponential = A^-1*b;
result = [{'Value'}; num2cell(exponential)];
