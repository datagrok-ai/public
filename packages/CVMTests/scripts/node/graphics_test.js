//name: NodeJS Graphics
//description: graphics output column input
//language: nodejs
//input: dataframe df
//input: column xName {type:numerical; allowNulls:false} [Column x]
//input: column yName {type:numerical; allowNulls:false} [Column y]
//output: graphics rect

$$.svg("<svg><rect width=180 height=190 style='fill: blue;'/></svg>");
