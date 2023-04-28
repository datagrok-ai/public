//name: NodeJS Params Test
//language: nodejs
//tags: test
//input: int i
//input: double d
//input: bool b
//input: string s
//input: datetime dt
//input_: map m
//input: dataframe df
//input_: column col
//output: int ri
//output: double rd
//output: bool rb
//output: string rs
//output: datetime rdt
//output_: map rm
//output: dataframe rdf
//output_: column rcol
ri = i / 2;
rd = d + 60;
rb = !b;
rs = s + s;
rdt = new Date(dt);
rdt.setDate(rdt.getDate() + 10);
rm = m;
rm['b'] = ri;
col = df.clone().col('height');
col.name = 'column';
df.columns.add(col);
rdf = df;
