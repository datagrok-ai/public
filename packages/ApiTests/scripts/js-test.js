//name: JavaScript Params Test
//language: javascript
//tags: test
//input: int i
//input: double d
//input: bool b
//input: string s
//input: datetime dt
//input: dataframe df
//input: column col
//output: int ri
//output: double rd
//output: bool rb
//output: string rs
//output: datetime rdt
//output_: dataframe rdf
//output_: column rcol
ri = i / 2;
rd = d + 60;
rb = !b;
rs = s + s;
rdt = new Date(dt);
rdt.setDate(rdt.getDate() + 10);
// rdf = df.col('height');
