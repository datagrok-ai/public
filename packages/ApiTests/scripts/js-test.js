//name: JavaScript Params Test
//language: javascript
//tags: test
//input: int i
//input: double d
//input: bool b
//input: string s
//input: datetime dt
//input: map m
//input: dataframe df
//input_: column col
//output: int ri
//output: double rd
//output: bool rb
//output: string rs
//output: datetime rdt
//output: map rm
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
df.columns.addNew('columnTest', 'int');
rdf = df;
