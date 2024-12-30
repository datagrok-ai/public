//name: LongScript
//language: javascript
//input: double a
//input: double b
//output: double res {category: Main result}
//output: double res2 {category: Additional data; precision: 2; units: cm}
//output: double res3 {category: Additional data; precision: 2; units: L}
//output: double res4 {category: Additional data}
//output: double res5 {category: Additional data; precision: 6; units: 1/F}
//output: dataframe df {viewer: ScatterPlot()}

const delay = (delayInms) => {
  return new Promise((resolve) => setTimeout(resolve, delayInms));
};
  
await delay(2000);

res = a+b;
df = grok.data.demo.demog(res);
res2 = Math.random();
res3 = Math.random()+100;
res4 = Math.round(Math.random()*100);
res5 = Math.random()*0.1;
