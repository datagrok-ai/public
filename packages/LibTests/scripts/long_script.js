//name: LongScript
//language: javascript
//input: double a
//input: double b
//output: double res
//output: dataframe df {viewer: ScatterPlot()}

const delay = (delayInms) => {
  return new Promise((resolve) => setTimeout(resolve, delayInms));
};
  
await delay(1000);

res = a+b;
df = grok.data.demo.demog(res);
