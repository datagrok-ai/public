1. SPGI_v2_infinity.csv - add viewers+layouts
2. Run the script and add viewers+layouts. For viewers, make sure 'height' and 'weight' columns are used


```js
let t = grok.data.demo.demog();
let view = grok.shell.addTableView(t);
let plot = view.scatterPlot({
x: 'height', y: 'weight',
size: 'age',
color: 'race',
});
t.col('height').set(0, NaN); 
t.col('weight').set(0, Infinity);
plot.setOptions({
showRegressionLine: true, markerType: 'square'
});
```

---
{
  "order": 18
}
