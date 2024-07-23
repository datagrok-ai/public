1. SPGI_v2_infinity.csv - add viewers+layouts
2. Run the script and add viewers+layouts


```js
let t = grok.data.demo.demog();
let view = grok.shell.addTableView(t);
let plot = view.scatterPlot({
x: 'height', y: 'weight',
size: 'age',
color: 'race',
});
t.col('height').getRawData()[0] = NaN; t.col('height').getRawData()[1] = Infinity;
plot.setOptions({
showRegressionLine: true, markerType: 'square'
});
```

---
{
  "order": 18
}
