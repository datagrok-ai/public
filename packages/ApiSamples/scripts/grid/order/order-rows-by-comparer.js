let t = grok.data.demo.demog();
let view = grok.shell.addTableView(t);

let height = t.col('height');
let weight = t.col('weight');
let bmi = (i) => height.get(i) / weight.get(i);

view.grid.sortIndexes((i, j) => bmi(i) - bmi(j));
