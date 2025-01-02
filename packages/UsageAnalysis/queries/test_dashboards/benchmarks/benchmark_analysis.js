//language: javascript
//input: dataframe result
//output: dataframe out  
let pivot = result
  .groupBy(['test'])
  .pivot('build_name')
  .min('duration')
  .aggregate(); 
let cols = [];
for(let col of pivot.columns){
	col.name = col.name.split(' ')[0]; 
  if(col.name != 'test')
  	cols.push(col);
}
 
for(let col of cols){
	pivot.columns.remove(col.name);
}

cols = cols.sort((a, b) => a.name.localeCompare(b.name));
for(let col of cols){
	pivot.columns.add(col);
}

out = pivot;
