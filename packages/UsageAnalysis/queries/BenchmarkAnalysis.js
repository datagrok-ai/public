//language: javascript
//input: dataframe result
//output: dataframe out  
let pivot = result
  .groupBy(['test'])
  .pivot('build_name')
  .min('duration')
  .aggregate(); 
for(let col of pivot.columns){
	col.name = col.name.split(' ')[0]; 
}

out = pivot;
