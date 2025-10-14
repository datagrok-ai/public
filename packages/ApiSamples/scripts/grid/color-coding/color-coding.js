const t = grok.data.demo.demog();

// Conditional color-coding for numerical columns
t.col('height').meta.colors.setConditional({'20-170': '#00FF00', '170-190': "#220505"});

// Linear color-coding for numerical columns
t.col('age').meta.colors.setLinear(['#ff0000', '#ffff00', '#00ff00'], {min: 19, belowMinColor: '#ac3400', max: 70});

// Linear color-coding with absolute values for numerical columns
t.col('weight').meta.colors.setLinearAbsolute({58.31: '#73aff5', 137.65: '#ffa500', 177.32: '#50af28', 217: '#9467bd'}, {aboveMaxColor: '#ff0000'});

// Categorical color-coding for string columns
t.col('race').meta.colors.setCategorical({'Asian': 4278190335, 'Black': 4286578816});

grok.shell.addTableView(t);

// To turn the color-coding off, use either approach:
// t.col('height').meta.colors.setDisabled();
