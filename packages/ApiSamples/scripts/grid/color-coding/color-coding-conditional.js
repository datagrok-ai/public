// Set color coding

let t = grok.data.demo.demog();

// Linear or conditional color coding for numerical or datetime columns
t.col('height').meta.colors.setConditional({"20-170": "#00FF00", "170-190": "#220505"});
t.col('age').meta.colors.setLinear([DG.Color.orange, DG.Color.green]);
t.col('started').meta.colors.setLinear([DG.Color.white, DG.Color.red]);
t.col('weight').meta.colors.setConditional({'<100': DG.Color.green, '100-200': '#ff0000'});

// Categorical for string columns
t.col('site').meta.colors.setCategorical({'New York': DG.Color.orange});

// Additional matching options for categorical
t.col('race').meta.colors.setCategorical(
  {'sian': DG.Color.orange},
  {fallbackColor: 'black', matchType: "regex"});

grok.shell.addTableView(t);

// Get color coding programmatically
grok.shell.info(`Column 'subj' color coding: ${t.col('subj').meta.colors.getType()}`);     // Off
grok.shell.info(`Column 'height' color coding: ${t.col('height').meta.colors.getType()}`); // Conditional
grok.shell.info(`Column 'age' color coding: ${t.col('age').meta.colors.getType()}`);       // Linear
