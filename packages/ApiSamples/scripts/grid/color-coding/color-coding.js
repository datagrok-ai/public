let t = grok.data.demo.demog();

// Conditional color-coding for numerical columns
t.col('height').tags[DG.TAGS.COLOR_CODING_TYPE] = 'Conditional';
t.col('height').tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = `{"20-170":"#00FF00","170-190":"#220505"}`;

// Linear color-coding for numerical columns
t.col('age').tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
t.col('age').tags[DG.TAGS.COLOR_CODING_LINEAR] = '["#ff0000", "#ffff00", "#00ff00"]';

// Categorical color-coding for string columns
t.col('race').tags[DG.TAGS.COLOR_CODING_TYPE] = 'Categorical';
t.col('race').tags[DG.TAGS.COLOR_CODING_CATEGORICAL] = `{"Asian":4278190335,"Black":4286578816}`;

grok.shell.addTableView(t);

// To turn the color-coding off, use either approach:
//t.col('height').tags[DG.TAGS.COLOR_CODING_TYPE] = 'Off';
//delete t.col('height').tags[DG.TAGS.COLOR_CODING_TYPE];