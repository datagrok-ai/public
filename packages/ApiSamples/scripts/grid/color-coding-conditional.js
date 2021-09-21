let t = grok.data.demo.demog();

t.col('height').tags[DG.TAGS.COLOR_CODING_TYPE] = 'Conditional';
t.col('height').tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = `{"20-170":"#00FF00","170-190":"#220505"}`;

t.col('age').tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
t.col('age').tags[DG.TAGS.COLOR_CODING_LINEAR] = `[${DG.Color.orange}, ${DG.Color.green}]`;

grok.shell.addTableView(t);