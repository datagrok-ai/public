let t = grok.data.demo.demog();
t.col('height').tags['color-coding-type'] = 'Conditional';
t.col('height').tags['color-coding-conditional'] = `{"20-170":"#00FF00","170-190":"#220505"}`;
grok.shell.addTableView(t);