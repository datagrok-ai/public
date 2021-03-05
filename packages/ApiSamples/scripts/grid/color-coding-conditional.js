let t = grok.data.demo.demog();
t.col('height').tags['color-coding-conditional'] = `{"20-30":"#00FF00","30-200":"#220505"}`;
grok.shell.addTableView(t);