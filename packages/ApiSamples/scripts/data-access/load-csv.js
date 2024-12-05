// Loading a csv file from the specified URL
let t = await grok.data.loadTable('https://public.datagrok.ai/demo/demog.csv')
grok.shell.addTableView(t);
