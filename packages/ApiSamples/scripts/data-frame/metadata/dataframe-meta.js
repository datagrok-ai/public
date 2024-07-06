// Different ways of associating metadata with dataframe and columns and manipulating it

const demog = grok.data.demo.demog();
demog.col('sex').tags.friendlyName = 'GENDER';
demog.col('weight').meta.format = '0.0';
grok.shell.addTableView(demog);