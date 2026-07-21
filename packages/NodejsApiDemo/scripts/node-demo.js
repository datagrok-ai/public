//name: NodeDemo
//description: Server-side js-api from a nodejs script - dapi, DataFrame/column ops, package python and JS functions, package files
//language: nodejs
//input: dataframe df
//input: string colName
//output: string user
//output: dataframe enriched
//output: double mean
//output: int jsSum
//output: string fileText

// grok/DG globals and the caller's auth come from the platform - no bootstrap needed
user = (await grok.dapi.users.current()).friendlyName;

// DataFrame + column ops
enriched = df.clone();
enriched.columns.addNewFloat('doubled').init((i) => df.col(colName).get(i) * 2);

// package python function - the dataframe crosses languages
mean = await grok.functions.call('NodejsApiDemo:PyMean', {df, colName});

// package JS function, loaded and executed in-process
await loadPackage('ApiTests');
jsSum = await grok.functions.call('ApiTests:dummyPackageFunction', {a: 20, b: 22});

// package AppData file via dapi
fileText = (await grok.dapi.files.readAsText('System:AppData/NodejsApiDemo/hello.txt')).trim();
