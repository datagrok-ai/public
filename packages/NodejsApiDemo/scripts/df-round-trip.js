//name: DfRoundTrip
//description: dataframe in/out + column-name param, computed with the js-api (DG)
//language: nodejs
//input: dataframe df
//input: string colName
//output: dataframe result
//output: double colMax
//output: int rows

let dgapi;
try { dgapi = require('datagrok-api/datagrok'); }
catch (_) { dgapi = require('/tmp/datagrok-api/datagrok.js'); }
if (!globalThis.__dgReady)
  globalThis.__dgReady = dgapi.startDatagrok({apiUrl: @DATAGROK_API_URL, apiToken: @USER_API_KEY, detached: true});
await globalThis.__dgReady;
grok.dapi.token = @USER_API_KEY;

// legacy JKG header delivers dataframe-js objects; the DG header delivers DG.DataFrame
const t = (typeof df.toCSV === 'function' && typeof df.toCsv !== 'function')
  ? DG.DataFrame.fromCsv(df.toCSV(true)) : df;
rows = t.rowCount;
colMax = t.col(colName).stats.max;
const out = t.clone();
out.columns.addNewFloat('doubled').init((i) => t.col(colName).get(i) * 2);
result = out;
result.toCSV = () => out.toCsv();  // legacy output path calls toCSV(true)
