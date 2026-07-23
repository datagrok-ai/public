//name: NodeJS Column List
//description: column list input
//language: nodejs
//input: dataframe df
//input: column_list cols
//output: dataframe result

// TODO(GROK-20443): drop the legacy dataframe-js branch once the Node js-api
// runtime (reddata side) is on master and deployed.
result = (typeof DG !== 'undefined')
  ? DG.DataFrame.fromColumns(cols.slice(1).map((name) => df.col(name)))
  : df.select(...cols.slice(1));
