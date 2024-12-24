
//tags: DataQuery
//help-url: https://datagrok.ai/help/access/databases#parameterized-queries
// An example of using parameterized query

let t = await grok.data
  .query("Samples:MSSQLByInt", { orderid: 10 })
grok.shell.info("Number of orders: " + t.rowCount);
