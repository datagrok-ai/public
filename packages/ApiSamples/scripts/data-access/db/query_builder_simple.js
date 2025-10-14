//tags: TableQueryBuilder, TableQuery, query

try {
  let df = await grok.data.db
  .buildQuery("Samples:PostgresNorthwind", "orders")
  .select(["shipcity", "shipcountry", "orderdate"])
  .where("orderdate", "before 1996-08-26", DG.TYPE.DATE_TIME)
  .sortBy("orderdate", false)
  .limit(10)
  .build()
  .executeTable()
  grok.shell.addTableView(df);
} catch (e) {
  grok.shell.error(e);
}
