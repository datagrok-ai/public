//tags: TableQueryBuilder, TableQuery

grok.data.db.buildQuery('Samples:PostgresNorthwind', 'orders')
        .select(['shipcity', 'shipcountry', 'orderdate'])
        .where('orderdate', 'before 1996-08-26', DG.TYPE.DATE_TIME)
        .sortBy('orderdate', false)
        .limit(10)
        .build()
        .executeTable()
    .then((df) => grok.shell.addTable(df))
    .catch((e) => grok.shell.error(e));
