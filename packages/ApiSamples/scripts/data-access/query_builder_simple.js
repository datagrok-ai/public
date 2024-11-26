//tags: TableQueryBuilder, TableQuery

grok.data.db.buildQuery('Samples:PostgresNorthwind', 'orders')
    .then((builder) => builder
        .select(['shipcity', 'shipcountry', 'orderdate']) // TableQueryBuilder
        .where('orderdate', 'before 1996-08-26', DG.TYPE.DATE_TIME) // TableQueryBuilder
        .sortBy('orderdate', false) // TableQueryBuilder
        .limit(10) // TableQueryBuilder
        .build() // TableQuery
        .executeTable()) // DataFrame
    .then((df) => grok.shell.addTable(df));
