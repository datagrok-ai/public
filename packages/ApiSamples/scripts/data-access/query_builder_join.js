//tags: TableQueryBuilder, TableQuery

grok.data.db.buildQuery('Samples:PostgresNorthwind', 'products')
    .then((builder) => builder
        .leftJoin('categories', ['categoryid'], ['categoryid']) // TableQueryBuilder
        .innerJoin('suppliers', ['supplierid'], ['supplierid']) // TableQueryBuilder
        .groupBy(['companyname', 'products.supplierid', 'country']) // TableQueryBuilder
        .having('country', 'in (USA, UK)') // TableQueryBuilder
        .limit(5) // TableQueryBuilder
        .build() // TableQuery
        .executeTable()) // DataFrame
    .then((df) => grok.shell.addTable(df));
