//tags: TableQueryBuilder, TableQuery

grok.data.db.buildQuery('Samples:PostgresNorthwind', 'products')
        .leftJoin('categories', ['categoryid'], ['categoryid'])
        .innerJoin('suppliers', ['supplierid'], ['supplierid'])
        .groupBy(['companyname', 'products.supplierid', 'country'])
        .having('country', 'in (USA, UK)')
        .limit(5)
        .build()
        .executeTable()
    .then((df) => grok.shell.addTable(df))
    .catch((e) => grok.shell.error(e));
