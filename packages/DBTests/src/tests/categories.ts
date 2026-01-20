import {before, category} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

const categories = {'Postgresql': 'PostgreSQLDBTests', 'Oracle': 'OracleDBTests', 'Neo4j': 'Neo4jDBTests', 'MySql': 'MySqlDBTests', 'MSSQL': 'MSSQLDBTests', 'MariaDb': 'MariaDbDBTests', 'ClickHouse': 'ClickHouseDBTests'};

for (const [cat, conn] of Object.entries(categories)) {
    category(cat, () => {
       before(async () => {
            const connection: DG.DataConnection = await grok.functions.eval(`DbTests:${conn}`);
            const res = await connection.test();
            if (res !== 'ok')
                throw new Error(`Connection ${conn} is not available`);
       });
    });
}