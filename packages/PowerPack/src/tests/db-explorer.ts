import * as DG from 'datagrok-api/dg';
import {category, expect, expectObject, test} from '@datagrok-libraries/test/src/test';
import {buildValueQuery, ValueQuery} from '@datagrok-libraries/db-explorer/src/query';
import {explicitCatalog, hasCatalogs, ownCatalog} from '@datagrok-libraries/db-explorer/src/utils';

const statement = (query: ValueQuery) =>
  query.sql.split('\n').map((line) => line.trim()).filter((line) => line && !line.startsWith('--')).join(' ');

const projectJoin = {fromSchema: 'chem', fromTable: 'compounds', columnName: 'project_id',
  onSchema: 'core', tableName: 'projects', onColumn: 'project_id', select: ['name as project_name', 'lead']};

const connection = (dataSource: string, db?: string) =>
  DG.DataConnection.create('db-explorer-test', {dataSource, server: 'localhost', ...(db ? {db} : {})});

// Development only: these stay skipped in CI. Take the reason off a test while working on the db-explorer library or
// on Configure Identifiers, and put it back before committing.
const DEV_ONLY = {skipReason: 'Development only: unskip while working on db-explorer or Configure Identifiers'};

category('DB Explorer', () => {
  test('A table of a database without catalogs is addressed as schema.table', async () => {
    expect(statement(buildValueQuery('Postgres', 'molecule_dictionary', 'molregno', 97, 'public')),
      'select a.* from "public"."molecule_dictionary" a where a.molregno = 97');
    expect(statement(buildValueQuery('Postgres', 'molecule_dictionary', 'chembl_id', 'CHEMBL25', 'public')),
      'select a.* from "public"."molecule_dictionary" a where a.chembl_id = \'CHEMBL25\'');
  }, DEV_ONLY);

  test('A catalog goes in front of the schema of the table and of every joined table', async () => {
    const query = buildValueQuery('Databricks', 'compounds', 'compound_id', 'DGC-000001', 'chem', [projectJoin],
      false, 'research');
    expect(statement(query), 'select a.* , b.`name` as `project_name`, b.`lead` ' +
      'from `research`.`chem`.`compounds` a ' +
      'left join `research`.`core`.`projects` b on a.`project_id` = b.`project_id` ' +
      'where a.compound_id = \'DGC-000001\'');
    expectObject(query.columnTables, {project_name: 'projects', lead: 'projects'});
    expectObject(query.columnSchemas, {project_name: 'core', lead: 'core'});
  }, DEV_ONLY);

  test('A join configured for another table is left out', async () => {
    const query = buildValueQuery('Snowflake', 'batches', 'batch_id', 'DGB-000001', 'chem', [projectJoin],
      false, 'RESEARCH');
    expect(statement(query), 'select a.* from "RESEARCH"."chem"."batches" a where a.batch_id = \'DGB-000001\'');
  }, DEV_ONLY);

  test('One row is asked for in the dialect of the database', async () => {
    expect(statement(buildValueQuery('MS SQL', 'batches', '', '', 'chem', [], false, 'research')),
      'select top 1 a.* from "research"."chem"."batches" a');
    expect(statement(buildValueQuery('Postgres', 'batches', '', '', 'chem')),
      'select a.* from "chem"."batches" a limit 1');
    expect(statement(buildValueQuery('MS SQL', 'batches', 'batch_id', 'DGB-000001', 'chem', [], true)),
      'select top 1 a.* from "chem"."batches" a where a.batch_id = \'DGB-000001\'');
  }, DEV_ONLY);

  test('A quote in the value does not end the string', async () => {
    expect(statement(buildValueQuery('Postgres', 'people', 'name', 'O\'Brien', 'public')),
      'select a.* from "public"."people" a where a.name = \'O\'\'Brien\'');
  }, DEV_ONLY);

  test('Catalogs are offered where one query can read across them', async () => {
    expect(hasCatalogs(connection('Databricks')), true);
    expect(hasCatalogs(connection('MS SQL', 'research')), true);
    expect(hasCatalogs(connection('Snowflake', 'TEST')), true);
    // a Postgres catalog is another database: a query cannot name it
    expect(hasCatalogs(connection('Postgres', 'chembl')), false);
    expect(hasCatalogs(connection('MySQL', 'world')), false);
    expect(ownCatalog(connection('MS SQL', 'research')), 'research');
    expect(ownCatalog(connection('Databricks')), null);
  }, DEV_ONLY);

  test('The connection\'s own catalog is not named, in whatever case it is written', async () => {
    const research = connection('Snowflake', 'research');
    expect(explicitCatalog(research, 'research') ?? null, null);
    expect(explicitCatalog(research, 'RESEARCH') ?? null, null);
    expect(explicitCatalog(research, null) ?? null, null);
    expect(explicitCatalog(research, 'ARCHIVE'), 'ARCHIVE');
    expect(explicitCatalog(connection('Databricks'), 'main'), 'main');
  }, DEV_ONLY);
});
