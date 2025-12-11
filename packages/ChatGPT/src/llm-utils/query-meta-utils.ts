/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBConnectionMeta, DBSchemaInfo} from './db-index-tools';

/** Gathers the information about connection metadata, tables, schemas, columns, relations, comments, embeddings, etc
 * @param connection DataConnection to gather the metadata for
 */
export async function getConnectionMeta(connection: DG.DataConnection) {
  // temporary filler for certain databases
//   let dbMeta: DBConnectionMeta | undefined = undefined;

  // kind of lazy load to save memory
  //   if (connection.name.toLowerCase() === 'chembl')
  //     dbMeta = (await import('./indexes/chembl-index')).chemblIndex;
  //   else if (connection.name.toLowerCase() === 'biologics')
  //     dbMeta = (await import('./indexes/biologics-index')).biologicsIndex;

  const dbInfo = await grok.data.db.getInfo(connection);

  return {
    getSchemas: async () => await dbInfo.getSchemas(),
    getSchemaNames: async () => (await dbInfo.getSchemas()).map((s) => s.name),
    getRelationsForTable: async (schemaName: string, tableName: string) => {
      const conInfo = await DBSchemaInfo.getSchemaInfo(connection.id, schemaName);
      const annotatedRelations = (await dbInfo.getRelations())
        .filter((rel) => (rel.fromSchema == schemaName && rel.fromTable == tableName) ||
                (rel.toSchema == schemaName && rel.toTable == tableName));
      const annotatedRelationsSet = new Set<string>();
      annotatedRelations.forEach((rel) => {
        annotatedRelationsSet.add(`${rel.fromSchema}.${rel.fromTable}.${(rel.fromColumns ?? []).join(',')}->${rel.toSchema}.${rel.toTable}.${(rel.toColumns ?? []).join(',')}`);
      });
      conInfo.descriptor.references
        .filter((obj) => obj.from.split('.')[0] === tableName || obj.to.split('.')[0] === tableName)
        .forEach((ref) => {
          const f = ref.from.split('.');
          const t = ref.to.split('.');
          const ft = f[0];
          const fc = f[1];
          const tt = t[0];
          const tc = t[1];
          if (annotatedRelationsSet.has(`${schemaName}.${ft}.${fc}->${schemaName}.${tt}.${tc}`))
            return;
          const meta = new BuiltinRelationMeta(
            schemaName,
            ft,
            [fc],
            schemaName,
            tt,
            [tc],
          );
          annotatedRelations.push(meta);
        });
    },
    getTablesInSchema: async (schemaName: string) => { 
        const schema = (await dbInfo.getSchemas()).find((s) => s.name === schemaName)!;
        return [await schema.getTables
    },


  };
}

class BuiltinRelationMeta implements DG.DbRelationInfo {
  get comment() {
    return undefined;
  }
  get llmComment() {
    return undefined;
  }
  get cardinality() {
    return undefined;
  }

  dart: any = null;
  get connection(): DG.DataConnection {
    throw new Error('Method not implemented.');
  }
  get isPrimaryPath(): boolean {
    return true;
  }
  constructor(
        public fromSchema: string,
        public fromTable: string,
        public fromColumns: string[],
        public toSchema: string,
        public toTable: string,
        public toColumns: string[],
  ) {
  }
}
