/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBSchemaInfo} from './db-index-tools';


/**
 * a subclass of DG.DbInfo that adds support for built-in DB relations that might not be indexed yet
 */
export class BuiltinDBInfoMeta extends DG.DbInfo {
  constructor(dart: any) {
    super(dart);
  }

  public static async fromConnection(connection: DG.DataConnection): Promise<BuiltinDBInfoMeta> {
    const dbInfo = await grok.data.db.getInfo(connection);
    return new BuiltinDBInfoMeta(dbInfo.dart);
  }

  private _cachedRelations: Map<string, DG.DbRelationInfo[]> = new Map();
  async getRelationsForSchema(schemaName: string): Promise<DG.DbRelationInfo[]> {
    if (this._cachedRelations.has(schemaName))
      return this._cachedRelations.get(schemaName)!;
    // for dbs without schemas, return all relations
    const annotatedRelations: DG.DbRelationInfo[] = (await super.getRelations())
      .filter((rel) => (!rel.fromSchema && !rel.toSchema) ||
              (rel.fromSchema == schemaName || rel.toSchema == schemaName));

    this._cachedRelations.set(schemaName, annotatedRelations);
    // this._relations = annotatedRelations;
    const annotatedRelationsSet = new Set<string>();
    annotatedRelations.forEach((rel) => {
      annotatedRelationsSet.add(`${rel.fromSchema}.${rel.fromTable}.${(rel.fromColumns ?? []).join(',')}->${rel.toSchema}.${rel.toTable}.${(rel.toColumns ?? []).join(',')}`);
    });
    // collect builtin relations here and add to _relations. remove this chunk in future once everything is indexed
    const conInfo = await DBSchemaInfo.getSchemaInfo(this.connection.id, schemaName);
    conInfo.descriptor.references
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
    return annotatedRelations;
  }
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

export function getDBTableMetaData(table: DG.TableInfo): DG.DbTableProperties {
  const props: (keyof DG.DbTableProperties)[] = ['comment', 'llmComment', 'domains', 'rowCount'];
  return Object.fromEntries(props.map((key) => ([key, table.tags.get(key) ?? undefined])));
}

export function getDBColumnMetaData(column: DG.ColumnInfo): DG.DbColumnProperties {
  const props: (keyof DG.DbColumnProperties)[] = ['comment', 'llmComment', 'isUnique', 'min', 'max', 'values', 'sampleValues', 'uniqueCount'];
  return Object.fromEntries(props.map((key) => ([key, column.tags.get(key) ?? undefined])));
}
