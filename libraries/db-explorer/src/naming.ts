/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {EntityOptions, LinkGroup, LinkRef, QueryJoinOptions, SchemaInfo} from './types';

/** Case-insensitive object key lookup */
export function ciGet<T>(obj: {[key: string]: T} | undefined | null, key: string): T | undefined {
  if (!obj || key == null) return undefined;
  if (obj[key] !== undefined) return obj[key];
  const lowerKey = key.toLowerCase();
  for (const k of Object.keys(obj)) {
    if (k.toLowerCase() === lowerKey)
      return obj[k];
  }
  return undefined;
}

export const ciEquals = (a?: string | null, b?: string | null) => (a ?? '').toLowerCase() === (b ?? '').toLowerCase();

/** Words that only say a column holds a key */
const KEY_WORDS = new Set(['id', 'ids', 'qualified', 'key', 'fk', 'identifier', 'ref', 'no', 'num']);
/** A table's name is cleaned of the words most of its schema's tables share (`stg`, `v1`) only in a schema with this many tables */
const MIN_TABLES_FOR_NOISE = 5;
const NOISE_SHARE = 0.6;
/** Card columns given for this many tables mean a curated configuration: a table it says nothing about is one to fold under More */
const CURATED_TABLES = 5;

export const tokens = (name: string) => (name ?? '').split(/[_\s\-.]+/).map((t) => t.toLowerCase()).filter((t) => t.length > 0);

/** The words most tables of a schema share, which say nothing about any one of them */
export function noiseTokens(tables: string[]): Set<string> {
  const noise = new Set<string>();
  if (tables.length < MIN_TABLES_FOR_NOISE)
    return noise;
  const counts = new Map<string, number>();
  for (const t of tables) {
    for (const tk of new Set(tokens(t)))
      counts.set(tk, (counts.get(tk) ?? 0) + 1);
  }
  counts.forEach((n, tk) => {
    if (n / tables.length >= NOISE_SHARE)
      noise.add(tk);
  });
  return noise;
}

const isCodeToken = (t: string, codes: Set<string>) => codes.has(t) || (t.length <= 5 && /^[a-z]+$/.test(t) && !/[aeiouy]/.test(t));

/** Words as a name: codes (identifier prefixes, words without vowels) in capitals, the first word capitalised */
export function asWords(tks: string[], codes: Set<string>): string {
  const ws = tks.map((t) => isCodeToken(t, codes) ? t.toUpperCase() : t);
  if (ws.length > 0)
    ws[0] = ws[0][0].toUpperCase() + ws[0].slice(1);
  return ws.join(' ');
}

const lastWord = (s: string) => {
  const m = s.match(/^(.*?)(\S+)$/);
  return m ? [m[1], m[2]] : ['', s];
};

export function pluralize(s: string): string {
  const [head, last] = lastWord(s);
  let p: string;
  if (/^[A-Z0-9]+$/.test(last)) p = last + 's';
  else if (/s$/i.test(last)) p = last;
  else if (/(ch|sh|x|z)$/i.test(last)) p = last + 'es';
  else if (/[^aeiou]y$/i.test(last)) p = last.slice(0, -1) + 'ies';
  else p = last + 's';
  return head + p;
}

export function singularize(s: string): string {
  const [head, last] = lastWord(s);
  let p = last;
  if (/^[A-Z0-9]+s$/.test(last)) p = last.slice(0, -1);
  else if (/ies$/i.test(last)) p = last.slice(0, -3) + 'y';
  else if (/(ch|sh|x|z|ss)es$/i.test(last)) p = last.slice(0, -2);
  else if (/[^su]s$/i.test(last) && !/is$/i.test(last)) p = last.slice(0, -1);
  return head + p;
}

/** Lower-cases the first letter, unless the first word is a code (ADC, CHEMBL) */
export const lcFirst = (s: string) => /^[A-Z0-9]{2,}\b/.test(s) || s.length === 0 ? s : s[0].toLowerCase() + s.slice(1);

/** The leading tokens every one of the semantic types shares (`DG_BIOLOGICS_`) */
function commonPrefix(lists: string[][]): string[] {
  if (lists.length < 2) return [];
  const out: string[] = [];
  for (let i = 0; ; i++) {
    const t = lists[0][i];
    if (t === undefined || lists.some((l) => l[i] !== t || l.length <= i + 1))
      return out;
    out.push(t);
  }
}

export type EntryPointInfo = {semType: string, table: string, column: string, code?: string};

export type TableName = {
  /** One record: `Drug`, `Purification batch` */
  singular: string;
  /** A list of them: `Drugs`, `Assay results` */
  plural: string;
  /** The identifier code shown after the name: `GROKMOL` */
  code?: string;
  /** Where the name came from */
  source: 'config' | 'semtype' | 'table' | 'physical';
};

/** Words that identify a table, used to tell what is left of a column name pointing at it ("parent") */
export type TableWords = Set<string>;

/** The code of an identifier pattern: `GROKMOL-` of GROKMOL-000001, `CHEMBL` of CHEMBL25 */
export function identifierCode(nonVariablePart?: string | null, matchRegexp?: string | null): string | undefined {
  const fromExample = (nonVariablePart ?? '').replace(/[^A-Za-z0-9]/g, '');
  if (fromExample) return fromExample;
  const m = (matchRegexp ?? '').replace(/^\^/, '').match(/^([A-Za-z][A-Za-z0-9]*)/);
  return m ? m[1] : undefined;
}

/**
 * The names the Links section gives a table, from what the explorer already knows: the configuration's `entities`,
 * the identifier registered for the table (its semantic type and pattern), and the table's name with the words most
 * tables of its schema share left out.
 */
export class TableNames {
  private entities: {[table: string]: EntityOptions} = {};
  private entryPoints: EntryPointInfo[] = [];
  private noise = new Map<string, Set<string>>();
  private comments: {[qualifiedTable: string]: string} = {};
  private schemaSeen: SchemaInfo | null = null;

  addEntities(entities: {[table: string]: EntityOptions}) {
    this.entities = {...this.entities, ...entities};
  }

  addEntryPoint(ep: EntryPointInfo) {
    this.entryPoints.push(ep);
  }

  /** Reads the schema's table lists and table comments, once per schema object */
  useSchema(schema: SchemaInfo) {
    if (this.schemaSeen === schema) return;
    this.schemaSeen = schema;
    this.noise.clear();
    for (const [schemaName, tables] of Object.entries(schema.references ?? {}))
      this.noise.set(schemaName.toLowerCase(), noiseTokens(Object.keys(tables)));
    this.comments = {};
    for (const [schemaName, tables] of Object.entries(schema.comments ?? {})) {
      for (const [table, comment] of Object.entries(tables))
        this.comments[`${schemaName}.${table}`.toLowerCase()] = comment;
    }
  }

  entity(schema: string, table: string): EntityOptions | undefined {
    return ciGet(this.entities, `${schema}.${table}`) ?? ciGet(this.entities, table);
  }

  entryPoint(table: string): EntryPointInfo | undefined {
    return this.entryPoints.find((ep) => ciEquals(ep.table, table));
  }

  description(schema: string, table: string): string | undefined {
    return this.entity(schema, table)?.description ?? this.comments[`${schema}.${table}`.toLowerCase()];
  }

  private get codes(): Set<string> {
    return new Set(this.entryPoints.map((ep) => ep.code?.toLowerCase()).filter((c): c is string => !!c));
  }

  /** The table's own words: its name without the schema's shared words */
  tableTokens(schema: string, table: string): string[] {
    const noise = this.noise.get(schema.toLowerCase()) ?? new Set<string>();
    const all = tokens(table);
    const kept = all.filter((t) => !noise.has(t));
    return kept.length > 0 ? kept : all;
  }

  /** The semantic type as a name, when it names the records (DG_BIOLOGICS_DRUG_ID) rather than the identifier (CHEMBL_ID) */
  private semTypeTokens(ep: EntryPointInfo): string[] | null {
    const prefix = commonPrefix(this.entryPoints.map((e) => tokens(e.semType)));
    const tk = tokens(ep.semType).slice(prefix.length);
    while (tk.length > 0 && KEY_WORDS.has(tk[tk.length - 1]))
      tk.pop();
    if (tk.length === 0) return null;
    if (ep.code && tk.join('') === ep.code.toLowerCase()) return null;
    if (ciEquals(ep.semType, ep.column)) return null;
    return tk;
  }

  name(schema: string, table: string): TableName {
    const cfg = this.entity(schema, table);
    const ep = this.entryPoint(table);
    const code = cfg?.short ?? ep?.code;
    // a short word of the table's own code is a code too: `adc` of GROKADC
    const own = (ep?.code ?? '').toLowerCase();
    const codes = new Set([...this.codes, ...[...tokens(table), ...(ep ? tokens(ep.semType) : [])]
      .filter((t) => own && t.length <= 4 && own.includes(t))]);
    const withCode = (n: Omit<TableName, 'code'>): TableName =>
      ({...n, code: code && !ciEquals(code, n.singular) && !ciEquals(code, n.plural) ? code : undefined});
    if (cfg?.label)
      return withCode({singular: cfg.label, plural: cfg.plural ?? pluralize(cfg.label), source: 'config'});
    const tableTk = this.tableTokens(schema, table);
    const tableName = asWords(tableTk, codes);
    const semTk = ep ? this.semTypeTokens(ep) : null;
    if (semTk) {
      // the table's name says more when it holds the type's words and others ("Purification batches" for PURIFICATION)
      const tableSingular = tokens(singularize(tableTk.join(' ')));
      const tableSaysMore = tableTk.length > semTk.length && semTk.every((t) => tableSingular.includes(t) || tableTk.includes(t));
      if (!tableSaysMore) {
        const singular = asWords(semTk, codes);
        return withCode({singular, plural: pluralize(singular), source: 'semtype'});
      }
    }
    return withCode({singular: singularize(tableName), plural: tableName, source: tableTk.length > 0 ? 'table' : 'physical'});
  }

  /** Words that identify the table: its code, its own words and their singulars, its type's words */
  words(schema: string, table: string, column?: string): TableWords {
    const w = new Set<string>();
    const add = (ts: string[]) => ts.forEach((t) => w.add(t));
    add(tokens(table));
    add(tokens(singularize(this.tableTokens(schema, table).join(' '))));
    const ep = this.entryPoint(table);
    if (ep) {
      add(tokens(ep.semType));
      if (ep.code) w.add(ep.code.toLowerCase());
    }
    const cfg = this.entity(schema, table);
    if (cfg?.short) w.add(cfg.short.toLowerCase());
    if (column) add(tokens(column));
    return w;
  }

  /** Whether the configuration curates its cards, so that a table it says nothing about is a view to fold away */
  isMore(schema: string, table: string, known: {unique: boolean, card: boolean, header: boolean}, cardTables: number): boolean {
    const cfg = this.entity(schema, table);
    if (cfg?.more !== undefined) return cfg.more;
    return cardTables >= CURATED_TABLES && !this.entryPoint(table) && !known.unique && !known.card && !known.header;
  }
}

/** What a column pointing at a table says beyond the table's name: `parent` of parent_molregno */
export function roleWords(column: string, target: TableWords): string[] {
  return tokens(column).filter((t) => !KEY_WORDS.has(t) && !target.has(t) && !target.has(singularize(t)));
}

/** A table name without the parent's name in front: `Regions` of `Sequence regions`, under a sequence */
export function withoutParent(name: string, parent: TableWords): string {
  const ws = name.split(' ');
  if (ws.length < 2 || !parent.has(ws[0].toLowerCase()))
    return name;
  ws.shift();
  ws[0] = ws[0][0].toUpperCase() + ws[0].slice(1);
  return ws.join(' ');
}

const tagOf = (col: DG.Column, tag: string): string | null => {
  const v = col.getTag(tag);
  return v && v !== 'null' ? v : null;
};

/**
 * The groups of the Links section of a row (the first row of [df], read from [curSchema].[curTable]): one for every
 * key of the row that other rows point at (the row itself, and each row a column of it points at). The value of a
 * group is the value of the column it is found from, never another column's.
 */
export function linkGroups(schema: SchemaInfo, curSchema: string, curTable: string, df: DG.DataFrame,
  joinOptions: QueryJoinOptions[] = []): LinkGroup[] {
  if (!df || df.rowCount === 0) return [];
  const own: LinkGroup[] = [];
  const others: LinkGroup[] = [];
  const seen = new Set<string>();
  const joins = joinOptions.filter((j) => ciEquals(j.fromTable, curTable) && (!j.fromSchema || ciEquals(j.fromSchema, curSchema)));
  const add = (into: LinkGroup[], isOwn: boolean, target: LinkRef, field: string, value: any, label?: string) => {
    if (value == null || value === '') return;
    const refBys = ciGet(ciGet(ciGet(schema.referencedBy, target.schema), target.table), target.column);
    if (!refBys || refBys.length === 0) return;
    const key = `${target.schema}.${target.table}.${target.column}=${value}`.toLowerCase();
    if (seen.has(key)) return;
    seen.add(key);
    into.push({own: isOwn, target, value, field, label,
      links: refBys.map((r) => ({schema: r.refSchema, table: r.refTable, column: r.refColumn, reverseLabel: r.reverseLabel}))});
  };
  for (const col of df.columns.toList()) {
    const table = tagOf(col, DG.Tags.DbTable) ?? curTable;
    const colSchema = tagOf(col, DG.Tags.DbSchema) ?? curSchema;
    const column = tagOf(col, DG.Tags.DbColumn) ?? col.name;
    const value = col.isNone(0) ? null : col.get(0);
    const isOwnTable = ciEquals(table, curTable) && ciEquals(colSchema, curSchema);
    // the column is a key other rows point at
    add(isOwnTable ? own : others, isOwnTable, {schema: colSchema, table, column}, col.name, value);
    // the column points at a key
    const ref = ciGet(ciGet(ciGet(schema.references, colSchema), table), column);
    if (ref) {
      add(others, false, {schema: ref.refSchema, table: ref.refTable, column: ref.refColumn}, col.name, value, ref.label);
      continue;
    }
    const join = isOwnTable ? joins.find((j) => ciEquals(j.columnName, column)) : undefined;
    if (join)
      add(others, false, {schema: join.onSchema ?? curSchema, table: join.tableName, column: join.onColumn}, col.name, value);
  }
  return [...own, ...others];
}
