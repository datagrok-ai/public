import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {plural, text} from '@datagrok-libraries/u2';
import type {IProperty, RowView} from '@datagrok-libraries/u2';
import {nounOf} from '../core/domain-backend.js';
import {DomainTable} from './facade.js';
import {SYSTEM_COLUMNS} from './backend.js';

// a measurement or a yes/no never tells one record from another
const MATCHABLE = ['string', 'int', 'bigint', 'datetime'];
const PAGE = 10000;
const DAY_MS = 24 * 60 * 60 * 1000;

/** A field a linked record is found by; `path` is what the query reads: `<ref>.<key>` for a link of its own. */
export interface LinkField {
  prop: IProperty;
  path: string;
}

/** A reference column of the imported table, offered in the mapping as `<column>.<field>`. */
export interface ImportLink {
  prop: IProperty;
  target: DomainTable;
  fields: LinkField[];
}

/** The record a file row names, or why it names none; both null for a row that names nothing. */
export interface LinkMatch {
  id: string | null;
  error: string | null;
}

/** Finds the records a file's columns name in the tables the imported one refers to: every row must find exactly
 * one, never a guess. */
export class ImportLinks {
  static isLink(prop: IProperty): boolean {
    return DomainTable.isReference(prop) && prop.semType !== 'User' && prop.semType !== 'Group';
  }

  /** The links among `targets`, each with the fields its records can be found by. */
  static async of(targets: IProperty[]): Promise<Map<string, ImportLink>> {
    const open = (p: IProperty) => DomainTable.shared(p.semType!).catch(() => null);
    const links = await Promise.all(targets.filter(ImportLinks.isLink).map(async (prop) => {
      const target = await open(prop);
      if (target === null)
        return null;
      const props = target.properties.filter((p) => !SYSTEM_COLUMNS.some(([name]) => name === p.name));
      const nested = await Promise.all(props.map((p) => ImportLinks.isLink(p) ? open(p) : null));
      const fields = props.flatMap((p, i): LinkField[] => nested[i] !== null ?
        [{prop: p, path: `${p.name}.${nested[i]!.info.businessKey[0] ?? nested[i]!.info.nameColumn ?? 'id'}`}] :
        !ImportLinks.isLink(p) && MATCHABLE.includes(p.propertyType ?? p.type ?? '') ? [{prop: p, path: p.name!}] : []);
      return {prop, target, fields};
    }));
    return new Map(links.flatMap((link) => link === null ? [] : [[link.prop.name!, link]]));
  }

  /** Per row of `values` (each row's values for `fields`, in order), the one record of `link` they find. */
  static async resolve(link: ImportLink, fields: LinkField[], values: unknown[][]): Promise<LinkMatch[]> {
    const keys = values.map((row) => row.map(ImportLinks._read));
    const named = keys.filter((key) => key.some((v) => v !== ''));
    const found = new Map<string, string[]>();
    if (named.length > 0) {
      const conditions = fields.map((f, i) =>
        ImportLinks._condition(f, [...new Set(named.map((key) => key[i]).filter((v) => v !== ''))]));
      const filter: DG.DomainConditionTree = conditions.flatMap((c, i) => i === 0 ? [c] : ['and' as const, c]);
      const expand = [...new Set(fields.filter((f) => f.path.includes('.')).map((f) => f.prop.name!))];
      for (let offset = 0; ; offset += PAGE) {
        const rows = await grok.dapi.domains.table(link.target.address)
          .query({filter, expand, sort: 'id', limit: PAGE, offset}) as RowView[];
        for (const row of rows) {
          const key = JSON.stringify(fields.map((f) => ImportLinks._read(row[f.path])));
          found.set(key, [...found.get(key) ?? [], String(row.id)]);
        }
        if (rows.length < PAGE)
          break;
      }
    }
    const one = nounOf(link.target.info, true);
    return keys.map((key): LinkMatch => {
      if (key.every((v) => v === ''))
        return {id: null, error: null};
      const ids = found.get(JSON.stringify(key)) ?? [];
      const said = fields.map((f, i) => `${DomainTable.captionOf(f.prop)} ${key[i]}`).join(', ');
      return ids.length === 1 ? {id: ids[0], error: null} : {id: null, error: ids.length === 0 ?
        `No ${one} with ${said}` : `${plural(ids.length, one, nounOf(link.target.info))} have ${said}: ` +
        'map one more column to tell them apart'};
    });
  }

  // a date is matched a day either side here and by its day after: the stored instant is the day in the zone it
  // was entered in
  private static _condition(field: LinkField, values: string[]): DG.DomainCondition | DG.DomainConditionTree {
    const type = field.prop.propertyType ?? field.prop.type;
    const days = values.map((v) => Date.parse(`${v}T00:00:00Z`)).filter((t) => !Number.isNaN(t));
    if (type === 'datetime' && days.length > 0) {
      return [{property: field.path, operator: '>=', value: new Date(Math.min(...days) - DAY_MS).toISOString()},
        'and', {property: field.path, operator: '<', value: new Date(Math.max(...days) + 2 * DAY_MS).toISOString()}];
    }
    const numeric = type === 'int' || type === 'bigint';
    return {property: field.path, operator: '=', value: numeric ?
      values.filter((v) => /^-?\d+$/.test(v)).map(Number) : values};
  }

  // a value as the file and the record both spell it: a date as its day
  private static _read(value: unknown): string {
    const date = value instanceof Date ? value : typeof (value as {toDate?: unknown})?.toDate === 'function' ?
      (value as {toDate(): Date}).toDate() : null;
    const iso = typeof value === 'string' ? /^(\d{4}-\d{2}-\d{2})T/.exec(value) : null;
    if (date === null)
      return iso?.[1] ?? text(value).trim();
    const pad = (n: number) => String(n).padStart(2, '0');
    return `${date.getFullYear()}-${pad(date.getMonth() + 1)}-${pad(date.getDate())}`;
  }
}
