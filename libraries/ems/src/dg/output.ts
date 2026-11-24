import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Form, Rows, component, div, notify, signal, untracked} from '@datagrok-libraries/u2';
import type {IProperty, ReadonlySignal, RowView} from '@datagrok-libraries/u2';
import {ColumnsInput, AggregatedColumnsInput} from '@datagrok-libraries/u2/dg';
import type {ColumnAggregation} from '@datagrok-libraries/u2/dg';
import type {DomainSource} from '../core/domain-source.js';
import {DomainTable} from './facade.js';
import {DomainApp} from './app.js';
import {SYSTEM_COLUMNS} from './backend.js';
import {DomainErrors} from './errors.js';
import {DomainDataTable} from './grid.js';

const LIMIT = 100000;
// the statistics in the field's own unit: a count, a variance or a shape has another one, or none
const IN_UNIT: string[] = [DG.AGG.FIRST, DG.AGG.MIN, DG.AGG.MAX, DG.AGG.SUM, DG.AGG.MED, DG.AGG.AVG, DG.AGG.STDEV,
  DG.AGG.Q1, DG.AGG.Q2, DG.AGG.Q3, 'geomean'];
// all else is computed from the values, so it means nothing across units
const COUNTS: string[] = [DG.AGG.TOTAL_COUNT, DG.AGG.VALUE_COUNT, DG.AGG.UNIQUE_COUNT, DG.AGG.MISSING_VALUE_COUNT];

/** Where an aggregated field lives: its table, the column naming the list's row (`batch_id.compound_id` one table
 * away), the field, and every table a write to which moves the value. */
interface Linked {
  address: string;
  path: string;
  column: string;
  tables: string[];
  measured: boolean;
}

const own = (p: IProperty) => !SYSTEM_COLUMNS.some(([name]) => name === p.name) && !Rows.isService(p.name!);

/** The aggregation's column header and the name its value is kept under: `avg(Batches.Purity (%))`. */
export function outputKey(a: ColumnAggregation): string {
  return `${a.aggregation}(${a.column})`;
}

/** The table view's Output and Aggregations under the filters: the platform's column pickers over the table's own
 * columns and the fields of the tables that refer to it, aggregated per row by the platform's `groupBy`. Built the
 * first time it is shown. */
export function domainOutput(source: DomainSource, table: DomainTable, exclude: string[],
  shown: ReadonlySignal<boolean>) {
  const props = table.properties.filter(own);
  const columns = signal(props.map((p) => p.name!).filter((name) => !exclude.includes(name)));
  const aggregations = signal<ColumnAggregation[]>([]);
  const values = signal(new Map<string, Record<string, unknown>>());
  const linked = new Map<string, Linked>();
  // a key names a row: typed as text, it is counted, never summed — and never the new row's first pick
  const keys = new Set<string>();
  let built = false;
  let gen = 0;

  // a domain property's type is a column type: the pickers offer each field what the platform offers its type
  const frame = (fields: [string, IProperty][]) => DG.DataFrame.fromColumns(fields.map(([name, p]) =>
    DG.Column.fromType((keys.has(name) ? DG.COLUMN_TYPE.STRING : p.propertyType ?? DG.COLUMN_TYPE.STRING) as
      DG.ColumnType, name, 0)));

  // the tables that refer to this one, and theirs: `Batches.Purity (%)`, `Assay Results.Reported CLint`; one that
  // cannot be opened is left out, and one reached twice is told apart by the link, as the children tabs do
  const linkedFields = async () => {
    type Ref = {schema: string, table: string, fkColumn: string, label: string};
    const open = (refs: Ref[]) => Promise.allSettled(refs.map((r) => DomainTable.shared(`${r.schema}.${r.table}`)));
    const reached: {from: DomainTable, path: string, via: string, tables: string[]}[] = [];
    const children = table.info.childTables;
    for (const [i, child] of (await open(children)).entries()) {
      if (child.status === 'rejected')
        continue;
      reached.push({from: child.value, path: children[i].fkColumn, via: children[i].label,
        tables: [child.value.address]});
      const grands = child.value.info.childTables.filter((g) => `${g.schema}.${g.table}` !== table.address);
      for (const [j, grand] of (await open(grands)).entries()) {
        if (grand.status === 'fulfilled') {
          reached.push({from: grand.value, path: `${grands[j].fkColumn}.${children[i].fkColumn}`,
            via: `${children[i].label} › ${grands[j].label}`, tables: [grand.value.address, child.value.address]});
        }
      }
    }
    const fields: [string, IProperty][] = [];
    for (const {from, path, via, tables} of reached) {
      const title = DomainApp.titleOf(from.info);
      const name = reached.filter((r) => r.from.address === from.address).length > 1 ? `${title} (${via})` : title;
      const units = from.properties.some((p) => p.name === 'unit');
      for (const p of from.properties.filter((p) => own(p) && !DomainTable.isReference(p))) {
        linked.set(`${name}.${DomainTable.captionOf(p)}`, {address: from.address, path, column: p.name!, tables,
          measured: units && /(^|_)quantity$/.test(p.name!)});
        fields.push([`${name}.${DomainTable.captionOf(p)}`, p]);
        if (from.info.businessKey.includes(p.name!))
          keys.add(`${name}.${DomainTable.captionOf(p)}`);
      }
    }
    return fields;
  };

  // the rows referring to the loaded ones, grouped by the row they refer to; a later change answering first wins
  const fetch = async (ids: string[], wanted: ColumnAggregation[]) => {
    const at = ++gen;
    const found = new Map<string, Record<string, unknown>>();
    const bySource = new Map<string, ColumnAggregation[]>();
    for (const a of wanted) {
      const {address, path} = linked.get(a.column)!;
      bySource.set(`${address}|${path}`, [...bySource.get(`${address}|${path}`) ?? [], a]);
    }
    for (const [key, list] of bySource) {
      const [address, path] = key.split('|');
      const rows = await grok.dapi.domains.table(address).queryDf({filter: {property: path, operator: '=', value: ids},
        expand: path.includes('.') ? [path.split('.')[0]] : [], limit: LIMIT});
      if (rows.rowCount === LIMIT)
        notify.warning(`Aggregations count the first ${LIMIT.toLocaleString()} rows of ${address} only`);
      const groups = rows.groupBy([path]);
      for (const a of list) {
        groups.add(a.aggregation as DG.AggregationType,
          a.aggregation === DG.AGG.TOTAL_COUNT ? null : linked.get(a.column)!.column, outputKey(a));
      }
      const measured = list.filter((a) => linked.get(a.column)!.measured && !COUNTS.includes(a.aggregation));
      if (measured.length > 0)
        groups.add(DG.AGG.UNIQUE_COUNT, 'unit', '#units').add(DG.AGG.FIRST, 'unit', '#unit');
      const result = groups.aggregate();
      for (let i = 0; i < result.rowCount; i++) {
        const row = found.get(String(result.get(path, i))) ?? {};
        for (const a of list) {
          const value = result.get(outputKey(a), i);
          row[outputKey(a)] = !measured.includes(a) ? value : result.get('#units', i) > 1 ? 'mixed units' :
            IN_UNIT.includes(a.aggregation) ? `${DomainDataTable.number(value)} ${result.get('#unit', i)}` : value;
        }
        found.set(String(result.get(path, i)), row);
      }
    }
    if (at === gen)
      values.value = found;
  };

  return component((self) => {
    const root = div([], 'u2-domain-output');
    const show = () => {
      if (built)
        return;
      built = true;
      self.task(async (t) => {
        const fields = await linkedFields();
        t.build(() => {
          const caption = (name: string) => DomainTable.captionOf(props.find((p) => p.name === name)!);
          const output = new ColumnsInput({label: 'Output', nullable: false, value: columns.peek().map(caption),
            table: frame(props.map((p) => [DomainTable.captionOf(p), p]))});
          self.effect(() => columns.value = (output.value ?? [])
            .map((c) => props.find((p) => DomainTable.captionOf(p) === c)!.name!));
          const form = new Form({layout: 'tall'}).add(output);
          if (fields.length > 0) {
            const picked = new AggregatedColumnsInput({label: 'Aggregations', table: frame(fields), value: []});
            self.effect(() => aggregations.value = picked.value ?? []);
            form.add(picked);
          }
          root.replaceChildren(form.root);
        });
        const ids = () => source.rows.items.value.filter((row: RowView) => !Rows.isDraft(row))
          .map((row) => String(row.id));
        self.effect(() => {
          if (shown.value)
            fetch(ids(), aggregations.value).catch((e) => DomainErrors.report(e));
        });
        const watching: Map<string, () => void> = new Map();
        const moved = () => {
          if (shown.peek())
            fetch(untracked(ids), aggregations.peek()).catch((e) => DomainErrors.report(e));
        };
        self.effect(() => {
          const tables = new Set(aggregations.value.flatMap((a) => linked.get(a.column)?.tables ?? []));
          for (const [address, unwatch] of [...watching]) {
            if (!tables.has(address)) {
              unwatch();
              watching.delete(address);
            }
          }
          for (const address of tables) {
            if (!watching.has(address))
              watching.set(address, source.watch(address, moved));
          }
        });
        self.own(() => {
          for (const unwatch of watching.values())
            unwatch();
        });
      });
    };
    return {root, columns, aggregations, values, show};
  });
}
