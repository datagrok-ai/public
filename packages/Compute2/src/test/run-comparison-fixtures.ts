import {ComparisonEntryNodes} from '../components/RunComparison/types';

export function makeEntry(
  entryId: string,
  scalars: {name: string, valueType?: string, units?: string, value?: number | null}[] = [],
  tables: {
    path: string, name?: string, nqName?: string, defaultIndexColumn?: string,
    columns: {name: string, type?: string, units?: string}[],
  }[] = [],
): ComparisonEntryNodes {
  return {
    entryId,
    entryName: entryId,
    scalars: scalars.map((s, i) => ({
      path: `io${i}`,
      name: s.name,
      valueType: s.valueType ?? 'double',
      units: s.units,
      value: s.value ?? 0,
    })),
    tables: tables.map((t) => ({
      path: t.path,
      name: t.name ?? t.path,
      nqName: t.nqName,
      defaultIndexColumn: t.defaultIndexColumn,
      columns: t.columns.map((c) => ({name: c.name, type: c.type ?? 'double', units: c.units})),
      rowCount: 10,
    })),
  };
}

export const indexMap = (spec: Record<string, Record<string, string>>) =>
  new Map(Object.entries(spec).map(([entryId, tables]) => [entryId, new Map(Object.entries(tables))]));
