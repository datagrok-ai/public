import * as DG from 'datagrok-api/dg';


export type Indexable = { [key: string]: any };

export const TAGS = ['Math', 'Logic', 'Stats', 'Date'];

export const DEFAULT_PARAM_VALUES: Indexable = {
  bool: true,
  int: 1,
  num: 0.5,
  string: '1',
  dynamic: 1,
};

export const SPECIAL_CASES: Indexable = {
  DateParse: { s: '2020-01-01' },
  RandBetween: { a: 0, b: 1 },
  TimeParse: { s: '10:05' },
};

export function findFuncsByTags(tags: string[]): DG.Func[] {
  const funcs = new Set<DG.Func>();
  tags.forEach((tag) => DG.Func.find({ tags: [tag] }).forEach((f) => funcs.add(f)));
  return [...funcs];
}
