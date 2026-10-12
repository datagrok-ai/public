import * as DG from 'datagrok-api/dg';

export type ResultKind = 'binary' | 'multiclass' | 'regression';
export type PanelKey = 'predictedVsActual' | 'rocCurve' | 'metrics' | 'confusions' | 'residuals' | 'distributions' |
  'insights' | 'wrongPredictions' | 'statistics' | 'features' | 'correlation' | 'classBalance' | `visualize:${number}`;

export interface ResultsRow {
  height: number;
  panels: PanelKey[];
  /** The panels' initial shares of the row's width, one per panel. */
  weights: number[];
}

/** The initial heights of rows 1-3, in px. */
const ROW_HEIGHTS = [350, 300, 450] as const;
const VISUALIZE_ROW_HEIGHT = 350;
const VISUALIZE_PER_ROW = 3;
/** The weight of **Class balance** when it joins row 2. */
const CLASS_BALANCE_ROW_2_WEIGHT = 1.5;

type Template = [PanelKey, number][];

/** The kind of the results of a model of [target]: a numerical target's, or a classifier's by its number of classes. */
export function resultKindOf(target: DG.Column): ResultKind {
  return target.matches(DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME) ? 'regression' :
    classesOf(target).length === 2 ? 'binary' : 'multiclass';
}

/** The classes of [target], its empty category aside. */
export function classesOf(target: DG.Column): string[] {
  return target.categories.filter((c) => c !== '');
}

const TEMPLATES: {[kind in ResultKind]: Template[]} = {
  binary: [[['predictedVsActual', 2], ['rocCurve', 2], ['metrics', 1]],
    [['insights', 1], ['confusions', 1.5], ['wrongPredictions', 2.5]],
    [['features', 1], ['correlation', 1], ['classBalance', 1]]],
  multiclass: [[['predictedVsActual', 1], ['confusions', 1]],
    [['insights', 1], ['metrics', 1.5], ['wrongPredictions', 2.5]],
    [['features', 1], ['correlation', 1], ['classBalance', 1]]],
  regression: [[['predictedVsActual', 1], ['residuals', 1], ['distributions', 1]],
    [['insights', 1], ['metrics', 1.5], ['statistics', 1.5]],
    [['features', 1], ['correlation', 1]]],
};

/** The rows of panels for [kind]: its template rows filtered to [available] (`metrics` is always available), an empty
 * row dropped, then the [visualizeCount] method panels in rows of at most {@link VISUALIZE_PER_ROW}. When row 2 is left
 * with fewer than three panels, **Class balance** moves from row 3 to its right end. */
export function resultsLayout(kind: ResultKind, available: ReadonlySet<PanelKey>, visualizeCount: number):
  ResultsRow[] {
  const templates = TEMPLATES[kind].map((template) =>
    template.filter(([key]) => key === 'metrics' || available.has(key)));
  const [, second, third] = templates;
  const balance = third.findIndex(([key]) => key === 'classBalance');
  if (balance >= 0 && second.length < 3) {
    third.splice(balance, 1);
    second.push(['classBalance', CLASS_BALANCE_ROW_2_WEIGHT]);
  }
  const rows: ResultsRow[] = [];
  for (const [i, template] of templates.entries()) {
    if (template.length > 0)
      rows.push({height: ROW_HEIGHTS[i], panels: template.map(([key]) => key), weights: template.map(([, w]) => w)});
  }
  for (let start = 0; start < visualizeCount; start += VISUALIZE_PER_ROW) {
    const panels: PanelKey[] = [];
    for (let i = start; i < Math.min(start + VISUALIZE_PER_ROW, visualizeCount); i++)
      panels.push(`visualize:${i}`);
    rows.push({height: VISUALIZE_ROW_HEIGHT, panels, weights: panels.map(() => 1)});
  }
  return rows;
}
