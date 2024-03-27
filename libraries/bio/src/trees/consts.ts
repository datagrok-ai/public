import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const NEWICK_EMPTY: string = ';';

export enum TAGS {
  NEWICK = '.newick',
  NEWICK_JSON = '.newick_json',
  OBJECTS = '.newick_objects',
  NEWICK_LEAF_COL_NAME = '.newickLeafColumn',
}

export enum TreeColorNames {
  Main = 'Main',
  Light = 'Light',
  Current = 'Current',
  MouseOver = 'MouseOver',
  Selection = 'Selection',
}

export enum DistanceMetric {
  Euclidean = 'euclidean',
  Manhattan = 'manhattan',
}

export enum LinkageMethod {
  Single = 'single',
  Complete = 'complete',
  Average = 'average',
  Weighted = 'weighted',
  Centroid = 'centroid',
  Median = 'median',
  Ward = 'ward',
}


export const TreeDefaultPalette: { [name: string]: number } = {
  [TreeColorNames.Main]: DG.Color.categoricalPalette[12],
  [TreeColorNames.Light]: DG.Color.categoricalPalette[13],
  [TreeColorNames.Current]: DG.Color.currentRow,
  [TreeColorNames.MouseOver]: DG.Color.mouseOverRows,
  [TreeColorNames.Selection]: DG.Color.selectedRows,
};
