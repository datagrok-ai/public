import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// mapping from internal Dart labels to JS DG.VIEWER labels
export const viewerTypesMapping: {[key: string]: string} = {
  ['barchart']: DG.VIEWER.BAR_CHART,
  ['boxplot']: DG.VIEWER.BOX_PLOT,
  ['calendar']: DG.VIEWER.CALENDAR,
  ['corrplot']: DG.VIEWER.CORR_PLOT,
  ['densityplot']: DG.VIEWER.DENSITY_PLOT,
  ['filters']: DG.VIEWER.FILTERS,
  ['form']: DG.VIEWER.FORM,
  ['globe']: DG.VIEWER.GLOBE,
  ['googlemap']: DG.VIEWER.GOOGLE_MAP,
  ['grid']: DG.VIEWER.GRID,
  ['heatmap']: DG.VIEWER.HEAT_MAP,
  ['histogram']: DG.VIEWER.HISTOGRAM,
  ['linechart']: DG.VIEWER.LINE_CHART,
  ['markup']: DG.VIEWER.MARKUP,
  ['matrixplot']: DG.VIEWER.MATRIX_PLOT,
  ['networkdiagram']: DG.VIEWER.NETWORK_DIAGRAM,
  ['pcplot']: DG.VIEWER.PC_PLOT,
  ['piechart']: DG.VIEWER.PIE_CHART,
  ['scatterplot']: DG.VIEWER.SCATTER_PLOT,
  ['3dscatterplot']: DG.VIEWER.SCATTER_PLOT_3D,
  ['shapemap']: DG.VIEWER.SHAPE_MAP,
  ['statistics']: DG.VIEWER.STATISTICS,
  ['tileviewer']: DG.VIEWER.TILE_VIEWER,
  ['timelines']: DG.VIEWER.TIMELINES,
  ['treemap']: DG.VIEWER.TREE_MAP,
  ['trellisplot']: DG.VIEWER.TRELLIS_PLOT,
  ['wordcloud']: DG.VIEWER.WORD_CLOUD,
} as const;

export const CARD_VIEW_TYPE = 'JsCardView' as const;
export const SCRIPTS_VIEW_TYPE = 'scripts' as const;
export const FUNCTIONS_VIEW_TYPE = 'functions' as const;
export const VIEWER_PATH = 'viewer' as const;
export const RESTRICTED_PATH = 'restrictedValues' as const;
export const EDIT_STATE_PATH = 'editState' as const;
export enum DIRECTION {
  INPUT = 'Input',
  OUTPUT = 'Output'
};

export const EXPERIMENTAL_TAG = 'experimental';
export const RUN_NAME_COL_LABEL = 'Run name' as const;
export const RUN_ID_COL_LABEL = 'RunId' as const;
export enum VISIBILITY_STATE {
  HIDDEN = 'hidden',
  VISIBLE = 'visible',
}

export enum ABILITY_STATE {
  ENABLED = 'enabled',
  DISABLED = 'disabled',
}

export type INPUT_STATE = 'disabled' | 'restricted' | 'restricted unlocked' | 'inconsistent' | 'user input';

export type VIEW_STATE = 'inconsistent' | 'consistent';

export const storageName = `ModelStorage`;

export const EXP_COLUMN_NAME = 'Source';
export const FAVORITE_COLUMN_NAME = 'Is favorite';
export const COMPLETE_COLUMN_NAME = 'Is complete';
export const ACTIONS_COLUMN_NAME = 'Delete';
export const AUTHOR_COLUMN_NAME = 'Author';
export const STARTED_COLUMN_NAME = 'Date';
export const TITLE_COLUMN_NAME = 'Title';
export const DESC_COLUMN_NAME = 'Desc.';
export const TAGS_COLUMN_NAME = 'Tags';
