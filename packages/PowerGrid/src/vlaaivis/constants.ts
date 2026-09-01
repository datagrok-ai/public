export const VLAAIVIS_METADATA_TAG = '.vlaaivis-metadata';

export const PANE_HEADER_SELECTOR = '.d4-accordion-pane-header';

export const DEFAULTS = {
  WEIGHT: 1,
  LOWER_BOUND: 0.8,
  UPPER_BOUND: 0.9,
  /// Columns put into their own sector when the editor opens on an unconfigured chart.
  AUTO_GROUP_COLUMNS: 3,
  PLOT_HEIGHT: 64,
};

export const LABELS = {
  UNASSIGNED: 'Unassigned columns',
  TIP: 'Drag a column name onto a sector to move it. Unassigned columns aren\'t drawn.',
  NO_UNASSIGNED: 'Every column is in a sector. Drag one here to leave it out.',
  LOWER_BOUND: 'Target min',
  UPPER_BOUND: 'Target max',
};

export const TOOLTIPS = {
  SECTOR: 'A group of columns sharing one color. Its share of the circle is the sum of its weights.',
  WEIGHT: 'How wide this column\'s wedge is. Relative — weights don\'t need to sum to 1.',
  LOWER_BOUND: 'Inner edge of the shaded ring drawn behind the wedges.',
  UPPER_BOUND: 'Outer edge of the shaded ring. Wedges reaching past it are meeting the target.',
  NEW_SECTOR: 'New sector',
  UNASSIGNED: 'Columns not in any sector aren\'t drawn. Drag one into a sector to include it.',
  DRAG: 'Drag into a sector to move it',
  COLOR: 'Sector color',
  RENAME: 'Rename sector',
  DELETE: 'Delete sector',
};
