/* The VlaaiVis editor in the Context Panel (`src/vlaaivis/editor.ts`): one accordion pane per sector,
   each holding a row per column (`createMpoRow` of the statistics library), and the pane of the
   columns no sector holds. Neither the panes nor the rows carry a name of their own beyond the
   sector's title, so they are reached by what a person reads: the sector's name, the column's name. */
import {kind} from '@datagrok-libraries/bdd';

kind('VlaaiVis sector', {
  selector: '.d4-accordion-pane:is(:has(> .d4-accordion-pane-header .power-grid-vlaaivis-swatch), [aria-label="Unassigned columns"])',
  match: ['aria'],
  parts: {color: '.power-grid-vlaaivis-swatch', header: '.d4-accordion-pane-header'},
  description: 'a sector of the VlaaiVis editor by its name ("Unassigned columns" is the pane of the columns no sector ' +
    'holds); "color of X VlaaiVis sector" is the swatch that opens its color picker',
});

kind('VlaaiVis column', {
  selector: '.statistics-mpo-row:has(.power-grid-vlaaivis-name)',
  match: ['label'],
  labelSelector: '.power-grid-vlaaivis-name',
  editorSelector: '.statistics-mpo-weight-cell input.ui-input-editor',
  parts: {name: '.power-grid-vlaaivis-name', weight: '.statistics-mpo-weight-cell .ui-input-root',
    curve: '.statistics-mpo-line-editor'},
  description: 'a column\'s row in a VlaaiVis sector by the column name: typing goes into its weight; "name of X ' +
    'VlaaiVis column" is the label dragged to another sector, "Settings icon in X VlaaiVis column" opens its score function',
});
