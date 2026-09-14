/* The names this package's features use. Platform names (toolbox, browse tab, context panel,
   console, status bar, open tableview, grid) are reserved; the app's own names live on a context,
   and apply after a step declared with `enters` (see steps.ts). Most of a u2 page needs no entry:
   "run button in toolbar" or "name input in dialog" resolve from the u2 contract alone. */
import {context} from '@datagrok-libraries/bdd';

export const app = context('Peptides app', {selector: '[data-u2-name="Peptides"]'});
// app.element('results', {selector: '[data-u2-name="results"]', aliases: ['results panel']});
