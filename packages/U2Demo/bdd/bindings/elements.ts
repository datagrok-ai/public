/* The MSA workbench vocabulary (src/pages/msa-workbench.ts) lives on a context: "user opens the
   MSA workbench" enters it, its names then resolve inside the workbench first, and platform names
   (toolbox, browse tab, grid …) keep their platform meaning. Most of the page needs no entry: the
   generic kinds resolve "run button in toolbar" or "sequence column input in MSA dialog" from the
   u2 contract alone. Entries exist for things worth a stable name, and for the override case. */
import {context} from '@datagrok-libraries/bdd';

export const workbench = context('MSA workbench', {selector: '[data-u2-name="msaWorkbench"]', aliases: ['workbench']});

workbench.element('alignment panel', {selector: '[data-u2-name="alignmentPanel"]', aliases: ['alignment form panel']});
workbench.element('results', {selector: '[data-u2-name="results"]', aliases: ['results panel', 'result panel']});
workbench.element('aligned sequences list', {selector: '[data-u2-name="alignedList"]', aliases: ['alignment list'], in: 'results'});
workbench.element('status line', {selector: '[data-u2-name="status"]', aliases: ['status'], in: 'results'});
workbench.element('log', {selector: '[data-u2-name="log"]', aliases: ['log panel'], in: 'results'});
workbench.element('settings panel', {selector: '[data-u2-name="settingsPanel"]'});

/* Override: composition would resolve "save button inside toolbar" as a button inside the toolbar
   (and would find it); this entry pins the phrase to one selector instead, which is what an author
   wants when the composed lookup is wrong or slow. The whole phrase wins before any splitting. */
workbench.element('save button inside toolbar', {selector: '[data-u2-name="toolbarSave"]',
  description: 'pinned as-is — never resolved by composition'});
