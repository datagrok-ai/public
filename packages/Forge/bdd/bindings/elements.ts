/* Forge's own names. Platform names (toolbox, context panel, grid …) are reserved; the inputs, buttons, dialogs and
   viewers of the Train view resolve from the platform's kinds ("Target input", "Save button", "grid viewer in
   "Performance" results panel"). What is here has no platform kind: Forge's collapsible groups, the panels of the
   results, the lines of the table pane, and a few elements with a stable class. The Train view publishes its
   readiness on `results area` (aria-busy / aria-invalid, ARCHITECTURE.md "Automation surface"). */
import {element, kind} from '@datagrok-libraries/bdd';

/* A collapsible group of the Train view or the Apply dialog ("Data group", "Columns group"): aria-expanded on its
   header, aria-invalid on its summary (`collapsible-group.ts`). */
kind('group', {
  aliases: ['collapsible group'],
  selector: '.forge-group',
  match: ['label'],
  labelSelector: '.forge-group-header label',
  parts: {header: '.forge-group-header', summary: '.forge-group-summary', body: '.forge-group-body'},
  description: 'a Forge collapsible group by its caption',
});

/* A panel of the results by its title ("Performance" results panel, "Runs" results panel). */
kind('results panel', {
  selector: '.forge-results-panel, .forge-runs',
  match: ['label'],
  labelSelector: ':scope > .forge-results-title',
  parts: {title: ':scope > .forge-results-title'},
  description: 'a panel of the Train or Evaluation view\'s results by its title',
});

/* The ROC Curve panel's title carries the AUC ("ROC Curve (AUC = 0.906)"), which the kind's whole-title match does not
   read as "ROC Curve". */
element('ROC Curve panel', {selector: '.forge-results-panel:has(> .forge-results-title:has-text("ROC Curve"))',
  description: 'the ROC Curve panel of the results, whatever AUC its title shows'});
element('results area', {selector: '.forge-train-results',
  description: 'the results of the Train view: aria-busy while it checks, trains or shows, aria-invalid for a problem'});
element('train form', {selector: '.forge-train-form', description: 'the inputs of the Train view'});
element('form splitter', {selector: '.ui-split-h:has(.forge-train-form) > .ui-split-h-divider',
  description: 'the divider between the inputs and the results of the Train view'});
element('seed copy icon', {selector: '.forge-results-metrics i[aria-label="Copy the seed"]',
  description: 'the copy icon after the seed in the Performance panel'});
/* LIB-PROMOTE: the "Select columns..." dialog's Search box is a bare input with no host or label; its grid reports
   only the rows it draws, so a column below them is searched for first. */
element('column search', {selector: 'input.d4-search-input[placeholder="Search"]',
  description: 'the Search box of the "Select columns..." dialog'});
element('tags box', {selector: 'input.d4-tags-selector-input',
  description: 'the text box of the Tags input, where a tag is typed and confirmed with Enter'});
/* LIB-PROMOTE: the platform's tags input shows a chip per tag in a container of its own; there is no kind for them. */
element('tag chips', {selector: '.d4-tags-selector-tags-container > *',
  description: 'the chips of a Tags input, one per tag'});

element('forge app', {selector: '.forge-app', description: 'the Forge view: Methods and the model catalog'});

/* LIB-PROMOTE: a line of the table pane Forge models (`table-models-panel.ts`) has the shape of the built-in Models
   pane's rows (`d4-entity-markup-row`): the run icon, the name as a link, the menu arrow of `ui.contextActions`. */
kind('model line', {
  selector: '.forge-model-line',
  match: ['label'],
  labelSelector: '.forge-model-name',
  parts: {'run icon': '.svg-model-run', name: '.forge-model-name', 'menu arrow': '.grok-context-arrow-down'},
  description: 'a model\'s line in the table pane Forge models by the model\'s name; its run icon, name and menu arrow',
});

/* LIB-PROMOTE: every entity's Details (and a package's Sharing) pane, and a model's tooltip, is a ui.tableFromMap, a
   <table> of caption and value rows; the library's "table row" kind has no value part. */
kind('detail', {
  aliases: ['detail row'],
  selector: '.d4-info-table tr',
  match: ['label'],
  labelSelector: ':scope > td:first-child',
  parts: {value: ':scope > td:last-child'},
  description: 'a row of a key-value table of the context panel or a tooltip (Details, Sharing, a model\'s tooltip) by ' +
    'its caption; its value part',
});