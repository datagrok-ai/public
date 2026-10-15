/* Forge's own names. Platform names (toolbox, context panel, grid …) are reserved; the inputs, buttons, dialogs and
   viewers of the Train view resolve from the platform's kinds ("Target input", "Save button", "grid viewer in
   "Performance" results panel"). What is here has no platform kind: Forge's collapsible groups, the panels of the
   results, and a few elements with a stable class. The Train view publishes its readiness on `results area`
   (aria-busy / aria-invalid, ARCHITECTURE.md "Automation surface"). */
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