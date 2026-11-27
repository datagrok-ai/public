/* Nouns the Select and filter how-to clips need and the library does not name. The hierarchical filter card
   draws a tree, not a grid, so it reports no hit areas: its rows are DOM nodes, named here by the value they
   show, with the parts a person clicks (the label keeps that branch alone, the box beside it toggles the branch,
   the triangle opens the next level). A row is found among the visible ones, so "Caucasian" is the one under
   the branch that was opened. */
import {kind} from '@datagrok-libraries/bdd';

kind('hierarchy node', {
  aliases: ['hierarchy row'],
  selector: '[name="viewer-Filters"] .d4-tree-view-node:has(> .d4-hierarchical-filter-caption-label)',
  match: ['label'],
  labelSelector: '.d4-hierarchical-filter-caption-value',
  parts: {
    label: ':scope > .d4-hierarchical-filter-caption-label',
    checkbox: ':scope > .d4-hierarchical-filter-checkbox-container',
    twistie: ':scope > .d4-tree-view-tri',
  },
  description: 'a row of a hierarchical filter card by the value it shows ("F" hierarchy node): its label keeps that ' +
    'branch alone, its checkbox toggles the branch, its twistie opens the next level',
});
