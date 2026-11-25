/* Flow's view as the features drive it: its start panel, its ribbon, its toolbox's node items and a Sketcher Input
   node's own editor, by their test ids (src/utils/test-ids.ts, `ff-…`); and the view itself as a widget the platform
   knows, whose status reports what its canvas did. */
import {element, kind} from '@datagrok-libraries/bdd';

element('flow editor widget', {selector: '.grok-view .funcflow-canvas-container',
  description: 'the canvas of the Flow view, a widget the platform knows (the view itself also holds the grid of a preview). Its readings: "parameter edits" (the node parameter edits its canvas reported since it opened), "nodes", "computed tables" (the output tabs that have a table), "unsaved changes" (whether closing it would ask)'});
kind('recent flow', {
  selector: '.funcflow-start-recent-item',
  match: ['text'],
  description: 'a saved flow among the start panel\'s "Recent flows", by its name; a click opens it',
});
element('flow run button', {selector: '[data-testid="ff-ribbon-run"]', description: 'Run on the Flow ribbon'});
element('flow save button', {selector: '[data-testid="ff-ribbon-save"]',
  description: 'Save on the Flow ribbon: opens the Save Flow dialog'});
element('flow blank canvas card', {selector: '[data-testid="ff-start-blank"]',
  description: 'the start panel\'s "Blank canvas" card: a new, empty flow'});
element('flow sketcher input item', {selector: '[data-testid="ff-browser-item-inputs-sketcher-input"]',
  description: 'Sketcher Input in the toolbox\'s Inputs section: a double click adds the node'});
element('flow sketcher preview', {selector: '[data-testid="ff-sketcher-compact"]',
  description: 'a Sketcher Input node\'s compact preview: a click expands it into the node\'s own sketcher'});
element('flow sketcher Done button', {selector: '[data-testid="ff-sketcher-done"]',
  description: 'Done above the node\'s sketcher: folds it back into the preview'});
