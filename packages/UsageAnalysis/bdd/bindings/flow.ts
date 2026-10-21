/* Flow's view as the features drive it: its start panel, its toolbox's node items and a Sketcher Input node's own
   editor, by their test ids (src/utils/test-ids.ts, `ff-…`); and the view itself as a widget the platform knows, whose
   status reports what its canvas did ("parameter edits": the node parameter edits reported since it opened). */
import {element} from '@datagrok-libraries/bdd';

element('flow editor widget', {selector: '.grok-view:has(.funcflow-canvas-container)',
  description: 'the Flow view, a widget the platform knows: its "parameter edits" reading counts the node parameter edits its canvas reported since it opened'});
element('flow blank canvas card', {selector: '[data-testid="ff-start-blank"]',
  description: 'the start panel\'s "Blank canvas" card: a new, empty flow'});
element('flow sketcher input item', {selector: '[data-testid="ff-browser-item-inputs-sketcher-input"]',
  description: 'Sketcher Input in the toolbox\'s Inputs section: a double click adds the node'});
element('flow sketcher preview', {selector: '[data-testid="ff-sketcher-compact"]',
  description: 'a Sketcher Input node\'s compact preview: a click expands it into the node\'s own sketcher'});
element('flow sketcher Done button', {selector: '[data-testid="ff-sketcher-done"]',
  description: 'Done above the node\'s sketcher: folds it back into the preview'});
