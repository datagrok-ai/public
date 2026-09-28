/* The Tutorials app is a panel docked to the left of the shell, not a view: the tracks with their
   cards, and the running tutorial — its header, progress and the step entries, which the engine
   names (`role="checkbox"`, `aria-label` = the instruction, `aria-checked`, `aria-current="step"`,
   `aria-invalid` for a step that could not complete). Platform names (toolbox, grid, context panel)
   are reserved; these apply everywhere, since the panel sits beside whatever view is current. */
import {element, kind} from '@datagrok-libraries/bdd';

element('Tutorials panel', {selector: '.tutorials-root', aliases: ['tutorials app panel'],
  description: 'the dock panel the Tutorials app opens: the tracks, their cards and the running tutorial'});

kind('tutorial step', {
  aliases: ['tutorial steps'],
  selector: '.grok-tutorial-entry[role="checkbox"]',
  match: ['aria'],
  description: 'a step entry of the running tutorial, by its instruction exactly as shown; checked once done, invalid when it could not complete',
});

kind('tutorial card', {
  aliases: ['tutorial cards'],
  selector: '.tutorials-card[role="button"]',
  match: ['aria'],
  description: 'a tutorial of a track in the Tutorials panel, by the tutorial\'s name',
});

element('tutorial title', {selector: '.tutorials-root-header h1', aliases: ['running tutorial title'],
  description: 'the name of the tutorial that is running'});

element('tutorial progress', {selector: '.tutorials-root-progress [role="progressbar"]',
  description: 'the running tutorial\'s progress bar (aria-valuenow of aria-valuemax), beside its "Step: N of M" caption'});
