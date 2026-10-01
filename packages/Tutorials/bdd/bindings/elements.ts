/* The Tutorials app is a panel docked to the left of the shell, not a view: the tracks with their
   cards, and the running tutorial — its header, progress and the step entries, which the engine
   names (`role="checkbox"`, `aria-label` = the instruction, `aria-checked`, `aria-current="step"`,
   `aria-invalid` for a step that could not complete). Platform names (toolbox, grid, context panel)
   are reserved; these apply everywhere, since the panel sits beside whatever view is current. */
import {element, kind} from '@datagrok-libraries/bdd';

element('Tutorials panel', {selector: '.tutorials-root', aliases: ['tutorials app panel'],
  description: 'the dock panel the Tutorials app opens: the tracks, their cards and the running tutorial'});

kind('tutorial card', {
  aliases: ['tutorial cards'],
  selector: '.tutorials-card[role="button"]',
  match: ['aria'],
  description: 'a tutorial of a track in the Tutorials panel, by the tutorial\'s name',
});

element('tutorial title', {selector: '.tutorials-root-header h1', aliases: ['running tutorial title'],
  description: 'the name of the tutorial that is running'});

kind('cliff molecule', {
  aliases: ['cliff molecules'],
  selector: '.chem-activity-cliffs-molecule[role="button"]',
  match: ['aria'],
  description: 'a molecule of the pair the Activity Cliffs pane shows (Chem), named "molecule of row N"; a click makes its row current',
});

kind('tutorial track', {
  aliases: ['tutorial tracks'],
  selector: '.tutorials-track[data-name]',
  match: ['label'],
  labelSelector: '.tutorials-track-title h1',
  description: 'a track of the Tutorials panel (its title, progress and cards), by the track\'s name',
});
