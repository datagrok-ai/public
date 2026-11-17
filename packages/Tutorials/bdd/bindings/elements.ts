/* The Tutorials app is a view: the tracks with their cards, a card starts its tutorial. A running
   tutorial shows its current step in a bar on top and all its steps, on demand, in a panel docked to
   the left — its header, progress and the step entries, which the engine names (`role="checkbox"`,
   `aria-label` = the instruction, `aria-checked`, `aria-current="step"`, `aria-invalid` for a step
   that could not complete). Platform names (toolbox, grid, context panel) are reserved; the bar and
   the panel apply everywhere, since they sit beside whatever view is current. */
import {element, kind} from '@datagrok-libraries/bdd';

element('Tutorials view', {selector: '.tutorials-view', aliases: ['tutorials app view'],
  description: 'the view the Tutorials app opens: the tracks with their progress and cards'});

element('Tutorials panel', {selector: '.tutorials-root', aliases: ['tutorial steps panel'],
  description: 'the panel with all steps of the running tutorial, docked to the left on demand'});

kind('tutorial card', {
  aliases: ['tutorial cards'],
  selector: '.tutorials-card[role="button"]',
  match: ['aria'],
  description: 'a tutorial of a track in the Tutorials view, by the tutorial\'s name',
});

element('tutorial title', {selector: '.tutorials-bar-title', aliases: ['running tutorial title'],
  description: 'the name of the tutorial that is running'});

element('tutorial bar', {selector: '.tutorials-bar', aliases: ['running tutorial bar'],
  description: 'the bar a running tutorial shows on top'});

element('all steps icon', {selector: '.tutorials-bar .fa-list', aliases: ['tutorial step list icon'],
  description: 'opens or closes the step list'});

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
  description: 'a track of the Tutorials view (its title, progress and cards), by the track\'s name',
});
