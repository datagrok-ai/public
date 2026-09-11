/* The names this package's features use. Platform names (toolbox, browse tab, context panel,
   console, status bar, open tableview, grid) are reserved; most of a model view needs no entry —
   "Fit ribbon item", "dose input" or "Multiaxis tab" resolve from the platform's kinds alone. */
import {element, kind} from '@datagrok-libraries/bdd';

/* The hub's cards: a template or a library model, named by the label in the card's header. */
kind('hub card', {selector: '.diff-studio-hub-card', match: ['label'], labelSelector: '.diff-studio-hub-card-header'});

/* The folder icon on the ribbon that opens the model menu — it has no name of its own, only the
   class the app puts on both of its ribbon dropdowns and the icon it is drawn with. */
element('open model button', {selector: '.diff-studio-ribbon-widget:has(.fa-folder-open)', aliases: ['open model icon']});

element('save to library icon', {selector: '.diff-studio-ribbon-save-to-model-catalog-icon',
  aliases: ['save to model hub icon']});

/* The Refresh icon of the Model Hub ribbon carries no label and no tooltip — only the name of the
   icon it is drawn with. */
element('model hub refresh icon', {selector: '[name="icon-sync"]', aliases: ['catalog refresh icon']});
