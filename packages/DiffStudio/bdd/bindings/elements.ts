/* The names this package's features use. Platform names (toolbox, browse tab, context panel,
   console, status bar, open tableview, grid) are reserved; most of a model view needs no entry —
   "Fit ribbon item", "dose input" or "Multiaxis tab" resolve from the platform's kinds alone. */
import {element} from '@datagrok-libraries/bdd';

/* The hub's model cards are the library's `model card` kind (alias `hub card`), shared with the Tutorials project. */
element('library section', {selector: '[name="section-Library"]', description: 'the Library cards of the hub'});

/* The folder icon on the ribbon that opens the model menu — it has no name of its own, only the
   class the app puts on both of its ribbon dropdowns and the icon it is drawn with. */
element('open model button', {selector: '.diff-studio-ribbon-widget:has(.fa-folder-open)', aliases: ['open model icon']});

element('save to library icon', {selector: '.diff-studio-ribbon-save-to-model-catalog-icon',
  aliases: ['save to model hub icon']});

/* Browse has its own Refresh action outside the current view's ribbon. */
element('model hub refresh icon', {selector: '.d4-ribbon [name="icon-sync"]', aliases: ['catalog refresh icon']});
