/* Chem's own screen parts: the calculators listed on the left of the Chemical Properties dialog
   (a checkbox and a name per registered calculator function), the MPO Profiles app — the
   editable title and description of a profile tab, its property rows and the rows of the list —
   the 3D view of a molecule, and the parts of the SAR Matrix Summary. */
import {element, kind} from '@datagrok-libraries/bdd';

kind('calculator', {
  selector: '.biochem-calc-nav-item',
  match: ['text'],
  editorSelector: 'input[type="checkbox"]',
  description: 'a calculator row of the Chemical Properties dialog, by its name ("Chemical Properties", "logP")',
});

kind('reaction', {
  aliases: ['reaction card'],
  selector: '.d4-dialog-contents .d4-flex-col:has(> canvas + label)',
  match: ['label'],
  labelSelector: ':scope > label:first-of-type',
  description: 'a reaction card of the Run Reaction and Two-Component Reaction dialogs, by its name ("Amide Coupling")',
});

element('MPO profile title', {selector: '.chem-profile-header',
  description: 'the editable name at the top of a profile tab; "Untitled Profile" until one is typed'});
element('MPO profile description', {selector: '.chem-profile-description',
  description: 'the editable description under the profile title'});

kind('MPO property', {
  selector: '.statistics-mpo-row',
  match: ['text'],
  editorSelector: '.statistics-mpo-property-cell input',
  description: 'a property row of the profile editor; its editor is the property-name field ("first MPO property")',
});

kind('MPO profile', {
  selector: '.chem-mpo-profiles-table tr:has(.chem-mpo-actions-button)',
  match: ['label'],
  labelSelector: ':scope > td:nth-child(2)',
  parts: {actions: '.chem-mpo-actions-button'},
  description: 'a row of the Manage Profiles list, by the profile name; "actions of X MPO profile" is its ⋮ button',
});

element('R-Groups settings icon', {selector: '.chem-rgroup-settings-icon',
  description: 'the gear of the R-Groups Analysis dialog: shows Matching strategy and Only match at R groups, which the dialog remembers per account'});

element('Key series card', {selector: '.chem-sar-sum-overview',
  description: 'the card of the SAR Matrix Summary with the series to open first, the fold-level chips in its title row'});
element('key series row', {selector: '.chem-sar-sum-overview .chem-sar-sum-list > .chem-sar-card',
  parts: {'open icon': '.fa-th', name: '.chem-sar-card-name'},
  description: 'a series of the Key series card; a click opens its best core under it ("first key series row"), "open icon of …" opens the series on the SAR Matrix tab at its best measured compound, "name of …" is the series name'});
element('key series detail', {selector: '.chem-sar-sum-overview .chem-sar-sum-expand',
  description: 'what an opened key series row shows under it: its best core, or its best R-group combination'});
element('Worth making card', {selector: '.chem-sar-sum-making-block',
  description: 'the card of the SAR Matrix Summary with the predicted analogs nobody has made, every one listed'});
element('analog worth making row', {selector: '.chem-sar-sum-making-block .chem-sar-card',
  parts: {cart: '.chem-sar-cart-icon', 'open icon': '.fa-th', name: '.chem-sar-card-name'},
  description: 'an analog of the Worth making card; a click opens how it was predicted under it, "cart of …" adds it to the Make list, "open icon of …" opens its cell on the SAR Matrix tab'});
element('analog worth making detail', {selector: '.chem-sar-sum-making-block .chem-sar-sum-expand',
  description: 'what an opened analog row shows under it'});
element('Never assayed card', {selector: '.chem-sar-sum-making > .chem-sar-sum-card:not(.chem-sar-sum-making-block)',
  description: 'the card of the SAR Matrix Summary above Worth making: compounds of the table with no activity value, ranked by their prediction; absent while there are none'});
element('never assayed row', {selector: '.chem-sar-sum-making > .chem-sar-sum-card:not(.chem-sar-sum-making-block) .chem-sar-card',
  description: 'a compound of the Never assayed card ("first never assayed row"): its id, its support ("n=3"), its place and its predicted value'});
element('Queue for testing button', {selector: '.chem-sar-sum-making > .chem-sar-sum-card:not(.chem-sar-sum-making-block) .chem-sar-sum-chips button',
  description: 'the button under the Never assayed card that puts every compound of it on the Make list; its caption counts them ("Queue these 3 for testing")'});
element('Add all to Make list button', {selector: '.chem-sar-sum-making-block .chem-sar-sum-chips button',
  description: 'the button under the Worth making list that collects every analog of it; its caption counts them ("Add all 24 to Make list")'});
element('Summary evidence', {selector: '.chem-sar-sum-detail:not(:empty)',
  description: 'the evidence an answer tile, or the best value or swap of a component, opens in place under it'});
element('evidence row', {selector: '.chem-sar-sum-detail .chem-sar-card',
  parts: {'open icon': '.fa-th', 'count badge': '.chem-sar-sum-click', value: '.chem-sar-card-title', name: '.chem-sar-card-name'},
  description: 'a row of the evidence open on the Summary ("first evidence row"); "open icon of …" opens it on the SAR Matrix tab, "count badge of …" ("12 cpd") selects the compounds it counts, "value of …" is the number on its right (an effect, or a swap\'s weakest pair)'});
element('setup line', {selector: '.chem-sar-sum-setup',
  description: 'the line of the SAR Matrix Summary saying how the matrices were built: which column runs across them'});
element('What to change title', {selector: '.chem-sar-sum-comp > .chem-sar-sum-title-row > .chem-sar-cp-section-title',
  description: 'the title of the What to change card when the components are columns, with their count ("What to change · 3 components")'});
element('more components fold', {selector: '.chem-sar-sum-comp > .chem-sar-sum-fold',
  description: 'the "N more" fold under the first rows of the What to change card; a click lists the other components'});
element('SAR Matrix note', {selector: '.chem-sar-empty-note',
  description: 'the note a tab of the SAR Matrix viewer shows in place of its content: why nothing was found, or what to do first'});
element('transfer title', {selector: '.chem-sar-xfer-panel .chem-sar-main-bar > .chem-sar-main-title',
  description: 'the title over the SAR Transfer pane: the source core, an arrow, and the target it is read against'});
element('Make list grid viewer', {selector: '.chem-sar-makelist-grid [name="viewer-Grid"]',
  description: 'the grid of the Make list tab, one collected compound per row'});

kind('level chip', {
  selector: '.chem-sar-sum-level',
  match: ['item'],
  description: 'a fold-level chip in the title row of the Key series card, by its level ("L2 level chip", "All level chip"); a click reads the whole Summary at that level, and the chip shows as selected',
});

kind('answer tile', {
  selector: '.chem-sar-sum-answer',
  match: ['label'],
  labelSelector: '.chem-sar-cp-section-title',
  description: 'an answer of the SAR Matrix Summary when the series are cut from the molecules, by its title ("What improves potency answer tile")',
});

kind('series card', {
  selector: '.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card',
  match: ['label'],
  labelSelector: '.chem-sar-card-name',
  parts: {name: '.chem-sar-card-name'},
  description: 'a card of the series list on the SAR Matrix tab, by its name ("Series 191 series card") or its place ("second series card", "name of second series card"); a click opens that series',
});

element('series filter icon', {selector: '.chem-sar-nav:not(.chem-sar-xfer-nav) .chem-sar-nav-controls .chem-sar-struct-icon',
  description: 'the funnel over the series list of the SAR Matrix tab; a click opens the filters of the series (Series, Core, Best, …) in a popup'});
element('series match count', {selector: '.chem-sar-nav:not(.chem-sar-xfer-nav) .chem-sar-nav-matches',
  description: 'the line under the series list\'s controls while a filter is on ("1 of 134 match"); empty without one'});
element('matrix cart', {selector: '.chem-sar-matrix .chem-sar-main-bar .chem-sar-cart-icon',
  description: 'the cart in the bar over the matrix (aria-label "Add to Make list"); a click adds the cell clicked last to the Make list'});
element('transfer filter icon', {selector: '.chem-sar-xfer-nav .chem-sar-filter-icon',
  description: 'the funnel over the transfer list of the SAR Transfer tab; a click opens the filters of the transfers (Source, Target, their cores) in a popup'});

element('Make list title', {selector: '.chem-sar-main-bar:has(+ .chem-sar-makelist-grid) > .chem-sar-main-title',
  description: 'the count over the collected compounds of the Make list tab ("1 compound")'});

kind('component row', {
  selector: '.chem-sar-sum-comp-row:has(> .chem-sar-sum-col)',
  match: ['label'],
  labelSelector: '.chem-sar-sum-col .chem-sar-card-title',
  parts: {value: '.chem-sar-sum-half:not(.chem-sar-sum-swap)', swap: '.chem-sar-sum-swap'},
  description: 'a row of the What to change card of the SAR Matrix Summary when the components are columns, by its name; "value of Warhead component row" is its best value, "swap of …" its best measured swap',
});
