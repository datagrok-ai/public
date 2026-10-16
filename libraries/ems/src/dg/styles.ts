/* Every sheet the domain stack paints with, in one import: an app built on `domains.*` (the
   three-line `table.app()` included) imports THIS and nothing else, and every control it renders
   — the app chrome, the form and its editors, the list and the grid, the filter box, the dialogs
   and balloons — is dressed. Side effects only; tokens first, since every other sheet reads them.
   A package rendering u2 controls of its own imports those sheets itself. */
import '@datagrok-libraries/u2/css/tokens.css';
import '@datagrok-libraries/u2/css/elements.css';
import '@datagrok-libraries/u2/css/buttons.css';
import '@datagrok-libraries/u2/css/inputs.css';
import '@datagrok-libraries/u2/css/number.css';
import '@datagrok-libraries/u2/css/date.css';
import '@datagrok-libraries/u2/css/choice.css';
import '@datagrok-libraries/u2/css/multi-select.css';
import '@datagrok-libraries/u2/css/radio.css';
import '@datagrok-libraries/u2/css/slider.css';
import '@datagrok-libraries/u2/css/color.css';
import '@datagrok-libraries/u2/css/font.css';
import '@datagrok-libraries/u2/css/icon-input.css';
import '@datagrok-libraries/u2/css/tags.css';
import '@datagrok-libraries/u2/css/typeahead.css';
import '@datagrok-libraries/u2/css/file.css';
import '@datagrok-libraries/u2/css/form.css';
import '@datagrok-libraries/u2/css/entity.css';
import '@datagrok-libraries/u2/css/list.css';
import '@datagrok-libraries/u2/css/tree.css';
import '@datagrok-libraries/u2/css/data-table.css';
import '@datagrok-libraries/u2/css/grid.css';
import '@datagrok-libraries/u2/css/menu.css';
import '@datagrok-libraries/u2/css/breadcrumbs.css';
import '@datagrok-libraries/u2/css/tabs.css';
import '@datagrok-libraries/u2/css/section.css';
import '@datagrok-libraries/u2/css/viewers.css';
import '@datagrok-libraries/u2/css/filter.css';
import '@datagrok-libraries/u2/css/filter-query.css';
import '@datagrok-libraries/u2/css/async.css';
import '@datagrok-libraries/u2/css/dialog.css';
import '@datagrok-libraries/u2/css/notify.css';
import '@datagrok-libraries/u2/css/tooltip.css';
import '@datagrok-libraries/u2/css/badge.css';
import '@datagrok-libraries/u2/css/chips.css';
import '@datagrok-libraries/u2/css/access-grid.css';
import '@datagrok-libraries/u2/css/splitter.css';
import '../../css/manifest-editor.css';
import '@datagrok-libraries/u2/css/wizard.css';
import '../../css/binding-dialog.css';
import '../../css/domain.css';
