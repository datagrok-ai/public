/* Every sheet the domain stack paints with, in one import: an app built on `domains.*` (the
   three-line `table.app()` included) imports THIS and nothing else, and every control it renders
   — the app chrome, the form and its editors, the list and the grid, the filter box, the dialogs
   and balloons — is dressed. Side effects only; tokens first, since every other sheet reads them.
   A package rendering u2 controls of its own imports those sheets itself. */
import '../../../css/tokens.css';
import '../../../css/elements.css';
import '../../../css/buttons.css';
import '../../../css/inputs.css';
import '../../../css/number.css';
import '../../../css/date.css';
import '../../../css/choice.css';
import '../../../css/multi-select.css';
import '../../../css/radio.css';
import '../../../css/slider.css';
import '../../../css/color.css';
import '../../../css/font.css';
import '../../../css/icon-input.css';
import '../../../css/tags.css';
import '../../../css/typeahead.css';
import '../../../css/file.css';
import '../../../css/form.css';
import '../../../css/entity.css';
import '../../../css/list.css';
import '../../../css/tree.css';
import '../../../css/data-table.css';
import '../../../css/grid.css';
import '../../../css/menu.css';
import '../../../css/breadcrumbs.css';
import '../../../css/tabs.css';
import '../../../css/section.css';
import '../../../css/viewers.css';
import '../../../css/filter.css';
import '../../../css/filter-query.css';
import '../../../css/async.css';
import '../../../css/dialog.css';
import '../../../css/notify.css';
import '../../../css/tooltip.css';
import '../../../css/badge.css';
import '../../../css/domain.css';
