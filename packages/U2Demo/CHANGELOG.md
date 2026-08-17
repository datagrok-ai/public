# U2Demo changelog

## v.next

* Reports Browser: rich list rows per the list-item-rendering recipe — hover-revealed copy/open-ticket actions with a right-click menu of all actions, compact timestamps with full date-time on hover, and adaptive hiding of the reporter as the pane narrows

* Both apps open from the Browse tree on single click (app functions return the view) and ride the shell chrome via u2's appView: search and refresh in the view ribbon, live "N of M reports" in the per-view status bar

* Introduced the U2 Demo app: a tabbed tour of @datagrok-libraries/u2 (inputs, combobox/async, virtual lists and trees, layout, popups, form and property grid, platform bridge)
* Added the Entities tab: user/group pickers over dapiSource + ObjectHandler rendering, entity chips with handler tooltips and current-object wiring
* Added the Reports Browser app: a master-detail browser over grok.dapi.reports on AsyncPager/dapiPager — virtualized incremental list, smart filter, and a detail pane with reporter/assignee chips
* Added the Spaces & dashboards tab: space and project pickers with a live ObjectHandler renderCard preview (entityCard), chip, and Open action
* Added the Molecules tab: the platform sketcher bridged as a u2 input (moleculeInput), semType-driven editor selection in propertyForm via the new Editors registry, and structure chips/typeahead rows via moleculeRenderer
