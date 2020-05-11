// Demonstrates handling of Grok-originated events

demog = grok.testData('demog', 5000);
view = grok.addTableView(demog);

function info(s) { grok.balloon.info(s); }

demog.onValuesChanged.subscribe((_) => info('ddt-values-changed'));
demog.onCurrentRowChanged.subscribe((_) => info('ddt-current-row-changed'));
demog.onMouseOverRowChanged.subscribe((_) => info('ddt-mouse-over-row-changed'));
demog.onCurrentColChanged.subscribe((_) => info('ddt-current-col-changed'));
demog.onMouseOverColChanged.subscribe((_) => info('ddt-mouse-over-col-changed'));
demog.onCurrentCellChanged.subscribe((_) => info('ddt-current-cell-changed'));
demog.onMouseOverRowGroupChanged.subscribe((_) => info('ddt-mouse-over-row-group-changed'));
demog.onNameChanged.subscribe((_) => info('ddt-table-name-changed'));
demog.onMetadataChanged.subscribe((_) => info('ddt-table-metadata-changed'));
demog.onColumnNameChanged.subscribe((_) => info('ddt-table-column-name-changed'));
demog.onColumnSelectionChanged.subscribe((_) => info('ddt-column-selection-changed'));
demog.onColumnsChanged.subscribe((_) => info('ddt-columns-changed'));
demog.onColumnsAdded.subscribe((_) => info('ddt-columns-added'));
demog.onColumnsRemoved.subscribe((_) => info('ddt-columns-removed'));
demog.onRowsRemoved.subscribe((_) => info('ddt-rows-removed'));
demog.onDataChanged.subscribe((_) => info('ddt-data-changed'));
demog.onFilterChanged.subscribe((_) => info('ddt-filter-changed'));
demog.onSelectionChanged.subscribe((_) => info('ddt-selection-changed'));