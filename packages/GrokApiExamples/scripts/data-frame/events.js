// Demonstrates handling of Grok-originated events

demog = grok.testData('demog', 5000);
view = grok.addTableView(demog);

function info(s) { grok.balloon.info(s); }

demog.onValuesChanged((_) => info('ddt-values-changed'));
demog.onCurrentRowChanged((_) => info('ddt-current-row-changed'));
demog.onMouseOverRowChanged((_) => info('ddt-mouse-over-row-changed'));
demog.onCurrentColChanged((_) => info('ddt-current-col-changed'));
demog.onMouseOverColChanged((_) => info('ddt-mouse-over-col-changed'));
demog.onCurrentCellChanged((_) => info('ddt-current-cell-changed'));
demog.onMouseOverRowGroupChanged((_) => info('ddt-mouse-over-row-group-changed'));
demog.onNameChanged((_) => info('ddt-table-name-changed'));
demog.onMetadataChanged((_) => info('ddt-table-metadata-changed'));
demog.onColumnNameChanged((_) => info('ddt-table-column-name-changed'));
demog.onColumnSelectionChanged((_) => info('ddt-column-selection-changed'));
demog.onColumnsChanged((_) => info('ddt-columns-changed'));
demog.onColumnsAdded((_) => info('ddt-columns-added'));
demog.onColumnsRemoved((_) => info('ddt-columns-removed'));
demog.onRowsRemoved((_) => info('ddt-rows-removed'));
demog.onDataChanged((_) => info('ddt-data-changed'));
demog.onFilterChanged((_) => info('ddt-filter-changed'));
demog.onSelectionChanged((_) => info('ddt-selection-changed'));