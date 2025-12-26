// Demonstrates handling of Datagrok-originated events

demog = grok.data.demo.demog();
view = grok.shell.addTableView(demog);

function info(s) {
  grok.shell.info(s);
}

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
demog.onRowsAdded.subscribe((_) => info('ddt-rows-added'));
demog.onRowsRemoved.subscribe((_) => info('ddt-rows-removed'));
demog.onRowsFiltering.subscribe((_) => info('ddt-rows-filtering'));
demog.onRowsFiltered.subscribe((_) => info('ddt-rows-filtered'));
demog.onDataChanged.subscribe((_) => info('ddt-data-changed'));
demog.onFilterChanged.subscribe((_) => info('ddt-filter-changed'));
demog.onSelectionChanged.subscribe((_) => info('ddt-selection-changed'));

// this will only get fired when a user changes conditional color-coding of a numerical column
demog.onMetadataChanged
  .pipe(rxjs.operators.filter(data => data.args.key === DG.TAGS.COLOR_CODING_CONDITIONAL))
  .subscribe((data) => info(`${data.args.change} - ${data.args.key} - ${data.args.value}`));
