/// this file was generated automatically from ddt classes declarations
import { toDart } from "../wrappers";
let api = <any>window;

export function histogram(col: any, bitset: any, flag: boolean, options?: {bins?: number, logScale?: boolean}): Int32List
  { return api.grok_histogram(toDart(col), toDart(bitset), toDart(flag), toDart(options?.bins), toDart(options?.logScale)); }

export class Tags {
  /// A user that created this entity
  static Description = 'description';

  static CreatedBy = 'createdBy';

  static SemanticDetectionDuration = '.semantic-detection-duration';

  static Tooltip = '.tooltip';

  static TooltipVisibility = '.tooltip-visibility';

  static TooltipForm = '.tooltip-form';

  static RowGroupTooltip = '.row-group-tooltip';

  static ValueFunction = '.value-function';

  static ColorCoding = '.color-coding';

  static ColorCodingType = '.color-coding-type';

  static ColorCodingText = '.color-coding-text';

  static ColorCodingLinear = '.color-coding-linear';

  static ColorCodingConditional = '.color-coding-conditional';

  static ColorCodingCategorical = '.color-coding-categorical';

  static ColorCodingValueColumn = '.color-coding-value-column';

  static ColorCodingSchemeMax = '.color-coding-scheme-max';

  static ColorCodingSchemeMin = '.color-coding-scheme-min';

  static DefaultAxisType = '.default-axis-type';

  static MarkerCoding = '.marker-coding';

  static HelpUrl = 'help.url';

  static OriginalTableName = '.orig.table.name';

  static OriginalViewName = '.orig.view.name';

  /// Comma-separated list of domains the entity belongs to
  static CellRenderer = 'cell.renderer';

  static Domains = 'domains';

  static Quality = 'quality';

  static Notation = 'notation';

  static LayoutColumnId = 'layout-column-id';

  static Units = 'units';

  static Format = 'format';

  static TooltipType = '.tooltip-type';

  static Markup = 'markup';

  /// JSON-encoded list of strings to be used in a cell editor.
  static SourcePrecision = '.source-precision';

  /// JSON-encoded list of strings to be used in a cell editor.
  /// Applicable for string columns only.
  static Formula = 'formula';

  static Choices = '.choices';

  static CategoryOrder = '.category-order';

  static DefaultFilter = '.default-filter';

  static IgnoreCustomFilter = '.ignore-custom-filter';

  /// When set to 'true', switches the cell editor to a combo box that only allows to choose values
  static CustomFilterType = '.custom-filter-type';

  /// When set to 'true', switches the cell editor to a combo box that only allows to choose values
  static Charts = '.charts';

  static AutoChoices = '.auto-choices';

  static ImportTime = 'import-time';

  static Source = 'source';

  static SourceUri = 'source.uri';

  static SourceFile = 'source.file';

  static Db = 'Db';

  static DbSchema = 'DbSchema';

  static DbTable = 'DbTable';

  static DbColumn = 'DbColumn';

  static DbPath = 'DbPath';

  static SparqlEndpoint = 'sparql.endpoint';

  static SparqlOntology = 'sparql.ontology';

  static AggregationPivotValues = 'aggr.pivot.values';

  static AggregationSourceTableName = 'aggr.src.table.name';

  static AggregationSourceColumnName = 'aggr.src.col.name';

  static AggregationType = 'aggr.type';

  static PredictiveModel = 'predictive.model';

  static Id = '.id';

  static TableInfo = '.table-info';

  static VariableName = '.VariableName';

  static DataConnectionId = '.data-connection-id';

  static DataConnection = '.data-connection';

  static DataSource = '.DataSource';

  static DataSchema = '.DataSchema';

  static DataQueryId = '.DataQuery.id';

  static DataQueryName = 'DataQuery.name';

  /// Expression that was used to derive the column.
  static DataQueryCall = '.DataQuery.query.call';

  static Presentation = '.presentation';

  static QueryJson = '.query-json';

  static Expression = 'expression';

  static TableSchema = 'table_schema';

  static ReferencesTable = 'references_table';

  static ReferencesColumn = 'references_column';

  static ChemDescriptor = 'chem-descriptor';

  static ChemFingerprinter = 'chem-fingerprinter';

  static DataHistory = '.history';

  static CreationScript = '.script';

  static DataSync = '.data-sync';

  static IsDbEntity = '.db.entity';

  static IsDbView = '.db.view';

  static ColumnConverter = '.column.converter';

  static ValueValidators = '.value-validators';

  static MultiValueSeparator = '.multi-value-separator';

  /// Applies to columns or dataframes.
  /// Comma-separated list of user or group names that are allowed to make changes to that column.
  static FriendlyName = 'friendlyName';

  static AllowRename = '.allow-rename';

  static EditableBy = 'editableBy';

  static PinIfEditable = 'pinIfEditable';

}
