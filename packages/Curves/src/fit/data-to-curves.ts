/* eslint-disable max-len */

import {FIT_FUNCTION_4PL_DOSE_RESPONSE, FIT_FUNCTION_4PL_REGRESSION, FitChartData, FitMarkerType, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as api from '../package-api';
import wu from 'wu';
import {FitFunctionType, FitFunctionTypes} from '@datagrok-libraries/statistics/src/fit/new-fit-API';

const groupColumnName = 'Compound|Assay|Target';

function setColumnDefaultValue(table: DG.DataFrame, input: DG.InputBase<DG.Column | null>, check: (col: DG.Column) => boolean, values: string[]) {
  const col = table.columns.toList().filter((c) => check(c)).sort((a, b) => a.name.length - b.name.length).find((c) => values.some((v) => c.name.toLowerCase().includes(v)));
  if (col)
    input.value = col;
}

function createFitParamsEditor(table: DG.DataFrame, type: FitFunctionType) {
  const inputs: DG.InputBase<DG.Column | null>[] = [];
  //['Top', 'Slope', 'EC50', 'Bottom'];
  switch (type) {
    case FitFunctionTypes.FOUR_PL_REGRESSION:
    case FitFunctionTypes.FOUR_PL_DOSE_RESPONSE:
    default:
      inputs.push(ui.input.column('Max', {table: table, nullable: false, tooltipText: 'Maximum value for the curve'}));
      inputs.push(ui.input.column('Hill', {table: table, nullable: false, tooltipText: 'Hill coefficient for the curve'}));
      inputs.push(ui.input.column('Inflection', {table: table, nullable: false, tooltipText: 'Inflection value for the curve'}));
      inputs.push(ui.input.column('Min', {table: table, nullable: false, tooltipText: 'Minimum value for the curve'}));
      setColumnDefaultValue(table, inputs[0], (c) => c.isNumerical && c.type !== DG.COLUMN_TYPE.DATE_TIME, ['max']);
      setColumnDefaultValue(table, inputs[1], (c) => c.isNumerical && c.type !== DG.COLUMN_TYPE.DATE_TIME, ['hill']);
      setColumnDefaultValue(table, inputs[2], (c) => c.isNumerical && c.type !== DG.COLUMN_TYPE.DATE_TIME, ['inflection']);
      setColumnDefaultValue(table, inputs[3], (c) => c.isNumerical && c.type !== DG.COLUMN_TYPE.DATE_TIME, ['min']);
      break;
  }
  return inputs;
}

export type FitFunctionEditorHistory = {
  fitParamCols: string[],
  fitFunction: FitFunctionType,
  useFitParamCols: boolean
}

function createFitFunctionEditor(table: DG.DataFrame) {
  const fitFunctions = [FitFunctionTypes.FOUR_PL_DOSE_RESPONSE, FitFunctionTypes.FOUR_PL_REGRESSION];
  const fitFunctionInput = ui.input.choice('Fit Function', {value: FitFunctionTypes.FOUR_PL_DOSE_RESPONSE, items: fitFunctions, nullable: false, tooltipText: 'Select the fit function to use for the curve fitting'});
  const form = ui.form([]);
  let fitParamInputs: DG.InputBase<DG.Column>[] = [];
  const useFitParamsInput = ui.input.bool('Use Prefit Parameters', {value: true, tooltipText: 'Use selected fit curve parameters instead of fitting the curve from scratch'});
  useFitParamsInput.onChanged.subscribe(() => {
    fitParamInputs.forEach((i) => i.enabled = useFitParamsInput.value);
  });
  // const checkContainer: HTMLDivElement | null = useFitParamsInput.root.querySelector('div.ui-input-editor');
  // if (checkContainer)
  //   checkContainer.style.setProperty('padding', '0px', 'important');
  fitFunctionInput.addOptions(useFitParamsInput.captionLabel);
  fitFunctionInput.addOptions(useFitParamsInput.input);
  useFitParamsInput.captionLabel.style.color = 'var(--grey-6)';
  const onChanged = () => {
    ui.empty(form);
    const inputs = createFitParamsEditor(table, fitFunctionInput.value!);
    fitParamInputs = inputs as DG.InputBase<DG.Column>[];
    ui.appendAll(form, [fitFunctionInput.root, ...inputs.map((i) => i.root)]);
    useFitParamsInput.fireChanged();
  };
  onChanged();
  fitFunctionInput.onChanged.subscribe(() => onChanged());
  return {
    form,
    getValues: () => useFitParamsInput.value ? fitParamInputs.map((i) => i.value) : [],
    getFitFunction: () => fitFunctionInput.value!,
    getHystory: (): FitFunctionEditorHistory => ({
      fitParamCols: fitParamInputs.map((i) => i.value?.name ?? ''),
      fitFunction: fitFunctionInput.value!,
      useFitParamCols: useFitParamsInput.value,
    }),
    applyHistory: (history: FitFunctionEditorHistory) => {
      // first, set the fit function
      fitFunctionInput.value = history.fitFunction as typeof fitFunctionInput.items[number];
      // then, set the fit parameters
      useFitParamsInput.value = history.useFitParamCols;
      if (history.fitParamCols?.length === fitParamInputs?.length) {
        for (let i = 0; i < fitParamInputs.length; i++) {
          const col = table.col(history.fitParamCols[i]);
          if (col)
            fitParamInputs[i].value = col;
        }
      }
    }
  };
}

export type ParentTableFormHistory = {
  tableName: string | null,
  reportedIC50Col: string | null,
  reportedQualifiedIC50Col: string | null,
  experimentIDCol: string | null,
  qualifierCol: string | null,
  additionalColumns: string[],
  fitParams: FitFunctionEditorHistory | null,
}

function createParentTableForm() {
  const tableInput = ui.input.table('Parent Table', {value: undefined, nullable: true, tooltipText: 'Table with parent data'});
  const form = ui.form([]);
  let fitFunctionEditor: ReturnType<typeof createFitFunctionEditor> | null = null;
  let reportedIC50Col: DG.InputBase<DG.Column | null>;
  let reportedQualifiedIC50Col: DG.InputBase<DG.Column | null>;
  let experimentIDCol: DG.InputBase<DG.Column | null>;
  let qualifierCol: DG.InputBase<DG.Column | null>;
  let additionalColumnsInput: DG.InputBase<DG.Column[] | null>;
  const onChanged = () => {
    ui.empty(form);
    form.appendChild(tableInput.root);
    if (!tableInput.value) {
      fitFunctionEditor = null;
      return;
    }
    reportedIC50Col = ui.input.column('Reported IC50', {table: tableInput.value!, nullable: true, tooltipText: 'Column with reported IC50 values'});
    reportedQualifiedIC50Col = ui.input.column('Reported Qualified IC50', {table: tableInput.value!, nullable: true, tooltipText: 'Column with reported qualified IC50 values'});
    experimentIDCol = ui.input.column('Experiment ID', {table: tableInput.value!, nullable: true, tooltipText: 'Column with experiment IDs'});
    qualifierCol = ui.input.column('Qualifier', {table: tableInput.value!, nullable: true, tooltipText: 'Column with qualifier values'});
    additionalColumnsInput = ui.input.columns('Additional Columns', {table: tableInput.value!, nullable: true, tooltipText: 'Additional columns to include in the curves data'});
    // the first two are sort of sketchy, because their columns (while expected to be numerical), can contain stuff like inconclusive, not applicable, etc.
    setColumnDefaultValue(tableInput.value!, reportedIC50Col, (c) => true && c.type !== DG.COLUMN_TYPE.DATE_TIME, ['reported ic50', 'reported_ic50', 'reported-ic50']);
    setColumnDefaultValue(tableInput.value!, reportedQualifiedIC50Col, (c) => true && c.type !== DG.COLUMN_TYPE.DATE_TIME, ['qualified_reported_ic50', 'qualified reported ic50', 'qualified-reported-ic50']);
    setColumnDefaultValue(tableInput.value!, experimentIDCol, (c) => c.isCategorical, ['experiment id', 'experiment_id', 'experiment-id']);
    setColumnDefaultValue(tableInput.value!, qualifierCol, (c) => c.isCategorical, ['qualifier']);
    const columnsForm = ui.form([reportedIC50Col, reportedQualifiedIC50Col, experimentIDCol, qualifierCol, additionalColumnsInput]);

    fitFunctionEditor = createFitFunctionEditor(tableInput.value!);
    ui.appendAll(form, [columnsForm, fitFunctionEditor.form]);
  };
  onChanged();
  tableInput.onChanged.subscribe(() => onChanged());

  return {
    getValue: () => ({
      fitParamCols: fitFunctionEditor?.getValues ? fitFunctionEditor.getValues() : null,
      reportedIC50Col: reportedIC50Col?.value,
      reportedQualifiedIC50Col: reportedQualifiedIC50Col?.value,
      experimentIDCol: experimentIDCol?.value,
      qualifierCol: qualifierCol?.value,
      additionalColumns: additionalColumnsInput?.value,
      fitFunction: fitFunctionEditor?.getFitFunction ? fitFunctionEditor.getFitFunction() : null,
    }),
    form,
    getTable: () => tableInput.value,
    getHistory: (): ParentTableFormHistory => ({
      tableName: tableInput.value?.name ?? null,
      reportedIC50Col: reportedIC50Col?.value?.name ?? null,
      reportedQualifiedIC50Col: reportedQualifiedIC50Col?.value?.name ?? null,
      experimentIDCol: experimentIDCol?.value?.name ?? null,
      qualifierCol: qualifierCol?.value?.name ?? null,
      additionalColumns: additionalColumnsInput?.value?.map((c) => c.name) ?? [],
      fitParams: fitFunctionEditor?.getHystory ? fitFunctionEditor.getHystory() : null,
    }),
    applyHistory: (history: ParentTableFormHistory) => {
      if (history.tableName) {
        const df = grok.shell.table(history.tableName);
        if (df)
          tableInput.value = df;
      }
      if (!tableInput.value)
        return;

      history.reportedIC50Col && tableInput.value.col(history.reportedIC50Col) && (reportedIC50Col.value = tableInput.value.col(history.reportedIC50Col));
      history.reportedQualifiedIC50Col && tableInput.value.col(history.reportedQualifiedIC50Col) && (reportedQualifiedIC50Col.value = tableInput.value.col(history.reportedQualifiedIC50Col));
      history.additionalColumns && (additionalColumnsInput.value = history.additionalColumns.map((c) => tableInput.value!.col(c)).filter((c) => c != null) as DG.Column[]);
      history.experimentIDCol && tableInput.value.col(history.experimentIDCol) && (experimentIDCol.value = tableInput.value.col(history.experimentIDCol));
      history.qualifierCol && tableInput.value.col(history.qualifierCol) && (qualifierCol.value = tableInput.value.col(history.qualifierCol));
      if (fitFunctionEditor && history.fitParams)
        fitFunctionEditor.applyHistory(history.fitParams);
    }
  };
}


function createWellTableForm() {
  const tableInput = ui.input.table('Well Table', {value: grok.shell.tv.dataFrame, nullable: false, tooltipText: 'Table with well level data'});
  let assayColInput: DG.InputBase<DG.Column | null>;
  let concentrationColInput: DG.InputBase<DG.Column | null>;
  let readoutColInput: DG.InputBase<DG.Column | null>;
  let excludeOutliersInput: DG.InputBase<DG.Column | null>;
  let compoundIDColInput: DG.InputBase<DG.Column | null>;
  let batchIDColInput: DG.InputBase<DG.Column | null>;
  let runIDColInput: DG.InputBase<DG.Column | null>;
  let targetEntityColInput: DG.InputBase<DG.Column | null>;
  const form = ui.form([]);

  const onTableInputChange = () => {
    const df = tableInput.value!;
    ui.empty(form);
    assayColInput = ui.input.column('Assay Name', {table: df, nullable: true, filter: (c) => c.isCategorical, tooltipText: 'Column with assay names'});
    concentrationColInput = ui.input.column('Concentration', {table: df, nullable: false, filter: (c) => c.isNumerical, tooltipText: 'Column with concentration values'});
    readoutColInput = ui.input.column('Readout', {table: df, nullable: false, filter: (c) => c.isNumerical, tooltipText: 'Column with readout values'});
    excludeOutliersInput = ui.input.column('Outliers', {table: df, nullable: true, filter: (c) => (c.isCategorical || c.type === DG.TYPE.BOOL) && c.categories.length === 2, tooltipText: 'Column with outliers to exclude from the analysis'});
    compoundIDColInput = ui.input.column('Compound ID', {table: df, nullable: true, tooltipText: 'Column with compound IDs'});
    batchIDColInput = ui.input.column('Batch ID', {table: df, nullable: false, tooltipText: 'Column with batch IDs. this column along with the Assay column will be used to group the data into curves'});
    runIDColInput = ui.input.column('Run ID', {table: df, nullable: true, tooltipText: 'Column with run IDs. this column along with the Assay column will be used to group the data into curves'});
    targetEntityColInput = ui.input.column('Target Entity', {table: df, nullable: false, filter: (c) => c.isCategorical, tooltipText: 'Column with target entity names. This column will be used to group the data into curves'});

    assayColInput.value = wu(df.columns.categorical).find((c) => c.name.toLowerCase().includes('assay') && c.name.toLowerCase().includes('name')) ?? null;
    concentrationColInput.value = df.columns.firstWhere((c) => c.isNumerical && c.type !== DG.COLUMN_TYPE.DATE_TIME && c.name.toLowerCase().includes('concentration')) ?? null;
    readoutColInput.value = df.columns.firstWhere((c) => c.isNumerical && c.type !== DG.COLUMN_TYPE.DATE_TIME && (c.name.toLowerCase().includes('inhibition'))) ?? null;
    excludeOutliersInput.value = df.columns.firstWhere((c) => (c.isCategorical || c.type === DG.TYPE.BOOL) && c.categories.length === 2 && (c.name.toLowerCase().includes('outlier') || (c.name.toLowerCase().includes('include') && c.name.toLowerCase().includes('exclude')))) ?? null;
    compoundIDColInput.value = df.columns.firstWhere((c) => c.isCategorical && c.name.toLowerCase().includes('compound_id') || c.name.toLowerCase().includes('compoundid') || c.name.toLowerCase().includes('compound id')) ?? null;
    batchIDColInput.value = df.columns.firstWhere((c) => c.name.toLowerCase().includes('batch') && c.name.toLowerCase().includes('id')) ?? null;
    runIDColInput.value = df.columns.firstWhere((c) => (c.name.toLowerCase().includes('run') || c.name.toLowerCase().includes('set')) && c.name.toLowerCase().includes('id')) ?? null;
    targetEntityColInput.value = df.columns.firstWhere((c) => c.isCategorical && (c.name.toLowerCase().includes('target') && c.name.toLowerCase().includes('entity'))) ?? null;
    ui.appendAll(form, [tableInput.root, assayColInput.root, batchIDColInput.root, runIDColInput.root, concentrationColInput.root, readoutColInput.root, compoundIDColInput.root, targetEntityColInput.root, excludeOutliersInput.root]);
  };
  onTableInputChange();
  tableInput.onChanged.subscribe(() => onTableInputChange());
  return {form, getValues: () => ({
    table: tableInput.value!,
    assayCol: assayColInput.value,
    batchIDCol: batchIDColInput.value,
    concentrationCol: concentrationColInput.value,
    readoutCol: readoutColInput.value,
    compoundIDCol: compoundIDColInput.value,
    excludeOutliersCol: excludeOutliersInput.value,
    targetEntityCol: targetEntityColInput.value,
    runIDCol: runIDColInput.value,
  }), applyHistory: (v: {[key in 'assayCol' | 'batchIDCol' | 'concentrationCol' | 'readoutCol' | 'compoundIDCol' | 'excludeOutliersCol' | 'runIDCol' | 'targetEntityCol'
  ]: string | null}
  ) => {
    const df = tableInput.value!;
    assayColInput.value = v.assayCol ? df.col(v.assayCol) : null;
    batchIDColInput.value = v.batchIDCol ? df.col(v.batchIDCol) : null;
    concentrationColInput.value = v.concentrationCol ? df.col(v.concentrationCol) : null;
    readoutColInput.value = v.readoutCol ? df.col(v.readoutCol) : null;
    compoundIDColInput.value = v.compoundIDCol ? df.col(v.compoundIDCol) : null;
    excludeOutliersInput.value = v.excludeOutliersCol ? df.col(v.excludeOutliersCol) : null;
    runIDColInput.value = v.runIDCol ? df.col(v.runIDCol) : null;
    targetEntityColInput.value = v.targetEntityCol ? df.col(v.targetEntityCol) : null;
  }
  };
}


export function dataToCurvesUI() {
  const tv = grok.shell.tv;
  if (!tv || !tv.dataFrame) {
    grok.shell.warning('No open tableview');
    return;
  }
  const wellLevelForm = createWellTableForm();
  const parentTableForm = createParentTableForm();
  const horzForm = ui.divH([wellLevelForm.form, parentTableForm.form], {style: {minWidth: '400px', gap: '20px'}});

  ui.dialog('Data to Curves')
    .add(horzForm)
    .onOK(async () => {
      const formRes = wellLevelForm.getValues();
      let parentLevelData: WellTableParentData | undefined = undefined;
      const parentTable = parentTableForm.getTable();
      const parentFormValues = parentTableForm.getValue();
      const fitParamValues = parentFormValues.fitParamCols;
      if (parentTable && fitParamValues && fitParamValues.length > 0) {
        parentLevelData = {
          table: parentTable,
          fitParamColumns: fitParamValues,
          fitFunction: parentFormValues.fitFunction!,
          reportedIC50Column: parentFormValues.reportedIC50Col!,
          reportedQualifiedIC50Column: parentFormValues.reportedQualifiedIC50Col!,
          experimentIDColumn: parentFormValues.experimentIDCol!,
          qualifierColumn: parentFormValues.qualifierCol!,
          additionalColumns: parentFormValues.additionalColumns ?? [],
        };
      }

      const {table, assayCol, batchIDCol, concentrationCol, readoutCol, compoundIDCol, excludeOutliersCol, runIDCol, targetEntityCol} = formRes;
      const func = DG.Func.find({name: 'dataToCurves'})[0];
      if (!func) {
        grok.shell.error('Function dataToCurves not found');
        return;
      }
      //parentTable?: DG.DataFrame, fitParamColumns?: DG.Column[], reportedIC50Column?: DG.Column, reportedQualifiedIC50Column?: DG.Column,
      //  experimentIDColumn?: DG.Column, qualifierColumn?: DG.Column, additionalColumns?: DG.Column[]
      const fc = func.prepare({
        df: table,
        concentrationCol: concentrationCol,
        readoutCol: readoutCol,
        batchIDCol: batchIDCol,
        assayCol: assayCol,
        runIDCol: runIDCol,
        compoundIDCol: compoundIDCol,
        targetEntityCol: targetEntityCol,
        excludeOutliersCol: excludeOutliersCol,
        // parentData:
        parentTable: parentLevelData?.table,
        fitParamColumns: parentLevelData?.fitParamColumns?.map((c) => c.name) ?? [],
        reportedIC50Column: parentLevelData?.reportedIC50Column?.name,
        reportedQualifiedIC50Column: parentLevelData?.reportedQualifiedIC50Column?.name,
        experimentIDColumn: parentLevelData?.experimentIDColumn?.name,
        qualifierColumn: parentLevelData?.qualifierColumn?.name,
        additionalColumns: parentLevelData?.additionalColumns?.map((c) => c.name) ?? [],
      });
      console.log(fc.toString());
      const needsCreationScript = table.tags[DG.Tags.CreationScript] && (!parentLevelData?.table || parentLevelData.table.tags[DG.Tags.CreationScript]);
      await fc.call(undefined, undefined, {processed: !needsCreationScript});
      //console.log(fc.getResultViews());
      //const resDf: DG.DataFrame | null = fc.getOutputParamValue();
      const tv = fc.getResultViews()[0];
      if (tv && tv instanceof DG.TableView && tv.dataFrame) {
        if (!needsCreationScript)
          grok.shell.addView(tv);

        tv.filters();
        const trellis = tv.trellisPlot({yColumnNames: [groupColumnName], xColumnNames: [], viewerType: 'MultiCurveViewer',
          showControlPanel: false, showXLabels: false, showYLabels: false, showXSelectors: false, showYSelectors: true, packCategories: true,
          onClick: 'Select'
        });
        tv.dockManager.dock(trellis, DG.DOCK_TYPE.TOP, tv.dockManager.findNode(tv.grid.root), 'Fitted Curves', 0.5);
        const selectionGrid = tv.addViewer(DG.VIEWER.GRID, {rowSource: 'Selected', selectedRowsColor: DG.Color.white});
        tv.dockManager.dock(selectionGrid, DG.DOCK_TYPE.FILL, tv.dockManager.findNode(tv.grid.root), 'Selected Curves');
        selectionGrid.props.title = 'Selected Curves';

        const multiCurveViewer = tv.addViewer('MultiCurveViewer', {showCurrentRowCurve: false, showMouseOverRowCurve: false, showSelectedRowsCurves: true});
        const md = tv.dockManager.dock(multiCurveViewer, DG.DOCK_TYPE.FILL, tv.dockManager.findNode(selectionGrid.root));
        tv.dockManager.dock(tv.grid, DG.DOCK_TYPE.FILL, md);
        tv.grid.props.title = 'Fitted Curves';
      }
    })
    .show()
    .history(() => {
      const formRes = wellLevelForm.getValues();
      const {assayCol, batchIDCol, concentrationCol, readoutCol, compoundIDCol, excludeOutliersCol, runIDCol, targetEntityCol} = formRes;
      const obj = {
        assayCol: assayCol?.name,
        batchIDCol: batchIDCol?.name,
        concentrationCol: concentrationCol?.name,
        readoutCol: readoutCol?.name,
        compoundIDCol: compoundIDCol?.name,
        excludeOutliersCol: excludeOutliersCol?.name,
        runIDCol: runIDCol?.name,
        targetEntityCol: targetEntityCol?.name,
      };
      return {wells: JSON.stringify(Object.fromEntries(Object.entries(obj).filter(([_, v]) => !!v)) ?? {}), parent: JSON.stringify(parentTableForm.getHistory() ?? {})};
    }, (v) => {
      const wells = JSON.parse(v?.wells ?? '{}');
      const parent = JSON.parse(v?.parent ?? '{}');
      wellLevelForm.applyHistory(wells);
      parentTableForm.applyHistory(parent);
    });
}
export type WellTableParentData = {
    table?: DG.DataFrame, fitParamColumns?: DG.Column[], fitFunction?: FitFunctionType, reportedIC50Column?: DG.Column, reportedQualifiedIC50Column?: DG.Column,
    experimentIDColumn?: DG.Column, qualifierColumn?: DG.Column, additionalColumns?: DG.Column[],
}

export async function convertDataToCurves(df: DG.DataFrame,
  concentrationCol: DG.Column, readoutCol: DG.Column, batchIDCol: DG.Column,
  assayCol: DG.Column, runIDCol: DG.Column, compoundIDCol: DG.Column, targetEntityCol: DG.Column, excludeOutliersCol?: DG.Column, parentData?: WellTableParentData
): Promise<DG.DataFrame> {
  if (!concentrationCol || !readoutCol || !batchIDCol) {
    grok.shell.warning('Please fill all required fields');
    throw new Error('Please fill all required fields');
  }

  const curvesObj: {[curveKey: string]: {x: number[], y: number[], outliers: boolean[], groupKey: string, runID: string, info: {[key: string]: any}}} = {};
  const assayColCats = assayCol?.categories;
  const assayColIndexes = assayCol?.getRawData();
  const consentrationRawData = concentrationCol.getRawData();
  const readoutRawData = readoutCol.getRawData();
  const batchIDCategories = batchIDCol?.categories;
  const batchIDRawData = batchIDCol?.getRawData();
  const compoundIDCategories = compoundIDCol?.categories;
  const compoundIDRawData = compoundIDCol?.getRawData();
  const runIDCategories = runIDCol?.categories;
  const runIDRawData = runIDCol?.getRawData();
  const targetEntityCategories = targetEntityCol?.categories;
  const targetEntityRawData = targetEntityCol?.getRawData();


  const getAssay = assayCol ? (row: number) => assayColCats![assayColIndexes![row]] : () => 'All';
  const otlierIndexes = excludeOutliersCol ? excludeOutliersCol.getRawData() : null;
  const isOutlier = excludeOutliersCol ? (excludeOutliersCol.type == DG.COLUMN_TYPE.BOOL ? (row: number) => !!otlierIndexes![row] : (row: number) => !otlierIndexes![row] ) : (_row: number) => false;
  const getConcentration = (row: number) => consentrationRawData[row];
  const getReadout = (row: number) => readoutRawData[row];
  const getBatchID = (row: number) => batchIDCategories[batchIDRawData[row]];
  const getCompoundID = compoundIDCol ? (row: number) => compoundIDCategories![compoundIDRawData![row]] : () => null;
  const getRunID = runIDCol ? (row: number) => runIDCategories![runIDRawData![row]] : () => null;
  const getTargetEntity = targetEntityCol ? (row: number) => targetEntityCategories![targetEntityRawData![row]] : () => null;
  const markerTypes: FitMarkerType[] = ['circle', 'triangle bottom', 'star', 'cross border', 'diamond', 'square', 'triangle left', 'triangle right', 'triangle top', 'asterisk'];
  let markerTypeIndex = 0; // will keep track of the markers
  const runIdToMarkerType: {[runID: string]: FitMarkerType} = {};
  const getMarkerType = (runID: string) => {
    if (runID in runIdToMarkerType) { return runIdToMarkerType[runID]; } else {
      const markerType = markerTypes[markerTypeIndex];
      runIdToMarkerType[runID] = markerType;
      markerTypeIndex = (markerTypeIndex + 1) % markerTypes.length;
      return markerType;
    }
  };

  const getGroupKey = (row: number) => `${getCompoundID(row)}||${getAssay(row)}||${getTargetEntity(row)}`;
  const getCurveKey = (row: number) => `${getGroupKey(row)}||${getRunID(row)}`;
  const getParentDataCurveKey = parentData && parentData.table?.col(compoundIDCol.name) && parentData.table.col(assayCol.name) && parentData.table.col(targetEntityCol.name) &&
  parentData.table.col(runIDCol.name) ? (row: number) =>
      `${parentData.table?.col(compoundIDCol.name)!.get(row)}||${parentData.table?.col(assayCol.name)!.get(row)}||${parentData.table?.col(targetEntityCol.name)!.get(row)}||${parentData.table?.col(runIDCol.name)!.get(row)}` :
    (_row: number) => null;
  // put the parent data into the parentObj
  const parentObj: {[curveKey: string]: {groupKey?: string, fitParams: number[], fitFunction: FitFunctionType, reportedIC50?: number, reportedQualifiedIC50?: number | string, experimentID?: string, qualifier?: string, additionalColumns?: {[colName: string]: any}}} = {};
  if (parentData) {
    for (let i = 0; i < (parentData.table?.rowCount ?? 0); i++) {
      const curveKey = getParentDataCurveKey ? getParentDataCurveKey(i) : null;
      if (!curveKey || parentObj[curveKey]) continue;
      parentObj[curveKey] = {
        fitParams: parentData.fitParamColumns?.map((c) => c.get(i)) ?? [],
        fitFunction: parentData.fitFunction ?? FIT_FUNCTION_4PL_DOSE_RESPONSE,
        reportedIC50: parentData.reportedIC50Column ? parentData.reportedIC50Column.get(i) : undefined,
        reportedQualifiedIC50: parentData.reportedQualifiedIC50Column ? parentData.reportedQualifiedIC50Column.get(i) : undefined,
        experimentID: parentData.experimentIDColumn ? parentData.experimentIDColumn.get(i) : undefined,
        qualifier: parentData.qualifierColumn ? parentData.qualifierColumn.get(i) : undefined,
        additionalColumns: parentData.additionalColumns ? Object.fromEntries(parentData.additionalColumns.map((c) => [c.name, c.get(i)])) : undefined,
      };
    }
  }

  for (let i = 0; i < df.rowCount; i++) {
    const assay = getAssay(i);
    const batchID = getBatchID(i);
    const concentration = getConcentration(i);
    const readout = getReadout(i);
    const outlier = isOutlier(i);
    const compoundID = getCompoundID(i);
    const targetEntity = getTargetEntity(i);
    const runID = getRunID(i) ?? '';

    const groupKey = getGroupKey(i);
    const curveKey = getCurveKey(i);

    if (!curvesObj[curveKey]) {
      curvesObj[curveKey] = {x: [], y: [], outliers: [], groupKey, runID: runID, info: {
        'Assay Name': assay,
        'Batch ID': batchID,
        'Compound ID': compoundID,
        'Target Entity': targetEntity,
        'Run ID': runID,
      }};
      if (parentData && parentObj[curveKey]) {
        curvesObj[curveKey].info['Reported IC50'] = parentObj[curveKey].reportedIC50;
        curvesObj[curveKey].info['Qualified Reported IC50'] = parentObj[curveKey].reportedQualifiedIC50;
        curvesObj[curveKey].info['Experiment ID'] = parentObj[curveKey].experimentID;
        curvesObj[curveKey].info['Qualifier'] = parentObj[curveKey].qualifier;
        if (parentObj[curveKey].additionalColumns)
          Object.assign(curvesObj[curveKey].info, parentObj[curveKey].additionalColumns);
        if (parentObj[curveKey].fitParams && parentObj[curveKey].fitParams.length > 0) {
          for (let j = 0; j < (parentData.fitParamColumns?.length ?? 0); j++) {
            const paramName = parentData.fitParamColumns![j].name;
            curvesObj[curveKey].info[paramName] = parentObj[curveKey].fitParams[j];
          }
        }
        parentObj[curveKey].groupKey = groupKey; // save the group key for the parent object
      }
    }
    const curve = curvesObj[curveKey];
    curve.x.push(concentration);
    curve.y.push(readout);
    curve.outliers.push(outlier);
  }

  // create fit series
  const colorHash: {[key: string]: string} = {};

  const getColor = (key: string): string => {
    if (colorHash[key]) return colorHash[key];
    const color = DG.Color.getCategoricalColor(Object.keys(colorHash).length);
    colorHash[key] = DG.Color.toHtml(color);
    return colorHash[key];
  };


  const curvesObjsArray = Object.entries(curvesObj);
  const tableColList: DG.Column[] = [];
  const curvesCol = DG.Column.string('Fitted Curve', curvesObjsArray.length);
  curvesCol.init((i) => {
    const [curveKey, curve] = curvesObjsArray[i];
    const x = curve.x;
    const y = curve.y;
    const outliers = curve.outliers;
    const runID = curve.runID;
    const markerType = getMarkerType(runID);
    const params = parentObj[curveKey] ? parentObj[curveKey].fitParams : undefined;
    const colorKey = (curve.info['Assay Name'] ?? '') + '||' + (curve.info['Target Entity'] ?? '');
    const color = getColor(colorKey);
    const fitFunctionName = parentObj?.[curveKey]?.fitFunction ? parentObj[curveKey].fitFunction : FIT_FUNCTION_4PL_DOSE_RESPONSE;
    const s: IFitSeries = {
      fit: undefined, fitFunction: fitFunctionName, clickToToggle: true, droplines: ['IC50'], name: `${runID}`,
      points: x.map((xv, i) => ({x: xv, y: y[i], outlier: outliers[i], marker: markerType, size: 5})).sort((a, b) => a.x - b.x), parameters: params, fitLineColor: color, pointColor: DG.Color.toHtml(DG.Color.gray)
    };
    const fitData: Omit<FitChartData, 'seriesOptions'> = {
      chartOptions: {
        xAxisName: concentrationCol.name,
        yAxisName: readoutCol.name,
        logX: true,
        title: `${curveKey}`,
      },
      series: [s],
    };
    return JSON.stringify(fitData);
  });
  tableColList.push(curvesCol);
  curvesCol.semType = 'fit';
  curvesCol.setTag('cell.renderer', 'fit');
  // add group column
  tableColList.push(DG.Column.fromStrings(groupColumnName, curvesObjsArray.map((c) => c[1].groupKey)));

  const otherColumns = curvesObjsArray.reduce((acc, c) => { Object.keys(c[1].info).forEach((k) => acc.add(k)); return acc; }, new Set<string>());

  //Object.keys(curvesObjsArray[0][1].info);

  for (const colName of otherColumns) {
    const col = DG.Column.fromStrings(colName, curvesObjsArray.map((c) =>
      c[1].info[colName]?.toString?.() ?? ''));
    tableColList.push(col);
  }


  const getNumericValue = (x: any) =>
    (x != undefined && !isNaN(x)) ?
      (typeof x === 'number' && x != DG.FLOAT_NULL && x != DG.INT_NULL) ? x :
        ((typeof x === 'string') ? ((Number.parseFloat(x) != undefined && !isNaN(Number.parseFloat(x))) ? Number.parseFloat(x) : undefined) : undefined) : undefined;

  // add additional calculated columns, based on the results
  // max percent inhibition, total results, all results and so on
  // all of them will be an aggregate of some kind, so we need to group the needed data
  const grouppedValues: {[groupKey: string]: {maxYs: number[], reportedIC50s: number[], qualifiedIC50s: string[],}} = {};
  for (const [curveKey, curve] of curvesObjsArray) {
    const groupKey = curve.groupKey;
    if (!groupKey) continue;
    if (!(groupKey in grouppedValues))
      grouppedValues[groupKey] = {maxYs: [], reportedIC50s: [], qualifiedIC50s: []};
    grouppedValues[groupKey].maxYs.push(Math.max(...curve.y.filter((_, i) => !curve.outliers[i])));
    if (parentData && parentObj[curveKey]) {
      if (parentObj[curveKey].reportedQualifiedIC50 ?? '' != '') {
        grouppedValues[groupKey].qualifiedIC50s.push(parentObj[curveKey].reportedQualifiedIC50 as string);
        let ic50 = getNumericValue(parentData.reportedQualifiedIC50Column ? parentObj[curveKey].reportedQualifiedIC50 : parentObj[curveKey].reportedIC50);
        if (ic50 == undefined) {
          try {
            ic50 = DG.Qnum.parse((parentData.reportedQualifiedIC50Column ? parentObj[curveKey].reportedQualifiedIC50 : parentObj[curveKey].reportedIC50) as string);
          } catch (_e) {
          }
        }
        if (ic50 != undefined)
          grouppedValues[groupKey].reportedIC50s.push(ic50);
      }
    }
  }
  const grouppedAggregations: {[groupKey: string]: {[key: string]: undefined | string | null}} = {};

  for (const groupKey in grouppedValues) {
    const values = grouppedValues[groupKey];
    grouppedAggregations[groupKey] = {
      'Max Percent Inhibition': values.maxYs.length ? Math.max(...values.maxYs)?.toString() : null,
      ...(parentData && parentData.table ? {
        'Number of total results': values.qualifiedIC50s.length?.toString(),
        'Number of reported results': values.reportedIC50s.length?.toString(),
        'Geomean IC50': values.reportedIC50s.length ? (Math.exp(values.reportedIC50s.reduce((a, b) => a + Math.log(b), 0) / values.reportedIC50s.length) || null)?.toString() : null, // log based, better for small values
        'Standard Deviation': values.reportedIC50s.length > 1 ? Math.sqrt(values.reportedIC50s.reduce((a, b) => a + (b - values.reportedIC50s.reduce((c, d) => c + d, 0) / values.reportedIC50s.length) ** 2, 0) / (values.reportedIC50s.length - 1))?.toString() : null,
        'All Results': values.qualifiedIC50s.join(', ')} : {}),
    };
  }
  const grouppedAggregationsColNames = Object.keys(grouppedAggregations[Object.keys(grouppedAggregations)[0]]);
  for (const colName of grouppedAggregationsColNames) {
    const col = DG.Column.fromStrings(colName, curvesObjsArray.map((c) => grouppedAggregations[c[1].groupKey]?.[colName] ?? ''));
    tableColList.push(col);
  }


  const resDF = DG.DataFrame.fromColumns(tableColList);
  resDF.name = 'Fitted Curves';

  if (!parentData || !parentData.table) {
    const actualStatNames = {
      'interceptX': 'IC50',
      'top': 'Max',
      'bottom': 'Min',
      'slope': 'Hill',
      'auc': 'AUC',
    };

    Object.entries(actualStatNames).forEach(([statName, alias]) => {
      const params = {table: resDF, colName: curvesCol.name, propName: statName, seriesNumber: 0};
      DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(params).callSync({processed: false}).getOutputParamValue().name = alias;
    });
    if (actualStatNames['interceptX'])
      resDF.col(actualStatNames['interceptX']) && (resDF.col(actualStatNames['interceptX'])!.meta.format = 'scientific');
  }
  return resDF;
};
