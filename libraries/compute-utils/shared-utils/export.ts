import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type ExcelJS from 'exceljs';
import type html2canvas from 'html2canvas';
import {interval, race, timer} from 'rxjs';
import {filter, switchMap, take} from 'rxjs/operators';
import {DEFAULT_FLOAT_FORMAT, G7_FLOAT_FORMAT} from '@datagrok-libraries/utils/src/format-version-utils';
import {viewerTypesMapping} from './consts';
import {getPropViewers, graphicsDataUrl, isSvgGraphics} from './utils';
import type {ValidationResult} from '../reactive-tree-driver/src/data/common-types';
import type {ConsistencyInfo} from '../reactive-tree-driver/src/runtime/StateTreeNodes';

// Exported chart raster size, and the per-chart row stride used to stack several charts for the
// same dataframe down a sheet. Excel's default row is 15pt (~20px), so a 720px image spans ~36
// rows; the stride must cover that (plus a 1-row gap) or consecutive charts overlap.
const EXPORT_VIEWER_WIDTH = 1280;
const EXPORT_VIEWER_HEIGHT = 720;
const EXCEL_ROW_HEIGHT_PX = 20;
const CHART_ROW_STRIDE = Math.ceil(EXPORT_VIEWER_HEIGHT / EXCEL_ROW_HEIGHT_PX) + 1;
const MAX_RENDER_WAIT_MS = 2000;

// isRenderPending clears shortly after the last render event, without another event. Viewers that send
// no render events (JS viewers, forms) or stay pending are captured after MAX_RENDER_WAIT_MS.
const viewerRendered = (viewer: DG.Viewer) => race(
  viewer.onViewerRendered.pipe(
    take(1),
    switchMap(() => interval(10)),
    filter(() => !viewer.isRenderPending),
  ),
  timer(MAX_RENDER_WAIT_MS),
).pipe(take(1)).toPromise();

const updateIndicatorWithText = (element: HTMLElement, updating: boolean, text?: string) => {
  ui.setUpdateIndicator(element, updating);
  const updatingLabel = element.querySelector('.d4-update-shadow .ui-label');
  if (updating && text && updatingLabel)
    updatingLabel.textContent = text;
};

export const richFunctionViewReport = async (
  format: string,
  func: DG.Func,
  lastCall: DG.FuncCall,
  dfToViewerMapping: {[key: string]: (DG.Viewer | undefined)[]},
  validationStates?: Record<string, ValidationResult>,
  consistencyStates?: Record<string, ConsistencyInfo>,
) => {
  const sheetNamesCache = {} as Record<string, string>;

  const getSheetName = (initialName: string, wb: ExcelJS.Workbook) => {
    if (sheetNamesCache[initialName]) return sheetNamesCache[initialName];

    let name = `${initialName}`;
    if (name.length > 31)
      name = `${name.slice(0, 31)}`;
    let i = 1;
    while (wb.worksheets.some((sheet) => sheet.name.toLowerCase() === name.toLowerCase())) {
      let truncatedName = `${initialName}`;
      if (truncatedName.length > (31 - `-${i}`.length))
        truncatedName = `${initialName.slice(0, 31 - `-${i}`.length)}`;
      name = `${truncatedName}-${i}`;
      i++;
    }

    sheetNamesCache[initialName] = name;

    return name;
  };

  if (format === 'Excel') {
    try {
      await DG.Utils.loadJsCss(['/js/common/exceljs.min.js', '/js/common/html2canvas.min.js']);
      //@ts-ignore
      const loadedExcelJS = window.ExcelJS as ExcelJS;
      //@ts-ignore
      const loadedHtml2canvas: typeof html2canvas = window.html2canvas;

      const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
      const exportWorkbook = new loadedExcelJS.Workbook() as ExcelJS.Workbook;

      ui.setDisplay(grok.shell.v.root, false);
      updateIndicatorWithText(grok.shell.v.root.parentElement!, true,
        'Generating report. Please do not switch the browser tab...');

      const plots: {sheet: ExcelJS.Worksheet, viewer: DG.Viewer | undefined, col: number, row: number}[] = [];

      // All charts render at once, then each is captured in turn
      const plotsToSheets = async () => {
        const boxes = plots.filter(({viewer}) => viewer?.dataFrame).map((plot) => {
          const newViewer = DG.Viewer.fromType(plot.viewer!.type, plot.viewer!.dataFrame.clone());
          newViewer.copyViewersLook(plot.viewer!);
          // the view container is a flex column, which would stretch the box past the image size
          const box = ui.div(newViewer.root, {style: {
            width: `${EXPORT_VIEWER_WIDTH}px`,
            height: `${EXPORT_VIEWER_HEIGHT}px`,
            flex: 'none',
          }});
          box.classList.add('ui-box');
          box.classList.remove('ui-div');
          return {...plot, newViewer, box, rendered: viewerRendered(newViewer)};
        });
        try {
          grok.shell.v.root.after(...boxes.map(({box}) => box));
          await Promise.all(boxes.map(({rendered}) => rendered));
          for (const {sheet, box, col, row} of boxes) {
            // html2canvas clones the whole document unless told to skip what is not around the box
            const canvas = await loadedHtml2canvas(box, {
              logging: false,
              ignoreElements: (el) => !el.contains(box) && !box.contains(el) && !document.head.contains(el),
            });
            sheet.addImage(exportWorkbook.addImage({base64: canvas.toDataURL(), extension: 'png'}), {
              tl: {col, row},
              ext: {width: EXPORT_VIEWER_WIDTH, height: EXPORT_VIEWER_HEIGHT},
            });
          }
        } finally {
          for (const {box, newViewer} of boxes) {
            box.remove();
            // removing the DOM node does not release the dart viewer or its cloned dataframe
            newViewer.detach();
          }
        }
      };

      const isScalarType = (type: DG.TYPE) => (DG.TYPES_SCALAR.has(type));

      const isDataFrame = (prop: DG.Property) => (prop.propertyType === DG.TYPE.DATA_FRAME);

      const dfInputs = func.inputs.filter((input) => isDataFrame(input));
      const scalarInputs = func.inputs.filter((input) => isScalarType(input.propertyType));
      const dfOutputs = func.outputs.filter((output) => isDataFrame(output));
      const scalarOutputs = func.outputs.filter((output) => isScalarType(output.propertyType));

      dfInputs.forEach((dfInput) => {
        const currentDf = lastCall.inputs[dfInput.name];
        if (isEmptyDf(currentDf))
          return;
        const visibleTitle = dfInput.options.caption || dfInput.name;
        const currentDfSheet =
      exportWorkbook.worksheets.find((ws) => ws.name === getSheetName(visibleTitle, exportWorkbook)) ??
      exportWorkbook.addWorksheet(getSheetName(visibleTitle, exportWorkbook));

        const validation = validationStates?.[dfInput.name];
        const consistency = consistencyStates?.[dfInput.name];
        dfToSheet({sheet: currentDfSheet, df: currentDf, validation, consistency});
      });

      if (scalarInputs.length) {
        const inputScalarsSheet = exportWorkbook.addWorksheet('Input scalars');
        scalarsToSheet(
          inputScalarsSheet,
          scalarInputs.map((scalarInput) => ({
            name: scalarInput.name,
            caption: scalarInput.options['caption'] ?? scalarInput.name,
            value: lastCall.inputs[scalarInput.name] ?? '',
            units: scalarInput.options['units'] ?? '',
            // the form's default, so the export matches it whether or not the form was shown
            format: scalarInput.format ||
              (scalarInput.propertyType === DG.TYPE.FLOAT ? G7_FLOAT_FORMAT : DEFAULT_FLOAT_FORMAT),
          })),
          validationStates,
          consistencyStates,
        );
      }

      dfOutputs.forEach((dfOutput) => {
        const currentDf = lastCall.outputs[dfOutput.name];
        if (isEmptyDf(currentDf))
          return;
        const visibleTitle = dfOutput.options.caption || dfOutput.name;
        const currentDfSheet =
      exportWorkbook.worksheets.find((ws) => ws.name === getSheetName(visibleTitle, exportWorkbook)) ??
      exportWorkbook.addWorksheet(getSheetName(visibleTitle, exportWorkbook));

        const validation = validationStates?.[dfOutput.name];
        const consistency = consistencyStates?.[dfOutput.name];
        dfToSheet({sheet: currentDfSheet, df: currentDf, validation, consistency});
      });


      if (scalarOutputs.length) {
        const outputScalarsSheet = exportWorkbook.addWorksheet('Output scalars');
        scalarsToSheet(
          outputScalarsSheet,
          scalarOutputs.map((scalarOutput) => ({
            name: scalarOutput.name,
            caption: scalarOutput.options['caption'] ?? scalarOutput.name,
            value: lastCall.outputs[scalarOutput.name] ?? '',
            units: scalarOutput.options['units'] ?? '',
            format: scalarOutput.format || DEFAULT_FLOAT_FORMAT,
          })),
          validationStates,
          consistencyStates,
        );
      }

      for (const graphicsOutput of func.outputs.filter((output) => output.propertyType === DG.TYPE.GRAPHICS)) {
        const value = lastCall.outputs[graphicsOutput.name];
        if (value == null || value === '')
          continue;
        const visibleTitle = graphicsOutput.options.caption || graphicsOutput.name;
        const sheet = exportWorkbook.addWorksheet(getSheetName(visibleTitle, exportWorkbook));
        const {dataUrl, width, height} = await graphicsToPng(value);
        sheet.addImage(exportWorkbook.addImage({base64: dataUrl, extension: 'png'}), {
          tl: {col: 0, row: 0},
          ext: {width, height},
        });
      }

      for (const inputProp of func.inputs.filter((prop) => isDataFrame(prop))) {
        const currentDf = lastCall.inputs[inputProp.name];
        if (isEmptyDf(currentDf)) continue;

        const nonGridViewers = (dfToViewerMapping[inputProp.name] ?? [])
          .filter((viewer) => viewer && viewer.type !== DG.VIEWER.GRID)
          .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer!.type));

        if (nonGridViewers.length === 0) continue;

        const visibleTitle = inputProp.options.caption || inputProp.name;

        for (const [index, viewer] of nonGridViewers.entries()) {
          plots.push({
            sheet: exportWorkbook.getWorksheet(getSheetName(visibleTitle, exportWorkbook))!,
            viewer,
            col: currentDf.columns.length + 2,
            row: index * CHART_ROW_STRIDE,
          });
        };
      }

      for (const outputProp of func.outputs.filter((prop) => isDataFrame(prop))) {
        const currentDf = lastCall.outputs[outputProp.name];
        if (isEmptyDf(currentDf)) continue;

        const nonGridViewers = (dfToViewerMapping[outputProp.name] ?? [])
          .filter((viewer) => viewer && viewer.type !== DG.VIEWER.GRID)
          .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer!.type));

        if (nonGridViewers.length === 0) continue;

        const visibleTitle = outputProp.options.caption || outputProp.name;

        for (const [index, viewer] of nonGridViewers.entries()) {
          if (!viewer)
            continue;
          if (viewer.type === DG.VIEWER.STATISTICS) {
            const length = currentDf.columns.length;
            const stats = DG.DataFrame.fromColumns([
              DG.Column.string('Name', length).init((i: number) => currentDf.columns.byIndex(i).name),
              DG.Column.int('Values', length).init((i: number) => currentDf.columns.byIndex(i).stats.valueCount),
              DG.Column.int('Nulls', length).init((i: number) => currentDf.columns.byIndex(i).stats.missingValueCount),
              DG.Column.float('Min', length).init((i: number) => currentDf.columns.byIndex(i).stats.min),
              DG.Column.float('Max', length).init((i: number) => currentDf.columns.byIndex(i).stats.max),
              DG.Column.float('Avg', length).init((i: number) => currentDf.columns.byIndex(i).stats.avg),
              DG.Column.float('Stdev', length).init((i: number) => currentDf.columns.byIndex(i).stats.stdev),
            ]);
            dfToSheet(
              {sheet: exportWorkbook.getWorksheet(getSheetName(visibleTitle, exportWorkbook))!, df: stats, column: currentDf.columns.length + 2, row: index * CHART_ROW_STRIDE},
            );
          } else {
            plots.push({
              sheet: exportWorkbook.getWorksheet(getSheetName(visibleTitle, exportWorkbook))!,
              viewer,
              col: currentDf.columns.length + 2,
              row: index * CHART_ROW_STRIDE,
            });
          }
        }
      }

      await plotsToSheets();

      const buffer = await exportWorkbook.xlsx.writeBuffer();

      const blob = new Blob([buffer], {type: BLOB_TYPE});

      return [blob, exportWorkbook] as const;
    } finally {
      ui.setDisplay(grok.shell.v.root, true);
      updateIndicatorWithText(grok.shell.v.root.parentElement!, false);
    }
  }

  throw new Error('Format is not supported');
};

const formatNumber = (val: any, format?: string) => {
  try {
    return (Number.isFinite(val) && format) ? DG.format(val, format) : val;
  } catch {
    return val;
  }
};

const scalarsToSheet = (sheet: ExcelJS.Worksheet,
  scalars: { name: string, caption: string, value: string, units: string, format?: string }[],
  validationStates?: Record<string, ValidationResult>,
  consistencyStates?: Record<string, ConsistencyInfo>,
) => {
  const cols = ['Parameter', 'Value', 'Units'];
  if (validationStates)
    cols.push('Validation');
  if (consistencyStates)
    cols.push('Consistency');
  sheet.addRow(cols).font = {bold: true};
  scalars.forEach((scalar) => {
    const vals = [scalar.caption, formatNumber(scalar.value, scalar.format), scalar.units];
    vals.push(getValidationString(validationStates?.[scalar.name]));
    vals.push(getConsistencyString(consistencyStates?.[scalar.name]));

    sheet.addRow(vals);
  });

  sheet.getColumn(1).width = Math.max(
    ...scalars.map((scalar) => scalar.caption.toString().length), 'Parameter'.length,
  ) * 1.2;
  sheet.getColumn(2).width = Math.max(
    ...scalars.map((scalar) => scalar.value.toString().length), 'Value'.length) * 1.2;
  sheet.getColumn(3).width = Math.max(
    ...scalars.map((scalar) => scalar.units.toString().length), 'Units'.length) * 1.2;

  sheet.getColumn(4).width = 100;
  sheet.getColumn(5).width = 100;
};

const getValidationString = (data?: ValidationResult) => {
  if (data == null)
    return '';
  const validations = (['errors', 'warnings'] as const)
    .filter((category) => !!data?.[category]?.length)
    .flatMap((category) => data![category]!.map((advice) =>({category, advice})))
    .map(({category, advice}) => `${category}: ${typeof advice === 'string' ? advice : advice.description}`)
    .join('\n');
  return validations;
};

const getConsistencyString = (data?: ConsistencyInfo) => {
  if (data == null)
    return '';
  if (data.inconsistent && (data.restriction === 'disabled' || data.restriction === 'restricted')) {
    const v = data.assignedValue;
    const formatted = (typeof v === 'number' && Number.isFinite(v)) ? DG.format(v, DEFAULT_FLOAT_FORMAT) : String(v);
    return `Inconsistent: value should be ${formatted}`;
  }
  return;
};

const dfToSheet = (
  {
    sheet,
    df,
    column,
    row,
    validation,
    consistency,
  } :
 { sheet: ExcelJS.Worksheet; df?: DG.DataFrame; column?: number; row?: number; validation?: ValidationResult; consistency?: ConsistencyInfo },
) => {
  if (!df)
    return;
  const columnKey = sheet.getColumn(column ?? 1).letter;
  const nCols = df.columns.length;
  const cols = df.columns.toList();
  const formats = cols.map((col) => col.type === 'double' ? col.tags?.['format'] ?? DEFAULT_FLOAT_FORMAT : undefined);
  const rows = [];
  for (let rowIdx = 0; rowIdx < df.rowCount; rowIdx++)
    rows.push(cols.map((col, colIdx) => formats[colIdx] ? formatNumber(col.get(rowIdx), formats[colIdx]) : col.get(rowIdx)));
  const tableConfig = {
    // Excel table names are unique in the workbook
    name: `ID_${sheet.workbook.worksheets.reduce((count, ws) => count + ws.getTables().length, 0)}`,
    ref: `${columnKey}${row ?? 1}`,
    columns: cols.map((col) => ({name: col.name, filterButton: false})),
    rows,
  };
  // Excel has no table without columns, so a table with none leaves the sheet empty
  if (nCols > 0) {
    sheet.addTable(tableConfig);
    sheet.columns.forEach((col) => {
      col.width = 25;
      col.alignment = {wrapText: true};
    });
  }
  const metaCol = nCols + 1;
  sheet.getColumn(metaCol).width = 100;
  sheet.getCell(1, metaCol).value = 'Validation and Consistency';
  if (validation)
    sheet.getCell(2, metaCol).value = getValidationString(validation);

  if (consistency)
    sheet.getCell(3, metaCol).value = getConsistencyString(consistency);

}

const isDataFrame = (prop: DG.Property) => (prop.propertyType === DG.TYPE.DATA_FRAME);

const loadImage = (src: string) => new Promise<HTMLImageElement>((resolve, reject) => {
  const img = new Image();
  img.onload = () => resolve(img);
  img.onerror = () => reject(new Error('Failed to load a graphics output'));
  img.src = src;
});

const svgSize = (markup: string) => {
  const svg = new DOMParser().parseFromString(markup, 'image/svg+xml').documentElement;
  const px = (attr: string) => /^\s*[\d.]+(px)?\s*$/.test(svg.getAttribute(attr) ?? '') ?
    parseFloat(svg.getAttribute(attr)!) : NaN;
  const [width, height] = [px('width'), px('height')];
  if (width > 0 && height > 0)
    return {width, height};
  const [, , vbWidth, vbHeight] = (svg.getAttribute('viewBox') ?? '').trim().split(/[\s,]+/).map(Number);
  if (vbWidth > 0 && vbHeight > 0)
    return {width: vbWidth, height: vbHeight};
  // no intrinsic proportions: square rather than stretched to the chart box
  return {width: 1, height: 1};
};

// Raster graphics keep their native pixel size; SVG is rasterized into the chart box, keeping its aspect ratio
const graphicsToPng = async (value: string) => {
  const dataUrl = graphicsDataUrl(value);
  const img = await loadImage(dataUrl);
  if (!isSvgGraphics(value))
    return {dataUrl, width: img.naturalWidth, height: img.naturalHeight};

  const size = svgSize(value);
  const scale = Math.min(EXPORT_VIEWER_WIDTH / size.width, EXPORT_VIEWER_HEIGHT / size.height);
  const width = Math.round(size.width * scale);
  const height = Math.round(size.height * scale);
  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  canvas.getContext('2d')!.drawImage(img, 0, 0, width, height);
  return {dataUrl: canvas.toDataURL('image/png'), width, height};
};

const isEmptyDf = (df?: DG.DataFrame) => !df || df.rowCount === 0;

const configToViewer = async (df: DG.DataFrame | undefined, config: Record<string, any>) => {
  if (!df)
    return;
  const type = config['type'];
  const viewer = await df.plot.fromType(type) as DG.Viewer;
  viewer.setOptions(config);

  return viewer;
};

export const dfToViewerMapping = async (funcCall: DG.FuncCall) => {
  const func = funcCall.func;

  const mapping = {} as Record<string, (DG.Viewer | undefined)[]>;
  await Promise.all(func.inputs
    .filter((output) => isDataFrame(output))
    .map(async (p) => {
      mapping[p.name] = await Promise.all(getPropViewers(p).config
        .map((config) => configToViewer(funcCall.inputs[p.name], config)));
    }));

  await Promise.all(func.outputs
    .filter((output) => isDataFrame(output))
    .map(async (p) => {
      mapping[p.name] = await Promise.all(getPropViewers(p).config
        .map((config) => configToViewer(funcCall.outputs[p.name], config)));
    }));

  return mapping;
};
