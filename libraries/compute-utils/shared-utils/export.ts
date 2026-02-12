import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type ExcelJS from 'exceljs';
import type html2canvas from 'html2canvas';
import {viewerTypesMapping} from './consts';
import {delay, getPropViewers} from './utils';

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

      let dfCounter = 0;

      const plotToSheet = async (
        sheet: ExcelJS.Worksheet,
        viewer: DG.Viewer | undefined,
        columnForImage: number, rowForImage = 0,
      ) => {
        if (!viewer || !viewer.dataFrame)
          return;
        const newViewer = DG.Viewer.fromType(viewer.type, viewer.dataFrame.clone());
        newViewer.copyViewersLook(viewer);

        let width = 1280;
        let height = 720;

        const viewerBox = ui.div(newViewer.root, {style: {
          width: `${width}px`,
          height: `${height}px`,
        }});
        viewerBox.classList.add('ui-box');
        viewerBox.classList.remove('ui-div');
        grok.shell.v.root.insertAdjacentElement('afterend', viewerBox);
        await delay(1000);
        const imageDataUrl = (await loadedHtml2canvas(viewerBox)).toDataURL();

        viewerBox.remove();

        const imageId = exportWorkbook.addImage({
          base64: imageDataUrl,
          extension: 'png',
        });

        sheet.addImage(imageId, {
          tl: {col: columnForImage, row: rowForImage},
          ext: {width, height},
        });
      };

      const isScalarType = (type: DG.TYPE) => (DG.TYPES_SCALAR.has(type));

      const isDataFrame = (prop: DG.Property) => (prop.propertyType === DG.TYPE.DATA_FRAME);

      const dfInputs = func.inputs.filter((input) => isDataFrame(input));
      const scalarInputs = func.inputs.filter((input) => isScalarType(input.propertyType));
      const dfOutputs = func.outputs.filter((output) => isDataFrame(output));
      const scalarOutputs = func.outputs.filter((output) => isScalarType(output.propertyType));

      dfInputs.forEach((dfInput) => {
        const visibleTitle = dfInput.options.caption || dfInput.name;
        const currentDfSheet =
      exportWorkbook.worksheets.find((ws) => ws.name === getSheetName(visibleTitle, exportWorkbook)) ??
      exportWorkbook.addWorksheet(getSheetName(visibleTitle, exportWorkbook));

        const currentDf = lastCall.inputs[dfInput.name];
        dfToSheet(currentDfSheet, dfCounter, currentDf);
        dfCounter++;
      });

      if (scalarInputs.length) {
        const inputScalarsSheet = exportWorkbook.addWorksheet('Input scalars');
        scalarsToSheet(inputScalarsSheet, scalarInputs.map((scalarInput) => ({
          caption: scalarInput.options['caption'] ?? scalarInput.name,
          value: lastCall.inputs[scalarInput.name] ?? '',
          units: scalarInput.options['units'] ?? '',
          format: scalarInput.format
        })));
      }

      dfOutputs.forEach((dfOutput) => {
        const visibleTitle = dfOutput.options.caption || dfOutput.name;
        const currentDfSheet =
      exportWorkbook.worksheets.find((ws) => ws.name === getSheetName(visibleTitle, exportWorkbook)) ??
      exportWorkbook.addWorksheet(getSheetName(visibleTitle, exportWorkbook));

        const currentDf = lastCall.outputs[dfOutput.name];
        dfToSheet(currentDfSheet, dfCounter, currentDf);
        dfCounter++;
      });


      if (scalarOutputs.length) {
        const outputScalarsSheet = exportWorkbook.addWorksheet('Output scalars');
        scalarsToSheet(outputScalarsSheet, scalarOutputs.map((scalarOutput) => ({
          caption: scalarOutput.options['caption'] ?? scalarOutput.name,
          value: lastCall.outputs[scalarOutput.name] ?? '',
          units: scalarOutput.options['units'] ?? '',
          format: scalarOutput.format,
        })));
      }

      for (const inputProp of func.inputs.filter((prop) => isDataFrame(prop))) {
        const nonGridViewers = (dfToViewerMapping[inputProp.name] ?? [])
          .filter((viewer) => viewer && viewer.type !== DG.VIEWER.GRID)
          .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer!.type));

        if (nonGridViewers.length === 0) continue;

        const visibleTitle = inputProp.options.caption || inputProp.name;
        const currentDf = lastCall.inputs[inputProp.name];

        for (const [index, viewer] of nonGridViewers.entries()) {
          await plotToSheet(
            exportWorkbook.getWorksheet(getSheetName(visibleTitle, exportWorkbook))!,
            viewer,
            currentDf.columns.length + 2,
            (index > 0) ? (index * 16) + 1 : 0,
          );
        };
      }

      for (const outputProp of func.outputs.filter((prop) => isDataFrame(prop))) {
        const nonGridViewers = (dfToViewerMapping[outputProp.name] ?? [])
          .filter((viewer) => viewer && viewer.type !== DG.VIEWER.GRID)
          .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer!.type));

        if (nonGridViewers.length === 0) continue;

        const visibleTitle = outputProp.options.caption || outputProp.name;
        const currentDf = lastCall.outputs[outputProp.name];

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
              exportWorkbook.getWorksheet(getSheetName(visibleTitle, exportWorkbook))!,
              dfCounter,
              stats,
              currentDf.columns.length + 2,
              (index > 0 && nonGridViewers[index-1]) ? Math.ceil(nonGridViewers[index-1]!.root.clientHeight / 20) + 1 : 0,
            );
            dfCounter++;
          } else {
            await plotToSheet(
              exportWorkbook.getWorksheet(getSheetName(visibleTitle, exportWorkbook))!,
              viewer,
              currentDf.columns.length + 2,
              (index > 0) ? (index * 16) + 1 : 0,
            );
          }
        }
      }

      const buffer = await exportWorkbook.xlsx.writeBuffer();

      const blob = new Blob([buffer], {type: BLOB_TYPE});

      return [blob, exportWorkbook] as const;
    } catch (e) {
      console.log(e);
    } finally {
      ui.setDisplay(grok.shell.v.root, true);
      updateIndicatorWithText(grok.shell.v.root.parentElement!, false);
    }
  }

  throw new Error('Format is not supported');
};

export const reportFuncCallExcel = async (funccall: DG.FuncCall) => {
  return richFunctionViewReport(
    'Excel',
    funccall.func,
    funccall,
    dfToViewerMapping(funccall),
  );
};

export const formatNumber = (val: any, format?: string) => {
  try {
    return (Number.isFinite(val) && format) ? DG.format(val, format) : val;
  } catch {
    return val;
  }
}

export const scalarsToSheet =
  (sheet: ExcelJS.Worksheet, scalars: { caption: string, value: string, units: string, format?: string }[]) => {
    sheet.addRow(['Parameter', 'Value', 'Units']).font = {bold: true};
    scalars.forEach((scalar) => {
      sheet.addRow([scalar.caption, formatNumber(scalar.value, scalar.format), scalar.units]);
    });

    sheet.getColumn(1).width = Math.max(
      ...scalars.map((scalar) => scalar.caption.toString().length), 'Parameter'.length,
    ) * 1.2;
    sheet.getColumn(2).width = Math.max(
      ...scalars.map((scalar) => scalar.value.toString().length), 'Value'.length) * 1.2;
    sheet.getColumn(3).width = Math.max(
      ...scalars.map((scalar) => scalar.units.toString().length), 'Units'.length) * 1.2;
  };


export const dfToSheet = (sheet: ExcelJS.Worksheet, dfCounter: number, df?: DG.DataFrame, column?: number, row?: number) => {
  if (!df)
    return;
  const columnKey = sheet.getColumn(column ?? 1).letter;
  const rows = [];
  for (let rowIdx = 0; rowIdx < df.rowCount; rowIdx++) {
    const row = [];
    for (let colIdx = 0; colIdx < df.columns.length; colIdx++) {
      const col = df.col(colIdx)!;
      const rawVal = df.get(col.name, rowIdx);
      if (col?.type === 'double') {
        const format = col?.tags?.['format'] ?? '.00';
        const val =  formatNumber(rawVal, format);
        row.push(val);
      } else {
        row.push(rawVal);
      }
    }
    rows.push(row);
  }
  const tableConfig = {
    name: `ID_${dfCounter.toString()}`,
    ref: `${columnKey}${row ?? 1}`,
    columns: df.columns.toList().map((col) => ({name: col.name, filterButton: false})),
    rows,
  };
  sheet.addTable(tableConfig);
  sheet.columns.forEach((col) => {
    col.width = 25;
    col.alignment = {wrapText: true};
  });
};

const isDataFrame = (prop: DG.Property) => (prop.propertyType === DG.TYPE.DATA_FRAME);

const configToViewer = async (df: DG.DataFrame | undefined, config: Record<string, any>) => {
  if (!df)
    return;
  const type = config['type'];
  const viewer = await df.plot.fromType(type) as DG.Viewer;
  viewer.setOptions(config);

  return viewer;
};

export const dfToViewerMapping = (funcCall: DG.FuncCall) => {
  const func = funcCall.func;

  const mapping = {} as Record<string, (DG.Viewer | undefined)[]>;
  Promise.all(func.inputs
    .filter((output) => isDataFrame(output))
    .map(async (p) => {
      mapping[p.name] = await Promise.all(getPropViewers(p).config
        .map((config) => configToViewer(funcCall.inputs[p.name], config)));

      return mapping[p.name];
    }));

  Promise.all(func.outputs
    .filter((output) => isDataFrame(output))
    .map(async (p) => {
      mapping[p.name] = await Promise.all(getPropViewers(p).config
        .map((config) => configToViewer(funcCall.outputs[p.name], config)));

      return mapping[p.name];
    }));

  return mapping;
};
