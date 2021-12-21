import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import ExcelJS from 'exceljs';
import {saveAs} from 'file-saver';


const getSheetName = (name: string, direction: DIRECTION) => {
  const idealName = `${direction} - ${name}`;
  return (idealName.length > 31) ? name.substring(0, 32) : idealName;
};

enum DIRECTION {
  INPUT = 'Input',
  OUTPUT = 'Output'
}

export async function exportFuncCall(call: DG.FuncCall) {
  //todo: check status
  // if (call.status != FuncCall.STATUS_COMPLETED) ...

  const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
  const exportWorkbook = new ExcelJS.Workbook();

  const isScalarType = (type: DG.TYPE) => {
    return DG.TYPES_SCALAR.has(type);
  };

  const isDataFrame = (type: DG.TYPE) => (type === DG.TYPE.DATA_FRAME);

  const dfInputs = call.func.inputs.filter((input) => isDataFrame(input.propertyType));
  const scalarInputs = call.func.inputs.filter((input) => isScalarType(input.propertyType));
  const dfOutputs = call.func.outputs.filter((output) => isDataFrame(output.propertyType));
  const scalarOutputs = call.func.outputs.filter((output) => isScalarType(output.propertyType));

  const realNames: Record<string, string> = {};

  dfInputs.forEach((dfInput) => {
    const visibleTitle = dfInput.options.caption || dfInput.name;
    const currentDfSheet = exportWorkbook.addWorksheet(getSheetName(visibleTitle, DIRECTION.INPUT));
    realNames[visibleTitle] = dfInput.name;
    const currentDf = (call.inputs[dfInput.name] as DG.DataFrame);
    currentDfSheet.addRow((currentDf.columns as DG.ColumnList).names());
    for (let i = 0; i < currentDf.rowCount; i++) {
      currentDfSheet.addRow([...currentDf.row(i).cells].map((cell: DG.Cell) => cell.value));
    }
  });

  if (scalarInputs.length) {
    const inputScalarsSheet = exportWorkbook.addWorksheet('Input scalars');
    scalarInputs.forEach((scalarInput) => {
      inputScalarsSheet.addRow([scalarInput.options['caption'] || scalarInput.name, call.inputs[scalarInput.name], scalarInput.options['units'] || '']);
    });
  }

  dfOutputs.forEach((dfOutput) => {
    const visibleTitle = dfOutput.options.caption || dfOutput.name;
    const currentDfSheet = exportWorkbook.addWorksheet(getSheetName(visibleTitle, DIRECTION.OUTPUT));
    realNames[visibleTitle] = dfOutput.name;
    const currentDf = (call.outputs[dfOutput.name] as DG.DataFrame);
    currentDfSheet.addRow((currentDf.columns as DG.ColumnList).names());
    for (let i = 0; i < currentDf.rowCount; i++) {
      currentDfSheet.addRow([...currentDf.row(i).cells].map((cell: DG.Cell) => cell.value));
    }
  });

  // UWAGA: very fragile solution
  const funcView = document.getElementsByClassName('grok-view grok-view-func ui-box')[0];
  const resultView = funcView.getElementsByClassName('d4-tab-host ui-box')[0];
  if (resultView) {
    const tabs = resultView.getElementsByClassName('d4-tab-header-stripe')[0].children;
    const selectedByUser = resultView.getElementsByClassName('d4-tab-header selected')[0];
    for (let i=0; i< tabs.length; i++) {
      (tabs[i] as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 150));

      const titleDivs = document.getElementsByClassName('grok-func-results-header');
      if (!titleDivs.length) continue;
      let skipNext = false;

      for (let i = 0; i < titleDivs.length; i++) {
        if (skipNext) {
          skipNext = false;
          continue;
        }
        const titleDiv = titleDivs[i];
        if (titleDiv.parentElement?.style.display === 'none') continue;

        const visibleTitle = titleDiv.firstChild?.textContent;
        if (visibleTitle === null || visibleTitle === undefined) continue;

        let imageId;
        let width = 0;
        let height = 0;

        const imageDiv = titleDiv.parentElement?.getElementsByClassName('grok-scripting-image-container')[0];
        if (imageDiv) {
          const regex = /background-image: url\(.*\)/;
          const urlAttr = imageDiv.outerHTML.match(regex)?.[0];
          if (!urlAttr) continue;
          const url = urlAttr.substring(23, urlAttr?.length - 2);
          const pic = (await (await fetch(url)).arrayBuffer());
          imageId = exportWorkbook.addImage({
            buffer: pic,
            extension: 'png',
          });
          width = imageDiv.clientWidth;
          height = imageDiv.clientHeight;
        } else {
          let plotDiv: ChildNode | null | undefined = titleDiv.nextSibling;

          if (plotDiv && !titleDiv.getElementsByClassName('svg-table').length && (plotDiv as HTMLElement).getAttribute('name') === 'viewer-Grid') continue;

          if (plotDiv && titleDiv.getElementsByClassName('svg-table').length) {
            if ((plotDiv as HTMLElement).getAttribute('name') === 'viewer-OutliersSelectionViewer') {
              plotDiv = titleDiv.parentElement?.nextSibling?.firstChild?.nextSibling;
              skipNext = true;
            }
            if ((plotDiv as HTMLElement).getAttribute('name') === 'viewer-Grid') {
              (titleDiv.getElementsByClassName('svg-table')[0] as HTMLElement).click();
              await new Promise((r) => setTimeout(r, 150));
              plotDiv = titleDiv.nextSibling;
            }
          }
          console.log(plotDiv);

          if (plotDiv) {
            const plot = (plotDiv as HTMLElement).getElementsByClassName('d4-layout-middle')[0];
            const canvas = await DG.HtmlUtils.renderToCanvas(plot as HTMLElement);
            const dataUrl = canvas.toDataURL('image/png');

            imageId = exportWorkbook.addImage({
              base64: dataUrl,
              extension: 'png',
            });
            width = canvas.width;
            height = canvas.height;
          }
        }

        if (imageId === null || imageId === undefined) continue;

        const worksheet = exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.INPUT)) ||
                          exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.OUTPUT));
        if (worksheet) {
          console.log(visibleTitle);
          console.log(realNames[visibleTitle]);
          const columnForImage = (((call
            .inputs[realNames[visibleTitle]] || call
            .outputs[realNames[visibleTitle]]) as DG.DataFrame)
            .columns as DG.ColumnList)
            .length + 1;
          worksheet.addImage(imageId, {
            tl: {col: columnForImage, row: 0},
            ext: {width, height},
          });
        } else {
          const newWorksheet = exportWorkbook.addWorksheet(`Output - Plot - ${visibleTitle}`);
          newWorksheet.addImage(imageId, {
            tl: {col: 0, row: 0},
            ext: {width, height},
          });
        }
      }
    }
    (selectedByUser as HTMLElement).click();
  } else {
    const resultView = document.getElementsByClassName('ui-panel grok-func-results')[0];
    const titleDivs = resultView.getElementsByClassName('grok-func-results-header');
    let skipNext = false;

    for (let i=0; i< titleDivs.length; i++) {
      if (skipNext) {
        skipNext = false;
        continue;
      }

      const titleDiv = titleDivs[i];
      const visibleTitle = titleDiv.firstChild?.textContent;
      if (visibleTitle === null || visibleTitle === undefined) continue;

      let imageId;
      let width = 0;
      let height = 0;

      const imageDiv = titleDiv.parentElement?.getElementsByClassName('grok-scripting-image-container')[0];
      if (imageDiv) {
        const regex = /background-image: url\(.*\)/;
        const urlAttr = imageDiv.outerHTML.match(regex)?.[0];
        if (!urlAttr) continue;
        const url = urlAttr.substring(23, urlAttr?.length - 2);
        const pic = (await (await fetch(url)).arrayBuffer());
        imageId = exportWorkbook.addImage({
          buffer: pic,
          extension: 'png',
        });
        width = imageDiv.clientWidth;
        height = imageDiv.clientHeight;
      } else {
        let plotDiv: ChildNode | null | undefined = titleDiv.nextSibling;
        if (plotDiv &&
            ((plotDiv as HTMLElement).getAttribute('name') === 'viewer-Grid' ||
            (plotDiv as HTMLElement).getAttribute('name') === 'viewer-OutliersSelectionViewer')) {
          if (titleDiv.getElementsByClassName('svg-table').length) {
            plotDiv = titleDiv.parentElement?.nextSibling?.firstChild?.nextSibling;
            skipNext = true;
          } else {
            continue;
          }
        }

        if (plotDiv) {
          const plot = (plotDiv as HTMLElement).getElementsByClassName('d4-layout-middle')[0];
          const canvas = await DG.HtmlUtils.renderToCanvas(plot as HTMLElement);
          const dataUrl = canvas.toDataURL('image/png');

          imageId = exportWorkbook.addImage({
            base64: dataUrl,
            extension: 'png',
          });
          width = canvas.width;
          height = canvas.height;
        }
      }

      if (imageId === null || imageId === undefined) continue;

      const worksheet = exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.OUTPUT));
      if (worksheet) {
        const columnForImage = ((call
          .outputs[realNames[visibleTitle]] as DG.DataFrame)
          .columns as DG.ColumnList)
          .length + 1;
        worksheet.addImage(imageId, {
          tl: {col: columnForImage, row: 0},
          ext: {width, height},
        });
      } else {
        const newWorksheet = exportWorkbook.addWorksheet(`Output - Plot - ${visibleTitle}`);
        newWorksheet.addImage(imageId, {
          tl: {col: 0, row: 0},
          ext: {width, height},
        });
      }
    }
  }

  if (scalarOutputs.length) {
    const outputScalarsSheet = exportWorkbook.addWorksheet('Output scalars');
    scalarOutputs.forEach((scalarOutput) => {
      outputScalarsSheet.addRow([scalarOutput.name, call.outputs[scalarOutput.name]]);
    });
  }

  exportWorkbook.xlsx.writeBuffer().then((data) => {
    const blob = new Blob([data], {type: BLOB_TYPE});
    saveAs(blob, call.func.name);
  });
}
