import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import ExcelJS from 'exceljs';
import {saveAs} from 'file-saver';

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

  dfInputs.forEach((dfInput) => {
    const currentDfSheet = exportWorkbook.addWorksheet(`Input - ${dfInput.name}`);
    const currentDf = (call.inputs[dfInput.name] as DG.DataFrame);
    currentDfSheet.addRow((currentDf.columns as DG.ColumnList).names());
    for (let i = 0; i < currentDf.rowCount; i++) {
      currentDfSheet.addRow([...currentDf.row(i).cells].map((cell: DG.Cell) => cell.value));
    }
  });

  const inputScalarsSheet = exportWorkbook.addWorksheet('Input scalars');
  scalarInputs.forEach((scalarInput) => {
    inputScalarsSheet.addRow([scalarInput.name, call.inputs[scalarInput.name]]);
  });

  dfOutputs.forEach((dfOutput) => {
    const currentDfSheet = exportWorkbook.addWorksheet(`Output - ${dfOutput.name}`);
    const currentDf = (call.outputs[dfOutput.name] as DG.DataFrame);
    currentDfSheet.addRow((currentDf.columns as DG.ColumnList).names());
    for (let i = 0; i < currentDf.rowCount; i++) {
      currentDfSheet.addRow([...currentDf.row(i).cells].map((cell: DG.Cell) => cell.value));
    }
  });

  // UWAGA: very fragile solution
  const funcView = document.getElementsByClassName('grok-view grok-view-func ui-box')[0];
  const resultView = funcView.getElementsByClassName('d4-tab-host ui-box')[0];
  const tabs = resultView.getElementsByClassName('d4-tab-header-stripe')[0].children;
  const selectedByUser = resultView.getElementsByClassName('d4-tab-header selected')[0];
  for (let i=0; i< tabs.length; i++) {
    (tabs[i] as HTMLElement).click();
    await new Promise((r) => setTimeout(r, 100));

    const titleDivs = document.getElementsByClassName('grok-func-results-header');
    if (!titleDivs.length) continue;

    for (let i = 0; i < titleDivs.length; i++) {
      const titleDiv = titleDivs[i];

      const title = titleDiv.firstChild?.textContent;
      if (!title) continue;

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
        const plot = titleDiv.nextSibling;
        // if plot does not exist or it is only grid, skip it
        if (!plot || (plot.firstChild as HTMLElement).getAttribute('name') === 'viewer-Grid') continue;

        if (plot) {
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

      const worksheet = exportWorkbook.getWorksheet(`Output - ${title}`);
      if (worksheet) {
        const columnForImage = ((call.outputs[title] as DG.DataFrame).columns as DG.ColumnList).length + 1;
        worksheet.addImage(imageId, {
          tl: {col: columnForImage, row: 0},
          ext: {width, height},
        });
      } else {
        const newWorksheet = exportWorkbook.addWorksheet(`Output - Plot - ${title}`);
        newWorksheet.addImage(imageId, {
          tl: {col: 0, row: 0},
          ext: {width, height},
        });
      }
    }
  }
  (selectedByUser as HTMLElement).click();

  const outputScalarsSheet = exportWorkbook.addWorksheet('Output scalars');
  scalarOutputs.forEach((scalarOutput) => {
    outputScalarsSheet.addRow([scalarOutput.name, call.outputs[scalarOutput.name]]);
  });

  exportWorkbook.xlsx.writeBuffer().then((data) => {
    const blob = new Blob([data], {type: BLOB_TYPE});
    saveAs(blob, call.func.name);
  });
}
