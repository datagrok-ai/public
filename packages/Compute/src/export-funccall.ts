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

  const addToFile = (imageId: number, visibleTitle: string, width: number, height: number, plotCount?: number) => {
    if (imageId === null || imageId === undefined) return;

    const worksheet = exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.INPUT)) ||
                      exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.OUTPUT));
    if (worksheet) {
      const columnForImage = (
        ((call.inputs[realNames[visibleTitle]] ||
          call.outputs[realNames[visibleTitle]]) as DG.DataFrame)
          .columns as DG.ColumnList)
        .length + 1 + (plotCount ? plotCount * (Math.ceil(width / 64.0) + 1): 0);
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
  };

  const exportImage = async (imageDiv: Element, visibleTitle: string) => {
    const regex = /background-image: url\(.*\)/;
    const urlAttr = imageDiv.outerHTML.match(regex)?.[0];
    if (!urlAttr) return;
    const url = urlAttr.substring(23, urlAttr?.length - 2);
    const pic = (await(await fetch(url)).arrayBuffer());
    const imageId = exportWorkbook.addImage({
      buffer: pic,
      extension: 'png',
    });
    addToFile(imageId, visibleTitle, imageDiv.clientWidth, imageDiv.clientHeight);
  };

  const exportPlot = async (plotDiv: ChildNode | null | undefined,
    blockDiv: Element, visibleTitle: string, plotCount? : number) => {
    const titleDiv = blockDiv.firstChild as HTMLElement;
    if (!titleDiv) return;

    if (plotDiv && !titleDiv.getElementsByClassName('svg-table').length &&
         (plotDiv as HTMLElement).getAttribute('name') === 'viewer-Grid') return;

    if (plotDiv && titleDiv.getElementsByClassName('svg-table').length) {
      if ((plotDiv as HTMLElement).getAttribute('name') === 'viewer-OutliersSelectionViewer') {
        return;
      }
      if ((plotDiv as HTMLElement).getAttribute('name') === 'viewer-Grid') {
        (titleDiv.getElementsByClassName('svg-table')[0] as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 150));
        plotDiv = titleDiv.nextSibling;
      }
    }

    if (plotDiv) {
      const plot = (plotDiv as HTMLElement).getElementsByClassName('d4-layout-middle')[0];
      const canvas = await DG.HtmlUtils.renderToCanvas(plot as HTMLElement);
      const dataUrl = canvas.toDataURL('image/png');

      const imageId = exportWorkbook.addImage({
        base64: dataUrl,
        extension: 'png',
      });
      addToFile(imageId, visibleTitle, canvas.width, canvas.height, plotCount);
    }
  };

  const exportTab = async (tab: HTMLElement) => {
    const blockDivs = document.getElementsByClassName('ui-div ui-block');
    let visibleTitle: string = '';
    let plotCount: number = 0;

    for (let i = 0; i < blockDivs.length; i++) {
      const blockDiv = blockDivs[i];

      if ((blockDiv as HTMLElement).style.display === 'none') continue;

      if (blockDiv.firstChild?.firstChild?.textContent) {
        visibleTitle = blockDiv.firstChild?.firstChild?.textContent;
        plotCount = 0;
      } else {
        plotCount++;
      }

      if ((blockDiv.firstChild?.nextSibling as HTMLElement).getAttribute('name') === 'viewer-OutliersSelectionViewer') {
        plotCount--;
      }

      const imageDiv = blockDiv.getElementsByClassName('grok-scripting-image-container')[0];
      if (imageDiv) {
        await exportImage(imageDiv, visibleTitle);
      } else {
        const plotDiv: ChildNode | null | undefined = blockDiv.firstChild?.nextSibling;
        await exportPlot(plotDiv, blockDiv, visibleTitle, plotCount);
      }
    };
  };

  // UWAGA: very fragile solution
  const funcView = document.getElementsByClassName('grok-view grok-view-func ui-box')[0];
  const resultView = funcView.getElementsByClassName('d4-tab-host ui-box')[0];
  if (resultView) {
    const tabs = resultView.getElementsByClassName('d4-tab-header-stripe')[0].children;
    const selectedByUser = resultView.getElementsByClassName('d4-tab-header selected')[0];
    for (let i=0; i< tabs.length; i++) {
      const tab = tabs[i] as HTMLElement;
      tab.click();
      await new Promise((r) => setTimeout(r, 150));
      await exportTab(tab);
    }
    (selectedByUser as HTMLElement).click();
  } else {
    const resultView = document.getElementsByClassName('ui-panel grok-func-results')[0];
    await exportTab(resultView as HTMLElement);
  }

  if (scalarOutputs.length) {
    const outputScalarsSheet = exportWorkbook.addWorksheet('Output scalars');
    scalarOutputs.forEach((scalarOutput) => {
      outputScalarsSheet.addRow([scalarOutput.name, call.outputs[scalarOutput.name]]);
    });
  }

  const metadataSheet = exportWorkbook.addWorksheet('Function run metadata');
  metadataSheet.addRow(['Author', `${call.func.author.lastName} ${call.func.author.firstName}`]);
  metadataSheet.addRow(['Login', call.func.author.toString()]);

  exportWorkbook.xlsx.writeBuffer().then((data) => {
    const blob = new Blob([data], {type: BLOB_TYPE});
    saveAs(blob, `${call.func.name} - ${new Date().toLocaleString()}`);
  });
}
