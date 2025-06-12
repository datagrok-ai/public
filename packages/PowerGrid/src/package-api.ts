import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function barCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:BarCellRenderer', {});
  }

  export async function sparklineCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:SparklineCellRenderer', {});
  }

  export async function barchartCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:BarchartCellRenderer', {});
  }

  export async function piechartCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:PiechartCellRenderer', {});
  }

  export async function radarCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:RadarCellRenderer', {});
  }

  export async function smartFormCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:SmartFormCellRenderer', {});
  }

  //Adds a sparkline column for the selected columns
  export async function summarizeColumns(columns: any): Promise<any> {
    return await grok.functions.call('PowerGrid:SummarizeColumns', { columns });
  }

  //Adds a 'form' column for the selected columns
  export async function addFormColumn(columns: any): Promise<any> {
    return await grok.functions.call('PowerGrid:AddFormColumn', { columns });
  }

  export async function testUnitsKgCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:TestUnitsKgCellRenderer', {});
  }

  export async function testUnitsTonCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:TestUnitsTonCellRenderer', {});
  }

  export async function addPinnedColumn(gridCol: any): Promise<any> {
    return await grok.functions.call('PowerGrid:AddPinnedColumn', { gridCol });
  }

  export async function demoTestUnitsCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:DemoTestUnitsCellRenderer', {});
  }

  export async function autoPowerGrid(): Promise<any> {
    return await grok.functions.call('PowerGrid:AutoPowerGrid', {});
  }

  //Forms viewer
  export async function formsViewer(): Promise<any> {
    return await grok.functions.call('PowerGrid:FormsViewer', {});
  }

  //Image content
  export async function imgContent(imageUrl: string): Promise<any> {
    return await grok.functions.call('PowerGrid:ImgContent', { imageUrl });
  }

  export async function demoCellTypes(): Promise<any> {
    return await grok.functions.call('PowerGrid:DemoCellTypes', {});
  }

  export async function scWebGPURender(sc: any, show: boolean): Promise<any> {
    return await grok.functions.call('PowerGrid:ScWebGPURender', { sc, show });
  }

  export async function scWebGPUPointHitTest(sc: any, pt: any): Promise<any> {
    return await grok.functions.call('PowerGrid:ScWebGPUPointHitTest', { sc, pt });
  }

  export async function isWebGPUAvailable(): Promise<any> {
    return await grok.functions.call('PowerGrid:IsWebGPUAvailable', {});
  }

  export async function isWebGPURenderValid(sc: any): Promise<any> {
    return await grok.functions.call('PowerGrid:IsWebGPURenderValid', { sc });
  }

  export async function binaryImageCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:BinaryImageCellRenderer', {});
  }

  export async function hyperlinkCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:HyperlinkCellRenderer', {});
  }

  export async function imageCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:ImageCellRenderer', {});
  }

  export async function multiChoiceCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:MultiChoiceCellRenderer', {});
  }

  export async function tagsCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:TagsCellRenderer', {});
  }

  export async function htmlTestCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:HtmlTestCellRenderer', {});
  }

  export async function scatterPlotCellRenderer(): Promise<any> {
    return await grok.functions.call('PowerGrid:ScatterPlotCellRenderer', {});
  }
}
