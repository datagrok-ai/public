import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function barCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:BarCellRenderer', {});
  }

  export async function sparklineCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:SparklineCellRenderer', {});
  }

  export async function barchartCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:BarchartCellRenderer', {});
  }

  export async function piechartCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:PiechartCellRenderer', {});
  }

  export async function radarCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:RadarCellRenderer', {});
  }

  export async function smartFormCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:SmartFormCellRenderer', {});
  }

  //Adds a sparkline column for the selected columns
  export async function summarizeColumns(columns: any): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:SummarizeColumns', { columns });
  }

  //Adds a 'form' column for the selected columns
  export async function addFormColumn(columns: any): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:AddFormColumn', { columns });
  }

  export async function testUnitsKgCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:TestUnitsKgCellRenderer', {});
  }

  export async function testUnitsTonCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:TestUnitsTonCellRenderer', {});
  }

  export async function addPinnedColumn(gridCol: any): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:AddPinnedColumn', { gridCol });
  }

  export async function demoTestUnitsCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:DemoTestUnitsCellRenderer', {});
  }

  export async function autoPowerGrid(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:AutoPowerGrid', {});
  }

  //Forms viewer
  export async function formsViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:FormsViewer', {});
  }

  //Image content
  export async function imgContent(imageUrl: string): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:ImgContent', { imageUrl });
  }

  export async function demoCellTypes(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:DemoCellTypes', {});
  }

  export async function scWebGPURender(sc: any, show: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:ScWebGPURender', { sc, show });
  }

  export async function scWebGPUPointHitTest(sc: any, pt: any): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:ScWebGPUPointHitTest', { sc, pt });
  }

  export async function isWebGPUAvailable(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:IsWebGPUAvailable', {});
  }

  export async function isWebGPURenderValid(sc: any): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:IsWebGPURenderValid', { sc });
  }

  export async function binaryImageCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:BinaryImageCellRenderer', {});
  }

  export async function hyperlinkCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:HyperlinkCellRenderer', {});
  }

  export async function imageCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:ImageCellRenderer', {});
  }

  export async function multiChoiceCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:MultiChoiceCellRenderer', {});
  }

  export async function tagsCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:TagsCellRenderer', {});
  }

  export async function htmlTestCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:HtmlTestCellRenderer', {});
  }

  export async function scatterPlotCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/power-grid:ScatterPlotCellRenderer', {});
  }
}
