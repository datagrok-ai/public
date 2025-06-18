import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace queries {
  export async function getPlates(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetPlates', {});
  }

  //Get all well level properties (either used in a well or specified in a template)
  export async function getWellLevelProperties(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetWellLevelProperties', {});
  }

  //Get all plate level properties (either used in a plate or specified in a template)
  export async function getPlateLevelProperties(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetPlateLevelProperties', {});
  }

  export async function getPropertyNames(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetPropertyNames', {});
  }

  export async function getPlateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetPlateTypes', {});
  }

  export async function getPlateTemplates(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetPlateTemplates', {});
  }

  export async function getWellRoles(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetWellRoles', {});
  }

  export async function getWellValuesByBarcode(barcode: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetWellValuesByBarcode', { barcode });
  }

  export async function getWellValuesById(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetWellValuesById', { id });
  }

  export async function getAllowedValues(propertyName: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetAllowedValues', { propertyName });
  }

  export async function getUniquePlatePropertyValues(propertyId: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetUniquePlatePropertyValues', { propertyId });
  }

  export async function getUniqueWellPropertyValues(propertyId: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetUniqueWellPropertyValues', { propertyId });
  }

  export async function createProperty(propertyName: string, valueType: string): Promise<number> {
    return await grok.data.query('@datagrok/curves:CreateProperty', { propertyName, valueType });
  }

  export async function createTemplate(name: string, description: string): Promise<number> {
    return await grok.data.query('@datagrok/curves:CreateTemplate', { name, description });
  }

  export async function getTemplateWellProperties(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetTemplateWellProperties', {});
  }

  export async function getTemplatePlateProperties(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/curves:GetTemplatePlateProperties', {});
  }
}

export namespace funcs {
  export async function platesAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('@datagrok/curves:PlatesAppTreeBrowser', { treeNode, browseView });
  }

  export async function fitChartCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/curves:FitChartCellRenderer', {});
  }

  //A viewer that superimposes multiple in-cell curves on one chart
  export async function multiCurveViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/curves:MultiCurveViewer', {});
  }

  //Curve fitting is the process of constructing a curve, or mathematical function, that has the best fit to a series of data points
  export async function curveFitDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/curves:CurveFitDemo', {});
  }

  //Assasy plates with concentration, layout and readout data
  export async function assayPlatesDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/curves:AssayPlatesDemo', {});
  }

  export async function initCurves(): Promise<any> {
    return await grok.functions.call('@datagrok/curves:InitCurves', {});
  }

  export async function addStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, seriesNumber: number): Promise<any> {
    return await grok.functions.call('@datagrok/curves:AddStatisticsColumn', { table, colName, propName, seriesNumber });
  }

  export async function addAggrStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, aggrType: string): Promise<any> {
    return await grok.functions.call('@datagrok/curves:AddAggrStatisticsColumn', { table, colName, propName, aggrType });
  }

  export async function platesFolderPreview(folder: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/curves:PlatesFolderPreview', { folder });
  }

  export async function previewPlate(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/curves:PreviewPlate', { file });
  }

  export async function importPlate(fileContent: string): Promise<any> {
    return await grok.functions.call('@datagrok/curves:ImportPlate', { fileContent });
  }

  export async function importPlateXlsx(fileContent: any): Promise<any> {
    return await grok.functions.call('@datagrok/curves:ImportPlateXlsx', { fileContent });
  }

  export async function previewPlateXlsx(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/curves:PreviewPlateXlsx', { file });
  }

  export async function checkExcelIsPlate(content: any): Promise<any> {
    return await grok.functions.call('@datagrok/curves:CheckExcelIsPlate', { content });
  }

  export async function checkFileIsPlate(content: string): Promise<any> {
    return await grok.functions.call('@datagrok/curves:CheckFileIsPlate', { content });
  }

  export async function platesApp(): Promise<any> {
    return await grok.functions.call('@datagrok/curves:PlatesApp', {});
  }

  export async function getPlateByBarcode(barcode: string): Promise<any> {
    return await grok.functions.call('@datagrok/curves:GetPlateByBarcode', { barcode });
  }

  export async function createDummyPlateData(): Promise<any> {
    return await grok.functions.call('@datagrok/curves:CreateDummyPlateData', {});
  }

  export async function plateGridCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/curves:PlateGridCellRenderer', {});
  }
}
