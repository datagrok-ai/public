import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Scripts {
  export async function fabArmExchange(meaConcentration: number, phOfFaeReaction: number, faeReactionTime: number, ionicStrenght: number, targetIgGConc: number, massOverload: number, P1MW: number, reactionVolume: number, partialPressureOfOxygen: number, doInitial: number, doBuffer: number, oxygenTransferRate: number, pressure: number, temperature: number, Lid: boolean, headSpace: number, dfVolume: number, filterateRateDF: number, phOfUfdfBuffer: number, phOfSecondUfdfBuffer: number, secondDfVolume: number, dfConcentration: number, filterateRateUf: number, holdTime: number, phDuringHoldTime: number, holdTime2: number, phDuringSecondHoldTime: number): Promise<number> {
    return await grok.functions.call('Compute:FabArmExchange', { meaConcentration, phOfFaeReaction, faeReactionTime, ionicStrenght, targetIgGConc, massOverload, P1MW, reactionVolume, partialPressureOfOxygen, doInitial, doBuffer, oxygenTransferRate, pressure, temperature, Lid, headSpace, dfVolume, filterateRateDF, phOfUfdfBuffer, phOfSecondUfdfBuffer, secondDfVolume, dfConcentration, filterateRateUf, holdTime, phDuringHoldTime, holdTime2, phDuringSecondHoldTime });
  }

  //Identify the beats in an ECG signal and compute the IBIs
  export async function intervalsFromECG(samplingFrequency: number, bpmMax: number, delta: number, k: number): Promise<any> {
    return await grok.functions.call('Compute:IntervalsFromECG', { samplingFrequency, bpmMax, delta, k });
  }

  //Predict the minimum filter size required for the separation of effluent from bioreactors, given the constraints of batch size and total batch time based on a training dataset of time vs filtrate volume for a given filter type, filter area and pressure
  export async function vmax(test_data: DG.DataFrame, test_area: number, vbatch: number, tbatch: number, sf: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Compute:Vmax', { test_data, test_area, vbatch, tbatch, sf });
  }
}

export namespace Funcs {
  export async function openModelFromFuncall(funccall: any): Promise<any> {
    return await grok.functions.call('Compute:OpenModelFromFuncall', { funccall });
  }

  //Creates an outliers selection viewer
  export async function outliersSelection(): Promise<any> {
    return await grok.functions.call('Compute:OutliersSelection', {});
  }

  export async function richFunctionViewEditor(call: any): Promise<any> {
    return await grok.functions.call('Compute:RichFunctionViewEditor', { call });
  }

  export async function pipelineStepEditor(call: any): Promise<any> {
    return await grok.functions.call('Compute:PipelineStepEditor', { call });
  }

  export async function renderPanel(func: any): Promise<any> {
    return await grok.functions.call('Compute:RenderPanel', { func });
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('Compute:Init', {});
  }

  export async function modelCatalog(): Promise<any> {
    return await grok.functions.call('Compute:ModelCatalog', {});
  }

  export async function modelCatalogTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('Compute:ModelCatalogTreeBrowser', { treeNode, browseView });
  }

  export async function customDataUploader(func: any): Promise<any> {
    return await grok.functions.call('Compute:CustomDataUploader', { func });
  }

  export async function customUploader(params: any): Promise<any> {
    return await grok.functions.call('Compute:CustomUploader', { params });
  }

  export async function customCustomizer(params: any): Promise<any> {
    return await grok.functions.call('Compute:CustomCustomizer', { params });
  }

  export async function simTimeValidator(params: any): Promise<any> {
    return await grok.functions.call('Compute:SimTimeValidator', { params });
  }

  export async function desiredTempValidator(params: any): Promise<any> {
    return await grok.functions.call('Compute:DesiredTempValidator', { params });
  }

  export async function initialTempValidator(params: any): Promise<any> {
    return await grok.functions.call('Compute:InitialTempValidator', { params });
  }

  export async function ambTempValidator(params: any): Promise<any> {
    return await grok.functions.call('Compute:AmbTempValidator', { params });
  }

  export async function heatCapValidator(params: any): Promise<any> {
    return await grok.functions.call('Compute:HeatCapValidator', { params });
  }

  export async function customStringInput(params: any): Promise<any> {
    return await grok.functions.call('Compute:CustomStringInput', { params });
  }

  export async function objectCoolingSelector(params: any): Promise<any> {
    return await grok.functions.call('Compute:ObjectCoolingSelector', { params });
  }
}
