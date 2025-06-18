import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  export async function fabArmExchange(meaConcentration: number, phOfFaeReaction: number, faeReactionTime: number, ionicStrenght: number, targetIgGConc: number, massOverload: number, P1MW: number, reactionVolume: number, partialPressureOfOxygen: number, doInitial: number, doBuffer: number, oxygenTransferRate: number, pressure: number, temperature: number, Lid: boolean, headSpace: number, dfVolume: number, filterateRateDF: number, phOfUfdfBuffer: number, phOfSecondUfdfBuffer: number, secondDfVolume: number, dfConcentration: number, filterateRateUf: number, holdTime: number, phDuringHoldTime: number, holdTime2: number, phDuringSecondHoldTime: number): Promise<number> {
    return await grok.functions.call('@datagrok/compute:FabArmExchange', { meaConcentration, phOfFaeReaction, faeReactionTime, ionicStrenght, targetIgGConc, massOverload, P1MW, reactionVolume, partialPressureOfOxygen, doInitial, doBuffer, oxygenTransferRate, pressure, temperature, Lid, headSpace, dfVolume, filterateRateDF, phOfUfdfBuffer, phOfSecondUfdfBuffer, secondDfVolume, dfConcentration, filterateRateUf, holdTime, phDuringHoldTime, holdTime2, phDuringSecondHoldTime });
  }

  //Identify the beats in an ECG signal and compute the IBIs
  export async function intervalsFromECG(samplingFrequency: number, bpmMax: number, delta: number, k: number): Promise<any> {
    return await grok.functions.call('@datagrok/compute:IntervalsFromECG', { samplingFrequency, bpmMax, delta, k });
  }

  //Predict the minimum filter size required for the separation of effluent from bioreactors, given the constraints of batch size and total batch time based on a training dataset of time vs filtrate volume for a given filter type, filter area and pressure
  export async function vmax(test_data: DG.DataFrame, test_area: number, vbatch: number, tbatch: number, sf: number): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/compute:Vmax', { test_data, test_area, vbatch, tbatch, sf });
  }
}

export namespace funcs {
  export async function openModelFromFuncall(funccall: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:OpenModelFromFuncall', { funccall });
  }

  //Creates an outliers selection viewer
  export async function outliersSelection(): Promise<any> {
    return await grok.functions.call('@datagrok/compute:OutliersSelection', {});
  }

  export async function richFunctionViewEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:RichFunctionViewEditor', { call });
  }

  export async function pipelineStepEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:PipelineStepEditor', { call });
  }

  export async function renderPanel(func: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:RenderPanel', { func });
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('@datagrok/compute:Init', {});
  }

  export async function modelCatalog(): Promise<any> {
    return await grok.functions.call('@datagrok/compute:ModelCatalog', {});
  }

  export async function modelCatalogTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('@datagrok/compute:ModelCatalogTreeBrowser', { treeNode, browseView });
  }

  export async function customDataUploader(func: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:CustomDataUploader', { func });
  }

  export async function customUploader(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:CustomUploader', { params });
  }

  export async function customCustomizer(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:CustomCustomizer', { params });
  }

  export async function simTimeValidator(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:SimTimeValidator', { params });
  }

  export async function desiredTempValidator(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:DesiredTempValidator', { params });
  }

  export async function initialTempValidator(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:InitialTempValidator', { params });
  }

  export async function ambTempValidator(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:AmbTempValidator', { params });
  }

  export async function heatCapValidator(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:HeatCapValidator', { params });
  }

  export async function customStringInput(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:CustomStringInput', { params });
  }

  export async function objectCoolingSelector(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute:ObjectCoolingSelector', { params });
  }
}
