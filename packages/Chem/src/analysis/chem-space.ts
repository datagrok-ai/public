import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// TODO: clean up this module
import {chemGetFingerprints} from '../chem-searches';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {Fingerprint} from '../utils/chem-common';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {IReduceDimensionalityResult} from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {ISequenceSpaceParams, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import { IDimReductionParam, ITSNEOptions, IUMAPOptions, TSNEOptions, T_SNE, UMAP, UMAPOptions } from '@datagrok-libraries/ml/src/reduce-dimensionality';


export async function chemSpace(spaceParams: ISequenceSpaceParams): Promise<ISequenceSpaceResult> {
  const fpColumn = await chemGetFingerprints(spaceParams.seqCol, Fingerprint.Morgan);
  const chemSpaceResult: IReduceDimensionalityResult = await reduceDimensinalityWithNormalization(
    fpColumn,
    spaceParams.methodName,
    spaceParams.similarityMetric as BitArrayMetrics,
    spaceParams.options);
  const cols: DG.Column[] = spaceParams.embedAxesNames.map((name: string, index: number) => DG.Column.fromFloat32Array(name, chemSpaceResult.embedding[index]));
  return {distance: chemSpaceResult.distance, coordinates: new DG.ColumnList(cols)};
}

export function getEmbeddingColsNames(df: DG.DataFrame) {
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it: string) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}

export class ChemSpaceFuncEditor {
  tableInput: DG.InputBase;
  molColInput: DG.InputBase;
  methodInput: DG.InputBase;
  similarityMetricInput: DG.InputBase;
  plotEmbeddingsInput: DG.InputBase;
  funcParamsDiv: HTMLDivElement;
  methodsParams: {[key: string]: UMAPOptions | TSNEOptions} = {
    [UMAP]: new UMAPOptions(),
    [T_SNE]: new TSNEOptions()
  };

  get algorithmOptions(): IUMAPOptions | ITSNEOptions {
    const algorithmParams: UMAPOptions | TSNEOptions = this.methodsParams[this.methodInput.value!];
    const options: any = {};
    Object.keys(algorithmParams).forEach((key: string) => {
      if ((algorithmParams as any)[key].value != null) 
        options[key] = (algorithmParams as any)[key].value;
    });
    return options;
  }

  get funcParams(): any {
    return {table: this.tableInput.value!, molecules: this.molColInput.value!, methodName: this.methodInput.value!,
      similarityMetric: this.similarityMetricInput.value!, plotEmbeddings: this.plotEmbeddingsInput.value!, 
      options: this.algorithmOptions}
  }

  get paramsUI(): HTMLDivElement{
    return this.funcParamsDiv;
  }

  constructor(){
    this.tableInput = ui.tableInput('Table', grok.shell.tv.dataFrame, undefined, () => {
      ui.empty(moleculesColDiv);
      this.molColInput = ui.columnInput('Molecules', this.tableInput.value!, this.tableInput.value!.columns.bySemType(DG.SEMTYPE.MOLECULE));
      this.molColInput.input.style.width = '160px';
      moleculesColDiv.append(this.molColInput.root);
    });

    this.molColInput = ui.columnInput('Molecules', this.tableInput.value!, this.tableInput.value!.columns.bySemType(DG.SEMTYPE.MOLECULE));
    this.molColInput.input.style.width = '160px';
    const moleculesColDiv: HTMLDivElement = ui.div();
    moleculesColDiv.append(this.molColInput.root);

    this.methodInput = ui.choiceInput('Method Name', UMAP, [UMAP, T_SNE], () => {
      this.createAlgorithmSettingsDiv(methodSettingsDiv, this.methodsParams[this.methodInput.value!]);
    });
    this.methodInput.input.style.width = '160px';
    this.methodInput.captionLabel.style.width = '122px';

    const methodSettingsIcon = ui.icons.settings(()=> {
      settingsOpened = !settingsOpened;
      if (!settingsOpened)
        ui.empty(methodSettingsDiv);
      else 
        this.createAlgorithmSettingsDiv(methodSettingsDiv, this.methodsParams[this.methodInput.value!]);
    }, 'Modify methods parameters');
    methodSettingsIcon.style.lineHeight = '1.7';
    methodSettingsIcon.style.fontSize = '18px';
    const methodSettingsDiv = ui.form([]);
    let settingsOpened = false;

    this.similarityMetricInput = ui.choiceInput('Similarity metric', 'Tanimoto', ['Tanimoto', 'Asymmetric', 'Cosine', 'Sokal']);
    
    this.plotEmbeddingsInput = ui.boolInput('Plot Embeddings', true);
    this.plotEmbeddingsInput.captionLabel.style.width = '130px';
    
    this.funcParamsDiv = ui.form([
      this.tableInput,
      //@ts-ignore
      moleculesColDiv,
      //@ts-ignore
      ui.divH([methodSettingsIcon, this.methodInput]),
      //@ts-ignore
      methodSettingsDiv,
      this.similarityMetricInput,
      this.plotEmbeddingsInput
    ])
  }

  createAlgorithmSettingsDiv(paramsForm: HTMLDivElement, params: UMAPOptions | TSNEOptions) {
    ui.empty(paramsForm);
    Object.keys(params).forEach((it: any) => {
      const param: IDimReductionParam = (params as any)[it];
      const input = ui.floatInput(param.uiName, param.value, () => {
        param.value = input.value;
      });
      ui.tooltip.bind(input.root, param.tooltip);
      paramsForm.append(input.root);
    });
    return paramsForm;
  }
}
