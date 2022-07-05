import * as DG from 'datagrok-api/dg';
import { AvailableMetrics } from '@datagrok-libraries/ml/src/typed-metrics';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {BitArrayMetrics, StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import { Matrix } from '@datagrok-libraries/utils/src/type-declarations';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

export interface ISequenceSpaceResult {
  distance: Matrix;
  coordinates: DG.ColumnList;
}

export async function sequenceSpace(molColumn: DG.Column, methodName: string, similarityMetric: string,
    axes: string[], options?: any): Promise<ISequenceSpaceResult> {
    let preparedData: any;
    if (!(molColumn!.tags[DG.TAGS.UNITS] === 'HELM')) {
      const sep = molColumn.getTag('separator');
      const sepFinal = sep ? sep === '.' ? '\\\.' : sep: '-';  
      var regex = new RegExp(sepFinal, "g");
      if (Object.keys(AvailableMetrics['String']).includes(similarityMetric)) {
          preparedData = molColumn.toList().map((v) => v.replace(regex, '')) as string[];
      } else {
          preparedData = molColumn.toList().map((v) => v.replace(regex, '')) as string[];
      }
    } else {
      preparedData = molColumn.toList();
    }
    
    const sequenceSpaceResult = await reduceDimensinalityWithNormalization(
      preparedData,
      methodName,
      similarityMetric as StringMetrics|BitArrayMetrics,
      options);
    const cols: DG.Column[] = axes.map((name, index) => DG.Column.fromFloat32Array(name, sequenceSpaceResult.embedding[index]))
    return {distance: sequenceSpaceResult.distance, coordinates: new DG.ColumnList(cols)};
  }


export function getEmbeddingColsNames(df: DG.DataFrame){
    const axes = ['Embed_X', 'Embed_Y'];
    const colNameInd = df.columns.names().filter((it) => it.includes(axes[0])).length + 1;
    return axes.map((it) => `${it}_${colNameInd}`);
  }