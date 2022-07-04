/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SequenceAlignment, Aligned} from './seq_align';

export const _package = new DG.Package();

import {WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {VdRegionsViewer} from './viewers/vd-regions-viewer';
import {runKalign, testMSAEnoughMemory} from './utils/multiple-sequence-alignment';

//name: sequenceAlignment
//input: string alignType {choices: ['Local alignment', 'Global alignment']}
// eslint-disable-next-line max-len
//input: string alignTable {choices: ['AUTO', 'NUCLEOTIDES', 'BLOSUM45', 'BLOSUM50','BLOSUM62','BLOSUM80','BLOSUM90','PAM30','PAM70','PAM250','SCHNEIDER','TRANS']}
//input: double gap
//input: string seq1
//input: string seq2
//output: object res
export function sequenceAlignment(alignType: string, alignTable: string, gap: number, seq1: string, seq2: string) {
  const toAlign = new SequenceAlignment(seq1, seq2, gap, alignTable);
  const res = alignType == 'Local alignment' ? toAlign.smithWaterman() : toAlign.needlemanWunch();
  return res;
}

//name: WebLogo
//description: WebLogo viewer
//tags: viewer, panel
//output: viewer result
export function webLogoViewer() {
  return new WebLogo();
}

//name: VdRegions
//description: V-Domain regions viewer
//tags: viewer, panel
//output: viewer result
export function vdRegionViewer() {
  return new VdRegionsViewer();
}

//top-menu: Bio | Activity Cliffs...
//name: Activity Cliffs
//description: detect activity cliffs
//input: dataframe df [Input data table]
//input: column smiles {type:categorical; semType: MacroMolecule}
//input: column activities
//input: double similarity = 80 [Similarity cutoff]
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
export async function activityCliffs(df: DG.DataFrame, smiles: DG.Column, activities: DG.Column,
  similarity: number, methodName: string) : Promise<void> {
}

//top-menu: Bio | Sequence Space...
//name: Sequence Space
//input: dataframe table
//input: column smiles { semType: MacroMolecule }
//input: string methodName { choices:["UMAP", "t-SNE", "SPE", "pSPE", "OriginalSPE"] }
//input: string similarityMetric { choices:["Tanimoto", "Asymmetric", "Cosine", "Sokal"] }
//input: bool plotEmbeddings = true
export async function chemSpaceTopMenu(table: DG.DataFrame, smiles: DG.Column, methodName: string,
  similarityMetric: string = 'Tanimoto', plotEmbeddings: boolean) : Promise<void> {
};

//top-menu: Bio | MSA...
//name: MSA
//input: dataframe table
//input: column sequence { semType: MacroMolecule }
export async function multipleSequenceAlignmentAny(table: DG.DataFrame, col: DG.Column): Promise<void> {
  const msaCol = await runKalign(col, false);
  table.columns.add(msaCol);
}