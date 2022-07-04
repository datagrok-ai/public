/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SequenceAlignment, Aligned} from './seq_align';

export const _package = new DG.Package();

import {WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {VdRegionsViewer} from './viewers/vd-regions-viewer';
import {runKalign, testMSAEnoughMemory} from './utils/multiple-sequence-alignment';
import {TableView} from 'datagrok-api/dg';

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
//input: column smiles {type:categorical; semType: Macromolecule}
//input: column activities
//input: double similarity = 80 [Similarity cutoff]
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
export async function activityCliffs(df: DG.DataFrame, smiles: DG.Column, activities: DG.Column,
  similarity: number, methodName: string): Promise<void> {
}

//top-menu: Bio | Sequence Space...
//name: Sequence Space
//input: dataframe table
//input: column smiles { semType: Macromolecule }
//input: string methodName { choices:["UMAP", "t-SNE", "SPE", "pSPE", "OriginalSPE"] }
//input: string similarityMetric { choices:["Tanimoto", "Asymmetric", "Cosine", "Sokal"] }
//input: bool plotEmbeddings = true
export async function chemSpaceTopMenu(table: DG.DataFrame, smiles: DG.Column, methodName: string,
  similarityMetric: string = 'Tanimoto', plotEmbeddings: boolean): Promise<void> {
};

//top-menu: Bio | MSA...
//name: MSA
//input: dataframe table
//input: column sequence { semType: Macromolecule }
export async function multipleSequenceAlignmentAny(table: DG.DataFrame, col: DG.Column): Promise<void> {
  const msaCol = await runKalign(col, false);
  table.columns.add(msaCol);
}

//name: Composition Analysis
//top-menu: Bio | Composition Analysis
//output: viewer result
export async function compositionAnalysis(): Promise<void> {
  const col = grok.shell.t.columns.bySemType('Macromolecule');//DG.SEMTYPE.MACROMOLECULE);
  if (col === null) {
    grok.shell.error('Current table does not contain sequences');
    return;
  }

  const wl = await col.dataFrame.plot.fromType('WebLogo', {});

  for (const v of grok.shell.views) {
    if (v instanceof TableView && (v as DG.TableView).dataFrame.name === col.dataFrame.name) {
      (v as DG.TableView).dockManager.dock(wl.root, 'down');
      break;
    }
  }
}

//name: importFasta
//description: Opens FASTA file
//tags: file-handler
//meta.ext: fasta, fna, ffn, faa, frn, fa
//input: string content
//output: list tables
export function importFasta(content: string): DG.DataFrame [] {
  const regex = /^>(.*)$/gm;
  let match;
  const descriptions = [];
  const sequences = [];
  let index = 0;
  while (match = regex.exec(content)) {
    descriptions.push(content.substring(match.index + 1, regex.lastIndex));
    if (index !== 0)
      sequences.push(content.substring(index, regex.lastIndex));
    index = regex.lastIndex + 1;
  }
  sequences.push(content.substring(index));
  return [DG.DataFrame.fromColumns([
    DG.Column.fromStrings('description', descriptions),
    DG.Column.fromStrings('sequence', sequences)
  ])];
}
