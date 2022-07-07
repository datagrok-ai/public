/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

import {WebLogo, SeqColStats} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {VdRegionsViewer} from './viewers/vd-regions-viewer';
import {runKalign, testMSAEnoughMemory} from './utils/multiple-sequence-alignment';
import {SequenceAlignment, Aligned} from './seq_align';
import {Nucleotides} from '@datagrok-libraries/bio/src/nucleotides';
import {Aminoacids} from '@datagrok-libraries/bio/src/aminoacids';
import {convert} from './utils/convert';
import {getEmbeddingColsNames, sequenceSpace} from './utils/sequence-space';
import {AvailableMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {getActivityCliffs} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {sequenceGetSimilarities, drawTooltip} from './utils/sequence-activity-cliffs';
import { getMolfilesFromSeq, HELM_CORE_LIB_FILENAME } from './utils/utils';

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

//top-menu: Bio | Sequence Activity Cliffs...
//name: Sequence Activity Cliffs
//description: detect activity cliffs
//input: dataframe df [Input data table]
//input: column sequence {semType: Macromolecule}
//input: column activities
//input: double similarity = 80 [Similarity cutoff]
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
export async function activityCliffs(df: DG.DataFrame, sequence: DG.Column, activities: DG.Column,
  similarity: number, methodName: string): Promise<void> {
  const axesNames = getEmbeddingColsNames(df);
  const options = {
    'SPE': {cycles: 2000, lambda: 1.0, dlambda: 0.0005},
  };
  const units = sequence!.tags[DG.TAGS.UNITS];
  await getActivityCliffs(
    df,
    sequence,
    axesNames,
    activities,
    similarity,
    'Levenshtein',
    methodName,
    DG.SEMTYPE.MACROMOLECULE,
    units,
    sequenceSpace,
    sequenceGetSimilarities,
    drawTooltip,
    (options as any)[methodName]);
}

//top-menu: Bio | Sequence Space...
//name: Sequence Space
//input: dataframe table
//input: column macroMolecule { semType: Macromolecule }
//input: string methodName { choices:["UMAP", "t-SNE", "SPE"] }
//input: string similarityMetric { choices:["Levenshtein", "Tanimoto"] }
//input: bool plotEmbeddings = true
export async function sequenceSpaceTopMenu(table: DG.DataFrame, macroMolecule: DG.Column, methodName: string,
  similarityMetric: string = 'Levenshtein', plotEmbeddings: boolean): Promise<void> {
  const embedColsNames = getEmbeddingColsNames(table);
  const chemSpaceParams = {
    seqCol: macroMolecule,
    methodName: methodName,
    similarityMetric: similarityMetric,
    embedAxesNames: embedColsNames
  };
  const sequenceSpaceRes = await sequenceSpace(chemSpaceParams);
  const embeddings = sequenceSpaceRes.coordinates;
  for (const col of embeddings)
    table.columns.add(col);
  if (plotEmbeddings) {
    for (const v of grok.shell.views) {
      if (v.name === table.name)
        (v as DG.TableView).scatterPlot({x: embedColsNames[0], y: embedColsNames[1]});
    }
  }
};

//top-menu: Bio | Molfiles From HELM...
//name: Molfiles From HELM
//description: returns molfiles for each monomer from HELM library
//input: dataframe df [Input data table]
//input: column sequence {semType: Macromolecule}
export async function molfilesFromHELM(df: DG.DataFrame, sequence: DG.Column): Promise<string[][] | null> {
  const monomersLibFile = await _package.files.readAsText(HELM_CORE_LIB_FILENAME);
  const monomersLibDf = DG.DataFrame.fromJson(monomersLibFile);
  const result = getMolfilesFromSeq(sequence, monomersLibDf);
  return result;
}


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
  // Higher priority for columns with MSA data to show with WebLogo.
  const tv = grok.shell.tv;
  const df = tv.dataFrame;
  const semTypeColList = df.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  let col: DG.Column | undefined = semTypeColList.find((col) => {
    const units = col.getTag(DG.TAGS.UNITS);
    return units ? units.indexOf('MSA') !== -1 : false;
  });
  if (!col)
    col = semTypeColList[0];

  if (!col) {
    grok.shell.error('Current table does not contain sequences');
    return;
  }

  tv.addViewer('WebLogo', {sequenceColumnName: col.name});
}

// helper function for importFasta
function parseMacromolecule(
  fileContent: string,
  startOfSequence: number,
  endOfSequence: number
): string {
  const seq = fileContent.slice(startOfSequence, endOfSequence);
  const seqArray = seq.split(/\s/);
  return seqArray.join('');
}

//name: importFasta
//description: Opens FASTA file
//tags: file-handler
//meta.ext: fasta, fna, ffn, faa, frn, fa
//input: string fileContent
//output: list tables
export function importFasta(fileContent: string): DG.DataFrame [] {
  const regex = /^>(.*)$/gm; // match lines starting with >
  const descriptionsArray = [];
  const sequencesArray: string[] = [];
  let startOfSequence = 0;
  let match; // match.index is the beginning of the matched line
  while (match = regex.exec(fileContent)) {
    const description = fileContent.substring(match.index + 1, regex.lastIndex);
    descriptionsArray.push(description);
    if (startOfSequence !== 0)
      sequencesArray.push(parseMacromolecule(fileContent, startOfSequence, match.index));
    startOfSequence = regex.lastIndex + 1;
  }
  sequencesArray.push(parseMacromolecule(fileContent, startOfSequence, -1));
  const descriptionsArrayCol = DG.Column.fromStrings('description', descriptionsArray);
  const sequenceCol = DG.Column.fromStrings('sequence', sequencesArray);
  sequenceCol.semType = 'Macromolecule';

  const stats: SeqColStats = WebLogo.getStats(sequenceCol, 5, WebLogo.splitterAsFasta);
  const seqType = stats.sameLength ? 'SEQ.MSA' : 'SEQ';
  const alphabetCandidates: [string, Set<string>][] = [
    ['NT', new Set(Object.keys(Nucleotides.Names))],
    ['PT', new Set(Object.keys(Aminoacids.Names))],
  ];
  // Calculate likelihoods for alphabet_candidates
  const alphabetCandidatesSim: number[] = alphabetCandidates.map(
    (c) => WebLogo.getAlphabetSimilarity(stats.freq, c[1]));
  const maxCos = Math.max(...alphabetCandidatesSim);
  const alphabet = maxCos > 0.65 ? alphabetCandidates[alphabetCandidatesSim.indexOf(maxCos)][0] : 'UN';
  sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
  const units: string = `fasta:${seqType}:${alphabet}`;
  sequenceCol.setTag(DG.TAGS.UNITS, units);

  return [DG.DataFrame.fromColumns([
    descriptionsArrayCol,
    sequenceCol,
  ])];
}

//name: Bio | Convert
//friendly-name: Bio | Convert
//tags: panel, bio
//input: column col {semType: Macromolecule}
export function convertPanel(col: DG.Column): void {
  convert(col);
}
