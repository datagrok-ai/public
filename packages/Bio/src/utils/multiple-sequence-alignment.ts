import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {FastaFileHandler} from '@datagrok-libraries/bio/src/utils/fasta-handler';
import {ALIGNMENT, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
//@ts-ignore: there are no types for this library
import Aioli from '@biowasm/aioli';

import {AlignedSequenceEncoder} from '@datagrok-libraries/bio/src/sequence-encoder';
import {kalignVersion} from './constants';

const fastaInputFilename = 'input.fa';
const fastaOutputFilename = 'result.fasta';

export class MsaWarning extends Error {
  constructor(
    public readonly element: HTMLElement, options?: ErrorOptions) {
    super(element.innerText, options);
  }
}

/**
 * Converts array of sequences into simple fasta string.
 *
 * @param {string[]} sequences Input list of sequences.
 * @return {string} Fasta-formatted string.
 */
function _stringsToFasta(sequences: string[]): string {
  return sequences.reduce((a, v, i) => a + `>sample${i + 1}\n${v}\n`, '');
}

/**
 * Runs Aioli environment with kalign tool.
 *
 * @param {DG.Column} srcCol Column with sequences.
 * @param {boolean} isAligned Whether the column is aligned.
 * @param {string | undefined} unUsedName
 * @param {DG.Column | null} clustersCol Column with clusters.
 * @param {number | undefined} gapOpen Gap open penalty.
 * @param {number | undefined} gapExtend Gap extend penalty.
 * @param {number | undefined} terminalGap Terminal gap penalty.
 * @return {Promise<DG.Column>} Aligned sequences.
 */
export async function runKalign(srcCol: DG.Column<string>, isAligned: boolean = false, unUsedName: string = '',
  clustersCol: DG.Column | null = null, gapOpen?: number, gapExtend?: number, terminalGap?: number,
): Promise<DG.Column> {
  let sequences: string[] = srcCol.toList();

  if (isAligned)
    sequences = sequences.map((v: string) => AlignedSequenceEncoder.clean(v).replace(/\-/g, ''));

  const sequencesLength = srcCol.length;
  clustersCol ??= DG.Column.string('Clusters', sequencesLength).init('0');
  if (clustersCol.type != DG.COLUMN_TYPE.STRING)
    clustersCol = clustersCol.convertTo(DG.TYPE.STRING);
  clustersCol.compact();

  //TODO: use fixed-size inner arrays, but first need to expose the method to get each category count
  const clustersColCategories = clustersCol.categories;
  const clustersColData = clustersCol.getRawData();
  const fastaSequences: string[][] = new Array(clustersColCategories.length);
  const clusterIndexes: number[][] = new Array(clustersColCategories.length);
  for (let rowIdx = 0; rowIdx < sequencesLength; ++rowIdx) {
    const clusterCategoryIdx = clustersColData[rowIdx];
    (fastaSequences[clusterCategoryIdx] ??= []).push(sequences[rowIdx]);
    (clusterIndexes[clusterCategoryIdx] ??= []).push(rowIdx);
  }
  checkForSingleSeqClusters(clusterIndexes, clustersColCategories);

  const CLI = await new Aioli([
    'base/1.0.0',
    {tool: 'kalign', version: kalignVersion, reinit: true},
  ]);
  const tgtCol = DG.Column.string(unUsedName, sequencesLength);

  for (let clusterIdx = 0; clusterIdx < clustersColCategories.length; ++clusterIdx) {
    const clusterSequences = fastaSequences[clusterIdx];
    const fasta = _stringsToFasta(clusterSequences);

    await CLI.fs.writeFile(fastaInputFilename, fasta);
    const gapOpenCommand = `${gapOpen !== undefined ? ` --gpo ${gapOpen}` : ''}`;
    const gapExtendCommand = `${gapExtend !== undefined ? ` --gpe ${gapExtend}` : ''}`;
    const terminalGapCommand = `${terminalGap !== undefined ? ` --tgpe ${terminalGap}` : ''}`;
    const extraParams = `${gapOpenCommand}${gapExtendCommand}${terminalGapCommand}`;

    const output = await CLI.exec(`kalign ${fastaInputFilename} -f fasta -o ${fastaOutputFilename}${extraParams}`);
    console.warn(output);

    const buf = await CLI.cat(fastaOutputFilename);
    if (!buf) {
      const errStr = parseKalignError(output, 1);
      throw new Error(errStr);
    }

    const ffh = new FastaFileHandler(buf);
    const aligned = ffh.sequencesArray; // array of sequences extracted from FASTA
    const clusterRowIds = clusterIndexes[clusterIdx];
    for (let clusterRowIdIdx = 0; clusterRowIdIdx < aligned.length; ++clusterRowIdIdx)
      tgtCol.set(clusterRowIds[clusterRowIdIdx], aligned[clusterRowIdIdx]);
  }

  // units
  const srcUnits = srcCol.meta.units;
  //aligned
  const tgtAligned = ALIGNMENT.SEQ_MSA;
  //alphabet
  const srcAlphabet = srcCol.getTag(bioTAGS.alphabet);

  tgtCol.meta.units = srcUnits;
  tgtCol.setTag(bioTAGS.aligned, tgtAligned);
  tgtCol.setTag(bioTAGS.alphabet, srcAlphabet);
  tgtCol.semType = DG.SEMTYPE.MACROMOLECULE;
  return tgtCol;
}

export async function testMSAEnoughMemory(col: DG.Column<string>): Promise<void> {
  const sequencesCount = col.length;
  const delta = sequencesCount / 100;

  for (let i = delta; i < sequencesCount; i += delta) {
    try {
      await runKalign(DG.Column.fromStrings(col.name, col.toList().slice(0, Math.round(i))));
      console.log(`runKalign succeeded on ${i}`);
    } catch (error) {
      console.log(`runKalign failed on ${i} with '${error}'`);
    }
  }
}

function parseKalignError(out: string, limit?: number): string {
  const errLineList: string[] = [];
  const errLineRe = /^.+ERROR : (.+)$/gm;
  let ma: RegExpExecArray | null;
  while ((ma = errLineRe.exec(out)) != null && (limit === undefined || errLineList.length < limit)) {
    //
    errLineList.push(ma[1]);
  }
  return errLineList.join('\n');
}

/** */
export function checkForSingleSeqClusters(clusterIndexes: number[][], clustersColCategories: string[]): void {
  const singleSeqClusterIdxList = clusterIndexes
    .map<[number[], number]>((idxs: number[], clusterI: number) => { return [idxs, clusterI]; })
    .filter(([idxs, _clusterIdx]) => idxs.length == 1)
    .map(([_idxs, clusterIdx]) => clusterIdx);
  if (singleSeqClusterIdxList.length > 0) {
    const errEl = ui.div([
      ui.divText(`MSA analysis is not available on single sequence clusters ` +
        `#${singleSeqClusterIdxList.length}:`),
      ...wu(singleSeqClusterIdxList).take(3)
        .map((clusterIdx) => {
          let clusterName = clustersColCategories[clusterIdx];
          if (clusterName.length > 25) clusterName = clusterName.slice(0, 25) + '...';
          return ui.divText(`"${clusterName}"${clusterIdx < singleSeqClusterIdxList.length - 1 ? ', ' : '.'}`);
        }).toArray(),
      ...singleSeqClusterIdxList.length > 3 ? [ui.divText('...')] : []
    ]);
    throw new MsaWarning(errEl);
  }
}
