import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {download} from '../../common/model/helpers';
import {SequenceToMolfileConverter} from './sequence-to-molfile';
import {MonomerNotFoundError} from './monomer-code-parser';
import {linkStrandsV3000} from './mol-transformations';
import {ITranslationHelper} from '../../../types';

export type StrandData = {
  strand: string,
  invert: boolean
}

/** Get a molfile for a single strand */
export function getMolfileForStrand(strand: string, invert: boolean, th: ITranslationHelper): string {
  if (strand === '')
    return '';
  const format = th.createFormatDetector(strand).getFormat();
  if (!format)
    throw new MonomerNotFoundError(`Unable to detect the format of the sequence '${strand}'.`);
  return (new SequenceToMolfileConverter(strand, invert, format)).convert();
}

/** Get molfile for single strand or linked strands */
export function getLinkedMolfile(
  ss: StrandData, as: StrandData, as2: StrandData, useChiral: boolean, th: ITranslationHelper
): string {
  const nonEmptyStrands = [ss, as, as2].filter((item) => item.strand !== '');
  if (nonEmptyStrands.length === 1) {
    return getMolfileForStrand(nonEmptyStrands[0].strand, nonEmptyStrands[0].invert, th);
  } else {
    const ssMol = getMolfileForStrand(ss.strand, ss.invert, th);
    const asMol = getMolfileForStrand(as.strand, as.invert, th);
    const as2Mol = getMolfileForStrand(as2.strand, as2.invert, th);

    // select only the non-empty strands
    const senseStrands = [ssMol].filter((item) => item !== '');
    const antiStrands = [asMol, as2Mol].filter((item) => item !== '');
    const resultingMolfile = linkStrandsV3000({senseStrands: senseStrands, antiStrands: antiStrands}, useChiral);

    return resultingMolfile;
  }
}

/** Save sdf in case ss and as (and optionally as2) strands entered */
export function saveSdf(
  ss: StrandData, as: StrandData, as2: StrandData, useChiral: boolean, oneEntity: boolean,
  th: ITranslationHelper
): void {
  if (ss.strand === '') {
    grok.shell.warning('Enter SENSE_STRAND and optionally ANTISENSE_STRAND/AS2 to save SDF');
  } else {
    let result: string;
    try {
      if (oneEntity) {
        result = getLinkedMolfile(ss, as, as2, useChiral, th) + '\n$$$$\n';
      } else {
        const ssMol = getMolfileForStrand(ss.strand, ss.invert, th);
        const asMol = getMolfileForStrand(as.strand, as.invert, th);
        const as2Mol = getMolfileForStrand(as2.strand, as2.invert, th);
        result = ssMol + '\n' +
          `> <Sequence>\nSense Strand\n$$$$\n`;
        if (asMol) {
          result += asMol + '\n' +
            `> <Sequence>\nAnti Sense\n$$$$\n`;
        }
        if (as2Mol) {
          result += as2Mol + '\n' +
            `> <Sequence>\nAnti Sense 2\n$$$$\n`;
        }
      }
    } catch (e: any) {
      grok.shell.warning('Unable to save SDF: ' + e.message);
      return;
    }

    // construct date-time in the form yyyy-mm-dd_hh-mm-ss
    const date = new Date();

    function pad(x: number): string {
      return (x >= 10) ? x.toString() : '0' + x.toString();
    }

    const dateString: string = date.getFullYear() + '-' + pad(date.getMonth() + 1) +
      '-' + pad(date.getDate()) + '_' + pad(date.getHours()) + '-' +
      pad(date.getMinutes()) + '-' + pad(date.getSeconds());

    download(`SequenceTranslator-${dateString}.sdf`, encodeURIComponent(result));
  }
}
