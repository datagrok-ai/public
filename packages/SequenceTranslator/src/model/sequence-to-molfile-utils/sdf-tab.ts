import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// datagrok libraries dependencies
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

// inner dependencies
import {download} from '../helpers';
import {sequenceToMolV3000} from './from-monomers';
import {linkStrandsV3000} from './mol-transformations';
import {isValidSequence} from '../code-converter/conversion-validation-tools';
import '../../view/css/sdf-tab.css';

export type StrandData = {
  strand: string,
  invert: boolean
}

/** Get a molfile for a single strand */
export function getMolfileForStrand(strand: string, invert: boolean): string {
  if (strand === '')
    return '';
  const validationOutput = isValidSequence(strand, null);
  const format = validationOutput.synthesizer![0];
  let molfile = '';
  try {
    molfile = sequenceToMolV3000(strand, invert, format!);
  } catch (err) {
    const errStr = errorToConsole(err);
    console.error(errStr);
  }
  return molfile;
}

/** Get molfile for single strand or linked strands */
export function getLinkedMolfile(
  ss: StrandData, as: StrandData, as2: StrandData, useChiral: boolean
): string {
  const nonEmptyStrands = [ss, as, as2].filter((item) => item.strand !== '');
  if (nonEmptyStrands.length === 1) {
    return getMolfileForStrand(nonEmptyStrands[0].strand, nonEmptyStrands[0].invert);
  } else {
    const ssMol = getMolfileForStrand(ss.strand, ss.invert);
    const asMol = getMolfileForStrand(as.strand, as.invert);
    const as2Mol = getMolfileForStrand(as2.strand, as2.invert);

    // select only the non-empty anti-strands
    const antiStrands = [asMol, as2Mol].filter((item) => item !== '');
    const resultingMolfile = linkStrandsV3000({senseStrands: [ssMol], antiStrands: antiStrands}, useChiral);

    return resultingMolfile;
  }
}

/** Save sdf in case ss and as (and optionally as2) strands entered */
export function saveSdf(
  ss: StrandData, as: StrandData, as2: StrandData, useChiral: boolean,
  oneEntity: boolean
): void {
  const nonEmptyStrands = [ss.strand, as.strand, as2.strand].filter((item) => item !== '');
  if (
    nonEmptyStrands.length === 0 ||
    nonEmptyStrands.length === 1 && ss.strand === ''
  ) {
    grok.shell.warning('Enter SS and AS/AS2 to save SDF');
  } else {
    let result: string;
    if (oneEntity) {
      result = getLinkedMolfile(ss, as, as2, useChiral) + '\n$$$$\n';
    } else {
      const ssMol = getMolfileForStrand(ss.strand, ss.invert);
      const asMol = getMolfileForStrand(as.strand, as.invert);
      const as2Mol = getMolfileForStrand(as2.strand, as2.invert);
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
