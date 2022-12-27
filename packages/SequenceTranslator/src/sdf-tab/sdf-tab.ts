import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

import * as rxjs from 'rxjs';

import {download} from '../utils/helpers';
import {sequenceToMolV3000} from '../utils/structures-works/from-monomers';
import {linkStrandsV3000} from '../utils/structures-works/mol-transformations';
import {getFormat} from './sequence-codes-tools';

import {drawMolecule} from '../utils/structures-works/draw-molecule';

/** Data associated with strands */
type StrandData = {
  strand: string,
  invert: boolean
}

/** Get a molfile for a single strand */
function getMolfileForStrand(strand: string, invert: boolean): string {
  const format = getFormat(strand);
  let molfile = '';
  try {
    molfile = sequenceToMolV3000(strand, invert, false, format!);
  } catch (err) {
    const errStr = errorToConsole(err);
    console.error(errStr);
  }
  return molfile;
}

/** Get molfile for single strand or linked strands */
function getLinkedMolfile(
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
  if (ss.strand === '' && (as.strand !== '' && as2.strand !== '')) {
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

/** UI of the SDF tab on the application's view */
export function getSdfTab(): HTMLDivElement {
  const onInput: rxjs.Subject<string> = new rxjs.Subject<string>();

  // label names
  const ssInputLabel = 'Sense Strand';
  const asInputLabel = 'Anti Sense';
  const as2InputLabel = 'Anti Sense 2';
  const ssDirectionLabel = 'SS direction';
  const asDirectionLabel = 'AS direction';
  const as2DirectionLabel = 'AS2 direction';
  // const inputColumnName = 'Sequences';
  // const directionColumnName = 'Direction';
  const inputColumnName = '';
  const directionColumnName = '';

  // default input values
  const straight = '5′ → 3′';
  const inverse = '3′ → 5′';
  let invertSS = false;
  let invertAS = false;
  const invertAS2 = false;

  // text inputs
  const ssInput = ui.textInput('', '', () => { onInput.next(); });
  const asInput = ui.textInput('', '', () => { onInput.next(); });
  const as2Input = ui.textInput('', '', () => { onInput.next(); });

  // bool inputs
  const saveEntity = ui.boolInput('Save as one entity', true);
  const useChiralInput = ui.boolInput('Use chiral', true);

  // choice inputs
  const ssDirection = ui.choiceInput(ssDirectionLabel, straight, [straight, inverse]);
  ssDirection.onChanged(() => { invertSS = ssDirection.value === inverse; });
  const asDirection = ui.choiceInput(asDirectionLabel, straight, [straight, inverse]);
  asDirection.onChanged(() => { invertAS = asDirection.value === inverse; });
  const as2Direction = ui.choiceInput(as2DirectionLabel, straight, [straight, inverse]);
  as2Direction.onChanged(() => { invertAS = asDirection.value === inverse; });

  // labels
  const ssLabel = ui.label(ssInputLabel);
  const asLabel = ui.label(asInputLabel);
  const as2Label = ui.label(as2InputLabel);

  // table layout
  const tableLayout = ui.table(
    ['ss', 'as1', 'as2'], (row, index) => {
      switch (row) {
      case 'ss':
        return [ssLabel, ssInput.root, ssDirection.root];
      case 'as1':
        return [asLabel, asInput.root, asDirection.root];
      case 'as2':
        return [as2Label, as2Input.root, as2Direction.root];
      }
    }, ['', inputColumnName, directionColumnName]
  );

  // text input label style
  for (const item of [ssLabel, asLabel, as2Label]) {
    const parentTd = item.parentElement!;
    parentTd.style.minWidth = '95px';
    parentTd.style.textAlign = 'right';
    parentTd.style.verticalAlign = 'top';
    parentTd.style.paddingTop = '3px';
  }

  // choice input label style
  for (const item of [ssDirection.root, asDirection.root, as2Direction.root]) {
    const parentTd = item.parentElement!;
    parentTd.style.minWidth = '100px';
    parentTd.style.textAlign = 'right';
    parentTd.style.verticalAlign = 'top';
    parentTd.style.float = 'right';
  }

  // text area style
  for (const item of [ssInput, asInput, as2Input]) {
    item.root.parentElement!.style.width = '100%';
    const textArea = item.root.children[1];
    //@ts-ignore
    textArea.style.width = '100%';
    //@ts-ignore
    textArea.style.resize = 'none';
  }

  // molecule image container
  const moleculeImgDiv = ui.block([]);

  DG.debounce<string>(onInput, 300).subscribe(async () => {
    let molfile = '';
    try {
      molfile = getLinkedMolfile(
        {strand: ssInput.value, invert: invertSS},
        {strand: asInput.value, invert: invertAS},
        {strand: as2Input.value, invert: invertAS2}, useChiralInput.value!
      );
    } catch (err) {
      const errStr = errorToConsole(err);
      console.error(errStr);
    }
    // todo: calculate relative numbers
    const canvasWidth = 650;
    const canvasHeight = 150;
    // todo: remove div with image if molfile empty
    await drawMolecule(moleculeImgDiv, canvasWidth, canvasHeight, molfile);
    // @ts-ignore
    moleculeImgDiv.children[0].style.float = 'right';
  });
  moleculeImgDiv.style.marginRight = '30px';
  moleculeImgDiv.style.float = 'right';

  const saveButton = ui.buttonsInput([
    ui.bigButton('Save SDF', () =>
      saveSdf(
        {strand: ssInput.value, invert: invertSS},
        {strand: asInput.value, invert: invertAS},
        {strand: as2Input.value, invert: invertAS2},
        useChiralInput.value!,
        saveEntity.value!
      )
    )
  ]);

  const boolInputsAndButtonArray = [saveEntity.root, useChiralInput.root, saveButton];
  const boolInputsAndButton = ui.divV(boolInputsAndButtonArray);
  for (const item of boolInputsAndButtonArray) {
    item.style.justifyContent = 'right';
    item.style.marginBottom = '10px';
  }

  const bottomDiv = ui.divH([boolInputsAndButton, moleculeImgDiv]);
  bottomDiv.style.flexDirection = 'row-reverse';
  bottomDiv.style.paddingTop = '20px';

  const body = ui.divV([tableLayout, bottomDiv]);
  body.style.paddingRight = '20px';

  return body;
}
