import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

// import $ from 'cash-dom';

import {download} from '../utils/helpers';

import {sequenceToMolV3000} from '../utils/structures-works/from-monomers';
import {linkStrandsV3000} from '../utils/structures-works/mol-transformations';
import {getFormat} from './sequence-codes-tools';

import {drawMolecule} from '../utils/structures-works/draw-molecule';

function getMolfileForImg(
  ss: string, as: string,
  as2: string | null = null,
  invertSS: boolean, invertAS: boolean, invertAS2: boolean,
  useChiral: boolean
): string {
  const formatAs = getFormat(as);
  const formatSs = getFormat(ss);
  let formatAs2: string | null = null;
  let molAS2: string | null = null;

  let molSS = '';
  let molAS = '';
  try {
    // a workaround to get the SS depicted in case AS is empty
    molSS = sequenceToMolV3000(ss, invertSS, false, formatSs!);
    molAS = sequenceToMolV3000(as, invertAS, false, formatAs!);
  } catch (err) {
    const errStr = errorToConsole(err);
    console.error(errStr);
  }

  if (as2 !== null && as2 !== '') {
    formatAs2 = getFormat(as2!);
    molAS2 = sequenceToMolV3000(as2, invertAS2!, false, formatAs2!);
  }

  const antiStrands = molAS2 === null ? [molAS] : [molAS, molAS2];
  const resultingMolfile = linkStrandsV3000({senseStrands: [molSS], antiStrands: antiStrands}, useChiral);
  return resultingMolfile;
}

export function saveSdf(as: string, ss: string,
  oneEntity: boolean, useChiral: boolean,
  invertSS: boolean, invertAS: boolean,
  as2: string | null = null, invertAS2: boolean | null
): void {
  if (ss === '' && as === '') {
    grok.shell.warning('Enter a sequence to save SDF');
  } else {
    const formatAs = getFormat(as);
    const formatSs = getFormat(ss);
    let formatAs2: string | null = null;
    let molAS2: string | null = null;

    const molSS = sequenceToMolV3000(ss, invertSS, false, formatSs!);
    const molAS = sequenceToMolV3000(as, invertAS, false, formatAs!);

    if (as2 !== null && as2 !== '') {
      formatAs2 = getFormat(as2!);
      molAS2 = sequenceToMolV3000(as2, invertAS2!, false, formatAs2!);
    }

    let result: string;
    if (oneEntity) {
      const antiStrands = molAS2 === null ? [molAS] : [molAS, molAS2];
      result = linkStrandsV3000({senseStrands: [molSS], antiStrands: antiStrands}, useChiral) + '\n$$$$\n';
    } else {
      result =
        molSS + '\n' +
        `> <Sequence>\nSense Strand\n$$$$\n` +
        molAS + '\n' +
        `> <Sequence>\nAnti Sense\n$$$$\n`;

      if (molAS2) {
        result += molAS2+ '\n' +
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

  // inputs
  const ssInput = ui.textInput('', '', () => { onInput.next();
    // hack
  });
  const asInput = ui.textInput('', '', () => { onInput.next(); });
  const as2Input = ui.textInput('', '', () => { onInput.next(); });
  const saveEntity = ui.boolInput('Save as one entity', true);
  const useChiralInput = ui.boolInput('Use chiral', true);

  // default values
  const straight = '5′ → 3′';
  const inverse = '3′ → 5′';
  let invertSS = false;
  let invertAS = false;
  const invertAS2 = false;

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

  // todo: port as much css as possible to external files
  // tableLayout.style.backgroundColor = 'var(--grey-1)';
  // tableLayout.style.borderBottom = 'solid';
  // tableLayout.style.borderBottomColor = 'var(--steel-5)';

  // table element styling
  for (const item of [ssLabel, asLabel, as2Label]) {
    const parentTd = item.parentElement!;
    parentTd.style.minWidth = '95px';
    parentTd.style.textAlign = 'right';
    parentTd.style.verticalAlign = 'top';
  }

  for (const item of [ssDirection.root, asDirection.root, as2Direction.root]) {
    const parentTd = item.parentElement!;
    parentTd.style.minWidth = '100px';
    parentTd.style.textAlign = 'right';
    parentTd.style.verticalAlign = 'top';
    parentTd.style.float = 'right';
  }

  for (const item of [ssInput, asInput, as2Input]) {
    item.root.parentElement!.style.width = '100%';
    const textArea = item.root.children[1];
    //@ts-ignore
    textArea.style.width = '100%';
    //@ts-ignore
    textArea.style.resize = 'none';
    //@ts-ignore
    //textArea.style.borderStyle = 'solid'; // for an unknown reason, borders disappear on input
    ////@ts-ignore
    //textArea.style.backgroundColor = 'pink';
    ////@ts-ignore
    //textArea.style.borderWidth = 'thin';
  }

  // molecule image
  const moleculeImgDiv = ui.block([]);

  DG.debounce<string>(onInput, 300).subscribe(async () => {
    let molfile = '';
    try {
      molfile = getMolfileForImg(
        ssInput.value, asInput.value, as2Input.value, invertSS, invertAS, invertAS2, useChiralInput.value!
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
    // moleculeImgDiv.style.borderStyle = 'solid';
    // moleculeImgDiv.style.borderColor = 'var(--blue-1)';
    // moleculeImgDiv.style.borderWidth = 'thin';
    // @ts-ignore
    moleculeImgDiv.children[0].style.float = 'right';
  });
  moleculeImgDiv.style.marginRight = '20px';
  moleculeImgDiv.style.float = 'right';

  const saveButton = ui.buttonsInput([
    ui.bigButton('Save SDF', () =>
      saveSdf(
        asInput.value, ssInput.value, saveEntity.value!,
        useChiralInput.value!, invertSS, invertAS, as2Input.value, invertAS2)
    )
  ]);

  // const lineBelowTable = ui.divH([saveEntity.root, useChiralInput.root, saveButton]);
  const boolInputsAndButtonArray = [saveEntity.root, useChiralInput.root, saveButton];
  const boolInputsAndButton = ui.divV(boolInputsAndButtonArray);
  for (const item of boolInputsAndButtonArray) {
    item.style.justifyContent = 'right';
    item.style.marginBottom = '10px';
  }
  // boolInputsAndButton.style.backgroundColor = 'var(--green-1)';
  const bottom = ui.divH([boolInputsAndButton, moleculeImgDiv]);
  bottom.style.flexDirection = 'row-reverse';
  bottom.style.paddingTop = '20px';
  const body = ui.divV([tableLayout, bottom]);
  body.style.paddingRight = '20px';

  return body;
}
