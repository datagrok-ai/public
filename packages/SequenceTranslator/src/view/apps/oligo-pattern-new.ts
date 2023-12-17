/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULT_PTO, DEFAULT_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH, USER_STORAGE_KEY, SS, AS, STRAND_NAME, STRANDS, TERMINAL, TERMINAL_KEYS, THREE_PRIME, FIVE_PRIME, JSON_FIELD as FIELD
} from '../../model/pattern-app/const';

import * as rxjs from 'rxjs';
import $ from 'cash-dom';

export class PatternLayoutHandler {
  get htmlDivElement(): HTMLDivElement {
    const leftSection = this.getLeftSection();
    // const rightSection = this.getRightSection();
    const layout = ui.splitH([
      leftSection,
      // rightSection
    ], {}, true);

    return layout;
  }

  private getLeftSection(): HTMLDivElement {
    return new LeftSection().getLayout();
  }

  // private getRightSection(): HTMLDivElement {
  // }
}

class LeftSection {
  private onCreateAsStrand = new rxjs.Subject<boolean>();
  private onStrandLengthChange = Object.fromEntries(
    STRANDS.map((strand) => [strand, new rxjs.Subject<number>()])
  );

  public getLayout(): HTMLDivElement {
    const patternControlsBlock = this.getPatternControlsBlock();
    const convertControlsBlock = this.getConvertControlsBlock();
    const leftSection = ui.box(
      ui.div([
        ...patternControlsBlock,
        ...convertControlsBlock,
      ], 'ui-form')
      , {style:{maxWidth:'450px'}});
    return leftSection;
  }

  private getPatternControlsBlock(): HTMLElement[] {
    const patternControlsBlock = [
      ui.h1('Pattern'),
      this.createAsStrandInput.root,
      this.strandLengthInput[SS].root,
      this.strandLengthInput[AS].root,
      // sequenceBase.root,
      // comment.root,
      // loadPatternDiv,
      // saveAs.root,
    ];

    return patternControlsBlock;
  }

  private getConvertControlsBlock(): HTMLDivElement[] {
    const convertControlsBlock = [
      ui.h1('Convert'),
      // tables.root,
      // strandColumnInput[SS],
      // strandColumnInput[AS],
      // inputIdColumn.root,
      // ui.buttonsInput([
      //   convertSequenceButton,
      // ]),
    ];

    return convertControlsBlock;
  }

  // todo: rename to input
  private get createAsStrandInput(): DG.InputBase<boolean> {
    const createAsStrand = ui.boolInput(
      `${STRAND_NAME.AS} strand`, true, (value: boolean) => {
        this.onCreateAsStrand.next(value)
      // modificationSection[AS].hidden = !v;
      // strandColumnInputDiv[AS].hidden = !v;
      // asLengthDiv.hidden = !v;
      // asModificationDiv.hidden = !v;
      // asExampleDiv.hidden = !v;
      // firstPto[AS].root.hidden = !v;
      // updateSvgScheme();
    });
    createAsStrand.setTooltip('Create antisense strand sections on SVG and table to the right');

    return createAsStrand as DG.InputBase<boolean>;
  }

  private get strandLengthInput(): {[strand: string]: DG.InputBase<number | null>} {
    const strandLengthInput = Object.fromEntries(STRANDS.map(
      (strand) => {
        const input = ui.intInput(
          `${STRAND_NAME[strand]} length`, DEFAULT_SEQUENCE_LENGTH,
          (value: number) => this.onStrandLengthChange[strand].next(value)
        );
        input.setTooltip(`Length of ${STRAND_NAME[strand].toLowerCase()}, including overhangs`);
        return [strand, input];
    }));

    this.onCreateAsStrand.subscribe((createAsCriterion: boolean) => {
      $(strandLengthInput[AS].root).toggle(createAsCriterion);
    })

    return strandLengthInput;
  }
}
