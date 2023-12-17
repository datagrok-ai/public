/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULT_PTO, DEFAULT_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH, USER_STORAGE_KEY, SS, AS, STRAND_NAME, STRANDS, TERMINAL, TERMINAL_KEYS, THREE_PRIME, FIVE_PRIME, JSON_FIELD as FIELD
} from '../../../model/pattern-app/const';
import {DataManager} from './utils';
import {PatternManager} from './pattern-manager';

import * as rxjs from 'rxjs';
import $ from 'cash-dom';

export class LeftSection {
  constructor(data: DataManager) {
    this.patternBlock = new PatternBlock(data);
  }

  private patternBlock: PatternBlock;

  public getLayout(): HTMLDivElement {
    const patternControlsBlock = this.patternBlock.getHtmlElements();
    const convertControlsBlock = this.getConvertControlsBlock();
    const leftSection = ui.box(
      ui.div([
        ...patternControlsBlock,
        ...convertControlsBlock,
      ], 'ui-form')
      , {style:{maxWidth:'450px'}});
    return leftSection;
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

}

class PatternBlock {
  constructor(private data: DataManager) {};

  private onCreateAsStrand = new rxjs.Subject<boolean>();
  private onStrandLengthChange = Object.fromEntries(
    STRANDS.map((strand) => [strand, new rxjs.Subject<number>()])
  );
  private onSequenceBaseChoice = new rxjs.Subject<string>();

  // todo: this should emit messages to external handlers that have all the
  // necessary information
  private triggerUpdateSvg = new rxjs.Subject<string>();

  getHtmlElements(): HTMLElement[] {
    const comment = ui.textInput('Comment', '', () => this.triggerUpdateSvg.next());

    // const saveAs = ui.textInput('Save as', 'Pattern name', () => this.triggerUpdateSvg.next());
    // saveAs.setTooltip('Name Of New Pattern');
    // saveAs.addOptions(savePatternButton);

    const loadPatternInputDiv = ui.div([]);

    const content = [
      this.createAsStrandInput,
      this.strandLengthInput[SS],
      this.strandLengthInput[AS],
      this.sequenceBaseInput,
      comment,
      loadPatternInputDiv,
      // saveAs,
    ].map((it) => it.root);
    const patternControlsBlock = [
      ui.h1('Pattern'),
      ...content
    ];

    return patternControlsBlock;
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

  private get sequenceBaseInput(): DG.InputBase<string> {
    const baseChoices = this.data.baseChoices;
    const defaultBase = baseChoices[0];
    const input = ui.choiceInput('Sequence basis', defaultBase, baseChoices,
      (value: string) => {
        this.onSequenceBaseChoice.next(value);
        // updateBases(v);
        // updateOutputExamples();
      });
    return input as DG.InputBase<string>;
  }

  private async getLoadPatternInput(): Promise<DG.InputBase<string>> {
    const patternManager = new PatternManager();
    const currentUserPatternList = patternManager.getCurrentUserPatternList();
    let loadPattern = ui.choiceInput(
      'Load pattern', '',
      currentUserPatternList,
      // todo: restore
      // (v: string) => parsePatternAndUpdateUi(v)
    );

    const currentUserName = (await grok.dapi.users.current()).friendlyName;
    const otherUsers = 'Other users';


    const patternListChoiceInput = ui.choiceInput(
      '', currentUserName, [currentUserName, otherUsers], (v: string) => {
      const currentList = v === currentUserName ? lstMy : lstOthers;
      loadPattern = ui.choiceInput('Load pattern', '', currentList, (v: string) => parsePatternAndUpdateUi(v));

      loadPattern.root.append(patternListChoiceInput.input);
      loadPattern.root.append(loadPattern.input);
      // @ts-ignore
      loadPattern.input.style.maxWidth = '120px';
      loadPattern.input.style.marginLeft = '12px';
      loadPattern.setTooltip('Apply Existing Pattern');

      loadPatternDiv.innerHTML = '';
      loadPatternDiv.append(loadPattern.root);
      loadPattern.root.append(
        ui.div([
          ui.button(ui.iconFA('trash-alt', () => {}), async () => {
            if (loadPattern.value === null)
              grok.shell.warning('Choose pattern to delete');
            else if (await isPatternCreatedByCurrentUser(saveAs.value))
              grok.shell.warning('Cannot delete pattern, created by other user');
            else {
              await grok.dapi.userDataStorage.remove(USER_STORAGE_KEY, loadPattern.value, false)
                .then(() => grok.shell.info('Pattern \'' + loadPattern.value + '\' deleted'));
            }
            await updatePatternsList();
          }),
        ], 'ui-input-options'),
      );
    });

    return loadPattern as DG.InputBase<string>;
  }
}

class ConvertBlock {

}
