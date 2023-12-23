/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULT_PHOSPHOROTHIOATE, DEFAULT_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH, USER_STORAGE_KEY, SENSE_STRAND, AS, STRAND_NAME, STRANDS, TERMINAL, TERMINAL_KEYS, THREE_PRIME, FIVE_PRIME, JSON_FIELD as FIELD
} from '../../../model/pattern-app/const';
import {generateExample, translateSequence, getShortName, isPatternCreatedByCurrentUser, findDuplicates, addColumnWithIds, addColumnWithTranslatedSequences} from '../../../model/pattern-app/oligo-pattern';
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
      // strandColumnInput[SENSE_STRAND],
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
    this.updatePatternsList(loadPatternInputDiv);

    const patternControlsBlock = [
      ui.h1('Pattern'),
      this.createAsStrandInput.root,
      this.strandLengthInput[SENSE_STRAND].root,
      this.strandLengthInput[AS].root,
      this.sequenceBaseInput.root,
      comment.root,
      loadPatternInputDiv,
      // saveAs,
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

  private async updatePatternsList(loadPatternDiv: HTMLDivElement) {
    const updatePatternsListFunction = async () => this.updatePatternsList(loadPatternDiv);

    function createPatternChoiceInput(patternList: string[], callback: (arg: string) => any) {
      return ui.choiceInput('Load pattern', '', patternList, (v: string) => callback(v));
    }

    function updateLoadPattern(): void {
      loadPattern.root.append(patternListChoiceInput.input);
      loadPattern.root.append(loadPattern.input);
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
            //todo: restore
            // else if (await isPatternCreatedByCurrentUser(saveAs.value))
            //   grok.shell.warning('Cannot delete pattern, created by other user');
            // else {
            //   await grok.dapi.userDataStorage.remove(USER_STORAGE_KEY, loadPattern.value, false)
            //     .then(() => grok.shell.info('Pattern \'' + loadPattern.value + '\' deleted'));
            // }
            await updatePatternsListFunction();
          }),
        ], 'ui-input-options'),
      );
    }

    const entities = await grok.dapi.userDataStorage.get(USER_STORAGE_KEY, false);

    const currentUser = await grok.dapi.users.current();
    const currentUserName = currentUser.friendlyName;

    const otherUsers = 'Other users';

    // todo: rename
    const lstMy: string[] = [];
    const lstOthers: string[] = [];

    // todo: rename
    for (const ent of Object.keys(entities)) {
      const isCreatedByCurrentUser = await isPatternCreatedByCurrentUser(ent);
      (isCreatedByCurrentUser ? lstOthers : lstMy).push(ent);
    }

    let loadPattern = createPatternChoiceInput(lstMy,
      (v: string) => {
        // todo: restore
        // parsePatternAndUpdateUi(v)
      }
    );

    const patternListChoiceInput = ui.choiceInput('', currentUserName, [currentUserName, otherUsers], (v: string) => {
      // todo: rename v to value
      const currentList = v === currentUserName ? lstMy : lstOthers;
      loadPattern = createPatternChoiceInput(
        currentList,
        (v: string) => {
          // todo: restore
          // parsePatternAndUpdateUi(v);
        }
      )
      updateLoadPattern();
    });
    patternListChoiceInput.input.style.maxWidth = '142px';

    updateLoadPattern();
  }
}

class ConvertBlock {

}
