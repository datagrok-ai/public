import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {SubscriptionManager} from '../../model/subscription-manager';

import '../style.css';

import {StringInput} from '../types';

import $ from 'cash-dom';
import {DataManager} from '../../model/data-manager';
import {EventBus} from '../../model/event-bus';

export class PatternLoadControlsManager {
  private subscriptions = new SubscriptionManager();
  private authorSelectedByUser = false;

  constructor(
    private eventBus: EventBus,
    private dataManager: DataManager
  ) {
    this.eventBus.patternLoadRequested$.subscribe((patternHash: string) => this.handlePatternChoice(patternHash));

    this.eventBus.patternDeletionRequested$.subscribe(async (patternName: string) => {
      await this.dataManager.deletePattern(patternName, this.eventBus);
    });
  }

  // private selectedPattern: string;

  private async handlePatternChoice(patternHash: string): Promise<void> {
    let patternConfiguration = await this.dataManager.getPatternConfig(patternHash);
    if (patternConfiguration === null)
      patternConfiguration = this.dataManager.getDefaultPatternConfig();

    this.eventBus.setPatternConfig(patternConfiguration);

    // const patternName = patternConfiguration.patternName;
    // this.selectedPattern = patternName;

    this.eventBus.updateControlsUponPatternLoaded(patternHash);
  }

  private isCurrentUserSelected(): boolean {
    return this.eventBus.getSelectedAuthor() !== this.dataManager.getOtherUsersAuthorshipCategory();
  }

  createControls(): HTMLElement[] {
    const inputsContainer = this.getPatternInputsContainer();
    return [
      ui.h1('Load'),
      inputsContainer,
    ];
  }

  private getPatternInputsContainer(): HTMLDivElement {
    const inputsContainer = ui.divH(this.createPatternInputs());

    this.eventBus.patternListUpdated$.subscribe(() => {
      this.subscriptions.unsubscribeAll();

      // todo: change values of selectedPattern and selectedUser here

      $(inputsContainer).empty();
      $(inputsContainer).append(this.createPatternInputs());
    });

    return inputsContainer;
  }

  private createPatternInputs(): HTMLElement[] {
    const authorChoiceInput = this.createAuthorChoiceInput();
    const patternChoiceInputContainer = this.createPatternChoiceInputContainer();
    const deletePatternButton = this.createDeletePatternButton();

    return [
      authorChoiceInput.root,
      patternChoiceInputContainer,
      deletePatternButton
    ];
  }

  private createPatternChoiceInputContainer(): HTMLDivElement {
    const patternChoiceInput = this.createPatternChoiceInput();
    const patternChoiceInputContainer = ui.div([patternChoiceInput.root]);

    const subscription = this.eventBus.userSelection$.subscribe(() => {
      $(patternChoiceInputContainer).empty();
      $(patternChoiceInputContainer).append(this.createPatternChoiceInput().root);
    });
    this.subscriptions.add(subscription);

    return patternChoiceInputContainer;
  }

  private createAuthorChoiceInput(): StringInput {
    const possibleValues = [this.dataManager.getCurrentUserAuthorshipCategory()];
    if (this.dataManager.getOtherUsersPatternNames().length > 0)
      possibleValues.push(this.dataManager.getOtherUsersAuthorshipCategory());

    const authorChoiceInput = ui.choiceInput(
      'Author',
      this.eventBus.getSelectedAuthor(),
      possibleValues,
      (userName: string) => {
        this.authorSelectedByUser = true;
        this.eventBus.selectAuthor(userName);
      }
    );
    this.setAuthorChoiceInputStyle(authorChoiceInput);
    authorChoiceInput.setTooltip('Select pattern author');

    return authorChoiceInput;
  }

  private setAuthorChoiceInputStyle(userChoiceInput: StringInput): void {
    $(userChoiceInput.input).css({
      'max-width': '100px',
      'min-width': '100px',
    });
    $(userChoiceInput.root).css({
      'padding-right': '30px',
      'padding-left': '30px',
    });
  }

  private createPatternChoiceInput(): StringInput {
    const patternList = this.isCurrentUserSelected() ?
      this.dataManager.getCurrentUserPatternNames() :
      this.dataManager.getOtherUsersPatternNames();

    if (this.authorSelectedByUser) {
      const patternHash = this.dataManager.getPatternHash(patternList[0], this.isCurrentUserSelected());
      this.eventBus.requestPatternLoad(patternHash);
      this.eventBus.updateUrlState(patternHash);
      this.authorSelectedByUser = false;
    }

    const defaultValue = this.getPatternName(patternList);
    const choiceInput = ui.choiceInput('Pattern', defaultValue, patternList);
    choiceInput.setTooltip('Select pattern to load');

    $(choiceInput.input).css({
      'max-width': '100px',
      'min-width': '100px',
    });

    this.subscriptions.add(
      choiceInput.onInput(
        () => {
          const patternHash = this.dataManager.getPatternHash(choiceInput.value!, this.isCurrentUserSelected());
          this.eventBus.requestPatternLoad(patternHash);
          this.eventBus.updateUrlState(patternHash);
        }
      )
    );

    this.subscriptions.add(
      this.eventBus.patternLoaded$.subscribe(() => {
        const patternName = this.eventBus.getPatternName();
        if (!choiceInput.value?.includes(patternName))
          choiceInput.value = this.getPatternName(patternList);
      })
    );

    return choiceInput;
  }

  private getPatternName(patternList: string[]): string {
    return patternList.find((patternName) => patternName.includes(this.eventBus.getPatternName())) ?? patternList[0];
  }

  private createDeletePatternButton(): HTMLButtonElement {
    const button = ui.button(
      ui.iconFA('trash-alt'),
      () => {
        if (this.eventBus.getPatternName() === this.dataManager.getDefaultPatternName()) {
          grok.shell.warning('Cannot delete example pattern');
          return;
        }
        this.showDeletePatternDialog();
      }
    );

    ui.tooltip.bind(button, 'Delete pattern from user storage');

    const subscription = this.eventBus.userSelection$.subscribe(() => {
      $(button).toggle(this.isCurrentUserSelected());
    });
    this.subscriptions.add(subscription);

    return button;
  }

  private showDeletePatternDialog(): void {
    const dialog = ui.dialog('Delete pattern');
    const patternName = this.eventBus.getPatternName();
    dialog.add(ui.divText(`Are you sure you want to delete pattern ${patternName}?`));
    dialog.onOK(() => this.eventBus.requestPatternDeletion(patternName));
    dialog.show();
  }
}
