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

  private selectedPattern: string;

  private async handlePatternChoice(patternHash: string): Promise<void> {
    let patternConfiguration = await this.dataManager.getPatternConfig(patternHash);
    if (patternConfiguration === null)
      patternConfiguration = this.dataManager.getDefaultPatternConfig();

    this.eventBus.setPatternConfig(patternConfiguration);

    const patternName = patternConfiguration.patternName;
    this.selectedPattern = patternName;

    this.eventBus.updateControlsUponPatternLoaded(patternHash);
  }

  private isCurrentUserSelected(): boolean {
    return this.eventBus.getSelectedUser() !== this.dataManager.getOtherUsersAuthorshipCategory();
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
    const userChoiceInput = this.createUserChoiceInput();
    const patternChoiceInputContainer = this.createPatternChoiceInputContainer();
    const deletePatternButton = this.createDeletePatternButton();

    return [
      userChoiceInput.root,
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

  private createUserChoiceInput(): StringInput {
    const possibleValues = [this.dataManager.getCurrentUserAuthorshipCategory()];
    if (this.dataManager.getOtherUsersPatternNames().length > 0)
      possibleValues.push(this.dataManager.getOtherUsersAuthorshipCategory());

    const userChoiceInput = ui.choiceInput(
      'Author',
      this.eventBus.getSelectedUser(),
      possibleValues,
      (userName: string) => {
        this.authorSelectedByUser = true;
        this.eventBus.selectUser(userName);
      }
    );
    this.setUserChoiceInputStyle(userChoiceInput);
    userChoiceInput.setTooltip('Select pattern author');

    return userChoiceInput;
  }

  private setUserChoiceInputStyle(userChoiceInput: StringInput): void {
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
      // [' '].concat(this.dataManager.getCurrentUserPatternNames()) :
      this.dataManager.getCurrentUserPatternNames() :
      this.dataManager.getOtherUsersPatternNames();
    this.selectedPattern = patternList[0];

    if (this.authorSelectedByUser) {
      const patternHash = this.dataManager.getPatternHash(this.selectedPattern, this.isCurrentUserSelected());
      this.eventBus.requestPatternLoad(patternHash);
      this.authorSelectedByUser = false;
    }

    const choiceInput = ui.choiceInput('Pattern', this.selectedPattern, patternList);
    choiceInput.setTooltip('Select pattern to load');

    $(choiceInput.input).css({
      'max-width': '100px',
      'min-width': '100px',
    });

    const subscription = choiceInput.onInput(
      () => {
        const patternHash = this.dataManager.getPatternHash(choiceInput.value!, this.isCurrentUserSelected());
        this.eventBus.requestPatternLoad(patternHash);
      }
    );
    this.subscriptions.add(subscription);

    return choiceInput;
  }

  private createDeletePatternButton(): HTMLButtonElement {
    const button = ui.button(
      ui.iconFA('trash-alt'),
      () => {
        if (this.selectedPattern === this.dataManager.getDefaultPatternName()) {
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
    dialog.add(ui.divText(`Are you sure you want to delete pattern ${this.selectedPattern}?`));
    dialog.onOK(() => this.eventBus.requestPatternDeletion(this.selectedPattern));
    dialog.show();
  }
}
