import * as ui from 'datagrok-api/ui';

import {OTHER_USERS} from '../../model/const';
import {SubscriptionManager} from '../../model/subscription-manager';

import '../style.css';

import {StringInput} from '../types';

import $ from 'cash-dom';
import {EventBus} from '../../model/event-bus';
import {PatternAppDataManager} from '../../model/external-data-manager';

export class PatternLoadControlsManager {
  private subscriptions = new SubscriptionManager();

  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager
  ) {
    this.eventBus.patternLoadRequested$.subscribe((value: string) => this.handlePatternChoice(value));

    const defaultUser = this.dataManager.getCurrentUserName();
    this.eventBus.selectUser(defaultUser);

    const defaultPattern = this.dataManager.getCurrentUserPatternNames()[0];
    this.selectedPattern = defaultPattern;
  }

  private selectedPattern: string;

  private async handlePatternChoice(patternName: string): Promise<void> {
    const patternConfiguration = await this.dataManager.getPatternConfig(patternName, this.isCurrentUserSelected());
    this.eventBus.setPatternConfig(patternConfiguration);
    this.selectedPattern = patternName;

    this.eventBus.updateControlsUponPatternLoaded();
  }

  private isCurrentUserSelected(): boolean {
    return this.eventBus.getSelectedUser() !== OTHER_USERS;
  }

  createUIComponents(): HTMLElement[] {
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
    const currentUser = this.dataManager.getCurrentUserName();
    const possibleValues = [currentUser + ' (me)', OTHER_USERS];

    const userChoiceInput = ui.choiceInput(
      'Author',
      this.eventBus.getSelectedUser(),
      possibleValues,
      (userName: string) => this.eventBus.selectUser(userName)
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
      [' '].concat(this.dataManager.getCurrentUserPatternNames()) :
      this.dataManager.getOtherUsersPatternNames();
    this.selectedPattern = patternList[0] || '<default>';
    this.eventBus.requestPatternLoad(this.selectedPattern);

    const choiceInput = ui.choiceInput('Pattern', this.selectedPattern, patternList);
    choiceInput.setTooltip('Select pattern to load');

    $(choiceInput.input).css({
      'max-width': '100px',
      'min-width': '100px',
    });

    const subscription = choiceInput.onInput(
      () => this.eventBus.requestPatternLoad(choiceInput.value!)
    );
    this.subscriptions.add(subscription);

    return choiceInput;
  }

  private createDeletePatternButton(): HTMLButtonElement {
    const button = ui.button(
      ui.iconFA('trash-alt'),
      () => this.showDeletePatternDialog()
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
