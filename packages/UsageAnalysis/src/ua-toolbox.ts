import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
// import * as grok from 'datagrok-api/grok';

import {BehaviorSubject} from 'rxjs';
import {UaFilter} from './filter';
import {ViewHandler} from './view-handler';
import {ChoiceInputGroups} from './elements/choice-input-groups';
import {ChoiceInputPackages} from './elements/choice-input-packages';
// import {UaView} from './tabs/ua';
import $ from 'cash-dom';

export class UaToolbox {
  rootAccordion: DG.Accordion;
  dateInput: DG.InputBase;
  groupsInput: ChoiceInputGroups;
  packagesInput: ChoiceInputPackages;
  filterStream: BehaviorSubject<UaFilter>;
  dateFromDD: DG.InputBase = ui.stringInput('From', '');
  dateToDD: DG.InputBase = ui.stringInput('To', '');
  usersDD: DG.InputBase = ui.stringInput('Users', '');
  packagesDD: DG.InputBase = ui.stringInput('Packages', '');
  formDD: HTMLDivElement;

  static async construct() {
    const date = 'this week';
    const packages = ['all'];
    const dateInput = ui.stringInput('Date', date);
    dateInput.addPatternMenu('datetime');
    dateInput.setTooltip('Set the date period');
    const groupsInput = await ChoiceInputGroups.construct();
    const packagesInput = await ChoiceInputPackages.construct();
    const groups = groupsInput.allUsers;
    const filterStream = new BehaviorSubject(new UaFilter({
      date: date,
      groups: groups,
      packages: packages,
    }));
    return new UaToolbox(dateInput, groupsInput, packagesInput, filterStream);
  }

  private constructor(dateInput: DG.InputBase, groupsInput: ChoiceInputGroups,
    packagesInput: ChoiceInputPackages, filterStream: BehaviorSubject<UaFilter>) {
    this.rootAccordion = ui.accordion();
    this.formDD = ui.div();
    this.dateInput = dateInput;
    this.groupsInput = groupsInput;
    this.packagesInput = packagesInput;
    this.filterStream = filterStream;

    const filters = this.rootAccordion.addPane('Filters', () => {
      const form = ui.narrowForm([
        dateInput,
        groupsInput.field,
        packagesInput.field,
      ]);
      const applyB = ui.bigButton('Apply', () => this.applyFilter());
      applyB.style.marginLeft = 'auto';
      $(form).append(applyB);
      this.dateFromDD.readOnly = true;
      this.dateToDD.readOnly = true;
      this.usersDD.readOnly = true;
      this.packagesDD.readOnly = true;
      this.formDD = ui.narrowForm([
        this.dateFromDD,
        this.dateToDD,
        this.usersDD,
        this.packagesDD,
      ]);
      this.formDD.style.display = 'none';
      const closeButton = ui.button('', () => {
        this.formDD.style.display = 'none';
        this.clearFormDD();
        ViewHandler.getCurrentView().getScatterPlot().reloadViewer();
      });
      closeButton.innerHTML = '<div class="tab-handle-close-button"><i class="fal fa-times"></i></div>';
      closeButton.style.margin = '0';
      closeButton.style.padding = '0';
      closeButton.style.marginLeft = 'auto';
      this.formDD.prepend(closeButton);
      this.formDD.style.border = '1px dashed #80949b';
      this.formDD.style.borderRadius = '5px';
      this.formDD.style.backgroundColor = '#F1FAFD';
      this.formDD.style.marginBottom = '24px';
      this.formDD.querySelectorAll('.ui-input-editor')
        .forEach((i) => (i as HTMLElement).style.backgroundColor = '#F1FAFD');
      return form;
    }, true);
    filters.root.before(this.formDD);

    ViewHandler.UA.tabs.onTabChanged.subscribe((tab) => {
      if (this.formDD.style.display === 'block') {
        this.formDD.style.display = 'none';
        this.applyFilter();
      }
      if (this.checkLabels()) this.formDD.style.display = 'block';
    });
  }

  checkLabels() {
    return [this.dateFromDD.value, this.dateToDD.value,
      this.usersDD.value, this.packagesDD.value].some((val) => val);
  }

  clearFormDD() {
    this.dateFromDD.value = '';
    this.dateToDD.value = '';
    this.usersDD.value = '';
    this.packagesDD.value = '';
  }

  getFilter() {
    return new UaFilter({
      date: this.dateInput.value,
      groups: this.groupsInput.getSelectedGroups(),
      packages: this.packagesInput.getSelectedPackages(),
    });
  }

  applyFilter() {
    this.clearFormDD();
    this.filterStream.next(this.getFilter());
    ViewHandler.getInstance().setUrlParam('date', this.dateInput.value, true);
    ViewHandler.getInstance().setUrlParam('users', this.groupsInput.getSelectedGroups().join(','), true);
    ViewHandler.getInstance().setUrlParam('packages', this.packagesInput.getSelectedPackages().join(','), true);
  }

  setDate(value: string) {
    this.dateInput.value = value;
  }

  setGroups(value: string) {
    this.groupsInput.addItems(value.split(','));
  }

  setPackages(value: string) {
    this.packagesInput.addItems(value.split(','));
  }
}
