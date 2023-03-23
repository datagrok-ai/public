import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
// import * as grok from 'datagrok-api/grok';
import {BehaviorSubject} from 'rxjs';
import {UaFilter} from './filter';
import {ViewHandler} from './view-handler';
import {ChoiceInputGroups} from './elements/choice-input-groups';
import {ChoiceInputPackages} from './elements/choice-input-packages';
import $ from 'cash-dom';

export class UaToolbox {
  rootAccordion: DG.Accordion;
  dateInput: DG.InputBase;
  groupsInput: ChoiceInputGroups;
  packagesInput: ChoiceInputPackages;
  filterStream: BehaviorSubject<UaFilter>;

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
    this.rootAccordion.addPane('Filters', () => {
      const form = ui.narrowForm([
        dateInput,
        groupsInput.field,
        packagesInput.field,
      ]);
      $(form).append(ui.bigButton('Apply', () => this.applyFilter()));
      return form;
    }, true);

    this.dateInput = dateInput;
    this.groupsInput = groupsInput;
    this.packagesInput = packagesInput;
    this.filterStream = filterStream;
  }

  getFilter() {
    return new UaFilter({
      date: this.dateInput.value,
      groups: this.groupsInput.getSelectedGroups(),
      packages: this.packagesInput.getSelectedPackages(),
    });
  }

  applyFilter() {
    this.filterStream.next(this.getFilter());
    ViewHandler.getInstance().setUrlParam('date', this.dateInput.value, true);
    ViewHandler.getInstance().setUrlParam('users', this.groupsInput.getSelectedGroups().join(','), true);
    ViewHandler.getInstance().setUrlParam('packages', this.packagesInput.getSelectedPackages().join(','), true);
  }

  setDate(value: string) {
    this.dateInput.value = value;
  }

  async setGroups(value: string) {
    this.groupsInput.addItems(value.split(','));
  }

  async setPackages(value: string) {
    this.packagesInput.addItems(value.split(','));
  }
}
