import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {BehaviorSubject} from 'rxjs';
import {UaFilter} from './filter';
import {ViewHandler} from './view-handler';
import {ChoiceInputGroups} from './elements/choice-input-groups';
import {ChoiceInputPackages} from './elements/choice-input-packages';
import $ from 'cash-dom';

export class UaToolbox {
  rootAccordion: DG.Accordion;
  dateInput: DG.InputBase;
  usersInput: ChoiceInputGroups;
  packagesInput: ChoiceInputPackages;
  filterStream: BehaviorSubject<UaFilter>;

  static async construct() {
    const date = 'this week';
    const users = ['all'];
    const packages = ['all'];

    const dateInput = ui.stringInput('Date', date);
    dateInput.addPatternMenu('datetime');
    dateInput.setTooltip('Set the date period');
    const usersInput = await ChoiceInputGroups.construct();
    const packagesInput = await ChoiceInputPackages.construct();

    const filterStream = new BehaviorSubject(new UaFilter({
      date: date,
      users: users,
      packages: packages,
    }));

    return new UaToolbox(dateInput, usersInput, packagesInput, filterStream);
  }

  private constructor(dateInput: DG.InputBase, usersInput: ChoiceInputGroups,
    packagesInput: ChoiceInputPackages, filterStream: BehaviorSubject<UaFilter>) {
    this.rootAccordion = ui.accordion();
    this.rootAccordion.addPane('Filters', () => {
      const form = ui.narrowForm([
        dateInput,
        usersInput.field,
        packagesInput.field,
      ]);
      $(form).append(ui.bigButton('Apply', () => this.applyFilter()));
      return form;
    }, true);

    this.dateInput = dateInput;
    this.usersInput = usersInput;
    this.packagesInput = packagesInput;
    this.filterStream = filterStream;
  }

  async getFilter() {
    return new UaFilter({
      date: this.dateInput.value,
      users: await UaToolbox.getUsersInGroups(this.usersInput.getSelectedGroups()),
      packages: this.packagesInput.getSelectedPackages(),
    });
  }

  static async getUsersInGroups(groups: string[]) {
    if (groups.length == 1 && groups[0] == 'all')
      return ['all'];
    const res = (await grok.data.query('UsageAnalysis:GetUsersInGroups', {groups: groups}))
      .columns.byName('login').toList();
    return res;
  }

  async applyFilter() {
    this.filterStream.next(await this.getFilter());
    ViewHandler.getInstance().setUrlParam('date', this.dateInput.value, true);
    ViewHandler.getInstance().setUrlParam('users', this.usersInput.getSelectedGroups().join(','), true);
    ViewHandler.getInstance().setUrlParam('packages', this.packagesInput.getSelectedPackages().join(','), true);
  }

  setDate(value: string) {
    this.dateInput.value = value;
  }

  async setUsers(value: string) {
    this.usersInput.addItems(value.split(','));
  }

  async setPackages(value: string) {
    this.packagesInput.addItems(value.split(','));
  }
}
