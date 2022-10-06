import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import { BehaviorSubject  } from "rxjs";
import {UaFilter} from "./filter2";
import {ViewHandler} from "./view-handler";
import {ChoiceInput} from "./elements/choice-input";

export class UaToolbox {
  rootAccordion: DG.Accordion;
  dateInput: DG.InputBase;
  usersInput: ChoiceInput;
  filterStream: BehaviorSubject<UaFilter>;

  static async construct() {
    let date = 'this month';
    let users = ['all'];

    let dateInput = ui.stringInput('Date', date);
    dateInput.addPatternMenu('datetime');
    dateInput.setTooltip('Set the date period');
    let usersInput = await ChoiceInput.construct();
    
    let filterStream = new BehaviorSubject(new UaFilter({
      date: date,
      users: users
    }));

    return new UaToolbox(dateInput, usersInput, filterStream)
  }

  private constructor(dateInput: DG.InputBase, usersInput: ChoiceInput, filterStream: BehaviorSubject<UaFilter>) {
    this.rootAccordion = ui.accordion();
    this.rootAccordion.addPane('Filters', () => ui.narrowForm([
      dateInput,
      // @ts-ignore
      usersInput.field,
      // @ts-ignore
      ui.buttonsInput([ui.bigButton('Apply', () => this.applyFilter())])
    ]), true);

    this.dateInput = dateInput;
    this.usersInput = usersInput;
    this.filterStream = filterStream;
  }

  async getFilter() {
    return new UaFilter({
      date: this.dateInput.value,
      users: await UaToolbox.getUsersInGroups(this.usersInput.getSelectedGroups())
    });
  }

  static async getUsersInGroups(groups: string[]) {
    if(groups.length == 1 && groups[0] == 'all')
      return ['all'];
    let res =  (await grok.data.query('UsageAnalysis:GetUsersInGroups', {groups: groups}))
        .columns.byName('login').toList();
    return res;
  }



  async applyFilter () {
    this.filterStream.next(await this.getFilter());
    ViewHandler.getInstance().setUrlParam('date', this.dateInput.value, true);
    ViewHandler.getInstance().setUrlParam('users', this.usersInput.getSelectedGroups().join(','), true);
  }

  setDate(value: string) {
    this.dateInput.value = value;
  }

  async setUsers(value: string) {
    this.usersInput.addItems(value.split(','));
  }
}