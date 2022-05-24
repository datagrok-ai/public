import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import {InputBase} from "datagrok-api/src/widgets";
import { BehaviorSubject  } from "rxjs";
import {UaFilter} from "./filter2";
import {ViewHandler} from "./view-handler";
import {ChoiceInput} from "./elements/choice-input";

export class UaToolbox {
  rootAccordion: DG.Accordion;
  dateInput: InputBase;
  usersInput: ChoiceInput;
  filterStream: BehaviorSubject<UaFilter>;

  getFilter() {
    return new UaFilter({
      date: this.dateInput.value,
      users: this.usersInput.getSelectedUsers()
    });
  }

  constructor() {
    this.dateInput = ui.stringInput('Date', 'this month');
    this.dateInput.addPatternMenu('datetime');
    this.dateInput.setTooltip('Set the date period');

    this.rootAccordion = ui.accordion();

    this.usersInput = new ChoiceInput(async () => (await grok.dapi.users.list()).map((u) => u.login));
    this.usersInput.ready.then(() => {
      this.rootAccordion.addPane('Filters', () => ui.narrowForm([
        this.dateInput,
        // @ts-ignore
        this.usersInput.root,
        // @ts-ignore
        ui.buttonsInput([ui.bigButton('Apply', () => this.applyFilter())])
      ]), true);
    });
    this.filterStream = new BehaviorSubject(this.getFilter());

  }

  applyFilter () {
    this.filterStream.next(this.getFilter());
    ViewHandler.getInstance().setUrlParam('date', this.dateInput.value, true);
    ViewHandler.getInstance().setUrlParam('users', this.usersInput.getSelectedUsers().join(','), true);
  }

  setDate(value: string) {
    this.dateInput.value = value;
  }

  setUsers(value: string) {
    this.usersInput.ready.then(() => {
      this.usersInput.addItems(value.split(','));
    });
  }
}