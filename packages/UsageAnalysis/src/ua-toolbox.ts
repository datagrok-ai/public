import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import {InputBase} from "datagrok-api/src/widgets";
import { BehaviorSubject  } from "rxjs";
import {UaFilter} from "./filter2";
import {ViewHandler} from "./view-handler";

export class UaToolbox {
  rootAccordion: DG.Accordion;
  dateInput: InputBase;
  usersInput: InputBase;
  filterStream: BehaviorSubject<UaFilter>;

  getFilter() {
    return new UaFilter({
      date: this.dateInput.value,
      users: [this.usersInput.value]
    });
  }

  constructor() {
    this.dateInput = ui.stringInput('Date', 'this month');
    this.dateInput.addPatternMenu('datetime');
    this.dateInput.setTooltip('Set the date period');

    this.usersInput = ui.stringInput('Users', 'all');
    this.usersInput.setTooltip('Comma-separated user logins');

    this.rootAccordion = ui.accordion();

    this.filterStream = new BehaviorSubject(this.getFilter());

    this.rootAccordion.addPane('Filters', () => ui.narrowForm([
      this.dateInput,
      this.usersInput,
      // @ts-ignore
      ui.buttonsInput([ui.bigButton('Apply', () => this.applyFilter())])
    ]), true);

  }

  applyFilter () {
    this.filterStream.next(this.getFilter());
    ViewHandler.getInstance().setUrlParam('date', this.dateInput.value, true);
    ViewHandler.getInstance().setUrlParam('users', this.usersInput.value, true);
  }

  setDate(value: string) {
    this.dateInput.value = value;
  }

  setUsers(value: string) {
    this.usersInput.value = value;
  }
}