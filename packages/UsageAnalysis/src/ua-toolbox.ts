import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import {InputBase} from "datagrok-api/src/widgets";
import { BehaviorSubject  } from "rxjs";
import {UaFilter} from "./filter2";

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
    this.dateInput = ui.stringInput('Date', 'this week');
    this.dateInput.addPatternMenu('datetime');
    this.dateInput.setTooltip('Set the date period');

    this.usersInput = ui.stringInput('User', 'all');
    this.usersInput.setTooltip('Enter users login');

    this.rootAccordion = ui.accordion();

    this.filterStream = new BehaviorSubject(this.getFilter());

    this.rootAccordion.addPane('Filters', () => ui.narrowForm([
      this.dateInput,
      this.usersInput,
      // @ts-ignore
      ui.buttonsInput([ui.bigButton('Apply', () => {
        this.filterStream.next(this.getFilter());
      })])
    ]), true);

  }

}