import "../../css/choice_input.css";

import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";

import Choices from "choices.js";
import {InputBase} from "datagrok-api/dg";

export class ChoiceInput {
  choices: Choices;
  field: InputBase;

  static async construct (funcToGetList: Function) {
   // @ts-ignore
    let field = ui.choiceInput('Groups',[], await funcToGetList());

    field.input.setAttribute('multiple', '');
    let choices = new Choices(field.input, {
      addItems: true,
      removeItems: true,
      removeItemButton: true,
      searchEnabled: true,
      searchChoices: true,
      itemSelectText: ''
    });

    return new ChoiceInput(choices, field);
  }

  private constructor(choices: Choices, field: InputBase) {
    this.choices = choices;
    this.field = field;
  }

  getSelectedGroups() : string[] {
    let users = [];
    // @ts-ignore
    users = this.choices?.getValue(true);
    if (users && users.length > 0)
      return users;
    else
      return ['all'];
  }

  addItems(items: string[]) {
    items.forEach((i) => this.choices?._addItem({value: i}));
  }
}
