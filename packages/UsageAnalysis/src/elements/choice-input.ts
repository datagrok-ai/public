import "../../css/choice_input.css";

import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";

import Choices from "choices.js";
import {InputBase} from "datagrok-api/dg";

export class ChoiceInput {

  choices: Choices | null = null;
  root: InputBase | null = null;
  ready: Promise<void>;

  constructor (funcToGetList: Function) {

    let init = async () => {
      // @ts-ignore
      let field = ui.choiceInput('Lookup field',[], await funcToGetList());

      field.input.setAttribute('multiple', '');
      this.choices = new Choices(field.input, {
        addItems: true,
        removeItems: true,
        removeItemButton: true,
        searchEnabled: true,
        searchChoices: true,
        itemSelectText: ''
      });

      this.root = field;
    };

    this.ready = init();
  }

  getSelectedUsers() : string[] {
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
