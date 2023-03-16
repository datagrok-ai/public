import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import '../../css/choice_input.css';
import Choices from 'choices.js';


export class ChoiceInputGroups {
  choices: Choices;
  field: DG.InputBase;

  static async construct() {
    const field = ui.choiceInput('Groups', [], []);

    field.input.setAttribute('multiple', '');
    const choices = new Choices(field.input, {
      addItems: true,
      removeItems: true,
      removeItemButton: true,
      searchEnabled: true,
      searchChoices: true,
      itemSelectText: '',
    });

    const groups = grok.dapi.groups;
    field.input.addEventListener('search', async (event) => {
      const newGroups: DG.Group[] = await groups.getGroupsLookup(choices.input.value);
      choices.clearChoices();
      choices.setChoices(() => newGroups.map((g: DG.Group) => {
        return {value: g.name, label: g.name};
      }));
    });

    return new ChoiceInputGroups(choices, field);
  }

  private constructor(choices: Choices, field: DG.InputBase) {
    this.choices = choices;
    this.field = field;
  }

  getSelectedGroups(): string[] {
    const users = this.choices?.getValue(true) as string[];
    if (users && users.length > 0)
      return users;
    else
      return ['all'];
  }

  addItems(items: string[]) {
    items.forEach((i) => this.choices?._addItem({value: i}));
  }
}
