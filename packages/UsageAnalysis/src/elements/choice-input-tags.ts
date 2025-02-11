import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import '../../css/choice_input.css';
import Choices from 'choices.js';


export class ChoiceInputTags {
    choices: Choices;
    field: DG.InputBase;

    static async construct() {
        const field = ui.input.choice('Tags');

        field.input.setAttribute('multiple', '');
        field.input.querySelectorAll('option').forEach((o) => {o.selected = false;});
        const choices = new Choices(field.input, {
            addItems: true,
            removeItems: true,
            removeItemButton: true,
            searchEnabled: true,
            searchChoices: true,
            itemSelectText: '',
        });
        field.input.addEventListener('change', (event) =>
            (document.querySelector('.ua-apply-button') as HTMLButtonElement).disabled = false);
        const tags: string[] = ((await grok.functions.call('UsageAnalysis:EntitiesTags')) as DG.DataFrame).getCol('tag').toList();

        choices.setChoices(() => tags.map((p: string) => {
            return {value: p, label: p};
        }));
        return new ChoiceInputTags(choices, field);
    }

    private constructor(choices: Choices, field: DG.InputBase) {
        this.choices = choices;
        this.field = field;
    }

    getSelectedTags(): string[] {
        const tags = this.choices?.getValue(true) as string[];
        if (tags && tags.length > 0 && tags[0] !== '[]')
            return tags;
        else
            return ['any'];
    }

    addItems(items: string[]) {
        items.forEach((i) => this.choices?._addItem({value: i}));
    }
}
