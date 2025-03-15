import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import '../../css/choice_input.css';
import Choices from 'choices.js';
import {ChoiceInputBase} from "./choice-input-base";


export class ChoiceInputProjects extends ChoiceInputBase {
    static async construct() {
        const field = ui.input.choice('Projects');

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
        field.input.addEventListener('change', (_) =>
            (document.querySelector('.ua-apply-button') as HTMLButtonElement).disabled = false);
        const projects: string[] = ((await grok.functions.call('UsageAnalysis:ProjectsList')) as DG.DataFrame).getCol('project_name').toList();

        choices.setChoices(() => projects.map((p: string) => {
            return {value: p, label: p};
        }));
        return new ChoiceInputProjects(choices, field);
    }

    private constructor(choices: Choices, field: DG.InputBase) {
        super(choices, field);
    }
}
