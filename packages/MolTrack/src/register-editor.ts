import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Scope } from './constants'
import dayjs from 'dayjs';

export class RegisterEditor {
    main_div: HTMLDivElement
    compound_registration_div: HTMLDivElement;
    batch_registration_div: HTMLDivElement;
    assay_run_registration_div: HTMLDivElement;
    assay_results_registration_div: HTMLDivElement;
    scope_input: DG.ChoiceInput<string | null>;
    structureInput: DG.InputBase;

    epa_batch_id_input: DG.InputBase;

    assay_name_input: DG.InputBase;
    assay_run_date_input: DG.DateInput;

    constructor() {
        this.main_div = ui.div('', 'moltrack-register-single-div');
        this.compound_registration_div = ui.div('', 'moltrack-register-single-div');
        this.batch_registration_div = ui.div('', 'moltrack-register-single-div');
        this.assay_run_registration_div = ui.div('', 'moltrack-register-single-div');
        this.assay_results_registration_div = ui.div('', 'moltrack-register-single-div');

        const scopeChoices = Object.values(Scope);
        this.scope_input = ui.input.choice('Scope', {
            value: scopeChoices[0], items: scopeChoices, onValueChanged: () => {
                ui.empty(this.main_div);
                this.main_div.append(this.getScopeDiv(this.scope_input.value!));
            }
        });

        // Compound registration
        this.structureInput = ui.input.molecule('Structure');
        this.structureInput.classList.add('moltrack-register-single-input');

        // Batch registration
        this.epa_batch_id_input = ui.input.string('EPA Batch ID');

        // Assay run registration
        this.assay_name_input = ui.input.string('Assay run name');
        this.assay_run_date_input = ui.input.date('Assay run date', { value: dayjs() });

        this.main_div.append(this.getScopeDiv(this.scope_input.value!));

    }

    private getScopeDiv(scope: string): HTMLDivElement {
        let ret_val = this.main_div;
        switch (scope) {
            case Scope.COMPOUNDS:
                ui.empty(this.compound_registration_div);
                this.compound_registration_div.append(this.scope_input.root);
                this.compound_registration_div.append(this.structureInput.root);
                ret_val = this.compound_registration_div;
                break;
            case Scope.BATCHES:
                ui.empty(this.batch_registration_div);
                this.batch_registration_div.append(this.scope_input.root);
                this.batch_registration_div.append(this.epa_batch_id_input.root);
                this.batch_registration_div.append(this.structureInput.root);
                ret_val = this.batch_registration_div;
                break;
            case Scope.ASSAY_RUNS:
                ui.empty(this.assay_run_registration_div);
                this.assay_run_registration_div.append(this.scope_input.root);
                this.assay_run_registration_div.append(this.assay_name_input.root);
                this.assay_run_registration_div.append(this.assay_run_date_input.root);
                ret_val = this.assay_run_registration_div;
                break;
            case Scope.ASSAY_RESULTS:
                ui.empty(this.assay_run_registration_div);
                this.assay_results_registration_div.append(this.scope_input.root);
                this.assay_results_registration_div.append(this.assay_name_input.root);
                this.assay_results_registration_div.append(this.assay_run_date_input.root);
                this.assay_results_registration_div.append(this.epa_batch_id_input.root);
                ret_val = this.assay_results_registration_div;
                break;
        }

        return ret_val;
    }

    public getEditor(): HTMLElement {
        return this.main_div;
    }

}