import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { AVAILABLE_FPS, FP_SETTINGS_KEY, FP_SETTINGS_STORAGE_KEY } from '../constants';
import { FpParams, _fpParamsMap, fpBuilderMap } from './fp-settings';
import { Fingerprint } from '../utils/chem-common';

export class FingerprintsSettingsEditor {

    fpTypeInput: DG.InputBase;
    fps = AVAILABLE_FPS;
    fpSettingsDiv = ui.inputs([]);


    constructor() {
        this.fpTypeInput = ui.choiceInput('FP', AVAILABLE_FPS[0], AVAILABLE_FPS, () => {
            this.createAlgorithmSettingsDiv(this.fpSettingsDiv, _fpParamsMap.get(this.fpTypeInput.value!) ?? undefined);
        });
        this.createAlgorithmSettingsDiv(this.fpSettingsDiv, _fpParamsMap.get(this.fpTypeInput.value!) ?? undefined);
    }

    getEditor(): HTMLElement {
        return this.fpSettingsDiv;
    }

    private createAlgorithmSettingsDiv(
        paramsForm: HTMLElement, params?: FpParams): HTMLElement {
        ui.empty(paramsForm);
        paramsForm.append(this.fpTypeInput.root);
        if (params) {
            Object.keys(params).forEach((it: any) => {
                const param = (params as any)[it];
                const input =
                    ui.intInput(param.uiName, param.value as any, () => {
                        param.value = input.value;
                    });
                ui.tooltip.bind(input.input ?? input.root, param.tooltip);
                paramsForm.append(input.root);
            });
        }
        return paramsForm;
    }
}


export async function getStoredFpSettings(): Promise<void> {
    let storedFpSettingsStr = await grok.dapi.userDataStorage.getValue(FP_SETTINGS_STORAGE_KEY, FP_SETTINGS_KEY, true);
    let storedFpSettingsObj: { [key: string]: { [key: string]: number } } = {};
    if (storedFpSettingsStr)
        storedFpSettingsObj = JSON.parse(storedFpSettingsStr);
    Object.keys(fpBuilderMap).forEach((fpType) => {
        _fpParamsMap.set(fpType as Fingerprint, fpBuilderMap[fpType](storedFpSettingsObj[fpType]));
    })
}


export function storeFpSettings(): void {
    let fpSettingsObj: { [key: string]: { [key: string]: number } } = {};
    _fpParamsMap.forEach((value: FpParams, fpType: Fingerprint) => {
        fpSettingsObj[fpType] = value.params;
    });
    grok.dapi.userDataStorage.postValue(FP_SETTINGS_STORAGE_KEY, FP_SETTINGS_KEY, JSON.stringify(fpSettingsObj), true);
}