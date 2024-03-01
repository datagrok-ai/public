import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { AVAILABLE_FPS, FP_SETTINGS_KEY, FP_SETTINGS_STORAGE_KEY } from '../constants';
import { FpParams, _fpParamsMap, fpBuilderMap } from './fp-settings';
import { Fingerprint } from '../utils/chem-common';
import { _package } from '../package';

const SKETCHER_FUNCS_FRIENDLY_NAMES: {[key: string]: string} = {
    OpenChemLib: 'OpenChemLib',
    Ketcher: 'Ketcher',
    Marvin: 'Marvin',
    ChemDraw: 'ChemDraw',
  };
  
  const PREVIOUS_SKETCHER_NAMES: {[key: string]: string} = {
    'Open Chem Sketcher': 'OpenChemLib',
    'ketcherSketcher': 'Ketcher',
    'Marvin JS': 'Marvin',
    'Chem Draw': 'ChemDraw',
  };

  const RENDERERS: string[] = [
    "RDKit",
    "OpenChemLib"
  ]

export class ChemSettingsEditor extends DG.Widget {

    fpTypeInput: DG.InputBase;
    rendererInput: DG.InputBase;
    fps = AVAILABLE_FPS;
    fpSettingsDiv = ui.inputs([]);
    sketcherInput: DG.InputBase;
    propList: DG.Property[];


    constructor(propList: DG.Property[]) {
        super(ui.div());
        this.propList = propList;

        //sketcher settings
        const sketcherNames = Object.values(SKETCHER_FUNCS_FRIENDLY_NAMES);
        const sketcherProp = this.propList.filter((it) => it.name === 'Sketcher');
        this.sketcherInput = ui.choiceInput('Sketcher', sketcherProp.length ? sketcherProp[0].get(null) : sketcherNames[0], sketcherNames, async () => {
            if (sketcherProp.length)
                sketcherProp[0].set(null, this.sketcherInput.value);
        });

        //renderer settings
        const rendererProp = this.propList.filter((it) => it.name === 'Renderer');
        this.rendererInput = ui.choiceInput('Renderer', rendererProp.length ? rendererProp[0].get(null) : RENDERERS[0], RENDERERS, async () => {
            if (rendererProp.length)
                rendererProp[0].set(null, this.sketcherInput.value);
        });

        //fingerprints settings
        const fps = AVAILABLE_FPS.filter((fp) => fp !== Fingerprint.MACCS);
        this.fpTypeInput = ui.choiceInput('FP', fps[0], fps, () => {
            this.createFpSettingsDiv(this.fpSettingsDiv, this.fpTypeInput.value!,
                _fpParamsMap.get(this.fpTypeInput.value!) ?? undefined);
        });
        const acc = ui.accordion();
        const pane = acc.addPane('Fingerprints settings',
            () => this.createFpSettingsDiv(this.fpSettingsDiv, this.fpTypeInput.value!,
                _fpParamsMap.get(this.fpTypeInput.value!) ?? undefined));
        pane.root.classList.add('chem-fingerprints-settings-acc-pane');
            this.root.append(//@ts-ignore
            ui.inputs([
                this.sketcherInput,
                this.rendererInput,
            ])
        );
        this.root.append(acc.root);
    }

    private createFpSettingsDiv(
        paramsForm: HTMLElement, fpName: string, params?: FpParams): HTMLElement {
        ui.empty(paramsForm);
        paramsForm.append(this.fpTypeInput.root);
        if (params) {
            Object.keys(params).forEach((it: any) => {
                const param = (params as any)[it];
                const input =
                    ui.intInput(param.uiName, param.value as any, () => {
                        const paramName = it;
                        const fpProp = this.propList.filter((it) => it.name === `${fpName}_${paramName}`);
                        if(fpProp.length)
                            fpProp[0].set(null, input.value);
                        param.value = input.value;
                    });
                ui.tooltip.bind(input.input, param.tooltip);
                paramsForm.append(input.root);
            });
        }
        return paramsForm;
    }
}

/* 
Updates current fingerprints settings. User storage takes precedence over package settings
 */
export async function loadFpSettings(): Promise<void> {
    let storedFpSettingsStr = await grok.dapi.userDataStorage.getValue(FP_SETTINGS_STORAGE_KEY, FP_SETTINGS_KEY, true);
    let storedFpSettingsObj: { [key: string]: { [key: string]: number } } = {};
    if (storedFpSettingsStr)
        storedFpSettingsObj = JSON.parse(storedFpSettingsStr);
    else {
        let regexStr = '^';
        AVAILABLE_FPS.forEach((fp, idx) => regexStr += idx === AVAILABLE_FPS.length - 1 ? `${fp}_` : `${fp}_|`);
        const regexp = new RegExp(regexStr);
        const properties = await _package.getProperties();
        Object.keys(properties).forEach((key) => {
            if(regexp.test(key)) {
                const underscoreIdx = key.indexOf('_');
                const fp = key.substring(0, underscoreIdx);
                if (!storedFpSettingsObj[fp]) {
                    storedFpSettingsObj[fp] = {};
                    storedFpSettingsObj[fp][key.substring(underscoreIdx+1)] = properties[key];
                }
            }
        })
    }
    Object.keys(fpBuilderMap).forEach((fpType) => {
        _fpParamsMap.set(fpType as Fingerprint, fpBuilderMap[fpType](storedFpSettingsObj[fpType]));
    })
}

export function storeFpSettingsToUserLocalStorage(): void {
    let fpSettingsObj: { [key: string]: { [key: string]: number } } = {};
    _fpParamsMap.forEach((value: FpParams, fpType: Fingerprint) => {
        fpSettingsObj[fpType] = value.params;
    });
    grok.dapi.userDataStorage.postValue(FP_SETTINGS_STORAGE_KEY, FP_SETTINGS_KEY, JSON.stringify(fpSettingsObj), true);
}


/* 
Updates current sketcher value. User storage takes precedence over package settings
 */
export async function loadSketcherSettings(): Promise<void> {
    const properties = await _package.getProperties();
    let storedSketcherType = await grok.dapi.userDataStorage.getValue(DG.chem.STORAGE_NAME, DG.chem.KEY, true);
    if (PREVIOUS_SKETCHER_NAMES[storedSketcherType])
      storedSketcherType = PREVIOUS_SKETCHER_NAMES[storedSketcherType];
    if (!storedSketcherType && properties.Sketcher)
      storedSketcherType = SKETCHER_FUNCS_FRIENDLY_NAMES[properties.Sketcher];
    const sketcherFunc = DG.Func.find({tags: ['moleculeSketcher']})
      .find((e) => e.name === storedSketcherType || e.friendlyName === storedSketcherType);
    if (sketcherFunc)
      DG.chem.currentSketcherType = sketcherFunc.friendlyName;
    else {
      if (!!storedSketcherType) {
        grok.shell.warning(
          `Package with ${storedSketcherType} function is not installed.Switching to ${DG.DEFAULT_SKETCHER}.`);
      }
  
      DG.chem.currentSketcherType = DG.DEFAULT_SKETCHER;
    }
}