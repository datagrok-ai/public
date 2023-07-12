import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';
import {ITemplate} from '../types';
import * as C from '../consts';

export async function createTemplateDialog() {
  return new Promise<ITemplate>(async (resolve) => {
    const functions = DG.Func.find({tags: [C.HitTriageComputeFunctionTag]});
    const availableTemplates = (await _package.files.list('templates')).map((file) => file.name.slice(0, -5));
    const availableSubmitFunctions = DG.Func.find({tags: [C.HitTriageSubmitTag]});
    const submitFunctionsMap: {[key: string]: DG.Func} = {};
    availableSubmitFunctions.forEach((func) => {
      submitFunctionsMap[func.friendlyName ?? func.name] = func;
    });
    const submitFunctionInput = ui.choiceInput('Submit function', null, Object.keys(submitFunctionsMap));
    submitFunctionInput.value = null;
    submitFunctionInput.setTooltip('Select function to be called upon submitting');
    const errorDiv = ui.divText('Template name is empty or already exists');
    errorDiv.style.color = 'red';
    errorDiv.style.textAlign = 'end';
    const TemplateNameInput = ui.stringInput('Template name', '', () => {
      if (TemplateNameInput.value === '' || availableTemplates.includes(TemplateNameInput.value)) {
        TemplateNameInput.root.style.borderBottom = '1px solid red';
        errorDiv.style.opacity = '100%';
      } else {
        TemplateNameInput.root.style.borderBottom = 'none';
        errorDiv.style.opacity = '0%';
      }
    });
    TemplateNameInput.root.style.borderBottom = 'none';
    errorDiv.style.opacity = '0%';
    const functionsMap: {[key: string]: string} = {};
    functionsMap['Descriptors'] = 'Descriptors';
    functions.forEach((func) => {
      functionsMap[func.friendlyName ?? func.name] = `${func.package.name}:${func.name}`;
    });
    const functionsInput = ui.multiChoiceInput('Select functions', [],
      Object.keys(functionsMap), () => {});
    functionsInput.setTooltip('Select functions to be applied to the data');
    functionsInput.root.style.maxHeight = '150px';
    functionsInput.root.style.overflowY = 'scroll';

    async function onOkProxy() {
      if (errorDiv.style.opacity === '100%') {
        grok.shell.error('Template name is empty or already exists');
        return;
      }
      const submitFunction = submitFunctionInput.value ? submitFunctionsMap[submitFunctionInput.value] : undefined;
      const out: ITemplate = {
        name: TemplateNameInput.value,
        compute: {
          descriptors: {
            enabled: functionsInput.value!.includes('Descriptors'),
          },
          functions: functionsInput.value!.filter((f) => f !== 'Descriptors')
            .map((f) => {
              const funcInfo = functionsMap[f].split(':');
              return {
                package: funcInfo[0],
                name: funcInfo[1],
              };
            }),
        },
        ...(submitFunction ? {submit: {fName: submitFunction.name, package: submitFunction.package.name}} : {}),
      };
      saveTemplate(out);
      resolve(out);
    }
    ui.dialog('Create template')
      .add(ui.divV([TemplateNameInput, errorDiv]))
      .add(functionsInput)
      .add(submitFunctionInput)
      .onOK(onOkProxy)
      .show();
  });
}

export function saveTemplate(template: ITemplate) {
  _package.files.writeAsText(`templates/${template.name}.json`, JSON.stringify(template));
}
