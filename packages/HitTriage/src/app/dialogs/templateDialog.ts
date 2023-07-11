import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';
import {UiUtils} from '@datagrok-libraries/compute-utils';
import {ITemplate} from '../types';
import * as C from '../consts';

export async function createTemplateDialog() {
  return new Promise<ITemplate>(async (resolve) => {
    const functions = DG.Func.find({tags: [C.HitTriageComputeFunctionTag]});
    const availableTemplates = (await _package.files.list('templates')).map((file) => file.name.slice(0, -5));
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
    let file: File | null = null;
    const fileInput = UiUtils.fileInput('', null, (f: File) => {file = f;}, undefined);
    const functionsInput = ui.multiChoiceInput('Select functions', [],
      Object.keys(functionsMap), () => {});

    async function onOkProxy() {
      if (!file) {
        grok.shell.error('File is not selected');
        return;
      }

      if (errorDiv.style.opacity === '100%') {
        grok.shell.error('Template name is empty or already exists');
        return;
      }
      const fileText = await file.text();
      const df = DG.DataFrame.fromCsv(fileText);
      await df.meta.detectSemanticTypes();
      const molCol = df.columns.bySemType(DG.SEMTYPE.MOLECULE);
      if (!molCol) {
        grok.shell.error('Molecule column is not found');
        return;
      }
      const molColName = molCol.name;
      const d = new Date();
      const time =
       `${d.getFullYear()}-${d.getMonth()}-${d.getDate()}_${d.getHours()}-${d.getMinutes()}-${d.getSeconds()}`;
      const folder = `${time}_${grok.shell.user.login}`;
      const fullName = `samples/${folder}/${file.name}`;
      await _package.files.writeAsText(`samples/${folder}/${file.name}`, fileText);
      const out: ITemplate = {
        name: TemplateNameInput.value,
        ingest: {
          type: 'File',
          query: 'System:AppData/HitTriage/' + fullName,
          molColName: molColName,
        },
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
      };
      saveTemplate(out);
      resolve(out);
    }
    ui.dialog('Create template')
      .add(ui.divV([TemplateNameInput, errorDiv]))
      .add(fileInput.root)
      .add(functionsInput)
      .onOK(onOkProxy)
      .show();
  });
}

export function saveTemplate(template: ITemplate) {
  _package.files.writeAsText(`templates/${template.name}.json`, JSON.stringify(template));
}
