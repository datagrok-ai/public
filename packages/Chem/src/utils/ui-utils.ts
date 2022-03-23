import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function createChoiceInputProperty(propName: string, obj: any, choiceInputsParams: any, choiceInputPostfix: string) {
    const choiceInputParams = choiceInputsParams[propName];
    obj[propName] = obj.string(propName, choiceInputParams.options[0], { choices: choiceInputParams.options });
    const InputPropName = `${propName}${choiceInputPostfix}`;
    obj[InputPropName] = ui.choiceInput(choiceInputParams.choiceInputName, obj[propName], choiceInputParams.options);
    obj[InputPropName].onChanged(() => {
        obj.setOptions({ [propName]: `${obj[InputPropName].value}` });
    });
}


export function createChoiceInputsDiv(obj: any, choiceInputsParams: any, choiceInputPostfix: string) {
    const div = ui.inputs([]);
    Object.keys(choiceInputsParams).forEach(key => {
        div.append(obj[`${key}${choiceInputPostfix}`].root);
    })
    return div;
}


export async function updateChoiceInputValue(obj: any, propertyName: string, choiceInputPostfix: string) {
    const choiceInput = obj[`${propertyName}${choiceInputPostfix}`];
    if (choiceInput) {
      let val = obj[propertyName];
      if (choiceInput.value !== val) {
        choiceInput.value = val;
      }
    }
  }