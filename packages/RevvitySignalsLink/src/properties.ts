import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { BaseConditionEditor, ConditionRegistry, Operators, SUGGESTIONS_FUNCTION } from './query-builder';
import { getUserIdByUserString, getUsersSuggestions, getUserStringIdById } from './users';

export const REVVITY_USER = 'revvity_user';

export const created = DG.Property.create('Created', DG.TYPE.DATE_TIME, (x: any) => x, (x: any, v) => x = v);
export const edited = DG.Property.create('Edited', DG.TYPE.DATE_TIME, (x: any) => x, (x: any, v) => x = v);
export const creator = DG.Property.create('Creator', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
creator.semType = REVVITY_USER;
creator.options[SUGGESTIONS_FUNCTION] = getUsersSuggestions;
export const editor = DG.Property.create('Editor', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
editor.semType = REVVITY_USER;
editor.options[SUGGESTIONS_FUNCTION] = getUsersSuggestions;
export const structure = DG.Property.create('Structure', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
//export const id = DG.Property.create('Id_float', DG.TYPE.FLOAT, (x: any) => x, (x: any, v) => x = v);
//export const idInt = DG.Property.create('Id_int', DG.TYPE.INT, (x: any) => x, (x: any, v) => x = v);
structure.semType = DG.SEMTYPE.MOLECULE;

export function getDefaultProperties(): DG.Property[] {
    return [created, edited, creator, editor, structure];
}

export const PROPERTY_NAMES_TO_QUERY_MAPPING = {
    'Created': 'createdAt',
    'Edited': 'editedAt',
    'Creator': 'createdBy',
    'Editor': 'editedBy',
}

export const NOT_IN_TAGS = [ 'Created', 'Edited', 'Creator', 'Editor', 'Structure', 'isMaterial', 'type', 'assetTypeEid'];

export const REVVITY_FIELD_TO_PROP_TYPE_MAPPING: {[key: string]: DG.TYPE} = {
    'double': DG.TYPE.FLOAT,
    'text': DG.TYPE.STRING,
    'date': DG.TYPE.DATE_TIME,
}

export class RevvityUserConditionEditor extends BaseConditionEditor<string> {
    override initializeEditor(prop: DG.Property): void {
        if (!this.condition.value)
            this.createUserInput('', prop);
        else
            getUserStringIdById(this.condition.value).then((userString: string) => {
                this.createUserInput(userString, prop);
            })
    }

    createUserInput(initValue: string, prop: DG.Property) {
        const userInput = ui.input.string('', {
            value: initValue,
            onValueChanged: async () => {
                this.condition.value = await getUserIdByUserString(userInput.value);
                this.onChanged.next(this.condition);
                this.addSuggestions(prop, userInput);
            }
        });
        this.root.append(userInput.root);
        this.initializeSuggestions(prop, userInput);
    }
}

//register operators for molecule semType, applicable for Revvity
ConditionRegistry.getInstance().registerSemTypeOperators(DG.SEMTYPE.MOLECULE, [Operators.CONTAINS, Operators.IS_SIMILAR]);
//register operators for user semtype
ConditionRegistry.getInstance().registerSemTypeOperators(REVVITY_USER, [Operators.EQ]);
ConditionRegistry.getInstance().registerEditor(DG.TYPE.STRING, REVVITY_USER, Operators.EQ, RevvityUserConditionEditor);
