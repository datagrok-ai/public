import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Subject } from 'rxjs';
import dayjs from 'dayjs';
import { IProperty } from 'datagrok-api/dg';

/**
 * Generates a color for nesting levels that progressively gets darker
 * @param level The nesting level (1-based)
 * @returns A hex color string
 */
function generateNestingColor(level: number): string {
    // Base color: light blue (#bbdefb)
    const baseColor = { r: 187, g: 222, b: 251 };
    
    // Calculate darkening factor (0.1 to 0.9)
    const darkeningFactor = Math.min(0.1 + (level - 1) * 0.08, 0.9);
    
    // Apply darkening
    const r = Math.round(baseColor.r * (1 - darkeningFactor));
    const g = Math.round(baseColor.g * (1 - darkeningFactor));
    const b = Math.round(baseColor.b * (1 - darkeningFactor));
    
    // Convert to hex
    return `#${r.toString(16).padStart(2, '0')}${g.toString(16).padStart(2, '0')}${b.toString(16).padStart(2, '0')}`;
}


// export interface SimpleConditionTest<T extends string = string> {
//     field?: string;
//     operator?: T;
//     value?: OperatorMapping[T];
// }

// type OperatorMapping = {
//    [Operators.CONTAINS]: string, 
//    [Operators.GT]: number,
//    [key: string]: any
// }

// const t: SimpleConditionTest = {
//     field: 'test',
//     operator: Operators.CONTAINS,
//     value: 5
// }

export type SimpleCondition<T = any> = {
    field: string;
    operator: string;
    value: T;
}

export type ComplexCondition = {
    logicalOperator: Operators.Logical;
    conditions: Expression[];
}

type Expression = SimpleCondition | ComplexCondition;

export class BaseConditionEditor<T = any> {
    root: HTMLDivElement = ui.div();
    condition: SimpleCondition<T>;
    onChanged: Subject<SimpleCondition<T>> = new Subject<SimpleCondition<T>>();

    constructor(prop: DG.Property, operator: string, initialCondition?: SimpleCondition<T>) {
        this.condition = initialCondition ?? {
            field: prop.name,
            operator: operator,
            value: undefined as T
        };
        this.initializeEditor(prop);
    }

    protected initializeEditor(prop: DG.Property): void {
        const input = ui.input.forProperty(prop, this.condition.value, {
            onValueChanged: () => {
                this.condition.value = input.value as T;
                this.onChanged.next(this.condition);
            }
        });
        input.addCaption('');
        this.root.append(input.root);
    }
}

export class DateBetweenConditionEditor extends BaseConditionEditor<Date> {
    // TODO: change UI
    override initializeEditor(prop: DG.Property): void {
        const input = ui.input.forProperty(prop, this.condition.value || new Date(), {
            onValueChanged: () => {
                this.condition.value = input.value;
                this.onChanged.next(this.condition);
            }
        });
        input.addCaption('');
        this.root.append(input.root);
    }
}

export class MultiValueConditionEditor extends BaseConditionEditor<any[]> {
    override initializeEditor(prop: DG.Property): void {
        const input = ui.input.string((this.condition.value || []).join(', '), {
            onValueChanged: () => {
                this.condition.value = input.value.split(',').map(s => s.trim()).filter(s => s.length > 0);
                this.onChanged.next(this.condition);
            }
        });
        input.addCaption('Comma-separated values');
        this.root.append(input.root);
    }
}

// Condition editor registry
export class ConditionEditorRegistry {
    private static instance: ConditionEditorRegistry;
    private editorRegistry: Map<string, new (prop: DG.Property, operator: string, initialCondition?: any) => BaseConditionEditor> = new Map();

    private constructor() {
        this.initializeDefaultEditors();
    }

    static getInstance(): ConditionEditorRegistry {
        if (!ConditionEditorRegistry.instance) {
            ConditionEditorRegistry.instance = new ConditionEditorRegistry();
        }
        return ConditionEditorRegistry.instance;
    }

    private initializeDefaultEditors(): void {
        // Register default editors by property type (with empty semType and operator)
        this.registerEditor(DG.TYPE.STRING, '', '', BaseConditionEditor);
        this.registerEditor(DG.TYPE.INT, '', '', BaseConditionEditor);
        this.registerEditor(DG.TYPE.FLOAT, '', '', BaseConditionEditor);
        this.registerEditor(DG.TYPE.DATE_TIME, '', '', BaseConditionEditor);

        // Register specialized editors for specific combinations
        this.registerEditor(DG.TYPE.STRING, '', 'in', MultiValueConditionEditor);
        this.registerEditor(DG.TYPE.DATE_TIME, '', 'between', DateBetweenConditionEditor);
    }

    // Register an editor for a specific combination
    registerEditor<T>(
        propertyType: string, 
        semType: string, 
        operator: string, 
        editorClass: new (prop: DG.Property, operator: string) => BaseConditionEditor<T>
    ): void {
        const key = this.createKey(propertyType, semType, operator);
        this.editorRegistry.set(key, editorClass);
    }

    // Create editor for a property and operator
    createEditor<T>(property: DG.Property, operator: string, initialCondition?: SimpleCondition<T>): BaseConditionEditor<T> {
        // Try to find a specialized editor first (propertyType + semType + operator)
        const specializedKey = this.createKey(property.propertyType, property.semType, operator);
        let EditorClass = this.editorRegistry.get(specializedKey);

        // If no specialized editor, try with empty semType (propertyType + operator)
        if (!EditorClass) {
            const typeOnlyKey = this.createKey(property.propertyType, '', operator);
            EditorClass = this.editorRegistry.get(typeOnlyKey);
        }

        // If still no specialized editor, use default editor for property type (propertyType only)
        if (!EditorClass) {
            const defaultKey = this.createKey(property.propertyType, '', '');
            EditorClass = this.editorRegistry.get(defaultKey);
        }

        // Fallback to base editor
        if (!EditorClass) {
            EditorClass = BaseConditionEditor;
        }

        return new EditorClass(property, operator, initialCondition);
    }

    private createKey(propertyType: string, semType: string, operator: string): string {
        return `${propertyType}:${semType}:${operator}`;
    }
}

export class OperatorRegistry {
    private static instance: OperatorRegistry;
    private typeOperators: Map<string, string[]> = new Map();
    private semTypeOperators: Map<string, string[]> = new Map();

    private constructor() {
        this.initializeDefaultOperators();
    }

    static getInstance(): OperatorRegistry {
        if (!OperatorRegistry.instance) {
            OperatorRegistry.instance = new OperatorRegistry();
        }
        return OperatorRegistry.instance;
    }

    private initializeDefaultOperators(): void {
        // Type operators
        this.registerTypeOperators(DG.TYPE.STRING, ['starts with', 'ends with', '=', '!=', 'in']);        
        this.registerTypeOperators(DG.TYPE.INT, ['>', '<', '>=', '<=', '=', '!=']);
        this.registerTypeOperators(DG.TYPE.FLOAT, ['>', '<', '>=', '<=', '=', '!=']);     
        this.registerTypeOperators(DG.TYPE.DATE_TIME, ['before', 'after', 'between']);

        // SemType operators (take precedence over type operators)
        this.registerSemTypeOperators(DG.SEMTYPE.MOLECULE, ['contains', 'is contained', '=', 'is similar']);
    }

    registerTypeOperators(propertyType: string, operators: string[]): void {
        this.typeOperators.set(propertyType, operators);
    }

    registerSemTypeOperators(semType: string, operators: string[]): void {
        this.semTypeOperators.set(semType, operators);
    }

    getOperatorsForProperty(property: DG.Property): string[] {
        // Check semType first, then fall back to type
        const semTypeOps = this.semTypeOperators.get(property.semType);
        if (semTypeOps) return semTypeOps;

        const typeOps = this.typeOperators.get(property.propertyType);
        return typeOps || [];
    }

    isOperatorSupported(property: DG.Property, operator: string): boolean {
        const operators = this.getOperatorsForProperty(property);
        return operators.includes(operator);
    }
}



export namespace Operators {
    export const CONTAINS = 'contains';
    export const IN = 'in';
    export const STARTS_WITH = 'starts with';
    export const ENDS_WITH = 'ends with';
    export const GT = '>';
    export const GTE = '>=';
    export const LT = '<';
    export const LTE = '<=';
    export const EQ = '=';
    export const NOT_EQ = '!=';
    export const BEFORE = 'before';
    export const AFTER = 'after';
    export const ON = 'on';
    export const BETWEEN = 'between';
    export const IS_SIMILAR = 'is similar';
    export const IS_CONTAINED = 'is contained';

    // Legacy properties - now use registry
    export const typeOperators: { [key: string]: string[] } = {
        [DG.TYPE.STRING]: [STARTS_WITH, ENDS_WITH, EQ, NOT_EQ, IN],
        [DG.TYPE.INT]: [GT, LT, GTE, LTE, EQ, NOT_EQ],
        [DG.TYPE.FLOAT]: [GT, LT, GTE, LTE, EQ, NOT_EQ],
        [DG.TYPE.DATE_TIME]: [BEFORE, AFTER, BETWEEN]
    };

    export const oneFieldOperators = [CONTAINS, STARTS_WITH, ENDS_WITH, GT, GTE, LT, LTE, EQ, NOT_EQ, BEFORE, AFTER, ON, IS_CONTAINED];

    export enum Logical {
        and = 'and',
        or = 'or',
    }

    export const semTypeOperators: { [key: string]: string[] } = {
        [DG.SEMTYPE.MOLECULE]: [CONTAINS, IS_CONTAINED, EQ, IS_SIMILAR]
    };
}

export function getConditionEditor<T = any>(property: DG.Property, operator: string, condition?: SimpleCondition<T>): BaseConditionEditor<T> {
    const operatorRegistry = OperatorRegistry.getInstance();
    const editorRegistry = ConditionEditorRegistry.getInstance();
    
    if (operatorRegistry.isOperatorSupported(property, operator)) {
        return editorRegistry.createEditor(property, operator, condition);
    } else {
        throw Error(`Operator '${operator}' is not supported for property type '${property.propertyType}'`);
    }
}

// Recursive UI builder for filter conditions
export function buildPropertyFilterUI(
    condition: ComplexCondition,
    properties: DG.Property[],
    onValueChange: () => void,
    onStructureChange: () => void,
    parentCondition?: ComplexCondition,
    nestingLevel: number = 0
): HTMLElement {

    if (!condition) {
        condition = {
            logicalOperator: Operators.Logical.and,
            conditions: []
        };
    }

    const removeOrReplaceCondition = (condToRemove: Expression, condArray?: Expression[], replaceCond?: ComplexCondition) => {
        if (condArray) {
            const condIdx = condArray.indexOf(condToRemove);
            if (condIdx !== -1)
                replaceCond ? condArray[condIdx] = replaceCond : condArray.splice(condIdx, 1)
        };
    }

    const createPropertyFilterContainer = (cond: SimpleCondition, parentCondition?: ComplexCondition) => {

        const createFilter = () => {
            ui.empty(operatorInputDiv);
            const property = properties.find(prop => prop.name === fieldChoiceInput.value!);
            if (property) {
                const registry = OperatorRegistry.getInstance();
                const operators = registry.getOperatorsForProperty(property);
                const initValue = cond.operator && cond.operator !== '' ? cond.operator : operators[0]
                const operatorsInput = ui.input.choice('', {
                    value: initValue,
                    items: operators,
                    onValueChanged: () => {
                        createEditor(property, operatorsInput.value!, cond);
                    },
                    nullable: false
                });
                operatorInputDiv.append(operatorsInput.root);
                createEditor(property, initValue, cond);
            }
        }

        const createEditor = (property: DG.Property, operator: string, condition: SimpleCondition) => {
            ui.empty(criteriaDiv);
            const editor = getConditionEditor(property, operator, condition);
            criteriaDiv.append(editor.root);
        }

        cond.field ??= (properties[0]?.name || '')

        const fieldChoiceInput = ui.input.choice('', {
            items: properties.map(prop => prop.name),
            value: cond.field,
            nullable: false,
            onValueChanged: () => {
                cond.field = fieldChoiceInput.value!;
                cond.value = undefined;
                cond.operator = '';
                createFilter();
                onValueChange();
            }
        });

        const deleteFieldIcon = ui.icons.delete(() => {
            container.removeChild(filterContainer);
            removeOrReplaceCondition(cond, parentCondition ? parentCondition.conditions : condition.conditions);
            onStructureChange();
        });

        const addNestedConditionIcon = ui.icons.add(() => {
            const parentCond = { logicalOperator: Operators.Logical.and, conditions: [cond] };
            removeOrReplaceCondition(cond, parentCondition ? parentCondition.conditions : condition.conditions, parentCond);
            onStructureChange();

        }, 'Add nested condition');

        const operatorInputDiv = ui.div();
        const criteriaDiv = ui.div();
        const filterContainer = ui.divH([fieldChoiceInput.root, operatorInputDiv, criteriaDiv,
            ui.divH([deleteFieldIcon, addNestedConditionIcon], 'add-delete-icons')], 'revvity-signals-filter-inputs');
        container.appendChild(filterContainer);
        createFilter();
        onValueChange();
    }

    // adding filter field
    const addFieldIcon = ui.icons.add(() => {
        const conditionForFilter: SimpleCondition = {
            field: '',
            value: '',
            operator: ''
        };
        condition.conditions?.push(conditionForFilter);
        condition.conditions.length < 2 ? container.classList.remove('property-query-builder-multipe-filelds') :
            container.classList.add('property-query-builder-multipe-filelds');
        onStructureChange();
    }, 'Add field to the condition');

    // AND/OR operators handling
    condition.logicalOperator ??= Operators.Logical.and;
    const logicalOperatorIcon = ui.button(condition.logicalOperator, () => {
        condition.logicalOperator = condition.logicalOperator === Operators.Logical.or ? Operators.Logical.and : Operators.Logical.or;
        logicalOperatorIcon.innerText = condition.logicalOperator.toUpperCase();
        onValueChange();
    }, 'Logical operator');

    logicalOperatorIcon.classList.add('property-query-builder-and-or-operator');

    const addDeleteIconsDiv = ui.divH([addFieldIcon], 'nested-add-delete-icons');
    const iconsDiv = ui.divH([logicalOperatorIcon, addDeleteIconsDiv], 'query-builder-nested-condition-operators');
    const container = ui.divV([iconsDiv], 'revvity-signals-search-query-operator');
    

    if (nestingLevel > 0) {
        container.setAttribute('data-nesting-level', nestingLevel.toString());
        const color = generateNestingColor(nestingLevel);
        container.style.setProperty('--nesting-color', color);
    }

    //delete icon for nested conditions to remove the whole nested condition at once
    if (parentCondition) {
        const deleteCondition = ui.icons.delete(() => {
            removeOrReplaceCondition(condition, parentCondition.conditions);
            onStructureChange();
        });
        addDeleteIconsDiv.append(deleteCondition);
        onValueChange();
    }

    for (const nestedCondition of condition.conditions!) {
        if ('field' in nestedCondition)
            createPropertyFilterContainer(nestedCondition);
        else {
            const nestedContainer = buildPropertyFilterUI(nestedCondition as ComplexCondition, properties, onValueChange, onStructureChange,
                condition, nestingLevel + 1);
            container.append(nestedContainer);
        }
    }

    condition.conditions.length < 2 ? container.classList.remove('property-query-builder-multipe-filelds') :
        container.classList.add('property-query-builder-multipe-filelds');

    return container;
}

// Main UI builder function for property filters
export function buildPropertyFilterForm(properties: DG.Property[], initialCondition?: ComplexCondition): HTMLElement {
    // Initialize with provided condition or create a default logical condition
    const rootCondition: ComplexCondition = initialCondition || {
            logicalOperator: Operators.Logical.and,
            conditions: []
    };

    // Preview area
    const jsonPreview = document.createElement('textarea');
    jsonPreview.readOnly = true;
    jsonPreview.classList.add('revvity-signals-search-query-preview');

    function updatePreview() {
        console.log(rootCondition);
        jsonPreview.value = JSON.stringify({ filter: rootCondition }, null, 2);
    }

    function rebuildUI() {
        updatePreview();
        builderDiv.replaceWith(builderDiv = buildPropertyFilterUI(rootCondition, properties, updatePreview, rebuildUI, undefined, 0));
    }

    let builderDiv = buildPropertyFilterUI(rootCondition, properties, updatePreview, rebuildUI, undefined, 0);

    const submitBtn = ui.button('Apply Filters', async () => {
        // Here you would apply the filters
        grok.shell.info('Filters applied: ' + JSON.stringify({ filter: rootCondition }));
    });

    const formContainer = ui.divV([
        ui.h3('Property Filters'),
        builderDiv,
        submitBtn,
    ], 'revvity-signals-search');

    const previewContainer = ui.divV([
        ui.h3('Filter Preview'),
        jsonPreview
    ]);

    const mainContainer = ui.splitV([
        formContainer,
        previewContainer
    ], {style: {width: '100%'}}, true);

    updatePreview();
    return mainContainer;
}