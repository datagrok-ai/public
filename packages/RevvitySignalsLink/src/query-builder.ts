import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Subject } from 'rxjs';
import dayjs from 'dayjs';
import { IProperty } from 'datagrok-api/dg';

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

const dict = {
    
}

export function getConditionEditor<T = any>(property: DG.Property, operator: string, condition?: SimpleCondition<T>): BaseConditionEditor<T> {
    if (Operators.oneFieldOperators.includes(operator))
        return new BaseConditionEditor<T>(property, operator, condition);
    else 
        throw Error(`Unknown operator`);
}

// Recursive UI builder for filter conditions
export function buildPropertyFilterUI(
    condition: ComplexCondition,
    properties: DG.Property[],
    onValueChange: () => void,
    onStructureChange: () => void,
    parentCondition?: ComplexCondition
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
                const operators = Operators.typeOperators[property.propertyType];
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
            deleteFieldIcon, addNestedConditionIcon], 'revvity-signals-filter-inputs');
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
    if (!parentCondition)
        logicalOperatorIcon.classList.add('property-query-builder-parent-level');

    const iconsDiv = ui.divH([addFieldIcon, logicalOperatorIcon], 'query-builder-nested-condition-operators');
    const container = ui.divV([iconsDiv], 'revvity-signals-search-query-operator');

    //delete icon for nested conditions to remove the whole nested condition at once
    if (parentCondition) {
        const deleteCondition = ui.icons.delete(() => {
            removeOrReplaceCondition(condition, parentCondition.conditions);
            onStructureChange();
        });
        iconsDiv.append(deleteCondition);
        onValueChange();
    }

    for (const nestedCondition of condition.conditions!) {
        if ('field' in nestedCondition)
            createPropertyFilterContainer(nestedCondition);
        else {
            const nestedContainer = buildPropertyFilterUI(nestedCondition as ComplexCondition, properties, onValueChange, onStructureChange,
                parentCondition ?? condition);
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
        builderDiv.replaceWith(builderDiv = buildPropertyFilterUI(rootCondition, properties, updatePreview, rebuildUI));
    }

    let builderDiv = buildPropertyFilterUI(rootCondition, properties, updatePreview, rebuildUI);

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