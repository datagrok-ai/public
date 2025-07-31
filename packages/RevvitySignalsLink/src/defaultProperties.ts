import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// Interface for filter conditions matching the query structure
export interface FilterCondition {
    field?: string;           // Field name
    operator?: string | 'and' | 'or';  // Operator like "=", ">", "IN", "IS SIMILAR", "AND", "OR", etc.
    values?: any[];             // Value for the condition
    threshold?: number;       // For similarity searches
    conditions?: FilterCondition[]; // For nested conditions
}

export const defaultFilters: { [key in DG.TYPE]?: (prop: DG.Property, cond: FilterCondition) => PropertyFilter } = {
 [DG.TYPE.DATE_TIME]: (prop: DG.Property, cond: FilterCondition) => new DateFilter(prop, cond),
 [DG.TYPE.STRING]: (prop: DG.Property, cond: FilterCondition) => new StringFilter(prop, cond),
 [DG.TYPE.FLOAT]: (prop: DG.Property, cond: FilterCondition) => new FloatFilter(prop, cond),
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

    export enum Logical {
        and = 'and',
        or = 'or',
    }

    export const semTypeOperators: { [key: string]: string[] } = {
        [DG.SEMTYPE.MOLECULE]: [CONTAINS, IS_CONTAINED, EQ, IS_SIMILAR]
    };

    export function getStringCondition(fieldName: string, op: string, values: any[]): string {
        if (!values || !values.length)
            return '';

        switch (op) {
            case STARTS_WITH: return `starts_with(${fieldName}, '${values[0]}')`;
            case ENDS_WITH: return `${fieldName} like '%${values[0]}'`;
            case LT: return `${fieldName} < ${values[0]}`;
            case GT: return `${fieldName} > ${values[0]}`;
            case LTE: return `${fieldName} <= ${values[0]}`;
            case GTE: return `${fieldName} >= ${values[0]}`;
            case EQ: return `${fieldName} = ${values[0]}`;
            case NOT_EQ: return `${fieldName} != ${values[0]}`;
            case IN: return `${fieldName} in(${values.join(',')})`;
            case BEFORE: return `${fieldName} < ${values[0]}`;
            case AFTER: return `${fieldName} > ${values[0]}`;
            case BETWEEN: return `${fieldName} between(${values[0]}, ${values[1]})`;
        }
        return '';
    }
}

export interface DateFilterParams {
    operator: string,
    date1: string;
    date2: string;
}

export interface FilterParams {
    operator: string,
    value: any
}


// Recursive UI builder for filter conditions
export function buildPropertyFilterUI(
    condition: FilterCondition,
    properties: DG.Property[],
    onValueChange: () => void,
    onStructureChange: () => void,
    parentCondition?: FilterCondition
): HTMLElement {

    if (!condition) {
        condition = {};
    }

    const removeOrReplaceCondition = (condToRemove: FilterCondition, condArray?: FilterCondition[], replaceCond?: FilterCondition) => {
        if (condArray) {
            const condIdx = condArray.indexOf(condToRemove);
            if (condIdx !== -1)
                replaceCond ? condArray[condIdx] = replaceCond : condArray.splice(condIdx, 1)
        };
    }

    const createPropertyFilterContainer = (cond: FilterCondition, parentCondition?: FilterCondition) => {

        const createFilter = () => {
            const property = properties.find(prop => prop.name === fieldChoiceInput.value!);
            if (property) {
                const filterFactory = defaultFilters[property.propertyType];
                if (filterFactory) {
                    // Clear existing field UI and add new filter UI
                    ui.empty(filterInputsContainer);
                    const fieldFilter = filterFactory(property, cond);
                    filterInputsContainer.appendChild(fieldFilter.root);
                }
            }
        }

        cond.field ??= (properties[0]?.name || '')

        const fieldChoiceInput = ui.input.choice('', {
            items: properties.map(prop => prop.name),
            value: cond.field,
            nullable: false,
            onValueChanged: () => {
                cond.field = fieldChoiceInput.value!;
                cond.values = undefined;
                cond.operator = undefined;
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
            const parentCond = { operator: Operators.Logical.and, conditions: [cond] };
            removeOrReplaceCondition(cond, parentCondition ? parentCondition.conditions : condition.conditions, parentCond);
            onStructureChange();

        }, 'Add nested condition');

        const filterInputsContainer = ui.div();
        const filterContainer = ui.divH([fieldChoiceInput.root, filterInputsContainer,
            deleteFieldIcon, addNestedConditionIcon], 'revvity-signals-filter-inputs');
        container.appendChild(filterContainer);
        createFilter();
    }

    // adding filter field
    const addFieldIcon = ui.icons.add(() => {
        const conditionForFilter = {};
        //if we are handling top level condition and adding second field - restructure so that we put conditions to array and add 'AND' condition by default
        if (!parentCondition && !condition.conditions?.length) {
            const firstCondition = Object.assign({}, condition);
            condition.operator = Operators.Logical.and;
            condition.field = undefined;
            condition.values = undefined;
            condition.conditions = [firstCondition, conditionForFilter];
            container.classList.add('property-query-builder-multipe-filelds');

        } else {
            condition.conditions?.push(conditionForFilter);
        }
        onStructureChange();
    }, 'Add field to the condition');

    // AND/OR operators handling
    let and = true;
    const logicalOperatorIcon = ui.button('AND', () => {
        and = !and;
        condition.operator = !and ? Operators.Logical.or : Operators.Logical.and;
        logicalOperatorIcon.innerText = !and ? 'OR' : 'AND';
        onValueChange();
    }, 'Logical operator');
    logicalOperatorIcon.classList.add('property-query-builder-and-or-operator');
    if (!parentCondition)
        logicalOperatorIcon.classList.add('property-query-builder-parent-level');


    const iconsDiv = ui.divH([addFieldIcon, logicalOperatorIcon]);
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

    if (condition.operator !== Operators.Logical.and && condition.operator !== Operators.Logical.or) {
        createPropertyFilterContainer(condition);
    }
    // If AND or OR is selected, create nested conditions
    else {
        if (condition.conditions) {
            for (const nestedCondition of condition.conditions!) {
                if (nestedCondition.operator !== Operators.Logical.and && nestedCondition.operator !== Operators.Logical.or)
                    createPropertyFilterContainer(nestedCondition);
                else {
                    const nestedContainer = buildPropertyFilterUI(nestedCondition, properties, onValueChange, onStructureChange,
                        parentCondition ?? condition);
                    container.append(nestedContainer);
                }
            }

            condition.conditions.length < 2 ? container.classList.remove('property-query-builder-multipe-filelds') :
                container.classList.add('property-query-builder-multipe-filelds');
        }
    }

    return container;
}

// Main UI builder function for property filters
export function buildPropertyFilterForm(properties: DG.Property[], initialCondition?: FilterCondition): HTMLElement {
    // Initialize with provided condition or create a default logical condition
    const rootCondition: FilterCondition = initialCondition || {};

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

export abstract class PropertyFilter {
  root: HTMLDivElement = ui.div();
  condition: FilterCondition;
  prop: DG.Property;

  abstract getStringCondition(): string;

  protected constructor(prop: DG.Property, condition: FilterCondition) {
    this.prop = prop;
    this.condition = condition;
  }
}

export class StringFilter extends PropertyFilter {
    operators = Operators.typeOperators[DG.TYPE.STRING];

    constructor(prop: DG.Property, condition: FilterCondition) {
        super(prop, condition);
        condition.operator ??= Operators.EQ;
        const operatorChoice = ui.input.choice('', {
            value: condition.operator,
            items: Object.values(this.operators),
            nullable: false,
            onValueChanged: () => {
                condition.operator = operatorChoice.value!;
            }
        });
        
        if (!condition.values || !Array.isArray(condition.values) || !condition.values.length)
            condition.values = [null];

        const stringInput = ui.input.forProperty(this.prop, null, {
            value: condition.values[0],
            onValueChanged: () => {
                condition.values![0] = stringInput.value!;
            }
        });
        stringInput.addCaption(''); // Hide the property name label

        this.root = ui.divH([
            operatorChoice.root,
            stringInput.root
        ]);

    }

    getStringCondition() {
        return Operators.getStringCondition(this.prop.name, this.condition.operator!, this.condition.values!);
    }
}

export class FloatFilter extends PropertyFilter {
    operators = Operators.typeOperators[DG.TYPE.FLOAT];

    constructor(prop: DG.Property, condition: FilterCondition) {
        super(prop, condition);
        condition.operator ??= Operators.EQ;
        const operatorChoice = ui.input.choice('', {
            value: condition.operator,
            items: Object.values(this.operators),
            nullable: false,
            onValueChanged: () => {
                condition.operator = operatorChoice.value!;
            }
        });

        if (!condition.values || !Array.isArray(condition.values) || !condition.values.length)
            condition.values = [null];

        const floatInput = ui.input.forProperty(this.prop, null, {
            value: condition.values[0],
            onValueChanged: () => {
                condition.values![0] = floatInput.value!;
            }
        });
        floatInput.addCaption(''); // Hide the property name label

        this.root = ui.divH([
            operatorChoice.root,
            floatInput.root
        ]);

    }

    getStringCondition() {
        return Operators.getStringCondition(this.prop.name, this.condition.operator!, this.condition.values!);
    }
}

export class DateFilter extends PropertyFilter {
    operators = Operators.typeOperators[DG.TYPE.DATE_TIME];
    datesDiv = ui.divH([]);

    constructor(prop: DG.Property, condition: FilterCondition) {
        super(prop, condition);
        condition.operator ??= Operators.BEFORE;
        const operatorChoice = ui.input.choice('', {
            value: condition.operator,
            items: Object.values(this.operators),
            nullable: false,
            onValueChanged: () => {
                condition.operator = operatorChoice.value!;
                this.createFilter();
            }
        });

        if (!condition.values || !Array.isArray(condition.values) || condition.values.length < 2)
            condition.values = [Date.now(), Date.now()];

        this.createFilter();
        this.root = ui.divH([
            operatorChoice.root,
            this.datesDiv,
        ]);

    }

    createFilter() {
        ui.empty(this.datesDiv);
        switch (this.condition.operator) {
            case Operators.BEFORE:
            case Operators.AFTER:
                const dateInput = ui.input.forProperty(this.prop, null, {
                    value: this.condition.values![0],
                    onValueChanged: () => {
                        this.condition.values![0] = dateInput.value!;
                    }
                });
                dateInput.addCaption(''); // Hide the property name label
                this.datesDiv.append(dateInput.root);
                return;
            case Operators.BETWEEN:
                const dateInput1 = ui.input.forProperty(this.prop, null, {
                    value: this.condition.values![0],
                    onValueChanged: () => {
                        this.condition.values![0] = dateInput1.value!;
                    }
                });
                dateInput1.addCaption(''); // Hide the property name label
                const dateInput2 = ui.input.forProperty(this.prop, null, {
                    value: this.condition.values![1],
                    onValueChanged: () => {
                        this.condition.values![1] = dateInput2.value!;
                    }
                });
                dateInput2.addCaption(''); // Hide the property name label
                this.datesDiv.append(dateInput1.root, dateInput2.root);
                return;
            default:
                this.datesDiv.append(ui.div('Unknown date operator'));

        }
    }

    getStringCondition() {
        return Operators.getStringCondition(this.prop.name, this.condition.operator!, this.condition.values!);
    }

}