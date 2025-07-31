import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// Interface for filter conditions matching the query structure
export interface FilterCondition {
    field?: string;           // Field name
    operator?: string | 'and' | 'or';  // Operator like "=", ">", "IN", "IS SIMILAR", "AND", "OR", etc.
    value?: any;             // Value for the condition
    threshold?: number;       // For similarity searches
    conditions?: FilterCondition[]; // For nested conditions
}

export const defaultFilters: { [key in DG.TYPE]?: (prop: DG.Property) => HTMLElement } = {
 [DG.TYPE.DATE_TIME]: (prop: DG.Property) => new DateFilter(prop).root,
 [DG.TYPE.STRING]: (prop: DG.Property) => new StringFilter(prop).root,
 [DG.TYPE.FLOAT]: (prop: DG.Property) => new FloatFilter(prop).root,
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
        [DG.TYPE.DATE_TIME]: [BEFORE, AFTER, ON, BETWEEN]
    };

    export enum ConditionOperators {
        property = 'property',
        and = 'and',
        or = 'or',
    }

    export const semTypeOperators: { [key: string]: string[] } = {
        [DG.SEMTYPE.MOLECULE]: [CONTAINS, IS_CONTAINED, EQ, IS_SIMILAR]
    };

    export function getCondition(fieldName: string, op: string, values: any[]): string | null {
        if (!values || !values.length)
            return null;

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
        return null;
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
    onStructureChange: () => void
): HTMLElement {
    
    // Handle empty or undefined condition by creating a default condition
    if (!condition) {
        condition = {};
    }
    
    // Create logical operator choice input
    const conditionInput = ui.input.choice('Condition', {
        items: Object.values(Operators.ConditionOperators),
        value: condition.operator || Operators.ConditionOperators.property,
        nullable: false,
        onValueChanged: () => {
            condition.operator = conditionInput.value!;
            onValueChange();
            onStructureChange();
        }
    });

    const container = ui.divV([conditionInput.root], 'revvity-signals-search-query-operator');

    // If PROPERTY is selected, show field selection
    if (conditionInput.value === Operators.ConditionOperators.property) {
        const fieldChoiceInput = ui.input.choice('', {
            items: properties.map(prop => prop.name),
            value: properties[0]?.name || '', // Set default to first property
            nullable: false, // Make it non-nullable
            onValueChanged: () => {
                const selectedField = fieldChoiceInput.value!;
                if (selectedField) {
                    condition.field = selectedField;
                    const property = properties.find(prop => prop.name === selectedField);
                    if (property) {
                        const filterFactory = defaultFilters[property.propertyType];
                        if (filterFactory) {
                            // Clear existing field UI and add new filter UI
                            const existingFieldUI = filterContainer.querySelector('.field-filter-ui');
                            if (existingFieldUI) {
                                existingFieldUI.remove();
                            }
                            
                            const fieldFilterUI = filterFactory(property);
                            fieldFilterUI.classList.add('field-filter-ui');
                            filterContainer.appendChild(fieldFilterUI);
                        }
                    }
                }
                onValueChange();
            }
        });
        
        // Create horizontal container for field choice and filter UI
        const filterContainer = ui.divH([fieldChoiceInput.root], 'revvity-signals-filter-inputs');
        container.appendChild(filterContainer);
        
        // Trigger the field selection to initialize the field filter UI
        if (fieldChoiceInput.value) {
            condition.field = fieldChoiceInput.value;
            const property = properties.find(prop => prop.name === fieldChoiceInput.value);
            if (property) {
                const filterFactory = defaultFilters[property.propertyType];
                if (filterFactory) {
                    const fieldFilterUI = filterFactory(property);
                    fieldFilterUI.classList.add('field-filter-ui');
                    filterContainer.appendChild(fieldFilterUI);
                }
            }
            onValueChange();
        }
    }
    
    // If AND or OR is selected, create nested conditions with light border
    else if (conditionInput.value === Operators.ConditionOperators.and || conditionInput.value === Operators.ConditionOperators.or) {
        // Ensure conditions array exists and is properly initialized
        if (!condition.conditions) {
            condition.conditions = [];
        }
        const conditions = condition.conditions;
        
        // Create container for nested conditions with light border
        const conditionsContainer = ui.divV([], 'revvity-signals-search-query-operator');
        
        // Function to rebuild the conditions UI for this specific level
        const rebuildConditionsUI = () => {
            ui.empty(conditionsContainer);
            
            // Add existing conditions
            conditions.forEach((childCondition, idx) => {
                const conditionUI = buildPropertyFilterUI(childCondition, properties, onValueChange, () => {
                    // When nested condition structure changes, rebuild this level and notify parent
                    rebuildConditionsUI();
                    onStructureChange();
                });
                const conditionContainer = ui.divH([
                    conditionUI,
                    ui.icons.delete(() => {
                        conditions.splice(idx, 1);
                        rebuildConditionsUI();
                        onStructureChange();
                    }, 'Remove condition')
                ]);
                conditionsContainer.appendChild(conditionContainer);
            });
            
            // Add plus icon to add new condition
            const addConditionIcon = ui.icons.add(() => {
                const newCondition: FilterCondition = {
                    operator: Operators.ConditionOperators.property // Initialize with NONE so it shows field selection
                };
                conditions.push(newCondition);
                // Add a temporary visual indicator
                rebuildConditionsUI();
                onStructureChange();
            });
            addConditionIcon.style.cssText = 'padding: 2px 6px; font-size: 12px; color: #2083D5;';
            conditionsContainer.appendChild(addConditionIcon);
        };
        
        // Initial build
        rebuildConditionsUI();
        
        container.appendChild(conditionsContainer);
    }

    return container;
}

// Main UI builder function for property filters
export function buildPropertyFilterForm(properties: DG.Property[], initialCondition?: FilterCondition): HTMLElement {
    // Initialize with provided condition or create a default logical condition
    const rootCondition: FilterCondition = initialCondition || {
        operator: Operators.ConditionOperators.property,
    };

    // Preview area
    const jsonPreview = document.createElement('textarea');
    jsonPreview.readOnly = true;
    jsonPreview.classList.add('revvity-signals-search-query-preview');

    function updatePreview() {
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

export class StringFilter {
    operators = Operators.typeOperators[DG.TYPE.STRING];
    root: HTMLElement;
    constructor(private prop: DG.Property, private value: string = '',
        private operator: string = Operators.EQ) {
        const operatorChoice = ui.input.choice('', {
            value: this.operator,
            items: Object.values(this.operators),
            nullable: false,
            onValueChanged: () => {
                this.operator = operatorChoice.value!;
            }
        });

        const stringInput = ui.input.forProperty(this.prop, null, {
            value: this.value,
            onValueChanged: () => {
                this.value = stringInput.value!;
            }
        });
        stringInput.addCaption(''); // Hide the property name label

        this.root = ui.divH([
            operatorChoice.root,
            stringInput.root
        ]);

    }

    getCondition() {
        return Operators.getCondition(this.prop.name, this.operator, [this.value]);
    }
}

export class FloatFilter {
    operators = Operators.typeOperators[DG.TYPE.FLOAT];
    root: HTMLElement;
    constructor(private prop: DG.Property, private value: number = 0,
        private operator: string = Operators.EQ) {
        const operatorChoice = ui.input.choice('', {
            value: this.operator,
            items: Object.values(this.operators),
            nullable: false,
            onValueChanged: () => {
                this.operator = operatorChoice.value!;
            }
        });

        const floatInput = ui.input.forProperty(this.prop, null, {
            value: this.value,
            onValueChanged: () => {
                this.value = floatInput.value!;
            }
        });
        floatInput.addCaption(''); // Hide the property name label

        this.root = ui.divH([
            operatorChoice.root,
            floatInput.root
        ]);

    }

    getCondition() {
        return Operators.getCondition(this.prop.name, this.operator, [this.value]);
    }
}

export class DateFilter {
    operators = Operators.typeOperators[DG.TYPE.DATE_TIME];
    root: HTMLElement;
    datesDiv = ui.divH([]);
    constructor(private prop: DG.Property, private date1: string = '', private date2: string = '',
        private operator: string = Operators.BEFORE) {
        const operatorChoice = ui.input.choice('', {
            value: this.operator,
            items: Object.values(this.operators),
            nullable: false,
            onValueChanged: () => {
                this.operator = operatorChoice.value!;
                this.createFilter();
            }
        });

        this.createFilter();
        this.root = ui.divH([
            operatorChoice.root,
            this.datesDiv,
        ]);

    }

    createFilter() {
        ui.empty(this.datesDiv);
        switch (this.operator) {
            case Operators.BEFORE:
            case Operators.AFTER:
                const dateInput = ui.input.forProperty(this.prop, null, {
                    value: this.date1,
                    onValueChanged: () => {
                        this.date1 = dateInput.value!;
                    }
                });
                dateInput.addCaption(''); // Hide the property name label
                this.datesDiv.append(dateInput.root);
                return;
            case Operators.BETWEEN:
                const dateInput1 = ui.input.forProperty(this.prop, null, {
                    value: this.date1,
                    onValueChanged: () => {
                        this.date1 = dateInput1.value!;
                    }
                });
                dateInput1.addCaption(''); // Hide the property name label
                const dateInput2 = ui.input.forProperty(this.prop, null, {
                    value: this.date1,
                    onValueChanged: () => {
                        this.date2 = dateInput2.value!;
                    }
                });
                dateInput2.addCaption(''); // Hide the property name label
                this.datesDiv.append(dateInput1.root, dateInput2.root);
                return;
            default:
                this.datesDiv.append(ui.div('Unknown date operator'));

        }
    }

    getCondition() {
        return Operators.getCondition(this.prop.name, this.operator, [this.date1, this.date2]);
    }
}
