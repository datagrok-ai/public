import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Subject } from 'rxjs';
import dayjs, { Dayjs } from 'dayjs';
import { IProperty } from 'datagrok-api/dg';


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

export const SUGGESTIONS_FUNCTION = 'suggestionsFunction';

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

    export enum Logical {
        and = 'and',
        or = 'or',
    }
}

export interface SimpleCondition<T = any> {
    field: string;
    operator: string;
    value: T;
}

export interface ComplexCondition {
    logicalOperator: Operators.Logical;
    conditions: Expression[];
}

export type Expression<T = any> = SimpleCondition<T> | ComplexCondition;


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
                if (showSuggestions && !menuItemClicked) {
                    addSuggestions();
                }
                menuItemClicked = false;
            }
        });
        input.addCaption('');
        this.root.append(input.root);
        const showSuggestions = prop.type === DG.TYPE.STRING && prop.options[SUGGESTIONS_FUNCTION];
        if (showSuggestions) {
           suggestionMenuKeyNavigation(this.root);
        }
        let menuItemClicked = false;
        const addSuggestions = async () => {
            if (input.value) {
                const suggestions: Array<string> = await prop.options[SUGGESTIONS_FUNCTION](input.value);
                if (suggestions.length === 1 && suggestions[0] === input.value)
                    return;
                const suggestionsMenu = DG.Menu.popup();
                suggestions.forEach((s) => {
                    suggestionsMenu!.item(s, () => {
                        input.value = s;
                        menuItemClicked = true;
                    });
                });
                suggestionsMenu.show({element: input.root, y: input.root.offsetHeight});
            }
        }
    }
}

export class DateBetweenConditionEditor extends BaseConditionEditor<Date[]> {
    override initializeEditor(prop: DG.Property): void {
        if (!this.condition.value || !Array.isArray(this.condition.value)) {
            this.condition.value = [new Date(), new Date()];
            this.onChanged.next(this.condition);
        }
        const date1 = ui.input.date('', {
            value: dayjs(this.condition.value[0]),
            onValueChanged: () => {
                this.condition.value[0] = date1.value!.toDate();
                this.onChanged.next(this.condition);
            }
        });
        const date2 = ui.input.date('', {
            value: dayjs(this.condition.value[1]),
            onValueChanged: () => {
                this.condition.value[1] = date2.value!.toDate();
                this.onChanged.next(this.condition);
            }
        });
        this.root.append(ui.divH([date1.root, date2.root], { style: { gap: '10px' } }));
    }
}

export class MoleculeConditionEditor extends BaseConditionEditor<string> {
    override initializeEditor(prop: DG.Property): void {
        //if we swith from similarity input to standard molecule input - need to modify value
        if (this.condition.value && typeof this.condition.value !== 'string') {
            if ('molecule' in this.condition.value)
                 this.condition.value = (this.condition.value as MoleculeSimilarity).molecule;
            else
                this.condition.value = '';
        }
        const input = ui.input.molecule('', {
            value: this.condition.value,
            onValueChanged: () => {
                this.condition.value = input.value;
                this.onChanged.next(this.condition);
            }
        });
        this.root.append(input.root);
        //to adjust the filter width considering molecule input
        this.root.style.paddingTop = '15px';
        this.root.style.paddingBottom = '15px';
    }
}

export type MoleculeSimilarity = {
    molecule: string,
    threshold: number
}
export class MoleculeSimilarityConditionEditor extends BaseConditionEditor<MoleculeSimilarity> {
    override initializeEditor(prop: DG.Property): void {
        //in case we switch from some other operator, where value is a string (contains, is contained)
        if (!this.condition.value || typeof this.condition.value === 'string') {
            this.condition.value = {molecule: this.condition.value ?? '', threshold: 0.7};
        }
        const moleculeInput = ui.input.molecule('', {
            onValueChanged: () => {
                this.condition.value.molecule = moleculeInput.value;
                this.onChanged.next(this.condition);
            }
        });
        const tresholdInput = ui.input.float('', {
            value: this.condition.value.threshold,
            min: 0,
            max: 1,
            step: 0.01
        });
        DG.debounce(tresholdInput.onChanged, 300).subscribe(() => {
            this.condition.value.threshold = tresholdInput.value!;
            this.onChanged.next(this.condition);
        })
        this.root.append(ui.divV([moleculeInput.root, tresholdInput.root]));
        this.root.style.paddingTop = '15px';
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

// Registry for both operators and editors
export class ConditionRegistry {
    private static instance: ConditionRegistry;
    private typeOperators: Map<string, string[]> = new Map();
    private semTypeOperators: Map<string, string[]> = new Map();
    private editorRegistry: Map<string, new (prop: DG.Property, operator: string, initialCondition?: any) => BaseConditionEditor> = new Map();

    private constructor() {
        this.initializeDefaultOperators();
        this.initializeDefaultEditors();
    }

    static getInstance(): ConditionRegistry {
        if (!ConditionRegistry.instance) {
            ConditionRegistry.instance = new ConditionRegistry();
        }
        return ConditionRegistry.instance;
    }

    private initializeDefaultOperators(): void {
        // Type operators
        this.registerTypeOperators(DG.TYPE.STRING, [Operators.STARTS_WITH, Operators.ENDS_WITH, Operators.EQ, Operators.NOT_EQ, Operators.IN]);
        this.registerTypeOperators(DG.TYPE.INT, [Operators.GT, Operators.LT, Operators.GTE, Operators.LTE, Operators.EQ, Operators.NOT_EQ]);
        this.registerTypeOperators(DG.TYPE.FLOAT, [Operators.GT, Operators.LT, Operators.GTE, Operators.LTE, Operators.EQ, Operators.NOT_EQ]);
        this.registerTypeOperators(DG.TYPE.DATE_TIME, [Operators.BEFORE, Operators.AFTER, Operators.BETWEEN]);

        // SemType operators (take precedence over type operators)
        this.registerSemTypeOperators(DG.SEMTYPE.MOLECULE, [Operators.CONTAINS, Operators.IS_CONTAINED, Operators.EQ, Operators.IS_SIMILAR]);
    }

    private initializeDefaultEditors(): void {
        // Register default editors by property type (with empty semType and operator)
        this.registerEditor(DG.TYPE.STRING, '', '', BaseConditionEditor);
        this.registerEditor(DG.TYPE.INT, '', '', BaseConditionEditor);
        this.registerEditor(DG.TYPE.FLOAT, '', '', BaseConditionEditor);
        this.registerEditor(DG.TYPE.DATE_TIME, '', '', BaseConditionEditor);

        // Register specialized editors for specific combinations
        this.registerEditor(DG.TYPE.STRING, '', Operators.IN, MultiValueConditionEditor);
        this.registerEditor(DG.TYPE.DATE_TIME, '', Operators.BETWEEN, DateBetweenConditionEditor);
        this.registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, Operators.CONTAINS, MoleculeConditionEditor);
        this.registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, Operators.IS_CONTAINED, MoleculeConditionEditor);
        this.registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, Operators.EQ, MoleculeConditionEditor);
        this.registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, Operators.IS_SIMILAR, MoleculeSimilarityConditionEditor);
    }

    // Operator registration methods
    registerTypeOperators(propertyType: string, operators: string[]): void {
        this.typeOperators.set(propertyType, operators);
    }

    registerSemTypeOperators(semType: string, operators: string[]): void {
        this.semTypeOperators.set(semType, operators);
    }

    // Get operators for a property
    getOperatorsForProperty(property: DG.Property): string[] {
        // Check semType first, then fall back to type
        const semTypeOps = this.semTypeOperators.get(property.semType);
        if (semTypeOps) return semTypeOps;

        const typeOps = this.typeOperators.get(property.propertyType);
        return typeOps || [];
    }

    // Check if operator is supported
    isOperatorSupported(property: DG.Property, operator: string): boolean {
        const operators = this.getOperatorsForProperty(property);
        return operators.includes(operator);
    }

    // Editor registration methods
    registerEditor<T>(
        propertyType: string,
        semType: string,
        operator: string,
        editorClass: new (prop: DG.Property, operator: string) => BaseConditionEditor<T>
    ): void {
        const key = this.createEditorKey(propertyType, semType, operator);
        this.editorRegistry.set(key, editorClass);
    }

    // Create editor for a property and operator
    createEditor<T>(property: DG.Property, operator: string, initialCondition?: SimpleCondition<T>): BaseConditionEditor<T> {
        // Try to find a specialized editor first (propertyType + semType + operator)
        const specializedKey = this.createEditorKey(property.propertyType, property.semType, operator);
        let EditorClass = this.editorRegistry.get(specializedKey);

        // If no specialized editor, try with empty semType (propertyType + operator)
        if (!EditorClass) {
            const typeOnlyKey = this.createEditorKey(property.propertyType, '', operator);
            EditorClass = this.editorRegistry.get(typeOnlyKey);
        }

        // If still no specialized editor, use default editor for property type (propertyType only)
        if (!EditorClass) {
            const defaultKey = this.createEditorKey(property.propertyType, '', '');
            EditorClass = this.editorRegistry.get(defaultKey);
        }

        // Fallback to base editor
        if (!EditorClass) {
            EditorClass = BaseConditionEditor;
        }

        return new EditorClass(property, operator, initialCondition);
    }

    // Get condition editor for a property and operator (with validation)
    getConditionEditor<T = any>(property: DG.Property, operator: string, condition?: SimpleCondition<T>): BaseConditionEditor<T> {
        if (this.isOperatorSupported(property, operator)) {
            return this.createEditor(property, operator, condition);
        } else {
            throw Error(`Operator '${operator}' is not supported for property type '${property.propertyType}'`);
        }
    }

    private createEditorKey(propertyType: string, semType: string, operator: string): string {
        return `${propertyType}:${semType}:${operator}`;
    }
}


export class QueryBuilder {
    root: HTMLElement;
    condition: ComplexCondition;
    properties: DG.Property[];
    structureChanged: Subject<ComplexCondition> = new Subject<ComplexCondition>();
    filterValueChanged: Subject<SimpleCondition> = new Subject<SimpleCondition>();

    constructor(properties: DG.Property[], initialCondition?: ComplexCondition) {
        // Initialize with provided condition or create a default logical condition
        this.condition = initialCondition || {
            logicalOperator: Operators.Logical.and,
            conditions: []
        };
        this.properties = properties;

        this.root = this.buildUI(this.condition, undefined, 0);
    }

    rebuildUI(): void {
        const newRoot = this.buildUI(this.condition, undefined, 0);
        this.root.replaceWith(newRoot);
        this.root = newRoot;
    }

    buildUI(
        condition: ComplexCondition,
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
                const property = this.properties.find((prop: DG.Property) => prop.name === fieldChoiceInput.value!);
                if (property) {
                    const registry = ConditionRegistry.getInstance();
                    const operators = registry.getOperatorsForProperty(property);
                    if (!cond.operator || cond.operator === '') {
                        cond.operator = operators[0];
                        cond.value = undefined;
                        this.filterValueChanged.next(cond);
                    }
                    const operatorsInput = ui.input.choice('', {
                        value: cond.operator,
                        items: operators,
                        onValueChanged: () => {
                            cond.operator = operatorsInput.value!;
                            this.rebuildUI();
                            this.structureChanged.next(this.condition);
                        },
                        nullable: false
                    });
                    operatorInputDiv.append(operatorsInput.root);
                    createEditor(property, cond.operator, cond);
                }
            }

            const createEditor = (property: DG.Property, operator: string, condition: SimpleCondition) => {
                ui.empty(criteriaDiv);
                const registry = ConditionRegistry.getInstance();
                const editor = registry.getConditionEditor(property, operator, condition);
                editor.onChanged.subscribe(() => {
                    this.filterValueChanged.next(editor.condition);
                })
                criteriaDiv.append(editor.root);
            }

            cond.field ??= (this.properties[0]?.name || '')

            const fieldChoiceInput = ui.input.choice('', {
                items: this.properties.map((prop: DG.Property) => prop.name),
                value: cond.field,
                nullable: false,
                onValueChanged: () => {
                    cond.field = fieldChoiceInput.value!;
                    cond.value = undefined;
                    cond.operator = '';
                    createFilter();
                    this.filterValueChanged.next(cond);
                }
            });

            const deleteFieldIcon = ui.icons.delete(() => {
                container.removeChild(filterContainer);
                removeOrReplaceCondition(cond, parentCondition ? parentCondition.conditions : condition.conditions);
                this.rebuildUI();
                this.structureChanged.next(this.condition);
            });

            const addNestedConditionIcon = ui.icons.add(() => {
                const parentCond = { logicalOperator: Operators.Logical.and, conditions: [cond] };
                removeOrReplaceCondition(cond, parentCondition ? parentCondition.conditions : condition.conditions, parentCond);
                this.rebuildUI();
                this.structureChanged.next(this.condition);

            }, 'Add nested condition');

            const operatorInputDiv = ui.div('', 'query-builder-filter-operator');
            const criteriaDiv = ui.div();
            const filterContainer = ui.divH([fieldChoiceInput.root, operatorInputDiv, criteriaDiv,
            ui.divH([deleteFieldIcon, addNestedConditionIcon], 'add-delete-icons')], 'revvity-signals-filter-inputs');
            container.appendChild(filterContainer);
            createFilter();
            this.filterValueChanged.next(cond);
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
            this.rebuildUI();
            this.structureChanged.next(this.condition);
        }, 'Add field to the condition');

        // AND/OR operators handling
        condition.logicalOperator ??= Operators.Logical.and;
        const logicalOperatorIcon = ui.button(condition.logicalOperator, () => {
            condition.logicalOperator = condition.logicalOperator === Operators.Logical.or ? Operators.Logical.and : Operators.Logical.or;
            logicalOperatorIcon.innerText = condition.logicalOperator.toUpperCase();
            this.rebuildUI();
            this.structureChanged.next(this.condition);
        }, 'Logical operator');

        logicalOperatorIcon.classList.add('property-query-builder-and-or-operator');

        const addDeleteIconsDiv = ui.divH([addFieldIcon], 'nested-add-delete-icons');
        const iconsDiv = ui.divH([logicalOperatorIcon, addDeleteIconsDiv], 'query-builder-nested-condition-operators');
        const container = ui.divV([iconsDiv], 'revvity-signals-search-query-operator');


        if (nestingLevel > 0) {
            container.setAttribute('data-nesting-level', nestingLevel.toString());
            const color = this.generateNestingColor(nestingLevel);
            container.style.setProperty('--nesting-color', color);
        }

        //delete icon for nested conditions to remove the whole nested condition at once
        if (parentCondition) {
            const deleteCondition = ui.icons.delete(() => {
                removeOrReplaceCondition(condition, parentCondition.conditions);
                this.rebuildUI();
                this.structureChanged.next(this.condition);
            });
            addDeleteIconsDiv.append(deleteCondition);
        }

        for (const nestedCondition of condition.conditions!) {
            if ('field' in nestedCondition)
                createPropertyFilterContainer(nestedCondition);
            else {
                const nestedContainer = this.buildUI(nestedCondition as ComplexCondition,
                    condition, nestingLevel + 1);
                container.append(nestedContainer);
            }
        }

        condition.conditions.length < 2 ? container.classList.remove('property-query-builder-multipe-filelds') :
            container.classList.add('property-query-builder-multipe-filelds');

        return container;
    }

    generateNestingColor(level: number): string {
        const baseColor = { r: 187, g: 222, b: 251 };

        const darkeningFactor = Math.min(0.1 + (level - 1) * 0.08, 0.9);

        const r = Math.round(baseColor.r * (1 - darkeningFactor));
        const g = Math.round(baseColor.g * (1 - darkeningFactor));
        const b = Math.round(baseColor.b * (1 - darkeningFactor));

        // Convert to hex
        return `#${r.toString(16).padStart(2, '0')}${g.toString(16).padStart(2, '0')}${b.toString(16).padStart(2, '0')}`;
    }
}


function suggestionMenuKeyNavigation(inputContainer: HTMLElement) {
  inputContainer.addEventListener('keydown', (e) => {
    if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp' && e.key !== 'Enter')
      return;
    const currentlySelected: HTMLElement | null = inputContainer.querySelector('.d4-menu-item-hover');
    const allItems: HTMLElement[] = Array.from(inputContainer.querySelectorAll('.d4-menu-item') ?? []);
    if (!allItems || allItems.length === 0)
      return;
    allItems.sort((a, b) => a.offsetTop - b.offsetTop); // sort by vertical position

    let currentIndex = currentlySelected ? allItems.indexOf(currentlySelected) : -1;
    if (e.key === 'ArrowDown')
      currentIndex = (currentIndex + 1) % allItems.length;
    else if (e.key === 'ArrowUp')
      currentIndex = (currentIndex - 1 + allItems.length) % allItems.length;
    else if (e.key === 'Enter' && currentlySelected) {
      currentlySelected.click();
      e.preventDefault();
      return;
    }
    currentlySelected?.classList.remove('d4-menu-item-hover');
    allItems[currentIndex].classList.add('d4-menu-item-hover');
  });
}

// Main UI builder function for property filters
export function buildPropertyBasedQueryBuilder(properties: DG.Property[], initialCondition?: ComplexCondition): HTMLElement {

    const builder = new QueryBuilder(properties, initialCondition);

    // Preview condition (will be further removed)
    const jsonPreview = document.createElement('textarea');
    jsonPreview.readOnly = true;
    jsonPreview.classList.add('revvity-signals-search-query-preview');

    const updatePreview = () => jsonPreview.value = JSON.stringify({ filter: builder.condition }, null, 2);

    builder.structureChanged.subscribe(() => updatePreview());
    builder.filterValueChanged.subscribe(() => updatePreview());

    const submitBtn = ui.button('Apply Filters', async () => {
        grok.shell.info('Filters applied: ' + JSON.stringify({ filter: builder.condition }));
    });

    const formContainer = ui.divV([
        ui.h3('Property Filters'),
        builder.root,
        submitBtn,
    ], 'revvity-signals-search');

    const previewContainer = ui.divV([
        ui.h3('Filter Preview'),
        jsonPreview
    ]);

    const mainContainer = ui.splitV([
        formContainer,
        previewContainer
    ], { style: { width: '100%' } }, true);

    updatePreview();
    return mainContainer;
}
