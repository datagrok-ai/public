import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Subject } from 'rxjs';
import '../../css/query-builder.css';
import dayjs from 'dayjs';

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

export enum QueryBuilderLayout {
    Standard = 'standard',
    Narrow = 'narrow'
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
    onValidationError: Subject<boolean> = new Subject<boolean>();
    showSuggestions = false;
    suggestionsMenuClicked = false;

    constructor(prop: DG.Property, operator: string, initialCondition?: SimpleCondition<T>) {
        this.condition = initialCondition ?? {
            field: prop.name,
            operator: operator,
            value: undefined as T
        };
        const inputs = this.initializeEditor(prop).then((inputs) => {
            inputs.forEach((input) => {
                input.onChanged.subscribe(() => {
                    this.onValidationError.next(!input.validate());
                });
            });
            if(inputs.some((it) => !it.validate()))
                this.onValidationError.next(true);;
        });
    }

    protected async initializeEditor(prop: DG.Property): Promise<DG.InputBase[]> {
        if (Array.isArray(this.condition.value))
            this.condition.value = undefined as T;
        const initVal = prop.type === DG.TYPE.DATE_TIME && typeof this.condition.value === 'string' ?
            dayjs(this.condition.value as string) : this.condition.value;
        const input = ui.input.forProperty(prop, undefined, {
            value: initVal,
            nullable: false,
            onValueChanged: () => {
                this.condition.value = input.value as T;
                this.onChanged.next(this.condition);
                this.addSuggestions(prop, input);
            }
        });
        input.addCaption('');
        this.root.append(input.root);
        this.initializeSuggestions(prop, input);
        return [input];
    }

    initializeSuggestions(prop: DG.Property, input: DG.InputBase) {
        this.showSuggestions = prop.type === DG.TYPE.STRING && prop.options[SUGGESTIONS_FUNCTION];
        if (this.showSuggestions) {
           suggestionMenuKeyNavigation(this.root);
        }
    }

    async addSuggestions(prop: DG.Property, input: DG.InputBase) {
        if (this.showSuggestions && !this.suggestionsMenuClicked) {
            if (input.value) {
            const suggestions: Array<string> = await prop.options[SUGGESTIONS_FUNCTION](input.value);
            if (suggestions.length === 1 && suggestions[0] === input.value)
                return;
            const suggestionsMenu = DG.Menu.popup();
            suggestions.forEach((s) => {
                suggestionsMenu!.item(s, () => {
                    input.value = s;
                    this.suggestionsMenuClicked = true;
                });
            });
            suggestionsMenu.show({ element: input.root, y: input.root.offsetHeight });
        }
        }
        this.suggestionsMenuClicked = false;
    }
}

export class BetweenConditionEditor extends BaseConditionEditor {
    override async initializeEditor(prop: DG.Property): Promise<DG.InputBase[]> {
        if (!this.condition.value || !Array.isArray(this.condition.value)) {
            this.condition.value = [this.condition.value, undefined];
            this.onChanged.next(this.condition);
        }
        const input1 = ui.input.forProperty(prop, undefined, {
            value: this.condition.value[0],
            nullable: false,
            onValueChanged: () => {
                this.condition.value[0] = input1.value!;
                this.onChanged.next(this.condition);
            }
        });
        input1.addCaption('');
        const input2 = ui.input.forProperty(prop, undefined, {
            value: this.condition.value[1],
            onValueChanged: () => {
                this.condition.value[1] = input2.value!;
                this.onChanged.next(this.condition);
            }
        });
        input2.addCaption('');
        this.root.append(ui.divH([input1.root, input2.root], { style: { gap: '10px' } }));
        return [input1, input2];
    }
}

export class BooleanConditionEditor extends BaseConditionEditor {
    override async initializeEditor(prop: DG.Property): Promise<DG.InputBase[]> {
        if (this.condition.value === undefined)
            this.condition.value = true;
        const options = ['true', 'false'];
        const input = ui.input.choice('', {
            items: options,
            nullable: false,
            value: this.condition.value ? 'true' : 'false',
            onValueChanged: () => {
                this.condition.value = input.value === 'true';
                this.onChanged.next(this.condition);
            }
        });
        this.root.append(input.root);
        return [input];
    }
}

export class MoleculeConditionEditor extends BaseConditionEditor<string> {
    override async initializeEditor(prop: DG.Property): Promise<DG.InputBase[]> {
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
        return [input];
    }
}

export type MoleculeSimilarity = {
    molecule: string,
    threshold: number
}
export class MoleculeSimilarityConditionEditor extends BaseConditionEditor<MoleculeSimilarity> {
    override async initializeEditor(prop: DG.Property): Promise<DG.InputBase[]> {
        //in case we switch from some other operator, where value is a string (contains, is contained)
        if (!this.condition.value || typeof this.condition.value === 'string') {
            this.condition.value = {molecule: this.condition.value ?? '', threshold: 0.7};
        }
        const moleculeInput = ui.input.molecule('', {
            value: this.condition.value.molecule,
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
        return [moleculeInput, tresholdInput];
    }
}


export class MultiValueConditionEditorString extends BaseConditionEditor<string[]> {
    override async initializeEditor(prop: DG.Property): Promise<DG.InputBase[]> {
        if (!Array.isArray(this.condition.value)) {
            this.condition.value = [];
        }
        const input = ui.input.list('', {
            nullable: false,
            value: this.condition.value,
            onValueChanged: () => {
                this.condition.value = input.value!;
                this.onChanged.next(this.condition);
            }
        });
        this.root.append(input.root);
        return [input];
    }
}

export class MultiValueConditionEditorInt extends BaseConditionEditor<number[]> {
    override async initializeEditor(prop: DG.Property): Promise<DG.InputBase[]> {
        if (!Array.isArray(this.condition.value)) {
            this.condition.value = [];
        }
        const input = ui.input.list('', {
            nullable: false,
            value: this.condition.value,
            onValueChanged: () => {
                this.condition.value = input.value!.map((it) => parseInt(it));
                this.onChanged.next(this.condition);
            }
        });
        this.root.append(input.root);
        return [input];
    }
}

export class MultiValueConditionEditorFloat extends BaseConditionEditor<number[]> {
    override async initializeEditor(prop: DG.Property): Promise<DG.InputBase[]> {
        if (!Array.isArray(this.condition.value)) {
            this.condition.value = [];
        }
        const input = ui.input.list('', {
            nullable: false,
            value: this.condition.value,
            onValueChanged: () => {
                this.condition.value = input.value!.map((it) => parseFloat(it));
                this.onChanged.next(this.condition);
            }
        });
        this.root.append(input.root);
        return [input];
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
        this.registerTypeOperators(DG.TYPE.BOOL, [Operators.EQ]);
        this.registerTypeOperators(DG.TYPE.STRING, [Operators.STARTS_WITH, Operators.CONTAINS, Operators.EQ, Operators.NOT_EQ, Operators.IN]);
        this.registerTypeOperators(DG.TYPE.INT, [Operators.GT, Operators.LT, Operators.GTE, Operators.LTE, Operators.EQ,
            Operators.NOT_EQ, Operators.BETWEEN, Operators.IN]);
        this.registerTypeOperators(DG.TYPE.FLOAT, [Operators.GT, Operators.LT, Operators.GTE, Operators.LTE, Operators.EQ,
            Operators.NOT_EQ, Operators.BETWEEN, Operators.IN]);
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
        this.registerEditor(DG.TYPE.BOOL, '', '', BaseConditionEditor);

        // Register specialized editors for specific combinations
        this.registerEditor(DG.TYPE.STRING, '', Operators.IN, MultiValueConditionEditorString);
        this.registerEditor(DG.TYPE.INT, '', Operators.IN, MultiValueConditionEditorInt);
        this.registerEditor(DG.TYPE.FLOAT, '', Operators.IN, MultiValueConditionEditorFloat);
        this.registerEditor(DG.TYPE.DATE_TIME, '', Operators.BETWEEN, BetweenConditionEditor);
        this.registerEditor(DG.TYPE.INT, '', Operators.BETWEEN, BetweenConditionEditor);
        this.registerEditor(DG.TYPE.FLOAT, '', Operators.BETWEEN, BetweenConditionEditor);
        this.registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, Operators.CONTAINS, MoleculeConditionEditor);
        this.registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, Operators.IS_CONTAINED, MoleculeConditionEditor);
        this.registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, Operators.EQ, MoleculeConditionEditor);
        this.registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, Operators.IS_SIMILAR, MoleculeSimilarityConditionEditor);
        this.registerEditor(DG.TYPE.BOOL, '', Operators.EQ, BooleanConditionEditor);
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


export type HistoryItem = {
    date: DG.TYPE.DATE_TIME,
    value: string,
}

export function saveHistory(key: string, value: string) {
    let collection: HistoryItem[] = JSON.parse(localStorage.getItem(key) ?? '[]');
    const newItem = { date: new Date(), value: value };
    if (!collection.length)
        localStorage.setItem(key, JSON.stringify([newItem]));
    else {
        collection = collection.filter((it) => it.value !== value)
        localStorage.setItem(key, JSON.stringify([newItem, ...collection.slice(0, 9)]));
    }
}

export function getHistory(key: string): HistoryItem[] {
    return JSON.parse(localStorage.getItem(key) ?? '[]');
}

export class QueryBuilder {
    root: HTMLElement;
    condition: ComplexCondition;
    properties: DG.Property[];
    structureChanged: Subject<ComplexCondition> = new Subject<ComplexCondition>();
    filterValueChanged: Subject<SimpleCondition> = new Subject<SimpleCondition>();
    validationError: Subject<boolean> = new Subject<boolean>();
    invalid = false;
    layout: QueryBuilderLayout;
    private _historyCacheKey: string;
    private historyIcon: HTMLElement | null = null;
    private historyIconOnTop = false;

    constructor(properties: DG.Property[], initialCondition?: ComplexCondition, layout: QueryBuilderLayout = QueryBuilderLayout.Standard,
        historyCacheKey: string = '', historyIconOnTop = false) {
        this.properties = properties;
        this.layout = layout;
        this._historyCacheKey = historyCacheKey;
        this.historyIconOnTop = historyIconOnTop;
        
        if (initialCondition) {
            this.condition = initialCondition;
        } else {
            const defaultCondition = this.createDefaultSimpleCondition();
            this.condition = {
                    logicalOperator: Operators.Logical.and,
                    conditions: [defaultCondition]
            };
        }

        this.root = this.buildUI(this.condition, undefined, 0);
        this.setupHistoryIcon();
    }

    set historyCache(key: string) {
        this._historyCacheKey = key;
        this.setupHistoryIcon(); // This will control visibility
    }

    saveConditionToHistory(): void {
        if (this._historyCacheKey)
            saveHistory(this._historyCacheKey, JSON.stringify(this.condition));
    }

    private setupHistoryIcon(): void {
        if (!this.historyIcon) {
            this.historyIcon = ui.iconFA('history', () => {
                this.showHistoryMenu();
            }, 'Query History');

            this.historyIcon.classList.add('query-builder-history-icon');
        }

        if (this.historyIconOnTop)
            this.root.prepend(this.historyIcon);
        else
            this.root.appendChild(this.historyIcon);
        
        // Set data attribute for CSS to control visibility
        if (this.root) {
            this.root.setAttribute('data-history', this._historyCacheKey  === '' ? 'false' : 'true');
        }
    }

    private showHistoryMenu(): void {
        const items: HistoryItem[] = getHistory(this._historyCacheKey);
        const menu = DG.Menu.popup();
        menu.items(items.map((it) => ui.tools.click(
            ui.divH([
                ui.divText(it.date.toString(), {style: {color: '#7990A5'}}),
                ui.divText(this.conditionToString(JSON.parse(it.value) as ComplexCondition)),
            ], {style: {gap: '5px'}}),
        () => {this.loadCondition(JSON.parse(it.value))})), () => {});
        menu.show();
    }

    loadCondition(condition: ComplexCondition): void {
        this.condition = JSON.parse(JSON.stringify(condition));
        this.rebuildUI();
        this.structureChanged.next(this.condition);
    }

    conditionToString(condition: ComplexCondition): string {
        if (!condition || !condition.conditions || condition.conditions.length === 0) {
            return 'No conditions';
        }

        // If there's only one condition, don't wrap it in logical operator
        if (condition.conditions.length === 1) {
            return this.expressionToString(condition.conditions[0]);
        }

        // For multiple conditions, wrap them in logical operator
        const logicalOp = condition.logicalOperator === Operators.Logical.and ? 'AND' : 'OR';
        const conditionStrings = condition.conditions.map(cond => this.expressionToString(cond));
        
        return `(${conditionStrings.join(` ${logicalOp} `)})`;
    }

    private expressionToString(expression: Expression): string {
        if ('field' in expression) {
            // This is a SimpleCondition
            return this.simpleConditionToString(expression);
        } else {
            // This is a ComplexCondition (recursive)
            return this.conditionToString(expression);
        }
    }

    private simpleConditionToString(condition: SimpleCondition): string {
        const field = condition.field || 'Unknown field';
        const operator = condition.operator || 'Unknown operator';
        const value = condition.value !== undefined && condition.value !== null ? 
            (Array.isArray(condition.value) ? `[${condition.value.join(', ')}]` : String(condition.value)) : 
            'No value';
        
        return `${field} ${operator} ${value}`;
    }

    createDefaultSimpleCondition(): SimpleCondition {
        const firstProperty = this.properties[0];
        if (firstProperty) {
            const registry = ConditionRegistry.getInstance();
            const operators = registry.getOperatorsForProperty(firstProperty);
            const firstOperator = operators[0];

            return {
                field: firstProperty.name,
                operator: firstOperator,
                value: undefined
            }
        } else
            return {
                field: '',
                value: undefined,
                operator: ''
            }
    }

    rebuildUI(): void {
        const newRoot = this.buildUI(this.condition, undefined, 0);
        this.root.replaceWith(newRoot);
        this.root = newRoot;
        
        // Recreate history icon after rebuilding UI
        this.setupHistoryIcon();
    }

    setLayout(layout: QueryBuilderLayout): void {
        this.layout = layout;

        if (layout === QueryBuilderLayout.Narrow)
            this.root.classList.add('query-builder-narrow');
        else
            this.root.classList.remove('query-builder-narrow');
        
        // Find all filter containers and apply/remove layout classes
        const filterContainers = this.root.querySelectorAll('.query-builder-filter-inputs, .query-builder-filter-inputs-narrow');
        
        filterContainers.forEach(container => {
            if (layout === QueryBuilderLayout.Narrow) {
                container.classList.remove('query-builder-filter-inputs');
                container.classList.add('query-builder-filter-inputs-narrow');
            } else {
                container.classList.remove('query-builder-filter-inputs-narrow');
                container.classList.add('query-builder-filter-inputs');
            }
        });
    }

    getLayout(): QueryBuilderLayout {
        return this.layout;
    }

    buildUI(
        condition: ComplexCondition,
        parentCondition?: ComplexCondition,
        nestingLevel: number = 0
    ): HTMLElement {

        //reset validation status
        this.invalid = false;
        this.validationError.next(false);

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
                const property = getPropByFriendlyName(fieldChoiceInput.value!);
                if (property) {
                    //in case of boolean input we need to hide operators input
                    if (property.type === DG.TYPE.BOOL)
                        filterContainer.classList.add('boolean-input');
                    else
                        filterContainer.classList.remove('boolean-input');
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
                });
                editor.onValidationError.subscribe((error) => {
                    this.invalid = error;
                    this.validationError.next(error);
                });
                criteriaDiv.append(editor.root);
            }

            const getPropByName = (name: string) => {
                const property = this.properties.find((prop: DG.Property) => prop.name === name);
                return property;
            }

            const getPropByFriendlyName = (friendlyName: string) => {
                let property = this.properties.find((prop: DG.Property) => prop.friendlyName === friendlyName);
                if (!property)
                    property = this.properties.find((prop: DG.Property) => prop.name === friendlyName);
                return property;
            }

            const getFieldChoiceInputVal = (name: string) => {
                const prop = getPropByName(name);
                return prop?.friendlyName ?? prop?.name ?? '';
            }

            const getPropNameByChoiceInputVal = (friendlyName: string) => {
                const prop = getPropByFriendlyName(friendlyName);
                return prop?.name ?? '';
            }

            cond.field ??= this.properties[0]?.name || '';

            const fieldChoiceInput = ui.input.choice('', {
                items: this.properties.map((prop: DG.Property) => prop.friendlyName ?? prop.name),
                value: getFieldChoiceInputVal(cond.field),
                nullable: false,
                onValueChanged: () => {
                    cond.field = getPropNameByChoiceInputVal(fieldChoiceInput.value!);
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
            }, 'Remove current condition');

            const addNestedConditionIcon = ui.icons.add(() => {
                const parentCond = { logicalOperator: Operators.Logical.and, conditions: [cond] };
                removeOrReplaceCondition(cond, parentCondition ? parentCondition.conditions : condition.conditions, parentCond);
                this.rebuildUI();
                this.structureChanged.next(this.condition);

            }, 'Add nested condition');

            const operatorInputDiv = ui.div('', 'query-builder-filter-operator');
            const criteriaDiv = ui.div();

            const filterContainer = ui.div([
                ui.div([fieldChoiceInput.root], 'query-builder-filter-field'),
                ui.div([operatorInputDiv], 'query-builder-filter-operator'),
                ui.div([criteriaDiv], 'query-builder-filter-value'),
                ui.div([addNestedConditionIcon, deleteFieldIcon], 'add-delete-icons')
            ], this.layout === QueryBuilderLayout.Narrow ? 'query-builder-filter-inputs-narrow' : 'query-builder-filter-inputs');
            
            container.appendChild(filterContainer);
            createFilter();
            this.filterValueChanged.next(cond);
        }

        // adding filter field
        const addFieldIcon = ui.icons.add(() => {
            const conditionForFilter: SimpleCondition = this.createDefaultSimpleCondition();
            condition.conditions?.push(conditionForFilter);
            condition.conditions.length < 2 ? container.classList.remove('property-query-builder-multipe-filelds') :
                container.classList.add('property-query-builder-multipe-filelds');
            this.rebuildUI();
            this.structureChanged.next(this.condition);
        }, 'Add field to the condition');

        // AND/OR operators handling
        condition.logicalOperator ??= Operators.Logical.and;

        const logicalOperatorChoice = ui.input.choice('', {
            value: condition.logicalOperator === Operators.Logical.and ? 'all' : 'any',
            items: ['all', 'any'],
            nullable: false,
            onValueChanged: () => {
                condition.logicalOperator = logicalOperatorChoice.value! === 'all' ? Operators.Logical.and : Operators.Logical.or;
                this.rebuildUI();
                this.structureChanged.next(this.condition);
            }
        })

        const logicalOperatorDiv = ui.divH([ui.divText('Match'), logicalOperatorChoice.root], 'query-builder-match-div');

        const addDeleteIconsDiv = ui.divH([addFieldIcon], 'nested-add-delete-icons');
        const iconsDiv = ui.divH([logicalOperatorDiv, addDeleteIconsDiv], 'query-builder-nested-condition-operators');
        const container = ui.divV([iconsDiv], 'query-builder-search-query-operator');


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

