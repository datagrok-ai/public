import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getUsers } from './package';
import { queryUsers } from './revvityApi';
import { getRevvityUsers } from './users';

export enum RevvityLogicalOperators {
    and = 'and',
    or = 'or',
}

export enum RevvityDateOperators {
    before = 'before',
    after = 'after',
    between = 'between',
}

export enum RevvityChemSearches {
    substructure = 'substructure',
    similar = 'similar',
    full = 'full',
    exact = 'exact',
    fullIncludeTautomers = 'full include tautomers'
}

export interface RevvityDateFilter {
    dates: string[];
    operator: RevvityDateOperators;
}

export interface RevvityChemFilter {
    structure: string;
    searchType: string;
    operator: RevvityChemSearches;
    logicalOperators: RevvityLogicalOperators;
}

export enum DefaultDateFilters {
    Created = 'Created',
    Edited = 'Edited',
}

export enum DefaultUserFilters {
    CreatedBy = 'Created by',
    EditedBy = 'Edited by',
}

export enum DefaultChemFilters {
    Chemical = 'Chemical'
}

export class RevvityFilters {
    dateFilters: {[key: string]: RevvityDateFilter[]} = {};
    userFilters: {[key: string]: DG.InputBase} = {};
    chemFilters: {[key: string]: RevvityChemFilter[]} = {};
    view: DG.ViewBase;
    usedFilters: Set<string> = new Set();
    filtersButton: HTMLButtonElement;
    filtersDiv = ui.div();
    isDocked = false;
    
    constructor(v: DG.View) {
        this.view = v;
        this.filtersButton = ui.button('Add filters', () => {
            this.showFiltersPopup();
        });
        this.view.setRibbonPanels( [[this.filtersButton]] );
    }

    private showFiltersPopup(): void {
        // Get all available filters including Chemical
        const allFilters = [
            ...Object.values(DefaultDateFilters),
            ...Object.values(DefaultUserFilters),
            ...Object.values(DefaultChemFilters)
        ];
        
        // Filter out already used filters
        const availableFilters = allFilters.filter(filter => !this.usedFilters.has(filter));
        
        if (availableFilters.length === 0) {
            grok.shell.info('No more filters available to add');
            return;
        }

        // Create tree view with filter categories
        const treeView = ui.tree();
        
        // Group filters by category
        const filterGroups = this.groupFiltersByCategory(availableFilters);

        // Add filter groups to tree
        Object.entries(filterGroups).forEach(([category, filters]) => {
            if (filters.length > 0) {
                const group = treeView.group(category);
                filters.forEach(filter => {
                    const item = group.item(filter);
                    item.onSelected.subscribe(() => {
                        this.selectFilter(filter, popup as HTMLElement);
                    });
                });
            }
        });
        
        const popupContent = ui.divV([
            ui.label('Select a filter to add:'),
            treeView
        ]);

        // Create and show popup
        const popup = ui.showPopup(popupContent, this.filtersButton, {
            vertical: false,
            dx: 0,
            dy: 0,
            smart: true
        });
    }

    private groupFiltersByCategory(filters: (DefaultDateFilters | DefaultUserFilters | DefaultChemFilters)[]): Record<string, (DefaultDateFilters | DefaultUserFilters | DefaultChemFilters)[]> {
        return {
            'Date': filters.filter(filter => Object.values(DefaultDateFilters).includes(filter as DefaultDateFilters)),
            'User': filters.filter(filter => Object.values(DefaultUserFilters).includes(filter as DefaultUserFilters)),
            'Chemical': filters.filter(filter => Object.values(DefaultChemFilters).includes(filter as DefaultChemFilters))
        };
    }

    private selectFilter(filterName: string, popup: HTMLElement): void {
        // Add to used filters
        this.usedFilters.add(filterName);
        
        // Show info message
        grok.shell.info(`Selected filter: ${filterName}`);
        
        // Close popup
        popup.remove();
        
        // Add filter to the filters div
        this.addFilterToDiv(filterName);
        
        // Dock the filters div if not already docked
        this.dockFiltersDiv();
    }

    private dockFiltersDiv(): void {
        if (!this.isDocked) {
            grok.shell.dockManager.dock(this.filtersDiv, 'left', null, 'Filters', 0.3);
            this.isDocked = true;
        }
    }

    private addFilterToDiv(filterName: string): void {
        // Clear the div if it's empty (first filter)
        if (this.filtersDiv.children.length === 0) {
            ui.empty(this.filtersDiv);
            this.filtersDiv.appendChild(ui.h3('Filters'));
        }

        // Create filter container
        const filterContainer = this.createFilterContainer(filterName);

        // Add close button for the filter
        const closeButton = this.createCloseButton(() => {
            this.removeFilter(filterName, filterContainer);
        });
        filterContainer.appendChild(closeButton);

        this.filtersDiv.appendChild(filterContainer);
    }

    private createFilterContainer(filterName: string): HTMLElement {
        const filterContainer = ui.divV([], 'revvity-signals-filter-row');
        
        // Add filter title
        const filterTitle = ui.divText(filterName, {style: {fontWeight: 'bold', marginBottom: '8px'}});
        filterContainer.appendChild(filterTitle);

        // Add appropriate filter controls based on filter type
        if (Object.values(DefaultDateFilters).includes(filterName as DefaultDateFilters)) {
            this.addDateFilterControls(filterContainer, filterName);
        } else if (Object.values(DefaultUserFilters).includes(filterName as DefaultUserFilters)) {
            this.addUserFilterControls(filterContainer, filterName);
        } else if (Object.values(DefaultChemFilters).includes(filterName as DefaultChemFilters)) {
            this.addUserChemFilterControls(filterContainer);
        }

        return filterContainer;
    }

    private createCloseButton(onClick: () => void): HTMLElement {
        const closeButton = ui.icons.close(onClick);
        closeButton.classList.add('revvity-signals-filter-close-button');
        return closeButton;
    }

    private addDateFilterControls(container: HTMLElement, filterName: string): void {
        this.createDateFilterRow(container, filterName, false);
    }

    private createDateFilterRow(container: HTMLElement, filterName: string, isAdditionalRow: boolean = false): void {
        const filterRow = this.createFilterRow(container);
        
        const { choiceInput, dateInput1, dateInput2, dateContainer } = this.createDateInputs();
        
        // Add inputs to row
        filterRow.appendChild(choiceInput.root);
        filterRow.appendChild(dateContainer);

        // Add 'or' button
        const orButton = this.createOrButton(() => {
            this.createDateFilterRow(container, filterName, true);
        });
        filterRow.appendChild(orButton);

        // Add close button for additional filter rows only
        if (isAdditionalRow) {
            const closeButton = this.createCloseButton(() => {
                // Find the index of this specific filter row in the container
                const filterRows = Array.from(container.children).filter(child => 
                    child.classList.contains('revvity-signals-filter-row')
                );
                if (filterRows.length >= 2) {
                    filterRows[filterRows.length - 2].classList.add('revvity-filter-last-row');
                }
                const rowIndex = filterRows.indexOf(filterRow);
                
                // Remove the filter at the correct index
                this.removeDateFilterByIndex(filterName, rowIndex);
                filterRow.remove();
            });
            filterRow.appendChild(closeButton);
        }

        // Save to appropriate date filter array
        this.saveDateFilterToArray(filterName, choiceInput, dateInput1, dateInput2);

        container.appendChild(filterRow);
    }

    private removeDateFilterByIndex(filterName: string, index: number): void {
        if (this.dateFilters[filterName] && index >= 0 && index < this.dateFilters[filterName].length) {
            // Remove the filter at the specific index
            this.dateFilters[filterName].splice(index, 1);
            
            // Remove key if array is empty
            if (this.dateFilters[filterName].length === 0) {
                delete this.dateFilters[filterName];
            }
        }
    }

    private createFilterRow(container: HTMLElement): HTMLElement {
        const filterRows = Array.from(container.children).filter(child =>
            child.classList.contains('revvity-signals-filter-row')
        );
        if (filterRows.length)
            filterRows[filterRows.length - 1].classList.remove('revvity-filter-last-row');
        const filterRow = ui.divH([]);
        filterRow.classList.add('revvity-signals-filter-row', 'revvity-filter-last-row');
        return filterRow;
    }

    private createDateInputs(): { choiceInput: DG.InputBase, dateInput1: DG.InputBase, dateInput2: DG.InputBase, dateContainer: HTMLElement } {
        const choiceInput = ui.input.choice('', {
            items: Object.keys(RevvityDateOperators),
            value: RevvityDateOperators.before,
            nullable: false,
        });

        const dateInput1 = ui.input.date('');
        const dateInput2 = ui.input.date('');

        const dateContainer = ui.divH([]);
        dateContainer.classList.add('revvity-signals-filter-inputs');
        dateContainer.appendChild(dateInput1.root);

        const updateDateInputs = () => {
            ui.empty(dateContainer);
            dateContainer.appendChild(dateInput1.root);
            
            if (choiceInput.value === RevvityDateOperators.between) {
                dateContainer.appendChild(dateInput2.root);
            }
        };

        choiceInput.onChanged.subscribe(() => {
            updateDateInputs();
        });

        updateDateInputs();

        return { choiceInput, dateInput1, dateInput2, dateContainer };
    }

    private createOrButton(onClick: () => void): HTMLElement {
        const orButton = ui.button('or', onClick);
        orButton.classList.add('revvity-signals-filters-or-button');
        return orButton;
    }

    private createButtonsContainer(onClick: () => void, showAndButton: boolean = false): HTMLElement {
        const buttonsContainer = ui.divH([]);
        buttonsContainer.style.gap = '8px';
        
        // Add 'or' button
        const orButton = this.createOrButton(onClick);
        buttonsContainer.appendChild(orButton);
        
        // Add 'and' button if requested
        if (showAndButton) {
            const andButton = ui.button('and', onClick);
            andButton.classList.add('revvity-signals-filters-and-button');
            buttonsContainer.appendChild(andButton);
        }
        
        return buttonsContainer;
    }

    private saveDateFilterToArray(filterName: string, choiceInput: DG.InputBase, dateInput1: DG.InputBase, dateInput2: DG.InputBase): void {
        // Extract date values from inputs
        const dates: string[] = [];
        if (dateInput1.value) {
            dates.push(dateInput1.value.toString());
        }
        if (choiceInput.value === RevvityDateOperators.between && dateInput2.value) {
            dates.push(dateInput2.value.toString());
        }

        const dateFilter: RevvityDateFilter = {
            dates: dates,
            operator: choiceInput.value as RevvityDateOperators
        };

        this.dateFilters[filterName] = [...(this.dateFilters[filterName] || []), dateFilter];
    }



    private addUserFilterControls(container: HTMLElement, filterName: string): void {
        getRevvityUsers().then((res) => {
            if (res) {
                const users = Object.values(res).map((user) => `${user.firstName} ${user.lastName} (${user.userName})`);
                const userChoice = ui.input.multiChoice('', {
                    items: users
                });
                container.appendChild(userChoice.root);
                
                // Save user filter to dictionary with the correct filter name as key
                this.userFilters[filterName] = userChoice;
            }
        });
    }

    private addUserChemFilterControls(container: HTMLElement): void {
        this.createChemFilterRow(container, false);
    }

    private createChemFilterRow(container: HTMLElement, isAdditionalRow: boolean = false): void {
        const filterRow = this.createFilterRow(container);
        
        const { moleculeInput, searchTypeInput, similaritySlider, inputsContainer } = this.createChemInputs();
        
        // Add inputs to row
        filterRow.appendChild(inputsContainer);
        
        // Create buttons container
        const buttonsContainer = this.createButtonsContainer(() => {
            this.createChemFilterRow(container, true);
        }, true); // Show AND button for chemical filters
        filterRow.appendChild(buttonsContainer);
        
        // Add close button for additional filter rows only
        if (isAdditionalRow) {
            const closeButton = this.createCloseButton(() => {
                // Find the index of this specific filter row in the container
                const filterRows = Array.from(container.children).filter(child => 
                    child.classList.contains('revvity-signals-filter-row')
                );
                if (filterRows.length >= 2) {
                    filterRows[filterRows.length - 2].classList.add('revvity-filter-last-row');
                }
                const rowIndex = filterRows.indexOf(filterRow);
                
                // Remove the filter at the correct index
                this.removeChemFilterByIndex(rowIndex);
                filterRow.remove();
            });
            filterRow.appendChild(closeButton);
        }
        
        // Save to chemFilters array
        this.saveChemFilterToArray(moleculeInput, searchTypeInput, similaritySlider);
        
        container.appendChild(filterRow);
    }

    private removeChemFilterByIndex(index: number): void {
        if (this.chemFilters[DefaultChemFilters.Chemical] && index >= 0 && index < this.chemFilters[DefaultChemFilters.Chemical].length) {
            // Remove the filter at the specific index
            this.chemFilters[DefaultChemFilters.Chemical].splice(index, 1);
            
            // Remove key if array is empty
            if (this.chemFilters[DefaultChemFilters.Chemical].length === 0) {
                delete this.chemFilters[DefaultChemFilters.Chemical];
            }
        }
    }

    private createChemInputs(): { moleculeInput: DG.InputBase, searchTypeInput: DG.InputBase, similaritySlider: DG.InputBase, inputsContainer: HTMLElement } {
        // Create molecule input
        const moleculeInput = ui.input.molecule('');
        
        // Create search type choice input
        const searchTypeInput = ui.input.choice('', {
            items: Object.values(RevvityChemSearches),
            value: RevvityChemSearches.substructure,
            nullable: false,
        });
        
        // Create similarity slider (initially hidden)
        const similaritySlider = ui.input.float('', {
            value: 0.8,
            min: 0,
            max: 1,
            step: 0.01,
            nullable: false,
        });
        similaritySlider.root.style.display = 'none';
        
        // Show/hide similarity slider based on search type
        const updateSimilaritySlider = () => {
            if (searchTypeInput.value === RevvityChemSearches.similar) {
                similaritySlider.root.style.display = 'block';
            } else {
                similaritySlider.root.style.display = 'none';
            }
        };
        
        searchTypeInput.onChanged.subscribe(() => {
            updateSimilaritySlider();
        });
        
        updateSimilaritySlider();
        
        // Create inputs container
        const inputsContainer = ui.divV([]);
        inputsContainer.classList.add('revvity-signals-filter-inputs');
        inputsContainer.appendChild(moleculeInput.root);
        inputsContainer.appendChild(searchTypeInput.root);
        inputsContainer.appendChild(similaritySlider.root);

        return { moleculeInput, searchTypeInput, similaritySlider, inputsContainer };
    }

    private saveChemFilterToArray(moleculeInput: DG.InputBase, searchTypeInput: DG.InputBase, similaritySlider: DG.InputBase): void {
        const chemFilter: RevvityChemFilter = {
            structure: moleculeInput.value?.toString() || '',
            searchType: searchTypeInput.value?.toString() || '',
            operator: searchTypeInput.value as RevvityChemSearches,
            logicalOperators: RevvityLogicalOperators.or // Default to 'or', can be enhanced later
        };
        this.chemFilters[DefaultChemFilters.Chemical] = [...(this.chemFilters[DefaultChemFilters.Chemical] || []), chemFilter];
    }



    private removeFilter(filterName: string, filterContainer: HTMLElement): void {
        // Remove from used filters - this will make it available again in the popup
        this.usedFilters.delete(filterName);
        
        // Remove filters from appropriate arrays
        this.removeFiltersFromArrays(filterName, filterContainer);
        
        // Remove from UI
        filterContainer.remove();
        
        // If no more filters, undock the div
        if (this.filtersDiv.children.length <= 1) { // Only the title remains
            if (this.isDocked) {
                grok.shell.dockManager.close(this.filtersDiv);
                this.isDocked = false;
            }
        }
        
        // Show info message that filter is available again
        grok.shell.info(`Filter "${filterName}" is now available to add again`);
    }

    private removeFiltersFromArrays(filterName: string, filterContainer: HTMLElement): void {
        // Remove date filters from arrays if it's a date filter
        if (Object.values(DefaultDateFilters).includes(filterName as DefaultDateFilters)) {
            this.removeDateFiltersFromContainer(filterName, filterContainer);
            // Remove key if array is empty
            if (this.dateFilters[filterName] && this.dateFilters[filterName].length === 0) {
                delete this.dateFilters[filterName];
            }
        }
        
        // Remove user filters from dictionary if it's a user filter
        if (Object.values(DefaultUserFilters).includes(filterName as DefaultUserFilters)) {
            delete this.userFilters[filterName];
        }
        
        // Remove chemical filters from chemFilters array if it's a chemical filter
        if (Object.values(DefaultChemFilters).includes(filterName as DefaultChemFilters)) {
            this.removeChemFiltersFromContainer(filterContainer);
            // Remove key if array is empty
            if (this.chemFilters[DefaultChemFilters.Chemical] && this.chemFilters[DefaultChemFilters.Chemical].length === 0) {
                delete this.chemFilters[DefaultChemFilters.Chemical];
            }
        }
    }

    private removeDateFiltersFromContainer(filterName: string, filterContainer: HTMLElement): void {
        if (this.dateFilters[filterName]) {
            // Remove all filters for this container
            delete this.dateFilters[filterName];
        }
    }

    private removeChemFiltersFromContainer(filterContainer: HTMLElement): void {
        if (this.chemFilters[DefaultChemFilters.Chemical]) {
            // Remove all filters for this container
            delete this.chemFilters[DefaultChemFilters.Chemical];
        }
    }
}
