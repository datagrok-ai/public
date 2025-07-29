import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getUsers } from './package';
import { queryUsers } from './revvityApi';
import { getRevvityUsers } from './users';

export enum RevvityDateOperators {
    before = 'before',
    after = 'after',
    between = 'between',
}

export interface RevvityDateFilter {
    filters: DG.InputBase[];
    operator: RevvityDateOperators;
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
    createdAt: RevvityDateFilter[] = [];
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
        // Get all available filters (excluding Chemical for now)
        const allFilters = [
            ...Object.values(DefaultDateFilters),
            ...Object.values(DefaultUserFilters)
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
        const dateFilters = availableFilters.filter(filter => 
            Object.values(DefaultDateFilters).includes(filter as DefaultDateFilters)
        );
        const userFilters = availableFilters.filter(filter => 
            Object.values(DefaultUserFilters).includes(filter as DefaultUserFilters)
        );

        const addFiltersGroup = (name: string, filters: (DefaultDateFilters | DefaultUserFilters)[]) => {
        if (filters.length > 0) {
            const group = treeView.group(name);
            filters.forEach(filter => {
                const item = group.item(filter);
                item.onSelected.subscribe(() => {
                    this.selectFilter(filter, popup as HTMLElement);
                });
            });
        }
        }

        // Add date filters group
        addFiltersGroup('Date', dateFilters);

        // Add user filters group
        addFiltersGroup('User', userFilters);
        
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
        const filterContainer = ui.divV([], 'revvity-signals-filter-row');
        
        // Add filter title
        const filterTitle = ui.divText(filterName, {style: {fontWeight: 'bold', marginBottom: '8px'}});
        filterContainer.appendChild(filterTitle);

        // Add appropriate filter controls based on filter type
        if (Object.values(DefaultDateFilters).includes(filterName as DefaultDateFilters)) {
            this.addDateFilterControls(filterContainer, filterName);
        } else if (Object.values(DefaultUserFilters).includes(filterName as DefaultUserFilters)) {
            this.addUserFilterControls(filterContainer);
        }

        // Add close button for the filter
        const closeButton = ui.icons.close(() => {
            this.removeFilter(filterName, filterContainer);
        });
        closeButton.classList.add('revvity-signals-filter-close-button');
        filterContainer.appendChild(closeButton);

        this.filtersDiv.appendChild(filterContainer);
    }

    private addDateFilterControls(container: HTMLElement, filterName: string): void {
        const filterRow = ui.divH([]);
        filterRow.classList.add('revvity-signals-filter-row');
        
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

        // Add inputs to row
        filterRow.appendChild(choiceInput.root);
        filterRow.appendChild(dateContainer);

        // Add 'or' button
        const orButton = ui.button('or', () => {
            this.addDateFilterToContainer(container, filterName);
        });
        orButton.classList.add('revvity-signals-filters-or-button');
        
        filterRow.appendChild(orButton);

        container.appendChild(filterRow);
    }

    private addDateFilterToContainer(container: HTMLElement, filterName: string): void {
        const filterRow = ui.divH([]);
        filterRow.classList.add('revvity-signals-filter-row');
        
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

        // Add inputs to row
        filterRow.appendChild(choiceInput.root);
        filterRow.appendChild(dateContainer);

        // Add 'or' button
        const orButton = ui.button('or', () => {
            this.addDateFilterToContainer(container, filterName);
        });
        orButton.classList.add('revvity-signals-filters-or-button');
        
        filterRow.appendChild(orButton);

        // Add close button for the additional filter row
        const closeButton = ui.icons.close(() => {
            filterRow.remove();
        });
        closeButton.classList.add('revvity-signals-filter-close-button');
        filterRow.appendChild(closeButton);

        container.appendChild(filterRow);
    }

    private addUserFilterControls(container: HTMLElement): void {
        getRevvityUsers().then((res) => {
            if (res) {
                const users = Object.values(res).map((user) => `${user.firstName} ${user.lastName} (${user.userName})`);
                const userChoice = ui.input.multiChoice('', {
                    items: users
                });
                container.appendChild(userChoice.root);
            }
        });
    }

    private removeFilter(filterName: string, filterContainer: HTMLElement): void {
        // Remove from used filters - this will make it available again in the popup
        this.usedFilters.delete(filterName);
        
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

}
