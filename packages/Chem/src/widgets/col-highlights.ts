import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { HIGHLIGHT_BY_SCAFFOLD_TAG } from '../constants';
import { IColoredScaffold } from '../rendering/rdkit-cell-renderer';
import { isSmarts } from '../package';
import { ItemType, ItemsGrid } from '@datagrok-libraries/utils/src/items-grid';
import { Subject, Subscription } from 'rxjs';

export function getmolColumnHighlights(col: DG.Column): DG.Widget {
    return new HighlightWidget(col);
}

export class HighlightWidget extends DG.Widget {

    props = [
        DG.Property.fromOptions({ name: 'molecule', type: DG.TYPE.STRING }),
        DG.Property.fromOptions({ name: 'color', inputType: 'Color', type: DG.TYPE.STRING })
    ];
    latestHighLightTag = '';
    itemsGridChanged = false;
    onMetaDataChangeSub: Subscription;
    col: DG.Column;

    constructor(col: DG.Column) {
        super(ui.div());
        this.col = col;
        this.createItemsGrid();
        this.onMetaDataChangeSub = this.col.dataFrame.onMetadataChanged.subscribe(() => {
            if ((this.latestHighLightTag && this.latestHighLightTag !== col.getTag(HIGHLIGHT_BY_SCAFFOLD_TAG)) && !this.itemsGridChanged) {
                this.createItemsGrid();
            }
            this.itemsGridChanged = false;
        })

    }

    createItemsGrid() {
        const scaffoldsString = this.col.getTag(HIGHLIGHT_BY_SCAFFOLD_TAG);
        const items: IColoredScaffold[] = scaffoldsString ? JSON.parse(scaffoldsString) : [];
        const itemsGrid = new ItemsGrid(this.props, items,
            {
                customInputs: {
                    'molecule': (item: ItemType) => {
                        return new CustomSketcherInput(item.molecule);
                    },
                },
                newItemFunction: () => ({ color: '#00FF00' }),
            });

        itemsGrid.onItemAdded.subscribe(() => {
            this.itemsGridChanged = true;
            this.latestHighLightTag = JSON.stringify(itemsGrid.items);
            this.col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, this.latestHighLightTag);
        });
        itemsGrid.onAddingItemChanged.subscribe((item) => {
            this.itemsGridChanged = true;
            if (itemsGrid.items.indexOf(item.item) === -1) {
                this.latestHighLightTag = JSON.stringify(itemsGrid.items.concat([item.item]));
                this.col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, this.latestHighLightTag);
            }
            if (itemsGrid.items.filter((it => item.item.molecule === it.molecule)).length > 0 && item.fieldName === 'molecule')
                grok.shell.warning(`Current structure has been added previously`);
        });
        itemsGrid.onItemRemoved.subscribe(() => {
            this.itemsGridChanged = true;
            this.latestHighLightTag = JSON.stringify(itemsGrid.items);
            this.col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, this.latestHighLightTag);
        });
        itemsGrid.onItemChanged.subscribe(() => {
            this.itemsGridChanged = true;
            this.latestHighLightTag = JSON.stringify(itemsGrid.items);
            this.col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, this.latestHighLightTag);
        });
        ui.empty(this.root);
        this.root.append(itemsGrid.root);
    };

    detach() {
        this.onMetaDataChangeSub.unsubscribe();
    }

}

export class CustomSketcherInput {
    sketcher: DG.chem.Sketcher;
    root: HTMLElement;
    onChangeFunc: Function = () => { };
    molecule: string | null = '';
    onChanged: Subject<any> = new Subject<any>();
    constructor(molecule: string) {
        this.sketcher = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.EXTERNAL);
        this.sketcher.syncCurrentObject = false;
        if (molecule) {
            this.molecule = molecule;
            isSmarts(molecule) ? this.sketcher.setSmarts(molecule) : this.sketcher.setMolFile(molecule);
        }
        this.root = ui.div(this.sketcher.root, 'chem-col-highlight-sketcher-input');
        this.root.addEventListener('click', e => {
            if (isSmarts(this.molecule!)) {
                e.stopPropagation();
                grok.shell.warning(`Smarts cannot be edited`);
            }
        }, true);
        this.sketcher.onChanged.subscribe(async () => {
            this.molecule = isSmarts(this.molecule!) ? await this.sketcher.getSmarts() : this.sketcher.getMolFile();
            this.onChanged.next();
        });
    }

    get value(): string | null {
        return this.molecule;
    }

    addValidator(f: Function) {}

}
