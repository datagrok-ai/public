import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { HIGHLIGHT_BY_SCAFFOLD_TAG } from '../constants';
import { IColoredScaffold } from '../rendering/rdkit-cell-renderer';
import { isSmarts } from '../package';
import { ItemType, ItemsGrid } from '@datagrok-libraries/utils/src/items-grid';

export function getmolColumnHighlights(col: DG.Column): DG.Widget {

    const props = [
        DG.Property.fromOptions({ name: 'molecule', type: DG.TYPE.STRING }),
        DG.Property.fromOptions({ name: 'color', inputType: 'Color', type: DG.TYPE.STRING })
    ];
    const scaffoldsString = col.getTag(HIGHLIGHT_BY_SCAFFOLD_TAG);
    const items: IColoredScaffold[] = scaffoldsString ? JSON.parse(scaffoldsString) : [];
    const itemsGrid = new ItemsGrid(props, items,
        {
            customInputs: {
                'molecule': (item: ItemType) => {
                    return new CustomSketcherInput(item.molecule);
                },
            },
            newItemFunction: () => ({color: '#00FF00'}),
        });

    itemsGrid.onItemAdded.subscribe(() => col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(itemsGrid.items)));
    itemsGrid.onItemRemoved.subscribe(() => col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(itemsGrid.items)));
    itemsGrid.onItemChanged.subscribe(() => col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(itemsGrid.items)));
    return new DG.Widget(itemsGrid.root);
}

export class CustomSketcherInput {
    sketcher: DG.chem.Sketcher;
    root: HTMLElement;
    onChangeFunc: Function = () => { };
    molecule: string | null = '';
    constructor(molecule: string) {
        this.sketcher = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.EXTERNAL);
        this.sketcher.syncCurrentObject = false;
        if (molecule) {
            this.molecule = molecule;
            isSmarts(molecule) ? this.sketcher.setSmarts(molecule) : this.sketcher.setSmiles(molecule);
        }
        this.root = ui.div(this.sketcher.root, 'chem-col-highlight-sketcher-input');
        this.root.addEventListener('click', e => {
            if (isSmarts(this.molecule!)) {
                e.stopPropagation();
                grok.shell.warning(`Smarts cannot be edited`);
            }
        }, true);
        this.sketcher.onChanged.subscribe(async () => {
            this.molecule = isSmarts(this.molecule!) ? await this.sketcher.getSmarts() : this.sketcher.getSmiles();
            this.onChangeFunc();
        });
    }

    get value(): string | null {
        return this.molecule;
    }

    onChanged(f: Function) {
        this.onChangeFunc = f;
    }

}
