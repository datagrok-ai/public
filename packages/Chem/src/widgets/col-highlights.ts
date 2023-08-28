import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { HIGHLIGHT_BY_SCAFFOLD_TAG } from '../constants';
import { IColoredScaffold } from '../rendering/rdkit-cell-renderer';


export function getmolColumnHighlights(col: DG.Column): DG.Widget {

    let editingMoleculeIdx: number | null = null;
    const scaffoldsString = col.getTag(HIGHLIGHT_BY_SCAFFOLD_TAG);

    const scaffolds: IColoredScaffold[] = scaffoldsString ? JSON.parse(scaffoldsString) : [];
    const scaffoldsDf = DG.DataFrame.create(scaffolds.length);
    const molCol = scaffoldsDf.columns.addNewString('Molecule').init((i) => scaffolds[i].molString);
    molCol.semType = DG.SEMTYPE.MOLECULE;
    scaffoldsDf.columns.addNewString('Color').init((i) => rgbPercentToHex(scaffolds[i].color!));
    scaffoldsDf.columns.addNewString('Edit');
    scaffoldsDf.columns.addNewString('Delete');

    const grid = scaffoldsDf.plot.grid();

    let molEditCol = grid.columns.byName('Edit')!;
    molEditCol.cellType = 'html';

    const resetSketcherAndFireDfChanges = () => {
        col.dataFrame.fireValuesChanged();
        sketcher.setMolFile(DG.WHITE_MOLBLOCK);
        sketcher.updateExtSketcherContent();
    }


    grid.onCellPrepare(function (gc) {
        if (gc.isTableCell && gc.gridColumn.name === 'Edit') {
            gc.style.element = ui.icons.edit(() => {
                colorPicker.value = gc.grid.dataFrame.get('Color', gc.tableRowIndex!)
                sketcher.setSmarts(gc.grid.dataFrame.get('Molecule', gc.tableRowIndex!));
                sketcher.updateExtSketcherContent();
                editingMoleculeIdx = gc.tableRowIndex!;
                ui.empty(buttonDiv);
                buttonDiv.append(saveEditedButton);
            })
        }
    });

    let molDeleteCol = grid.columns.byName('Delete')!;
    molDeleteCol.cellType = 'html';

    grid.onCellPrepare(function (gc) {
        if (gc.isTableCell && gc.gridColumn.name === 'Delete') {
            gc.style.element = ui.icons.delete(() => {
                scaffoldsDf.rows.removeAt(gc.tableRowIndex!, 1);
                const scaffolds = scaffoldsDf.rowCount ? convertDfToObj(scaffoldsDf) : [];
                col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(scaffolds));
                resetSketcherAndFireDfChanges();
            });
        }
    });

    const addMoleculeButton = ui.button('Add', async () => {
        const currentSmarts = await sketcher.getSmarts();
        if (currentSmarts) {
            const scaffolds = convertDfToObj(scaffoldsDf);
            const substrIdx = scaffolds.findIndex((it) => it.molString === currentSmarts);
            if (scaffoldsDf.rowCount === 0 || substrIdx == -1) {
                scaffoldsDf.rows.addNew([currentSmarts, colorPicker.value]);
                scaffolds.push({molString: currentSmarts, color: hexToPercentRgb(colorPicker.value)!});
            } else {
                scaffoldsDf.set('Color', substrIdx, colorPicker.value);
                scaffolds[substrIdx].color = hexToPercentRgb(colorPicker.value)!;
            }
            col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(scaffolds));
            resetSketcherAndFireDfChanges();
        }
    })

    const saveEditedButton = ui.button('Save', async () => {
        const currentSmarts = await sketcher.getSmarts();
        if (currentSmarts) {
            scaffoldsDf.col('Molecule')?.set(editingMoleculeIdx!, currentSmarts);
            scaffoldsDf.col('Color')?.set(editingMoleculeIdx!, colorPicker.value);
            editingMoleculeIdx = null;
            const scaffolds = convertDfToObj(scaffoldsDf);
            col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(scaffolds));
            resetSketcherAndFireDfChanges();
            ui.empty(buttonDiv);
            buttonDiv.append(addMoleculeButton);
        }
    })

    const buttonDiv = ui.div(addMoleculeButton);

  const sketcher = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.EXTERNAL);
  sketcher.syncCurrentObject = false;
  sketcher.setMolFile(DG.WHITE_MOLBLOCK);
  sketcher.root.classList.add('ui-input-editor');
  sketcher.root.style.marginTop = '3px';

  const colorPicker = ui.colorInput('Color', '#ff0000');
   
    const host = ui.divV([
        sketcher,
        colorPicker,
        buttonDiv,
        grid
    ])

    return new DG.Widget(host);
}

function convertDfToObj(df: DG.DataFrame): IColoredScaffold[] {
    const arr: IColoredScaffold[] = [];

    for (let i = 0; i < df.rowCount; i++) {
        arr.push({molString: df.get('Molecule', i), color: hexToPercentRgb(df.get('Color', i))!})
    }

    return arr;
}

export function hexToPercentRgb(hex: string) {
    const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);   
    return result ? [
        parseInt(result[1], 16)/256,
        parseInt(result[2], 16)/256,
        parseInt(result[3], 16)/256,
    ] : null;
  }

  function componentToHex(c: number): string {
    var hex = (c*256).toString(16);
    return hex.length == 1 ? "0" + hex : hex;
  }
  
  function rgbPercentToHex(color: number[]): string {
    return "#" + componentToHex(color[0]) + componentToHex(color[1]) + componentToHex(color[2]);
  }
