import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function getPTCombineDialog() {

    const createEditor = () => {
        const editorDiv = ui.divH([]);
        const tableEditor = ui.input.table('Table', {value: undefined, tooltipText: 'Table with sequences'});
        const columnEditor = ui.input.choice('Column', {items: [] as string[], value: undefined, tooltipText: 'Sequence column'});

        tableEditor.onChanged.subscribe(async () => {
            const table = tableEditor.value;
            if (!table) {
                columnEditor.items = [];
                columnEditor.value = null;
                return;
            }
            await table.meta.detectSemanticTypes();
            await grok.data.detectSemanticTypes(table);
            const stringColumns = Array.from(table.columns.categorical);
            const stringColumnsNames = stringColumns.map((c) => c.name);
            columnEditor.items = stringColumnsNames;
            columnEditor.value = stringColumns.find((c) => c.semType === DG.SEMTYPE.MACROMOLECULE)?.name ?? stringColumnsNames.find((c) => {
                const lc = c.toLowerCase();
                return lc.includes('seq') || lc.includes('pep');
            }) ?? stringColumnsNames[0];
        });
        editorDiv.appendChild(tableEditor.root);
        editorDiv.appendChild(columnEditor.root);
        editorDiv.style.alignItems = 'center';
        return {root: editorDiv, getValue: () => ({table: tableEditor.value, column: columnEditor.value})};
    }
    type Editor = ReturnType<typeof createEditor>;

    const editors: Editor[] = [];
    const editorsDiv = ui.divV([]);

    const insertNewEditor = (index: number) => {
        const editorDiv = createEditor();
        const removeButton = ui.icons.delete(() => {
            if (!editorDiv.root.parentElement || editors.length < 2)
                return;
            editorDiv.root.remove();
            const editorIndex = editors.indexOf(editorDiv);
            if (editorIndex !== -1) {
                editors.splice(editorIndex, 1);
            }
        }, 'Remove');
        const addEditor = ui.icons.add(() => {
            let editorIndex = editors.indexOf(editorDiv);
            if (editorIndex === -1) {
                editorIndex = editors.length;
            }
            insertNewEditor(editorIndex + 1);
        }, 'Add');
        editorDiv.root.appendChild(removeButton);
        editorDiv.root.appendChild(addEditor);
        removeButton.style.marginLeft = '8px';
        removeButton.style.marginRight = '8px';
        removeButton.style.color = 'var(--blue-1)';
        addEditor.style.color = 'var(--blue-1)';

        
        const prevEditor = editors[index];
        if (prevEditor) {
            editorsDiv.insertBefore(editorDiv.root, prevEditor.root);
        } else {
            editorsDiv.appendChild(editorDiv.root);
        }
        editors.splice(index, 0, editorDiv);
    }
    insertNewEditor(0);

    function validate() {
        const values = editors.map((e) => e.getValue());
        const tables = values.map((v) => v.table);
        const columns = values.map((v) => v.column);
        return tables.every((t) => !!t) && columns.every((c) => !!c);
    }

    const separatorInput = ui.input.string('Separator', {value: '-', tooltipText: 'Separator for sequences', nullable: false});

    ui.dialog('Combine Sequences')
        .add(editorsDiv)
        .add(separatorInput.root)
        .onOK(async () => {
            if (!validate()) {
                grok.shell.error('Please fill all the fields');
                return;
            }
            const values = editors.map((e) => e.getValue());
            const tables = values.map((v) => v.table);
            const columns = values.map((v) => v.column);
            const separator = separatorInput.value;
            const cols = columns.map((c, i) => tables[i]!.col(c!)!.toList().filter((v) => !!v));
            // do a DFS to generate all combinations
            let curIndex = 0;
            const totalLength = cols.reduce((acc, c) => acc * c.length, 1);
            if (totalLength > 10_000_000) {
                grok.shell.error('Too many combinations. Maximum allowed is 10M');
                return;
            }
            const result: string[] = new Array(totalLength).fill(null);

            const dfs = (cur: string, depth: number) => {
                if (depth === cols.length) {
                    result[curIndex++] = cur;
                    return;
                }
                const actCur = `${cur}${cur ? separator : ''}`;
                const col = cols[depth];
                for (let i = 0; i < cols[depth].length; i++) {
                    dfs(actCur + col[i], depth + 1);
                }
            }
            dfs('', 0);
            const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('Combined', result)]);
            df.name = 'Combined Sequences';
            await df.meta.detectSemanticTypes();
            await grok.data.detectSemanticTypes(df);
            grok.shell.addTableView(df);
        })
        .show({resizable: true});

   




}