/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import SmartLabel from 'fusioncharts-smartlabel';


export const _package = new DG.Package();

//name: info
export function info() {
    grok.shell.info(_package.webRoot);
}

//name: complement
//tags: panel, widgets
//input: string nucleotides {semType: dna_nucleotide}
//output: widget result
//condition: true
export function complement(nucleotides) {
    let complement = new Map([
        ['A', 'T'],
        ['C', 'G'],
        ['G', 'C'],
        ['T', 'A'],
    ]);
    let result = '';
    for (let i = 0; i < nucleotides.length; i++) {
        let complementNucleotide = complement.get(nucleotides[i]);
        if (complementNucleotide) {
            result += complementNucleotide;
        } else {
            result += nucleotides[i];
        }
    }
    return new DG.Widget(ui.divText(result));
}

//*helper method for fuzzyJoin
function makeTable(str) {


    let table = new Array(str.length);
    let maxPrefix = 0;

    table[0] = 0;


    for (let i = 1; i < str.length; i++) {
        while (maxPrefix > 0 && str.charAt(i) !== str.charAt(maxPrefix)) {

            maxPrefix = table[maxPrefix - 1];
        }

        if (str.charAt(maxPrefix) === str.charAt(i)) {
            maxPrefix++;
        }
        table[i] = maxPrefix;
    }
    return table;
}

function kmpMatching(str, word) {
    let count = 0;
    let prefixes = makeTable(word);

    let j = 0;
    let i = 0;
    while (i < str.length) {

        if (str.charAt(i) === word.charAt(j)) {
            i++;
            j++;
        }

        if (j === word.length) {
            count += 1;

            j = prefixes[j - 1];
        } else if (str.charAt(i) !== word.charAt(j)) {

            if (j !== 0) {
                j = prefixes[j - 1];
            } else {

                i++;
            }
        }
    }

    return count;
}

function matchAllSequnces(word, df) {
    let count = 0;
    let iteratorRow = 0;
    let column = df.columns.bySemType('dna_nucleotide')
    let row = column.get(iteratorRow);
    while ((row != null)) {
        // row to string
        count += kmpMatching(row, word);
        iteratorRow += 1;
        row = column.get(iteratorRow);
    }
    return count;
}


function matchAllSequncesForSubstring(str, N, df, preCounted) {
    let countAll = 0;
    for (let i = 0; i + N < str.length; i++) {
        let word = str.substring(i, i + N);
        if (!preCounted[word]) {
            let count = matchAllSequnces(word, df);
            preCounted[word] = count;
            countAll += count;
        } else {
            countAll += preCounted[word];
        }
    }
    return countAll;
}

function addCountToDataFrame(df1, df2, N) {
    let preCounted = {};
    let counts = df1.columns.byName('Count');
    if (counts == null) 
        counts = df1.columns.addNew('Count', 'int');
    
    let iteratorRow = 0;
    let column = df1.columns.bySemType('dna_nucleotide')
    let row1 = column.get(iteratorRow);
    while (row1 != null) {
        let count = matchAllSequncesForSubstring(row1, N, df2, preCounted);
        counts.set(iteratorRow, count);
        iteratorRow += 1;
        row1 = column.get(iteratorRow);
    }
}

//name: fuzzyJoin
//input: dataframe df1
//input: dataframe df2
//input: int N
export function fuzzyJoin(df1, df2, N) {
    addCountToDataFrame(df1, df2, N);
    addCountToDataFrame(df2, df1, N);
    grok.shell.addTableView(df1.append(df2));
}

class NucleotideBoxCellRenderer extends DG.GridCellRenderer {
    get name() {
        return 'Nucleotide cell renderer';
    }

    get cellType() {
        return 'dna_nucleotide';
    }

    render(g, x, y, w, h, gridCell, cellStyle) {
        let seq = gridCell.cell.value;
        const sl = new SmartLabel('id', true);
        console.log(ce)
        sl.setStyle(cellStyle);
        let ctx = g.canvas.getContext("2d");
        ctx.font = '11px courier';
        ctx.fillStyle = "black";
        let labelObj = sl.getSmartText(seq, w, h);
        labelObj = SmartLabel.textToLines(labelObj);
        let lines = labelObj.lines;
        for (let i = 0; i < lines.length; i++)
            ctx.fillText(lines[i], x, y + (i * 11) + 11);
    }
}

//name: nucleotideBoxCellRenderer
//tags: cellRenderer
//meta.cellType: dna_nucleotide
//output: grid_cell_renderer result
export function nucleotideBoxCellRenderer() {
    return new NucleotideBoxCellRenderer();
}

//name: testENASwagger
export async function testENASwagger() {
    let data = await grok.data.query('ArtemSequense:PerformATextSearchAndDownloadDataInXMLFormat', {
        'query': 'coronavirus',
        'result': 'assembly'
    });
    grok.shell.addTableView(data);
}

//name: enaSequence
//tags: panel, widgets
//input: string cellText {semType: dna_nucleotide}
//output: widget result
//condition: isPotentialENAId(cellText)
export async function enaSequence(cellText) {
    const url = `https://www.ebi.ac.uk/ena/browser/api/fasta/${cellText}`;
    const fasta = await (await grok.dapi.fetchProxy(url)).text();
    const lines = fasta.split('\n').slice(1);
    const sequence = lines.reduce((acc, line) => acc + line, '');
    const IdUi = ui.stringInput('ID: ', cellText);
    const sequenceUi = ui.divText(sequence);
    return new DG.Widget(ui.box(ui.splitV([ui.div([IdUi]), sequenceUi])));
}

async function _fetchENASequence(url) {
    const fasta = await (await grok.dapi.fetchProxy(url)).text();
    let sequences = [];
    let sequence = '';
    let lines = fasta.split('\n');
    let started = false;
    for (let i = 0; i < lines.length; i++) {
        let line = lines[i];
        if (line.startsWith('SQ')) {
            sequence = '';
            started = true;
        } else if (line.startsWith('//') && started) {
            sequences.push(sequence);
            sequence = '';
            started = false;
        } else {
            // remove numbers from line
            line = line.replace(/\d+/g, '');
            // remove double spaces
            line = line.replace(/\s+/g, ' ');
            // Remove all line breaks from a string
            line = line.replace(/(\r\n|\n|\r)/gm, "");
            sequence += line;
        }
    }
    let ids = [];
    for (let i = 0; i < lines.length; i++) {
        let line = lines[i];
        if (line.startsWith('ID')) {
            let id = line.split(';')[0].split('   ')[1];
            ids.push(id);
        }
    }
    // shrink ids if ids.length > sequences.length
    if (sequences.length < ids.length) {
        ids = ids.slice(0, sequences.length);
    }
    let df = DG.DataFrame.fromColumns([
        DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', ids),
        DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sequence', sequences)
    ]);
    return df;
}

//name: formENADataTable
export async function formENADataTable() {
    let url = 'https://www.ebi.ac.uk/ena/browser/api/embl/textsearch?result=sequence&query=coronavirus&limit=10';
    let df = await _fetchENASequence(url);
    let grid = DG.Viewer.grid(df);
    let limitInput = ui.intInput('How many rows: ', 100);
    let queryInput = ui.stringInput('Query: ', 'coronavirus');
    let button = ui.button('Preview');
    ui.dialog('Create sequences table')
        .add(ui.splitV([
            ui.splitH([
                ui.span([queryInput.root]),
                button
            ]),
            ui.div([grid]),
            ui.div([limitInput])
        ]))
        .onOK(async () => {
            /* Handle table creation */
            // Display the resulting table
            url = `https://www.ebi.ac.uk/ena/browser/api/embl/textsearch?result=sequence&query=${queryInput.value}&limit=${limitInput.value}`;
            let df = await _fetchENASequence(url);
            console.log(limitInput);
            console.log(queryInput);
            console.log(url);

            grok.shell.addTableView(df);
        })
        .show();
}
