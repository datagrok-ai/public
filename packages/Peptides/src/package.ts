/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import React from 'react';
import ReactDOM from 'react-dom';
import {ProteinLogo} from 'logojs-react';
import {SARViewer} from './peptide-sar-viewer/sar-viewer';
import {AlignedSequenceCellRenderer, AminoAcidsCellRenderer} from './utils/cell-renderer';
import {DataFrame} from 'datagrok-api/dg';
import {splitAlignedPeptides} from './split-aligned';
import {StackedBarChart} from './stacked-barchart/stacked-barchart-viewer';

export const _package = new DG.Package();

async function main(chosenFile : string) {
  const pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  //let peptides =
  //  await grok.data.loadTable('https://datagrok.jnj.com/p/ejaeger.il23peptideidp5562/il-23_peptide_idp-5562');
  const path = _package.webRoot + 'files/' + chosenFile;
  const peptides = (await grok.data.loadTable(path));
  peptides.name = "Peptides";
  peptides.setTag("dataType", "peptides");

  const view = grok.shell.addTableView(peptides);
  view.name = 'PeptidesView';

  view.grid.onCellRender.subscribe(function (args) {
    if(args.cell.isColHeader){
      let textSize = args.g.measureText(args.cell.gridColumn.name);
      args.g.fillText(args.cell.gridColumn.name, args.bounds.x + (args.bounds.width - textSize.width)/2, 
      args.bounds.y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent));
      args.g.fillStyle = '#4b4b4a';
      args.preventDefault();
    }
  });

  pi.close();
}

//name: Peptides
//tags: app
export function Peptides() {
  const appDescription = ui.info(
    [
      ui.span(['For more details see LINK ']),
      ui.divText('\n To start the application :', {style: {'font-weight': 'bolder'}}),
      ui.divText('Select the corresponding .csv table with peptide sequences'),
    ], 'Transform peptide sequence data to research insights',
  );
  const annotationViewerDiv = ui.div();

  const windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showHelp = false;
  windows.showProperties = false;

  const mainDiv = ui.div();
  grok.shell.newView('Peptides', [
    appDescription,
    ui.h2('Choose .csv file'),
    ui.div([
      ui.block25([
        ui.button('Open simple case demo', () => main('aligned.csv'), ''),
        ui.button('Open complex case demo', () => main('aligned_2.csv'), '')
      ]),
      ui.block75([annotationViewerDiv]),
    ]),
    mainDiv,
  ]);
}

//name: SAR Viewer Help
//tags: panel, widget
//input: column _ {semType: alignedSequence}
//output: widget result
export function SARViewerHelp(_: DG.Column): DG.Widget {
  const helpStr =
  'Circle size in the viewer is based on ratio and color is based on MAD.\n' +
  '\n' +
  'Statistics:\n' +
  'Count - number of peptides containing AAR at position\n' +
  'MAD - Median absolute deviation\n' +
  'Median - median of activity\n' +
  'IQR - interquartile range\n' +
  'CQV - coefficient of quartile variation (quartile coefficient of dispersion)\n' +
  'Ratio - share of peptides containing AAR at position\n';
  const div = ui.divV(helpStr.split('\n').map((line) => {
    return ui.divText(line);
  }));
  return new DG.Widget(div);
}

//name: Analyze Peptides
//tags: panel, widgets
//input: column col {semType: alignedSequence}
//output: widget result
export function analyzePeptides(col: DG.Column): DG.Widget {
  // let defaultColumn: DG.Column | null = col;
  let tempCol = null;
  for (const column of col.dataFrame.columns.numerical) {
    tempCol = column;
    break;
  }
  const defaultColumn: DG.Column = col.dataFrame.col('activity') || col.dataFrame.col('IC50') || tempCol;

  const activityColumnChoice = ui.columnInput('Activity column', col.dataFrame, defaultColumn);
  const activityScalingMethod = ui.choiceInput('Activity scaling', 'none', ['none', 'lg', '-lg']);
  // const showHistogram = ui.boolInput('Show histogram', false);

  const startBtn = ui.button('Start', async () => {
    if (activityColumnChoice.value.type === DG.TYPE.FLOAT) {
      const options = {
        'activityColumnColumnName': activityColumnChoice.value.name,
        'activityScalingMethod': activityScalingMethod.value,
        // 'showHistogram': showHistogram.value
      };
      (grok.shell.v as DG.TableView).addViewer('peptide-sar-viewer', options);
    } else {
      grok.shell.error('The activity column must be of floating point number type!');
    }
  });
  //showHistogram removed
  return new DG.Widget(ui.divV([ui.inputs([activityColumnChoice, activityScalingMethod]), startBtn]));
}

//name: peptide-sar-viewer
//description: Peptides SAR Viewer
//tags: viewer
//output: viewer result
export function sar(): SARViewer {
  return new SARViewer();
}

//name: StackedBarchart Widget
//tags: panel, widgets
//input: column col {semType: aminoAcids}
//output: widget result

export async function stackedBarchartWidget(col:DG.Column):Promise<DG.Widget> {
  const viewer = await col.dataFrame.plot.fromType('StackedBarChartAA');
  const panel = ui.divH([viewer.root]);
  return new DG.Widget(panel);
}

//name: StackedBarChartAA
//description: Creates an awesome viewer
//tags: viewer
//output: viewer result
export function stackedBarChart():DG.JsViewer {
  return new StackedBarChart();
}

//name: alignedSequenceCellRenderer
//tags: cellRenderer, cellRenderer-alignedSequence
//meta-cell-renderer-sem-type: alignedSequence
//output: grid_cell_renderer result
export function alignedSequenceCellRenderer() {
  return new AlignedSequenceCellRenderer();
}

//name: aminoAcidsCellRenderer
//tags: cellRenderer, cellRenderer-aminoAcids
//meta-cell-renderer-sem-type: aminoAcids
//output: grid_cell_renderer result
export function aminoAcidsCellRenderer() {
  return new AminoAcidsCellRenderer();
}



class Logo extends DG.JsViewer {
  initialized: boolean;
  option: any;
  colSemType: string;
  splitted: DataFrame | null;
  ppm: Array<Array<number>>;
  reactHost: HTMLDivElement | null;
  PROT_NUMS: { [id: string]: number };
  LET_COLORS: Array<any>;
  target: DataFrame | undefined | null;

  constructor() {
    super();
    this.initialized = false;
    this.colSemType = this.string('colSemType', 'alignedSequence');

    this.splitted = null;
    this.ppm = [];
    this.reactHost = null;
    this.target = null;
    this.PROT_NUMS = {
      'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,
      'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'U': 18, 'V': 19, 'W': 20, 'Y': 21, 'Z': 22,
    };
    this.LET_COLORS = [
      {color: 'rgb(44,160,44)', regex: 'A'},
      {color: 'rgb(44,160,44)', regex: 'B'},
      {color: 'rgb(188,189,34)', regex: 'C'},
      {color: 'rgb(31,119,180)', regex: 'D'},
      {color: 'rgb(30,110,96)', regex: 'E'},
      {color: 'rgb(24,110,79)', regex: 'F'},
      {color: 'rgb(214,39,40)', regex: 'G'},
      {color: 'rgb(158,218,229)', regex: 'H'},
      {color: 'rgb(122, 102, 189)', regex: 'I'},
      {color: 'rgb(108, 218, 229)', regex: 'K'},
      {color: 'rgb(158,218,229)', regex: 'L'},
      {color: 'rgb(141, 124, 217)', regex: 'M'},
      {color: 'rgb(235,137,70)', regex: 'N'},
      {color: 'rgb(255,152,150)', regex: 'P'},
      {color: 'rgb(205, 111, 71)', regex: 'Q'},
      {color: 'rgb(23,190,207)', regex: 'R'},
      {color: 'rgb(255,187,120)', regex: 'S'},
      {color: 'rgb(245,167,100)', regex: 'T'},
      {color: 'rgb(188,189,34)', regex: 'U'},
      {color: 'rgb(23,190,207', regex: 'V'},
      {color: 'rgb(182, 223, 138)', regex: 'W'},
      {color: 'rgb(152,223,138)', regex: 'Y'},
      {color: 'rgb(205, 111, 71)', regex: 'Z'},
    ];
  }

  init() {
    this.initialized = true;
    this.reactHost = ui.div([], {style: {height: '500px', width: '500px'}});
    console.log('INIT');
    this.target = this.dataFrame;
    this.splitted = splitAlignedPeptides(this.dataFrame!.columns.bySemType(this.colSemType));
  }

  onTableAttached() {
    if (typeof this.dataFrame !== 'undefined') {
      if (!this.initialized) {
        this.init();
      }

      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_: any) => this.render()));
      this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_: any) => this.render()));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_: any) => this.render()));
    }

    this.render();
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);

    this.render();
  }

  async render() {
    const bits = this.dataFrame!.selection;
    let selected = false;
    if (bits.trueCount > 0) {
      selected = true;
      this.target = this.dataFrame!
        .groupBy([this.dataFrame!.columns.bySemType(this.colSemType).name])
        .whereRowMask(this.dataFrame!.selection)
        .aggregate();
    }
    if (selected) {
      this.splitted = splitAlignedPeptides(this.target!.columns.bySemType(this.colSemType));
    } else this.splitted = splitAlignedPeptides(this.dataFrame!.columns.bySemType(this.colSemType));
    $(this.root).empty();

    if (typeof this.dataFrame !== 'undefined') {
      this.findLogo();

      if (this.reactHost !== null) {
        this.root.appendChild(this.reactHost);
      }
    }
  }

  async findLogo() {
    this.getInfoFromDf();
    ReactDOM.render(React.createElement(ProteinLogo, {alphabet: this.LET_COLORS, ppm: this.ppm}, null), this.reactHost);
  }

  getInfoFromDf() {
    let index: number = 0;
    this.ppm = [];

    for (const col of this.splitted!.columns) {
      const size = col.length;
      this.ppm.push(new Array(22).fill(0));
      for (let i = 0; i < col.length; i++) {
        if (col.get(i) != '-') {
          this.ppm[index][this.PROT_NUMS[col.get(i)]] += 1 / size;
        }
      }
      index++;
    }
  }
}

//name: Show logo
//tags: panel, widget
//input: column col {semType: alignedSequence}
//output: widget result
export function showl(col: DG.Column): DG.Widget {
  const startBtn = ui.button('Start', async () => {
    const peptidesView = grok.shell.v;
    (<DG.TableView>peptidesView).addViewer('peptide-logo-viewer');
  });
  return new DG.Widget(ui.div(startBtn));
}

//name: peptide-logo-viewer
//tags: viewer, panel
//output: viewer result
export function logov() {
  return new Logo();
}
