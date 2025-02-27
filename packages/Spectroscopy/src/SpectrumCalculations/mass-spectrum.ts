/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as MassSpectrum from 'ms-spectrum';
import * as jcampParser from 'jcampconverter';


export function buildMassSpecViewer(name: string, data: string, molForm?: string) {
  let panel = ui.accordion('testData');
  var result = jcampParser.convert(data);
  //@ts-ignore
  let xColumn = DG.Column.fromFloat32Array('x', (result.flatten[0].spectra[0].data.x ?? []) as Float32Array);
  //@ts-ignore
  let yColumn = DG.Column.fromFloat32Array('y', (result.flatten[0].spectra[0].data.y ?? []) as Float32Array);
  let df = DG.DataFrame.fromColumns([xColumn, yColumn]);
  let lineChart = DG.Viewer.lineChart(df);
  let view = DG.View.fromRoot(lineChart.root);
  view.name = name;
  lineChart.root.style.width = '100%';
  lineChart.root.style.height = '100%';
  if (result.flatten[0].dataType)
    df.setTag('dataType', result.flatten[0].dataType);
  let spec = new MassSpectrum.Spectrum({ x: result.flatten[0].spectra[0].data.x, y: result.flatten[0].spectra[0].data.y });
  let onCurrentViewChangedSub = grok.events.onCurrentViewChanged.subscribe((e: any) => {
    //@ts-ignore
    if (grok.shell.v?.id ?? '' === view.id) {
      buildMassSpecContextPanel(view, lineChart, spec, molForm);
    }
  });

  let onCurrentViewClosedSub = grok.shell.dockManager.onClosed.subscribe((e: any) => {
    if (e.children[0] === view.root) {
      onCurrentViewChangedSub.unsubscribe();
      onCurrentViewClosedSub.unsubscribe();
      if(grok.shell.o === panel.root)
        grok.shell.o = null;
    }
  });
  grok.shell.addView(view);


  function buildMassSpecContextPanel(view: DG.View, lineChartView: DG.Viewer<DG.ILineChartSettings>, data: MassSpectrum.Spectrum, baseMol?: string) {
    const peaksInput = ui.input.choice('Show Peaks', { items: ['None', 'All', 'Best', 'Fragmented', 'Charged By Similarity'], value: 'None' });

    peaksInput.onChanged.subscribe(async (e: any) => {
      peaksOptions.innerHTML = '';
      lineChartView.props.segmentColumnName = '';
      if (peaksInput.value === 'All') {
        if (!lineChartView.dataFrame.columns.names().includes('IsPeak')) {
          const peaks = data.getPeaks(data);
          const peaksMap = new Map(peaks.map(obj => [obj.x, obj.y]));
          addColumnByFunction('IsPeak', lineChartView.dataFrame, (row) => peaksMap.get(row['x']) === row['y'], DG.COLUMN_TYPE.BOOL);
        }
        lineChartView.props.segmentColumnName = 'IsPeak';
      }
      else if (peaksInput.value === 'Best') {
        if (!lineChartView.dataFrame.columns.names().includes('IsBestPeak')) {
          const peaks = data.getBestPeaks(data);
          const peaksMap = new Map(peaks.map(obj => [obj.x, obj.y]));
          addColumnByFunction('IsBestPeak', lineChartView.dataFrame, (row) => peaksMap.get(row['x']) === row['y'], DG.COLUMN_TYPE.BOOL);
        }
        lineChartView.props.segmentColumnName = 'IsBestPeak';
      }
      else if (peaksInput.value === 'Fragmented') {
        peaksOptions.append(buildFragmenetedPeaksForm(lineChartView, data, baseMol))
      }
      else if (peaksInput.value === 'Charged By Similarity') {
        peaksOptions.append(buildPeackChargeBySimilarityForm(lineChartView, data))
      }
    });

    let info: any = data.info;
    info.isContinuous = data.isContinuous();
    const card = ui.tableFromMap(info);

    let peaksOptions = ui.div();

    panel.addTitle(ui.h1(`Mass Spectrum ${view.name}`));
    panel.addPane('Info', () => card, true);
    panel.addPane('Peaks', () => ui.divV([peaksInput, peaksOptions]), false);
    panel.addPane('MARA', () => buildMARAInput(lineChartView, data), false);
    panel.root.classList.add("mass-spectrum-context-panel");
    grok.shell.o = panel.root;
  }

}

function buildMARAInput(lineChartView: DG.Viewer<DG.ILineChartSettings>, data: MassSpectrum.Spectrum): HTMLDivElement {
  const maraMassInput = ui.input.float('Mass', { onValueChanged: (e) => updateLineChart() });
  const maraDeltaInput = ui.input.float('Delta', { min: 0.001, max: 100, value: 0.001, onValueChanged: (e) => updateLineChart() });

  async function updateLineChart() {
    if (maraMassInput.value) {
      let mols = await data.getMassRemainderFct(maraMassInput.value, { delta: maraDeltaInput.value });
      console.log(mols);
      addColumnByFunction('FragnmentedPeak', lineChartView.dataFrame, (df) => {
        //TODO
        return '';
      })
    }
  }
  updateLineChart();
  return ui.divV([maraMassInput, maraDeltaInput]);
}

function buildFragmenetedPeaksForm(lineChartView: DG.Viewer<DG.ILineChartSettings>, data: MassSpectrum.Spectrum, baseMol?: string): HTMLDivElement {
  let baseMolLabel = ui.label(`Molecule Formula: none`);
  baseMolLabel.style.padding = '6px 0';
  if (baseMol)
    baseMolLabel.textContent = (`Molecule Formula: ${baseMol}`);

  let mol = baseMol;
  // const molInput = ui.input.molecule('Molecule', {});
  const sketcher = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.EXTERNAL);
  sketcher.syncCurrentObject = false;
  const from = ui.input.qNum('From', { onValueChanged: (e) => { if (mol) updateLineChart() } });
  const to = ui.input.qNum('To', { onValueChanged: (e) => { if (mol) updateLineChart() } });
  const threshold = ui.input.float('Threshold', { onValueChanged: (e) => { if (mol) updateLineChart() } });
  const limit = ui.input.qNum('Limit', { onValueChanged: (e) => { if (mol) updateLineChart() } });

  sketcher.onChanged.subscribe(async (e) => {
    let smiles = sketcher.getSmiles();
    if (smiles) {
      mol = await grok.functions.call('Chem:getMolecularFormula', { molecule: smiles });
      baseMolLabel.textContent = (`Molecule Formula: ${mol}`);
    }
    updateLineChart();
  });

  async function updateLineChart() {
    if (mol) {
      let mols = await data.getFragmentPeaksFct(mol, { from: from.value, to: to.value, threshold: threshold.value, limit: limit.value });
      console.log(mols);
      addColumnByFunction('FragnmentedPeak', lineChartView.dataFrame, (df) => {
        //TODO
        return '';
      })
    }
  }
  updateLineChart();
  return ui.divV([baseMolLabel, ui.div(sketcher.root, { style: { position: 'relative' } }), from, to, threshold]);
}

//TODO
//  data.getSelectedPeaksWithCharge()
// function buildSelectedPeaksWithChargeForm(lineChartView: DG.Viewer<DG.ILineChartSettings>, data: Spectrum): HTMLDivElement{
//   const mass = ui.input.qNum('Mass', { onValueChanged: (e) => { updateLineChart() }, value: 1 });
//   const maxCharge = ui.input.qNum('Min Charge', { onValueChanged: (e) => { updateLineChart() }, value: 1 });
//   const minCharge = ui.input.qNum('Max Charge', { onValueChanged: (e) => { updateLineChart() }, value: 10 });
//   // const widthBottom = ui.input.qNum('mass', { onValueChanged: (e) => { updateLineChart() } , value: 1});
//   // const widthTop = ui.input.qNum('mass', { onValueChanged: (e) => { updateLineChart() } , value: 1});
//   // const similarityZoneLow = ui.input.qNum('mass', { onValueChanged: (e) => { updateLineChart() } , value: 1});
//   // const similarityZoneHigh = ui.input.qNum('mass', { onValueChanged: (e) => { updateLineChart() } , value: 1});
//   async function updateLineChart() {
//     if (mass.value) {
//       let mols = await data.getSelectedPeaksWithCharge(mass.value, { maxCharge: maxCharge.value, minCharge: minCharge.value });
//       console.log(mols);
//       addColumnByFunction('FragnmentedPeak', lineChartView.dataFrame, (df) => {
//         return '';
//       })
//     }
//   }
//   updateLineChart();
//   return ui.divV([mass, maxCharge, minCharge]);
// }

function buildPeackChargeBySimilarityForm(lineChartView: DG.Viewer<DG.ILineChartSettings>, data: MassSpectrum.Spectrum) {
  const mass = ui.input.qNum('Mass', { onValueChanged: (e) => { updateLineChart() }, value: 1 });
  const maxCharge = ui.input.qNum('Min Charge', { onValueChanged: (e) => { updateLineChart() }, value: 1 });
  const minCharge = ui.input.qNum('Max Charge', { onValueChanged: (e) => { updateLineChart() }, value: 10 });
  // const widthBottom = ui.input.qNum('mass', { onValueChanged: (e) => { updateLineChart() } , value: 1});
  // const widthTop = ui.input.qNum('mass', { onValueChanged: (e) => { updateLineChart() } , value: 1});
  // const similarityZoneLow = ui.input.qNum('mass', { onValueChanged: (e) => { updateLineChart() } , value: 1});
  // const similarityZoneHigh = ui.input.qNum('mass', { onValueChanged: (e) => { updateLineChart() } , value: 1});
  //TO DO
  // * @param {number}   [options.minCharge=1]
  // * @param {number}   [options.maxCharge=10]
  // * @param {object}   [options.similarity={}]
  // * @param {number}   [options.similarity.widthBottom]
  // * @param {number}   [options.similarity.widthTop]
  // * @param {object}   [options.similarity.widthFunction] - function called with mass that should return an object width containing top and bottom
  // * @param {object}   [options.similarity.zone={}]
  // * @param {number}   [options.similarity.zone.low=-0.5] - window shift based on observed monoisotopic mass
  // * @param {number}   [options.similarity.zone.high=2.5] - 

  async function updateLineChart() {
    if (mass.value) {
      let mols = await data.getPeakChargeBySimilarity(mass.value, { maxCharge: maxCharge.value, minCharge: minCharge.value });
      console.log(mols);
      addColumnByFunction('FragnmentedPeak', lineChartView.dataFrame, (df) => {
        //TODO
        return '';
      })
    }
  }
  updateLineChart();
  return ui.divV([mass, maxCharge, minCharge]);
}

function addColumnByFunction(name: string, df: DG.DataFrame, func: ((df: DG.Row) => any), colType: DG.ColumnType = DG.COLUMN_TYPE.STRING) {
  if (df.columns.names().includes(name))
    df.columns.remove(name);
  let resultList = [];
  for (let row of df.rows) {
    resultList.push(func(row));
  }

  let peaksCol = DG.Column.fromList(colType, name, resultList);
  df.columns.add(peaksCol);
}
