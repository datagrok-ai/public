import map from '../docs/_static/example_data/S5_iJO1366.Glycolysis_PPP_AA_Nucleotides.json'
import model from '../docs/_static/example_data/iJO1366.json'
import { Builder, libs } from '../src/ts/main'

const sampleMap = new Map();
sampleMap.set('TALA', [1,2,3,4,5]);
window.builder = new Builder( // eslint-disable-line no-new
  map,
  model,
  null,
  libs.d3_select('#root'),
  {
    fill_screen: true,
    never_ask_before_quit: true,
    scroll_behavior: 'zoom',
    samplingFunction: (r) => new Promise((res) => {console.log(r.reactions.find(a => a.id === 'TALA')); res({upper_bound: 1, lower_bound: 0, data: sampleMap})}),
  }
)
