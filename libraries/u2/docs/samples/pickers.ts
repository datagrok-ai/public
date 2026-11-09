import {ChipsInput, Control, Form, SuggestInput, TagsInput, component} from '@datagrok-libraries/u2';

const ASSAYS = ['Binding (Kd)', 'Cell viability', 'Enzyme inhibition (IC50)', 'hERG patch clamp', 'Solubility'];

export function assayFilters(searchTargets: (query: string, abort: AbortSignal) => Promise<string[]>): Control {
  return component(() => new Form().addAll([
    new TagsInput<string>({label: 'Targets', source: searchTargets, value: ['EGFR'], placeholder: 'Add a target'}),
    new SuggestInput({label: 'Assay', openOnFocus: true,
      source: async (q) => ASSAYS.filter((a) => a.toLowerCase().includes(q.toLowerCase()))}),
    new ChipsInput({label: 'Species', items: ['Human', 'Mouse', 'Rat', 'Dog'], value: ['Human']}),
  ]), {name: 'assay-filters'});
}

// Live demo
const TARGETS = ['ABL1', 'BRAF', 'BTK', 'CDK4', 'EGFR', 'ERBB2', 'JAK2', 'KIT', 'MET', 'PIK3CA', 'SRC'];

async function searchTargets(query: string, abort: AbortSignal): Promise<string[]> {
  await new Promise((resolve) => setTimeout(resolve, 200));
  abort.throwIfAborted();
  return TARGETS.filter((t) => t.toLowerCase().includes(query.toLowerCase()));
}

export function render(host: HTMLElement): Control {
  const demo = assayFilters(searchTargets);
  host.append(demo.root);
  return demo;
}
