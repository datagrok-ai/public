import type {CobraModelData} from '../../escher_src/src/ts/types';

export class Model {
  id: string = '';
  name: string = '';
  description: string = '';
  notes: any;
  metabolites: CobraModelData['metabolites'] = [];
  reactions: CobraModelData['reactions'] = [];
  genes: CobraModelData['genes'] = [];
}

export class Solution {
  fluxes: { [key: string]: number; };
  objectiveValue: number;
  constructor (objectiveValue: number, fluxes: { [key: string]: number }) {
    this.objectiveValue = objectiveValue;
    this.fluxes = fluxes;
  }
}

export function modelFromWorkerData (data: CobraModelData | null) {
  const model = new Model();
  model.reactions = data!.reactions;
  model.metabolites = data!.metabolites;
  model.genes = data!.genes;
  model.id = data!.id;
  model.notes = data!.notes;
  model.description = data!.description;
  return model;
}

export function modelFromJsonData (data: CobraModelData | null) {
  if (data === null) return null;

  const model = new Model();

  model.reactions = data.reactions.map(r => {
    let coeff = 0;
    if (r.objective_coefficient && r.objective_coefficient !== 0) {
      if (r.objective_coefficient < 0) coeff = -1;
      else coeff = 1;
    }
    return ({ ...r, objective_coefficient: coeff });
  });

  model.metabolites = data.metabolites.map(x => ({...x}));
  model.genes = data.genes.map(x => ({...x}));
  model.id = data.id;
  model.notes = data.notes;
  model.description = data.description;
  return model;
}
