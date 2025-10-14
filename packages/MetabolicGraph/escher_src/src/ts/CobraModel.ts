/* eslint-disable camelcase */
/**
 * CobraModel
 */

import * as utils from './utils';
import * as dataStyles from './dataStyles';
import {CobraGene, CobraModelData} from './types';

/**
 * Return a reaction string for the given stoichiometries. Adapted from
 * cobra.core.Reaction.build_reaction_string().
 * @param {Object} stoichiometries - An object with metabolites as keys and
 * stoichiometries as values.
 * @param {Boolean} is_reversible - Whether the reaction is reversible.
 */
export function build_reaction_string(stoichiometries: { [key: string]: number }, is_reversible: boolean) {
  const format = function(number: number) {
    if (number == 1)
      return '';

    return String(number) + ' ';
  };
  const reactant_bits = [];
  const product_bits = [];
  for (const the_metabolite in stoichiometries) {
    const coefficient = stoichiometries[the_metabolite];
    if (coefficient > 0)
      product_bits.push(format(coefficient) + the_metabolite);
    else
      reactant_bits.push(format(Math.abs(coefficient)) + the_metabolite);
  }
  let reaction_string = reactant_bits.join(' + ');
  if (is_reversible)
    reaction_string += ' ↔ ';
  else
    reaction_string += ' → ';

  reaction_string += product_bits.join(' + ');
  return reaction_string;
}

export function from_cobra_json(model_data: CobraModelData) {
  /** Use a JSON Cobra model exported by COBRApy to make a new CobraModel
      object.

      The COBRA "id" becomes a "bigg_id", and "upper_bound" and "lower_bound"
      bounds become "reversibility".

      Fills out a "genes" list.

  */
  // reactions and metabolites
  if (!(model_data.reactions && model_data.metabolites))
    throw new Error('Bad model data.');


  // make a gene dictionary
  const genes: { [key: string]: CobraGene } = {};
  for (let i = 0, l = model_data.genes.length; i < l; i++) {
    const r = model_data.genes[i];
    const the_id = r.id;
    genes[the_id] = r;
  }

  const model = new CobraModel();

  model.reactions = {};
  for (let i = 0, l = model_data.reactions.length; i < l; i++) {
    const r = model_data.reactions[i];
    const the_id = r.id;
    const reaction = utils.clone(r);
    delete reaction.id;
    reaction.bigg_id = the_id;
    reaction.data_string = '';
    // add the appropriate genes
    reaction.genes = [];

    // set reversibility
    reaction.reversibility = (reaction.lower_bound < 0 && reaction.upper_bound > 0);
    if (reaction.upper_bound <= 0 && reaction.lower_bound < 0) {
      // reverse stoichiometries
      for (const met_id in reaction.metabolites)
        reaction.metabolites[met_id] = -reaction.metabolites[met_id];
    }
    delete reaction.lower_bound;
    delete reaction.upper_bound;

    if ('gene_reaction_rule' in reaction) {
      const gene_ids = dataStyles.genes_for_gene_reaction_rule(reaction.gene_reaction_rule!);
      gene_ids.forEach(function(gene_id) {
        if (gene_id in genes) {
          const gene = utils.clone(genes[gene_id]);
          // rename id to bigg_id
          gene.bigg_id = gene.id;
          delete gene.id;
          reaction.genes!.push(gene);
        } else
          console.warn('Could not find gene for gene_id ' + gene_id);
      });
    }
    model.reactions[the_id] = reaction;
  }
  model.metabolites = {};
  for (let i = 0, l = model_data.metabolites.length; i < l; i++) {
    const r = model_data.metabolites[i];
    const the_id = r.id;
    const met = utils.clone(r);
    delete met.id;
    met.bigg_id = the_id;
    model.metabolites[the_id] = met;
  }
  return model;
}


export default class CobraModel {
  reactions: { [key: string]: any };
  metabolites: { [key: string]: any };


  // instance methods
  constructor() {
    this.reactions = {};
    this.metabolites = {};
  }

  /**
 * Apply data to model. This is only used to display options in
 * BuildInput. apply_reaction_data overrides apply_gene_data.
 */
  apply_reaction_data(reaction_data: { [key: string]: (number | string | null)[]}, styles: string[], compare_style: string) {
    dataStyles.apply_reaction_data_to_reactions(this.reactions, reaction_data,
      styles, compare_style);
  }

  /**
 * Apply data to model. This is only used to display options in BuildInput.
 */
  apply_metabolite_data(metabolite_data: { [key: string]: (number | string | null)[]}, styles: string[], compare_style: string) {
    dataStyles.apply_metabolite_data_to_nodes(this.metabolites, metabolite_data,
      styles, compare_style);
  }

  /**
 * Apply data to model. This is only used to display options in
 * BuildInput. Overrides apply_reaction_data.
 */
  apply_gene_data(gene_data_obj: { [key: string]: { [key: string]: (number | string | null)[] } }, styles: string[], identifiers_on_map: string,
    compare_style: string, and_method_in_gene_reaction_rule: string) {
    dataStyles.apply_gene_data_to_reactions(this.reactions, gene_data_obj,
      styles, identifiers_on_map,
      compare_style,
      and_method_in_gene_reaction_rule);
  }

  static from_cobra_json(model_data: CobraModelData) {
    return from_cobra_json(model_data);
  }

  static build_reaction_string(stoichiometries: { [key: string]: number }, is_reversible: boolean) {
    return build_reaction_string(stoichiometries, is_reversible);
  }
}

// class methods
