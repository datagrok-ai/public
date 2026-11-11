/**
 * Treated-arm archetypes per kind of therapy: the same declared effect is reached by different curves.
 * The drug acts as a kill rate kappa (1/day) on the control's Gompertz growth from `onsetDays` after the
 * start of treatment, d log V / dt = alpha - beta (log V - log v_ref) - kappa; kappa varies between animals
 * (log-normal, sd `kappaSd`). `responders` < 1 makes the arm bimodal (the others grow as controls).
 *
 * Defaults estimated offline from the atlas's treated arms (2026-10-21, the 3 original lines):
 * 1,670 treated blocks classified by modality and mechanism; mono arms with a base-pool control, TGI 30-80%.
 *   kappaSd   matched to the arms' spread of the animals' log growth against their control's (directly acting
 *             1.0; immune small molecules 1.4; immune antibodies / proteins 2.4; cell therapy 1.0-1.2; oligo / RNA
 *             1.6). How a spread became a kappaSd was not recorded and cannot be reconstructed.
 * Checked on 2026-11-11 against the 12 atlas lines (v3 regeneration; preclinical bundle analysis/archetype_check.py):
 * simulating each of 322 single-agent arms (105 papers) with these constants, the engine's 90% band holds 89-95% of
 * the observed animals in every archetype but vaccine / peptide (83%, 7 papers), papers as the unit. Constants
 * re-derived from the 12 lines did not do better without losing coverage in some line, so these stay. They are
 * conservative: the engine's arms are wider than the observed ones (observed / simulated spread 0.81 for small
 * molecules, 0.72 for immune antibodies).
 *   onset     the day the arm separates from its control (mean log V 0.2 below) is reproduced with no onset
 *             (2-4 days come from the kill rate building up); cell therapy separates ~2 days later
 *   responders  no bimodality: the animals' growth is one continuous distribution at 30-80% and at > 80% TGI
 *             (12 lines: MC38's immune antibodies have a tail of regressions, 13% of animals, not a second mode)
 * Published arms are mostly the ones that worked: the shape (onset, spread) is what these defaults carry.
 */
export type Therapy = 'small_molecule' | 'biologic' | 'immunotherapy' | 'immune_agonist' | 'cell_therapy' |
  'vaccine_peptide' | 'oligo_rna';

export interface Archetype {
  label: string;
  onsetDays: number;
  responders: number;
  kappaSd: number;
  /** Translate (later): how a mouse effect becomes a human dose for this kind of therapy. */
  translation: string;
  /** where the defaults come from */
  basis: string;
}

export const ARCHETYPES: {[t in Therapy]: Archetype} = {
  small_molecule: {label: 'Small molecule (cytotoxic or targeted)', onsetDays: 0, responders: 1, kappaSd: 0.3,
    translation: 'tumour static concentration (TSC = alpha / k2) against human exposure',
    basis: 'atlas: 59 directly acting small-molecule arms, 37 papers (spread 1.0)'},
  biologic: {label: 'Biologic (antibody, ADC, fusion protein)', onsetDays: 0, responders: 1, kappaSd: 0.3,
    translation: 'PK scaling and target occupancy',
    basis: 'assumed equal to directly acting small molecules: 11 directly acting biologic arms from 4 papers ' +
      'are too few alone (with small molecules: 70 arms, spread 1.1)'},
  immunotherapy: {label: 'Immunotherapy biologic (checkpoint antibody, agonist antibody, cytokine)', onsetDays: 0,
    responders: 1, kappaSd: 0.85, translation: 'receptor occupancy / pharmacodynamics; no TSC',
    basis: 'atlas: 51 immune-mediated antibody and protein arms, 33 papers (spread 1.3; MC38 the widest)'},
  immune_agonist: {label: 'Immune agonist, small molecule (STING, TLR, innate)', onsetDays: 0, responders: 1,
    kappaSd: 0.55,
    translation: 'exposure at the target / pharmacodynamics; no TSC',
    basis: 'atlas: 56 immune-mediated small-molecule arms, 17 papers (spread 1.5)'},
  cell_therapy: {label: 'Cell therapy', onsetDays: 2, responders: 1, kappaSd: 0.35,
    translation: 'cell dose and expansion; no TSC',
    basis: 'atlas: 6 cell-therapy arms from 5 papers (spread 1.2, separates ~2 days later), as immune antibodies'},
  vaccine_peptide: {label: 'Vaccine or peptide', onsetDays: 3, responders: 1, kappaSd: 0.5,
    translation: 'immunogenicity; no TSC',
    basis: 'atlas: 10 arms from 5 papers: weak estimate'},
  oligo_rna: {label: 'Oligonucleotide or RNA', onsetDays: 1, responders: 1, kappaSd: 0.6,
    translation: 'PK scaling and target knockdown',
    basis: 'atlas: 5 oligo / RNA arms from 3 papers (spread 1.6): weak estimate'},
};

export const THERAPIES = Object.keys(ARCHETYPES) as Therapy[];
