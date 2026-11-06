/**
 * Treated-arm archetypes per kind of therapy: the same declared effect is reached by different curves.
 * The drug acts as a kill rate kappa (1/day) on the control's Gompertz growth from `onsetDays` after the
 * start of treatment, d log V / dt = alpha - beta (log V - log v_ref) - kappa; kappa varies between animals
 * (log-normal, sd `kappaSd`). `responders` < 1 makes the arm bimodal (the others grow as controls).
 *
 * Defaults estimated offline from the atlas's treated arms (2026-10-21):
 * 1,670 treated blocks classified by modality and mechanism; mono arms with a base-pool control, TGI 30-80%.
 *   kappaSd   matched to the arms' spread of the animals' log growth against their control's (directly acting
 *             1.0; immune small molecules 1.4; immune antibodies / proteins 2.4; cell therapy 1.0-1.2)
 *   onset     the day the arm separates from its control (mean log V 0.2 below) is reproduced with no onset
 *             (2-4 days come from the kill rate building up); cell therapy separates ~2 days later
 *   responders  no bimodality: the animals' growth is one continuous distribution at 30-80% and at > 80% TGI
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
    basis: 'atlas: 11 directly acting small-molecule arms (spread 1.0)'},
  biologic: {label: 'Biologic (antibody, ADC, fusion protein)', onsetDays: 0, responders: 1, kappaSd: 0.3,
    translation: 'PK scaling and target occupancy',
    basis: 'assumed equal to directly acting small molecules: 1 directly acting biologic arm in the atlas'},
  immunotherapy: {label: 'Immunotherapy biologic (checkpoint antibody, agonist antibody, cytokine)', onsetDays: 0,
    responders: 1, kappaSd: 0.85, translation: 'receptor occupancy / pharmacodynamics; no TSC',
    basis: 'atlas: 30 immune-mediated antibody and protein arms (spread 2.4: the most heterogeneous)'},
  immune_agonist: {label: 'Immune agonist, small molecule (STING, TLR, innate)', onsetDays: 0, responders: 1,
    kappaSd: 0.55,
    translation: 'exposure at the target / pharmacodynamics; no TSC',
    basis: 'atlas: 48 immune-mediated small-molecule arms (spread 1.4)'},
  cell_therapy: {label: 'Cell therapy', onsetDays: 2, responders: 1, kappaSd: 0.35,
    translation: 'cell dose and expansion; no TSC',
    basis: 'atlas: 6-21 cell-therapy arms (spread 1.0-1.2, separates ~2 days later)'},
  vaccine_peptide: {label: 'Vaccine or peptide', onsetDays: 3, responders: 1, kappaSd: 0.5,
    translation: 'immunogenicity; no TSC',
    basis: 'atlas: few arms (6-10): weak estimate'},
  oligo_rna: {label: 'Oligonucleotide or RNA', onsetDays: 1, responders: 1, kappaSd: 0.6,
    translation: 'PK scaling and target knockdown',
    basis: 'atlas: 15 oligo / RNA arms (spread 1.6)'},
};

export const THERAPIES = Object.keys(ARCHETYPES) as Therapy[];
