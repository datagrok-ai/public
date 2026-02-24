export interface DomainInfo {
  description: string;
  category: string;
}

export const DOMAINS_DESCRIPTIONS: {[key: string]: DomainInfo} = {
  ag: {description: 'Additional Genetic Tests', category: 'Genetics / Genomics'},
  bg: {description: 'Background Genetics', category: 'Genetics / Genomics'},
  bw: {description: 'Body Weights', category: 'Findings'},
  cl: {description: 'Clinical Observations', category: 'Findings'},
  co: {description: 'Comments', category: 'Core / Special Purpose'},
  dm: {description: 'Demographics', category: 'Core / Special Purpose'},
  ds: {description: 'Disposition', category: 'Events'},
  eg: {description: 'Electrocardiograms', category: 'Findings'},
  ey: {description: 'Eye Examinations', category: 'Findings'},
  fw: {description: 'Food and Water Consumption', category: 'Findings'},
  ig: {description: 'In-Life Gestational Findings', category: 'Findings'},
  is: {description: 'Immunogenicity', category: 'Findings'},
  lb: {description: 'Laboratory Test Results', category: 'Findings'},
  ma: {description: 'Macroscopic Findings', category: 'Findings'},
  mb: {description: 'Microbiology', category: 'Findings'},
  mi: {description: 'Microscopic Findings', category: 'Findings'},
  ml: {description: 'Morphology / Localization', category: 'Genetics / Genomics'},
  nc: {description: 'Neurobehavioral Findings', category: 'Findings'},
  om: {description: 'Organ Measurements', category: 'Findings'},
  pc: {description: 'Pharmacokinetic Concentrations', category: 'Findings'},
  pm: {description: 'Palpable Masses', category: 'Findings'},
  po: {description: 'Physical Observations', category: 'Findings'},
  pp: {description: 'Pharmacokinetic Parameters', category: 'Pharmacokinetics / Modeling'},
  re: {description: 'Reproductive/Developmental', category: 'Findings'},
  relrec: {description: 'Related Records', category: 'Core / Special Purpose'},
  rp: {description: 'Respiratory', category: 'Findings'},
  sc: {description: 'Subject Characteristics', category: 'Core / Special Purpose'},
  se: {description: 'Subject Elements', category: 'Core / Special Purpose'},
  sm: {description: 'Subject Milestones', category: 'Core / Special Purpose'},
  sr: {description: 'Short-Term Routine Observations', category: 'Findings'},
  sv: {description: 'Subject Visits', category: 'Core / Special Purpose'},
  ta: {description: 'Trial arms', category: 'Trial Design'},
  te: {description: 'Trial elements', category: 'Trial Design'},
  tf: {description: 'Tumor Findings', category: 'Findings'},
  tx: {description: 'Trial Sets', category: 'Trial Design'},
  ts: {description: 'Trial Summary', category: 'Trial Design'},
  vs: {description: 'Vital Signs', category: 'Findings'},
};

export const SUPP_DOMAIN_CATEGORY = 'Supplemental';

export const DOMAINS_CATEGORIES_LIST = ['Core / Special Purpose', 'Events', 'Findings',
  'Genetics / Genomics', 'Pharmacokinetics / Modeling', 'Trial Design', 'Supplemental', 'Other'];
