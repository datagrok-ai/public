export interface DomainInfo {
  description: string;
  category: string;
}

export const DOMAINS_DESCRIPTIONS: {[key: string]: DomainInfo} = {
  // Core / Special Purpose
  dm: {description: 'Demographics', category: 'Core / Special Purpose'},
  co: {description: 'Comments', category: 'Core / Special Purpose'},
  relrec: {description: 'Related Records', category: 'Core / Special Purpose'},
  sc: {description: 'Subject Characteristics', category: 'Core / Special Purpose'},
  se: {description: 'Subject Elements', category: 'Core / Special Purpose'},
  sm: {description: 'Subject Milestones', category: 'Core / Special Purpose'},
  sv: {description: 'Subject Visits', category: 'Core / Special Purpose'},
  // Trial Design
  ta: {description: 'Trial Arms', category: 'Trial Design'},
  te: {description: 'Trial Elements', category: 'Trial Design'},
  ti: {description: 'Trial Inclusion/Exclusion Criteria', category: 'Trial Design'},
  td: {description: 'Trial Disease Assessments', category: 'Trial Design'},
  tm: {description: 'Trial Disease Milestones', category: 'Trial Design'},
  ts: {description: 'Trial Summary', category: 'Trial Design'},
  tv: {description: 'Trial Visits', category: 'Trial Design'},
  tx: {description: 'SEND Trial Summary', category: 'Trial Design'},
  // Interventions
  cm: {description: 'Concomitant Medications', category: 'Interventions'},
  ec: {description: 'Exposure as Collected', category: 'Interventions'},
  ex: {description: 'Exposure', category: 'Interventions'},
  pr: {description: 'Procedures', category: 'Interventions'},
  su: {description: 'Substance Use', category: 'Interventions'},
  da: {description: 'Drug Accountability', category: 'Interventions'},
  // Events
  ae: {description: 'Adverse Events', category: 'Events'},
  ce: {description: 'Clinical Events', category: 'Events'},
  dd: {description: 'Death Details', category: 'Events'},
  ds: {description: 'Disposition', category: 'Events'},
  dv: {description: 'Protocol Deviations', category: 'Events'},
  ho: {description: 'Healthcare Encounters', category: 'Events'},
  ie: {description: 'Inclusion/Exclusion Criteria Not Met', category: 'Events'},
  // Findings - General / Clinical Findings
  lb: {description: 'Laboratory Test Results', category: 'Findings'},
  vs: {description: 'Vital Signs', category: 'Findings'},
  pe: {description: 'Physical Examination', category: 'Findings'},
  eg: {description: 'Electrocardiograms', category: 'Findings'},
  ft: {description: 'Functional Tests', category: 'Findings'},
  qs: {description: 'Questionnaires', category: 'Findings'},
  oe: {description: 'Objective Evidence', category: 'Findings'},
  cl: {description: 'Clinical Observations', category: 'Findings'},
  fa: {description: 'Findings About', category: 'Findings'},
  // Findings - SEND In-Life / Observational Findings
  bw: {description: 'Body Weights', category: 'Findings'},
  fw: {description: 'Food and Water Consumption', category: 'Findings'},
  po: {description: 'Physical Observations', category: 'Findings'},
  sr: {description: 'Short-Term Routine Observations', category: 'Findings'},
  ey: {description: 'Eye Examinations', category: 'Findings'},
  ig: {description: 'In-Life Gestational Findings', category: 'Findings'},
  // Findings - SEND System-Specific Findings
  cv: {description: 'Cardiovascular System Findings', category: 'Findings'},
  cg: {description: 'Cardiovascular System Findings', category: 'Findings'},
  hp: {description: 'Hematology Findings', category: 'Findings'},
  la: {description: 'Laboratory Findings', category: 'Findings'},
  nc: {description: 'Neurobehavioral/Cognitive Findings', category: 'Findings'},
  nv: {description: 'Nervous System Findings', category: 'Findings'},
  re: {description: 'Reproductive System Findings', category: 'Findings'},
  ur: {description: 'Urinalysis Findings', category: 'Findings'},
  // Findings - Pathology Findings (SEND)
  ma: {description: 'Macropathology Findings', category: 'Findings'},
  mk: {description: 'Microscopic Findings', category: 'Findings'},
  mo: {description: 'Morphology Findings', category: 'Findings'},
  om: {description: 'Organ Measurements', category: 'Findings'},
  // Findings - Microbiology / Immunology
  mi: {description: 'Microbiology Findings (SEND)', category: 'Findings'},
  mb: {description: 'Microbiology Findings (Clinical)', category: 'Findings'},
  ms: {description: 'Microbiology Susceptibility', category: 'Findings'},
  is: {description: 'Immunogenicity Specimen Assessments', category: 'Findings'},
  // Genetics / Genomics
  bg: {description: 'Background Genetics', category: 'Genetics / Genomics'},
  gf: {description: 'Genetic Features', category: 'Genetics / Genomics'},
  pc: {description: 'Pharmacogenetics/Genomics', category: 'Genetics / Genomics'},
  ag: {description: 'Additional Genetic Tests', category: 'Genetics / Genomics'},
  ml: {description: 'Morphology/Localization', category: 'Genetics / Genomics'},
  // Pharmacokinetics / Modeling
  pp: {description: 'Pharmacokinetic Parameters', category: 'Pharmacokinetics / Modeling'},
  pm: {description: 'Pharmacometrics', category: 'Pharmacokinetics / Modeling'},
  // Associated Persons
  apce: {description: 'Associated Persons Clinical Events', category: 'Associated Persons'},
  apcm: {description: 'Associated Persons Concomitant Medications', category: 'Associated Persons'},
  apeg: {description: 'Associated Persons ECG', category: 'Associated Persons'},
  aplb: {description: 'Associated Persons Laboratory Tests', category: 'Associated Persons'},
  apmh: {description: 'Associated Persons Medical History', category: 'Associated Persons'},
  apvs: {description: 'Associated Persons Vital Signs', category: 'Associated Persons'},
  // Medical History
  mh: {description: 'Medical History', category: 'Medical History'},
  // Tumor / Disease Response
  rs: {description: 'Disease Response / Clinical Classification', category: 'Tumor / Disease Response'},
  tr: {description: 'Tumor/Lesion Results', category: 'Tumor / Disease Response'},
  tu: {description: 'Tumor/Lesion Identification', category: 'Tumor / Disease Response'},
};

export const SUPP_DOMAIN_CATEGORY = 'Supplemental';

export const DOMAINS_CATEGORIES_LIST = ['Core / Special Purpose', 'Interventions', 'Events', 'Findings',
  'Tumor / Disease Response', 'Medical History', 'Genetics / Genomics', 'Pharmacokinetics / Modeling',
  'Associated Persons', 'Trial Design', 'Supplemental', 'Other'];
