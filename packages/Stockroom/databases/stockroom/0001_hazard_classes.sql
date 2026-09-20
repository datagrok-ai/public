-- GHS hazard classes (UNECE GHS Rev.10, parts 2-4). Applied once via package_db_ups after the
-- schema deploys; ids are md5-derived from the chapter number so re-running is a no-op.

INSERT INTO stockroom.hazard_class (id, code, name, hazard_group) VALUES
  (md5('stockroom.hazard_class:2.1')::uuid,  '2.1',  'Explosives', 'physical'),
  (md5('stockroom.hazard_class:2.2')::uuid,  '2.2',  'Flammable gases', 'physical'),
  (md5('stockroom.hazard_class:2.3')::uuid,  '2.3',  'Aerosols and chemicals under pressure', 'physical'),
  (md5('stockroom.hazard_class:2.4')::uuid,  '2.4',  'Oxidizing gases', 'physical'),
  (md5('stockroom.hazard_class:2.5')::uuid,  '2.5',  'Gases under pressure', 'physical'),
  (md5('stockroom.hazard_class:2.6')::uuid,  '2.6',  'Flammable liquids', 'physical'),
  (md5('stockroom.hazard_class:2.7')::uuid,  '2.7',  'Flammable solids', 'physical'),
  (md5('stockroom.hazard_class:2.8')::uuid,  '2.8',  'Self-reactive substances and mixtures', 'physical'),
  (md5('stockroom.hazard_class:2.9')::uuid,  '2.9',  'Pyrophoric liquids', 'physical'),
  (md5('stockroom.hazard_class:2.10')::uuid, '2.10', 'Pyrophoric solids', 'physical'),
  (md5('stockroom.hazard_class:2.11')::uuid, '2.11', 'Self-heating substances and mixtures', 'physical'),
  (md5('stockroom.hazard_class:2.12')::uuid, '2.12', 'Substances and mixtures which, in contact with water, emit flammable gases', 'physical'),
  (md5('stockroom.hazard_class:2.13')::uuid, '2.13', 'Oxidizing liquids', 'physical'),
  (md5('stockroom.hazard_class:2.14')::uuid, '2.14', 'Oxidizing solids', 'physical'),
  (md5('stockroom.hazard_class:2.15')::uuid, '2.15', 'Organic peroxides', 'physical'),
  (md5('stockroom.hazard_class:2.16')::uuid, '2.16', 'Corrosive to metals', 'physical'),
  (md5('stockroom.hazard_class:2.17')::uuid, '2.17', 'Desensitized explosives', 'physical'),
  (md5('stockroom.hazard_class:3.1')::uuid,  '3.1',  'Acute toxicity', 'health'),
  (md5('stockroom.hazard_class:3.2')::uuid,  '3.2',  'Skin corrosion/irritation', 'health'),
  (md5('stockroom.hazard_class:3.3')::uuid,  '3.3',  'Serious eye damage/eye irritation', 'health'),
  (md5('stockroom.hazard_class:3.4')::uuid,  '3.4',  'Respiratory or skin sensitization', 'health'),
  (md5('stockroom.hazard_class:3.5')::uuid,  '3.5',  'Germ cell mutagenicity', 'health'),
  (md5('stockroom.hazard_class:3.6')::uuid,  '3.6',  'Carcinogenicity', 'health'),
  (md5('stockroom.hazard_class:3.7')::uuid,  '3.7',  'Reproductive toxicity', 'health'),
  (md5('stockroom.hazard_class:3.8')::uuid,  '3.8',  'Specific target organ toxicity - single exposure', 'health'),
  (md5('stockroom.hazard_class:3.9')::uuid,  '3.9',  'Specific target organ toxicity - repeated exposure', 'health'),
  (md5('stockroom.hazard_class:3.10')::uuid, '3.10', 'Aspiration hazard', 'health'),
  (md5('stockroom.hazard_class:4.1')::uuid,  '4.1',  'Hazardous to the aquatic environment', 'environmental'),
  (md5('stockroom.hazard_class:4.2')::uuid,  '4.2',  'Hazardous to the ozone layer', 'environmental')
ON CONFLICT (id) DO NOTHING;
