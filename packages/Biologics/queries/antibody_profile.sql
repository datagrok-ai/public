--name: antibodyProfileByOrganismAndTarget
--connection: Biologics:biologics
--input: string organism {choices: query("SELECT DISTINCT name FROM biologics.target_organisms")}
--input: string target {choices: query("SELECT DISTINCT name FROM biologics.targets")}
--description: "Get comprehensive antibody profile filtered by organism and target, including structure, HC/LC regions, and all assay data."
--meta.searchPattern: "biologics antibody profile for ${organism} targeting ${target}"
SELECT
    adc.glyph,
    adc.identifier AS adc_identifier,
    adc.name AS adc_name,
    seq.identifier AS sequence_identifier,
    seq.sequence,
    SUBSTRING(seq.sequence FROM hc.start_pos FOR (hc.end_pos - hc.start_pos + 1)) AS heavy_chain,
    SUBSTRING(seq.sequence FROM lc.start_pos FOR (lc.end_pos - lc.start_pos + 1)) AS light_chain,
    at.name AS assay_type,
    at.category AS assay_category,
    at.fit_model,
    ar.result_value,
    ar.units,
    tgt.name AS target_name,
    org.name AS organism_name,
    d.name AS drug_name,
    d.smiles AS drug_smiles,
    l.linker_type
FROM biologics.assay_results ar
JOIN biologics.assay_types at ON ar.assay_id = at.id
JOIN biologics.target_organisms org ON ar.target_organism_id = org.id
JOIN biologics.targets tgt ON ar.target_id = tgt.id
JOIN biologics.adc adc ON ar.adc_id = adc.id
JOIN biologics.sequences seq ON adc.antibody_id = seq.id
LEFT JOIN biologics.drugs d ON adc.drug_id = d.id
LEFT JOIN biologics.linkers l ON adc.linker_id = l.id
LEFT JOIN biologics.sequence_regions hc ON hc.sequence_id = seq.id AND hc.region_type = 'HC'
LEFT JOIN biologics.sequence_regions lc ON lc.sequence_id = seq.id AND lc.region_type = 'LC'
WHERE org.name ILIKE @organism
  AND tgt.name ILIKE @target
ORDER BY adc.identifier, at.category, at.name;
-- end
