--name: _PK data from 'Curated Drug Pharmacokinetic Data' source for @drug
--friendlyName: Browse | PK for @drug
--connection: Chembl
--meta.searchPattern: "Pharmacokinetic Data for ${drug}" 
--input: string drug = "LEVOFLOXACIN"
SELECT DISTINCT
  d.title,
  min(case when ap.standard_type = 'DATASET' then coalesce(ap.standard_value::text, ap.standard_text_value) else null end)         dataset,
  a.assay_id,
  a.description,
  min(case when actp.standard_type = 'DOSED_COMPOUND_NAME' then coalesce(actp.standard_value::text, actp.standard_text_value) || ' ' || actp.standard_units  else null end)  dosed_compound_name,
  min(case when actp.standard_type = 'DOSE'                then coalesce(actp.standard_value::text, actp.standard_text_value) || ' ' || actp.standard_units  else null end)  dose,
  min(case when actp.standard_type = 'DOSAGE_FORM'         then coalesce(actp.standard_value::text, actp.standard_text_value) || ' ' || actp.standard_units  else null end)  dosage_form,
  min(case when actp.standard_type = 'REGIMEN'             then coalesce(actp.standard_value::text, actp.standard_text_value) || ' ' || actp.standard_units  else null end)  regimen,
  min(case when actp.standard_type = 'ROUTE'               then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) route,
  min(case when actp.standard_type = 'GENDER'              then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) gender,
  min(case when actp.standard_type = 'AGE_RANGE'           then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) age_range,
  min(case when actp.standard_type = 'HEALTH_STATUS'       then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) health_status,
  min(case when actp.standard_type = 'TISSUE'              then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) tissue,
  cr.molregno,
  cr.compound_name,
  act.activity_id,
  act.toid,
  act.standard_type,
  act.standard_relation,
  act.standard_value,
  act.standard_units,
  act.activity_comment
FROM source s
  JOIN compound_records cr ON s.src_id = cr.src_id
  JOIN docs d ON d.doc_id = cr.doc_id
  JOIN activities act ON cr.record_id = act.record_id AND cr.doc_id = act.doc_id
  JOIN activity_properties actp ON act.activity_id = actp.activity_id
  JOIN assays a ON act.assay_id = a.assay_id
  JOIN assay_parameters ap ON a.assay_id = ap.assay_id
                              AND s.src_description = 'Curated Drug Pharmacokinetic Data'
                              AND cr.compound_name ILIKE @drug
GROUP BY d.title, a.assay_id, a.description, cr.molregno, cr.compound_name, act.activity_id, act.toid,
  act.standard_type, act.standard_relation, act.standard_value, act.standard_units, act.activity_comment
ORDER BY cr.compound_name, act.toid, act.standard_type;
--end