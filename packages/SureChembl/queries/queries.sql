--name: searchPatentBySubstructure
--connection: Vhlushchenko:SureCHEMBL
--input: string pattern {semType: Substructure}
--input: int maxMols = 10
with molecules as (select m, schembl_chem_id from rdk.mols mols where m@>@pattern::qmol limit @maxMols)
 select m as smiles, text as title, sd.id as doc_surechembl_id, lang as language, assign_applic, scpn as doc_id, published,
 string_agg(CASE WHEN sdc.field=1 THEN 'DESCRIPTION'
             WHEN sdc.field=2 THEN 'CLAIMS'
             WHEN sdc.field=3 THEN 'ABSTRACT'
             WHEN sdc.field=4 THEN 'TITLE'
             WHEN sdc.field=5 THEN 'IMAGES'
             WHEN sdc.field=6 THEN 'ATTACHMENTS'
             ELSE 'other'
             END, ',') as fields
 from molecules mols
 join schembl_document_chemistry sdc on sdc.schembl_chem_id = mols.schembl_chem_id
 join schembl_document sd on sdc.schembl_doc_id = sd.id
 join schembl_document_title sdt on sdt.schembl_doc_id = sd.id
 where sdc.frequency > 0
 group by m, text, sd.id, lang, assign_applic, scpn, published
--end
