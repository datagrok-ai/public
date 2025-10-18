-- name: getPlates
-- connection: Curves:plates
SELECT * FROM plates.plates;
-- end


-- name: getWellLevelProperties
-- connection: Curves:plates
SELECT DISTINCT p.*
FROM plates.properties p
WHERE p.scope = 'well'
ORDER BY name;
-- end


-- name: getPlateLevelProperties
-- connection: Curves:plates
SELECT DISTINCT p.*
FROM plates.properties p
WHERE p.scope = 'plate'
ORDER BY name;
-- end

-- name: getPropertyNames
-- connection: Curves:plates
SELECT name FROM plates.properties;
-- end


-- name: getPlateTypes
-- connection: Curves:plates
SELECT * FROM plates.plate_types;
-- end


-- name: getPlateTemplates
-- connection: Curves:plates
SELECT * FROM plates.templates;
-- end

-- name: getProperties
-- connection: Curves:plates
-- output: dataframe result
SELECT * FROM plates.properties;
-- end


-- name: getWellRoles
-- connection: Curves:plates
SELECT pav.id, pav.value_string AS name
FROM plates.property_allowed_values pav
JOIN plates.properties p ON pav.property_id = p.id
WHERE p.name = 'Well Role';
-- end


-- name: getWellValuesByBarcode
-- connection: Curves:plates
-- input: string barcode
select v.row, v.col, v.value_num, v.value_string, v.value_bool, v.property_id from plates.plate_well_values v
join plates.plates p on v.plate_id = p.id
where p.barcode = @barcode
order by property_id, row, col
-- end


-- name: getWellValuesById
-- connection: Curves:plates
-- input: int id
select v.row, v.col, v.value_num, v.value_string, v.value_bool, v.property_id from plates.plate_well_values v
join plates.plates p on v.plate_id = p.id
where p.id = @id
order by property_id, row, col
-- end


-- name: getAllowedValues
-- connection: Curves:plates
-- input: string propertyName { choices: getPropertyNames() }
SELECT pav.id, pav.value_string AS name
FROM plates.property_allowed_values pav
JOIN plates.properties p ON pav.property_id = p.id
WHERE p.name = @propertyName;
-- end


-- name: getUniquePlatePropertyValues
-- connection: Curves:plates
SELECT DISTINCT p.name, pd.value_string
FROM plates.plate_details pd
JOIN plates.properties p ON pd.property_id = p.id
WHERE p.type = 'string';
-- end


-- name: getUniqueWellPropertyValues
-- connection: Curves:plates
SELECT DISTINCT p.name, pwv.value_string
FROM plates.plate_well_values pwv
JOIN plates.properties p ON pwv.property_id = p.id
WHERE p.type = 'string';
-- end


-- name: createProperty
-- connection: Curves:plates
-- input: string propertyName
-- input: string valueType
-- input: string scope
-- input: string choices { nullable: true }
-- input: double min { nullable: true }
-- input: double max { nullable: true }
-- output: int propertyId
INSERT INTO plates.properties(name, type, scope, choices, min, max)
VALUES(@propertyName, @valueType, @scope, @choices, @min, @max)
RETURNING id;
-- end


-- name: createTemplate
-- connection: Curves:plates
-- input: string name
-- input: string description
-- output: int templateId
INSERT INTO plates.templates(name, description)
VALUES(@name, @description)
RETURNING id;
-- end

-- name: getTemplateProperties
-- connection: Curves:plates
-- input: int templateId
-- output: dataframe result
SELECT * FROM plates.properties
WHERE template_id = @templateId;
-- end


-- name: getPlateByBarcode
-- connection: Curves:Plates
-- input: string barcode
-- output: dataframe result
SELECT * FROM plates.plates WHERE barcode = @barcode;
-- end


-- name: createAnalysisRun
-- connection: Curves:plates
-- input: int plateId
-- input: string analysisType
-- input: list<string> groups
-- output: int runId
INSERT INTO plates.analysis_runs(plate_id, analysis_type, groups)
VALUES (@plateId, @analysisType, @groups)
RETURNING id;
-- end


-- name: saveAnalysisRunParameter
-- connection: Curves:plates
-- input: int analysisRunId
-- input: int propertyId
-- input: string valueString {nullable: true}
-- input: double valueNum {nullable: true}
-- input: bool valueBool {nullable: true}
-- input: string valueJsonb {nullable: true}
INSERT INTO plates.analysis_run_parameters(
    analysis_run_id, property_id, value_string, value_num, value_bool, value_jsonb
) VALUES (
    @analysisRunId, @propertyId, @valueString, @valueNum, @valueBool, CAST(@valueJsonb AS jsonb)
);
-- end


-- name: saveAnalysisResult
-- connection: Curves:plates
-- input: int analysisRunId
-- input: list<string> groupCombination
-- input: int propertyId
-- input: string valueString {nullable: true}
-- input: double valueNum {nullable: true}
-- input: bool valueBool {nullable: true}
-- input: string valueJsonb {nullable: true}
INSERT INTO plates.analysis_results(
    analysis_run_id, group_combination, property_id,
    value_string, value_num, value_bool, value_jsonb
) VALUES (
    @analysisRunId, @groupCombination, @propertyId,
    @valueString, @valueNum, @valueBool, CAST(@valueJsonb AS jsonb)
);
-- end


-- name: getAnalysisRunGroups
-- connection: Curves:plates
-- input: string analysisType
-- output: dataframe result
SELECT DISTINCT unnest(ar.groups) AS "group"
FROM plates.analysis_runs ar
WHERE ar.analysis_type = @analysisType
ORDER BY "group";
-- end


-- name: queryAnalyses
-- connection: Curves:plates
-- input: string fullQuery
-- output: dataframe result
@fullQuery
-- end

-- name: queryAnalysesTemplate
SELECT
    ar.id as run_id,
    p.id as plate_id,
    p.barcode,
    res_pivot.group_combination,
    res_pivot.properties
FROM
    plates.analysis_runs ar
JOIN
    plates.plates p ON ar.plate_id = p.id
CROSS JOIN LATERAL (
    SELECT
        res.group_combination,
        jsonb_object_agg(
            prop.name,
            CASE
                WHEN res.value_string IS NOT NULL THEN to_jsonb(res.value_string)
                WHEN res.value_num IS NOT NULL THEN to_jsonb(res.value_num)
                WHEN res.value_bool IS NOT NULL THEN to_jsonb(res.value_bool)
                WHEN res.value_jsonb IS NOT NULL THEN res.value_jsonb
                ELSE 'null'::jsonb
            END
        )::text as properties
    FROM
        plates.analysis_results res
    JOIN
        plates.properties prop ON res.property_id = prop.id
    WHERE
        res.analysis_run_id = ar.id
    GROUP BY
        res.group_combination
) as res_pivot
WHERE
    ar.analysis_type = '${analysisType}'
    ${finalWhereClause}
-- end



-- name: addTemplatePlateProperty
-- connection: Curves:plates
-- input: int templateId
-- input: int propertyId
-- input: bool isRequired
-- input: string defaultValue { nullable: true }
INSERT INTO plates.template_plate_properties (template_id, property_id, is_required, default_value)
VALUES (@templateId, @propertyId, @isRequired, @defaultValue)
ON CONFLICT (template_id, property_id) 
DO UPDATE SET is_required = EXCLUDED.is_required, default_value = EXCLUDED.default_value;
-- end

-- name: addTemplateWellProperty
-- connection: Curves:plates
-- input: int templateId
-- input: int propertyId
-- input: bool isRequired
-- input: string defaultValue { nullable: true }
INSERT INTO plates.template_well_properties (template_id, property_id, is_required, default_value)
VALUES (@templateId, @propertyId, @isRequired, @defaultValue)
ON CONFLICT (template_id, property_id) 
DO UPDATE SET is_required = EXCLUDED.is_required, default_value = EXCLUDED.default_value;
-- end

-- name: getTemplatePlateProperties
-- connection: Curves:plates
-- input: int templateId
SELECT p.*, tpp.is_required, tpp.default_value
FROM plates.template_plate_properties tpp
JOIN plates.properties p ON tpp.property_id = p.id
WHERE tpp.template_id = @templateId
ORDER BY p.name;
-- end

-- name: getTemplateWellProperties
-- connection: Curves:plates
-- input: int templateId
SELECT p.*, twp.is_required, twp.default_value
FROM plates.template_well_properties twp
JOIN plates.properties p ON twp.property_id = p.id
WHERE twp.template_id = @templateId
ORDER BY p.name;
-- end

-- name: removeTemplatePlateProperty
-- connection: Curves:plates
-- input: int templateId
-- input: int propertyId
DELETE FROM plates.template_plate_properties 
WHERE template_id = @templateId AND property_id = @propertyId;
-- end

-- name: removeTemplateWellProperty
-- connection: Curves:plates
-- input: int templateId
-- input: int propertyId
DELETE FROM plates.template_well_properties 
WHERE template_id = @templateId AND property_id = @propertyId;
-- end