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
-- input: int templateId { nullable: true }
-- input: string scope
-- input: string choices { nullable: true }
-- input: double min { nullable: true }
-- input: double max { nullable: true }
-- input: int originPlateId { nullable: true }
-- output: int propertyId
INSERT INTO plates.properties(name, type, template_id, scope, choices, min, max, origin_plate_id)
VALUES(@propertyName, @valueType, @templateId, @scope, @choices, @min, @max, @originPlateId)
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

-- name: getTemplateWellProperties
-- connection: Curves:plates
SELECT template_id, property_id
FROM plates.template_well_properties
-- end


-- name: getTemplatePlateProperties
-- connection: Curves:plates
SELECT template_id, property_id
FROM plates.template_plate_properties
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
-- input: string analysisName
-- input: string parameters
-- output: int runId
INSERT INTO plates.analysis_runs(plate_id, analysis_name, parameters)
VALUES (@plateId, @analysisName, CAST(@parameters AS jsonb))
RETURNING id;
-- end


-- name: saveCurveResult
-- connection: Curves:plates
-- input: int runId
-- input: string seriesName {nullable: true}
-- input: string curveJson
-- input: double ic50 {nullable: true}
-- input: double hillSlope {nullable: true}
-- input: double rSquared {nullable: true}
-- input: double minValue {nullable: true}
-- input: double maxValue {nullable: true}
-- input: double auc {nullable: true}
INSERT INTO plates.analysis_results_curves(
    analysis_run_id, series_name, curve_json, ic50, hill_slope, r_squared, min_value, max_value, auc
) VALUES (
    @runId, @seriesName, CAST(@curveJson AS jsonb), @ic50, @hillSlope, @rSquared, @minValue, @maxValue, @auc
);
-- end
