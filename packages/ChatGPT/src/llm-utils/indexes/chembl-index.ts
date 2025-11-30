import {DBConnectionMeta} from '../db-index-tools';

/* eslint-disable max-len */
export const chemblIndex: DBConnectionMeta = {
  'name': 'Chembl',
  'schemas': [
    {
      'name': 'public',
      'tables': [
        {
          'name': 'action_type',
          'schema': 'public',
          'columns': [
            {
              'name': 'action_type',
              'table': 'action_type',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'AGONIST',
                'POSITIVE MODULATOR',
                'BLOCKER',
                'NEGATIVE MODULATOR',
                'ACTIVATOR',
                'OTHER',
                'BINDING AGENT',
                'OXIDATIVE ENZYME',
                'ANTISENSE INHIBITOR',
                'EXOGENOUS PROTEIN',
                'POSITIVE ALLOSTERIC MODULATOR',
                'INVERSE AGONIST',
                'PARTIAL AGONIST',
                'EXOGENOUS GENE',
                'CHELATING AGENT',
                'NEGATIVE ALLOSTERIC MODULATOR',
                'STABILISER',
                'REDUCING AGENT',
                'VACCINE ANTIGEN',
                'SEQUESTERING AGENT',
                'DISRUPTING AGENT',
                'DEGRADER',
                'CROSS-LINKING AGENT',
                'ALLOSTERIC ANTAGONIST',
                'INHIBITOR',
                'OPENER',
                'SUBSTRATE',
                'MODULATOR',
                'METHYLATING AGENT',
                'PROTEOLYTIC ENZYME',
                'HYDROLYTIC ENZYME',
                'RELEASING AGENT',
                'RNAI INHIBITOR',
                'ANTAGONIST'
              ],
              'isUnique': true,
              'comment': 'Type of action, indicating the pharmacological effect.',
              'LLMComment': 'Represents the pharmacological action type, such as AGONIST or BLOCKER.'
            },
            {
              'name': 'description',
              'table': 'action_type',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'semanticType': 'Text',
              'comment': 'Detailed description of the action type.',
              'LLMComment': 'A unique textual description providing more context about the action type.'
            },
            {
              'name': 'parent_type',
              'table': 'action_type',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'POSITIVE MODULATOR',
                'NEGATIVE MODULATOR',
                'OTHER'
              ],
              'isUnique': false,
              'comment': 'Indicates the broader category of the action type, such as POSITIVE MODULATOR or NEGATIVE MODULATOR.',
              'LLMComment': 'Categorizes the action type under a parent type, helping to understand its relationship with other action types.'
            }
          ],
          'rowCount': 34,
          'comment': 'This table defines various types of actions that can be associated with activities in the database, including their descriptions and hierarchical relationships.',
          'LLMComment': 'The \'public.action_type\' table serves as a foundational reference for categorizing different actions related to activities within the database. It contains 34 distinct action types, each with a descriptive text and a potential parent type, allowing for a structured hierarchy of actions. This is crucial for understanding how different activities are classified and related to one another in the context of biological and chemical research, facilitating better data organization and retrieval.'
        },
        {
          'name': 'activities',
          'schema': 'public',
          'columns': [
            {
              'name': 'activity_id',
              'table': 'activities',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Unique identifier for each activity.'
            },
            {
              'name': 'assay_id',
              'table': 'activities',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the assay associated with the activity.'
            },
            {
              'name': 'doc_id',
              'table': 'activities',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the document that contains the activity data.'
            },
            {
              'name': 'record_id',
              'table': 'activities',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Unique identifier for the record in the database.'
            },
            {
              'name': 'molregno',
              'table': 'activities',
              'schema': 'public',
              'type': 'bigint',
              'semanticType': 'molregno',
              'comment': 'Molecular registration number, used to identify chemical substances.',
              'LLMComment': 'A unique identifier for a chemical compound in the database.'
            },
            {
              'name': 'standard_relation',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Describes the relationship of the standard value to the measured value.',
              'LLMComment': 'Indicates how the standard value relates to the activity measurement.'
            },
            {
              'name': 'standard_value',
              'table': 'activities',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'The standard value associated with the activity measurement.',
              'LLMComment': 'Numerical representation of the standard measurement.'
            },
            {
              'name': 'standard_units',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Units of measurement for the standard value.',
              'LLMComment': 'Specifies the units in which the standard value is expressed.'
            },
            {
              'name': 'standard_flag',
              'table': 'activities',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Flag indicating the status of the standard value (e.g., validated, provisional).',
              'LLMComment': 'Indicates the reliability or validation status of the standard value.'
            },
            {
              'name': 'standard_type',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of standard value (e.g., mean, median).',
              'LLMComment': 'Describes the statistical type of the standard value.'
            },
            {
              'name': 'activity_comment',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Comments or notes related to the activity.',
              'LLMComment': 'Additional information or observations regarding the activity.'
            },
            {
              'name': 'data_validity_comment',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Comments on the validity of the data associated with the activity.',
              'LLMComment': 'Notes on the reliability or quality of the data.'
            },
            {
              'name': 'potential_duplicate',
              'table': 'activities',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the activity may be a duplicate of another entry.',
              'LLMComment': 'Flag to identify possible duplicate activities.'
            },
            {
              'name': 'pchembl_value',
              'table': 'activities',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Calculated pChEMBL value for the activity.',
              'LLMComment': 'A logarithmic measure of the potency of a compound.'
            },
            {
              'name': 'bao_endpoint',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Biological Activity Ontology endpoint related to the activity.',
              'LLMComment': 'Specifies the biological endpoint for the activity measurement.'
            },
            {
              'name': 'uo_units',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Units of measurement for the unstandardized value.',
              'LLMComment': 'Specifies the units for the unstandardized activity measurement.'
            },
            {
              'name': 'qudt_units',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'semanticType': 'URL',
              'comment': 'Units of measurement defined by the QUDT standard.',
              'LLMComment': 'Units based on the Quantities, Units, Dimensions and Types (QUDT) framework.'
            },
            {
              'name': 'toid',
              'table': 'activities',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Identifier for the target of the activity.',
              'LLMComment': 'Unique identifier for the biological target associated with the activity.'
            },
            {
              'name': 'upper_value',
              'table': 'activities',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Upper limit of the activity measurement range.',
              'LLMComment': 'Specifies the maximum value for the activity measurement.'
            },
            {
              'name': 'standard_upper_value',
              'table': 'activities',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Upper limit of the standard value range.',
              'LLMComment': 'Indicates the maximum standard value for comparison.'
            },
            {
              'name': 'src_id',
              'table': 'activities',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Identifier for the source of the activity data.',
              'LLMComment': 'Indicates the origin or source of the activity information.'
            },
            {
              'name': 'type',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of activity (e.g., inhibition, activation).',
              'LLMComment': 'Categorizes the nature of the biological activity.'
            },
            {
              'name': 'relation',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Describes the relationship between the activity and other measurements.',
              'LLMComment': 'Indicates how this activity relates to other data points.'
            },
            {
              'name': 'value',
              'table': 'activities',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Measured value of the activity.',
              'LLMComment': 'Numerical representation of the observed activity.'
            },
            {
              'name': 'units',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Units of measurement for the activity value.',
              'LLMComment': 'Specifies the units in which the activity value is expressed.'
            },
            {
              'name': 'text_value',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Textual representation of the activity value.',
              'LLMComment': 'Describes the activity value in a non-numeric format.'
            },
            {
              'name': 'standard_text_value',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Textual representation of the standard value.',
              'LLMComment': 'Describes the standard value in a non-numeric format.'
            },
            {
              'name': 'action_type',
              'table': 'activities',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of action associated with the activity (e.g., experimental, predicted).',
              'LLMComment': 'Indicates whether the activity was experimentally measured or predicted.'
            }
          ],
          'rowCount': 20772701,
          'comment': 'This table stores detailed information about various biological activities associated with assays, including measurements, units, and comments related to the activities.',
          'LLMComment': 'The \'public.activities\' table is a crucial component of a biological database, containing over 20 million records of activities linked to assays. Each record includes identifiers for the activity and assay, measurement values, units, and various flags and comments that provide context about the validity and nature of the data. This table is essential for researchers analyzing biological interactions, as it aggregates quantitative data that can be used to assess the efficacy and characteristics of compounds in drug discovery and development.'
        },
        {
          'name': 'activity_properties',
          'schema': 'public',
          'columns': [
            {
              'name': 'ap_id',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'activity_id',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'type',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'The category or classification of the activity property, indicating its nature or purpose.',
              'LLMComment': 'Describes the specific type of property associated with the activity, such as \'duration\', \'intensity\', etc.'
            },
            {
              'name': 'relation',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'The relationship of this property to other properties or entities, which may define how it interacts or correlates with them.',
              'LLMComment': 'Indicates how this property relates to other data points or standards, such as \'greater than\', \'less than\', etc.'
            },
            {
              'name': 'value',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'The numerical value representing the measurement or quantity of the activity property.',
              'LLMComment': 'Holds the primary numeric measurement for the activity property, such as a count, duration, or score.'
            },
            {
              'name': 'units',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'The units of measurement for the value, specifying the scale or dimension (e.g., seconds, meters, etc.).',
              'LLMComment': 'Defines the measurement units for the \'value\' column, providing context for the numeric data.'
            },
            {
              'name': 'text_value',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'A textual representation of the value, which may provide additional context or clarification.',
              'LLMComment': 'Contains a string version of the value, useful for human-readable formats or descriptions.'
            },
            {
              'name': 'standard_type',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'The standardized category or classification for the property, used for comparison or benchmarking.',
              'LLMComment': 'Indicates the standard classification of the property, facilitating consistency across datasets.'
            },
            {
              'name': 'standard_relation',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'The standardized relationship of this property, which may be used for comparative analysis.',
              'LLMComment': 'Describes the standard relational context for the property, aiding in data normalization.'
            },
            {
              'name': 'standard_value',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'The standardized numeric value for the property, allowing for uniform comparisons across different datasets.',
              'LLMComment': 'Represents the normalized measurement of the property, useful for standardization purposes.'
            },
            {
              'name': 'standard_units',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'The standardized units of measurement for the standard value, ensuring consistency in interpretation.',
              'LLMComment': 'Specifies the units for the \'standard_value\', ensuring clarity in data comparisons.'
            },
            {
              'name': 'standard_text_value',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'A textual representation of the standard value, providing additional context or clarification.',
              'LLMComment': 'Contains a string version of the \'standard_value\', enhancing readability and understanding.'
            },
            {
              'name': 'comments',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Additional notes or remarks regarding the activity property, which may provide context or insights.',
              'LLMComment': 'Allows for supplementary information about the property, useful for documentation or clarification.'
            },
            {
              'name': 'result_flag',
              'table': 'activity_properties',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'A flag indicating the result status of the activity property, which may denote success, failure, or other states.',
              'LLMComment': 'Indicates the outcome status of the property, often used for validation or error checking.'
            }
          ],
          'rowCount': 10314826,
          'comment': 'This table stores properties associated with various activities, including their types, values, and relationships to standards.',
          'LLMComment': 'The `public.activity_properties` table is a crucial component in the database that captures detailed attributes of activities within a biological or chemical research context. Each record links an activity (identified by `activity_id`) to specific properties such as type, value, and units, which can be numeric or textual. This table supports the analysis of activities by providing standardized comparisons through fields like `standard_type` and `standard_value`. With over 10 million entries, it plays a significant role in enabling researchers to query and analyze the characteristics of various activities, facilitating insights into biological assays, drug mechanisms, and other scientific investigations.'
        },
        {
          'name': 'activity_smid',
          'schema': 'public',
          'columns': [
            {
              'name': 'smid',
              'table': 'activity_smid',
              'schema': 'public',
              'type': 'bigint'
            }
          ],
          'rowCount': 1732478,
          'comment': 'This table stores unique identifiers for activities, represented by the \'smid\' column, which is a bigint type.',
          'LLMComment': 'The \'public.activity_smid\' table serves as a key reference point within the broader database schema, containing a large number of unique activity identifiers (smid). These identifiers are crucial for linking various activity-related data across multiple tables, such as \'activities\', \'action_type\', and \'activity_properties\'. This structure supports complex queries and data retrieval processes, enabling researchers and analysts to efficiently access and analyze activity-related information in the context of biotherapeutics, drug mechanisms, and biological assays.'
        },
        {
          'name': 'activity_stds_lookup',
          'schema': 'public',
          'columns': [
            {
              'name': 'std_act_id',
              'table': 'activity_stds_lookup',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 161,
              'isUnique': true,
              'comment': 'Standard Activity ID'
            },
            {
              'name': 'standard_type',
              'table': 'activity_stds_lookup',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Type of standard activity, indicating the category or classification.',
              'LLMComment': 'Describes the classification of the standard activity, which may include types like \'biological\', \'chemical\', etc.'
            },
            {
              'name': 'definition',
              'table': 'activity_stds_lookup',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Detailed description of the standard activity.',
              'LLMComment': 'Provides a comprehensive explanation of what the standard activity entails.'
            },
            {
              'name': 'standard_units',
              'table': 'activity_stds_lookup',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                '-',
                'ug.g-1',
                '%',
                'ng.eq.g-1',
                'hr',
                'mL.min-1.kg-1',
                'g',
                'cells.uL-1',
                'uM.hr',
                'pg',
                'degrees C',
                'mg.kg-1',
                'L.kg-1',
                'mm',
                'nM',
                'M-1.s-1',
                'umol.kg-1',
                'mEq.L-1',
                'mL.min-1.g-1',
                'ug.mL-1',
                'ng.eq.mL-1',
                's-1',
                'fL',
                'ng.hr.g-1',
                's',
                'U.L-1',
                'uL.min-1.(10^6cells)-1',
                '/100WBC',
                'ng.hr.mL-1',
                'IU.L-1'
              ],
              'isUnique': false,
              'comment': 'Units of measurement for the standard activity, indicating how values are expressed.',
              'LLMComment': 'Specifies the measurement units used for the standard, such as micrograms per gram (ug.g-1) or percentage (%).'
            },
            {
              'name': 'normal_range_min',
              'table': 'activity_stds_lookup',
              'schema': 'public',
              'type': 'numeric',
              'min': -10,
              'max': 100,
              'isUnique': false,
              'comment': 'Minimum value of the normal range for the standard activity.',
              'LLMComment': 'Indicates the lower limit of acceptable values for the standard activity.'
            },
            {
              'name': 'normal_range_max',
              'table': 'activity_stds_lookup',
              'schema': 'public',
              'type': 'numeric',
              'min': 1,
              'max': 1000000000,
              'isUnique': false,
              'comment': 'Maximum value of the normal range for the standard activity.',
              'LLMComment': 'Indicates the upper limit of acceptable values for the standard activity.'
            }
          ],
          'rowCount': 150,
          'comment': 'This table contains a lookup for standard activities related to various assays, including their definitions, types, and normal ranges.',
          'LLMComment': 'The `activity_stds_lookup` table serves as a reference for standard activities in the context of biological assays. It includes key information such as the unique identifier for each standard activity (`std_act_id`), the type of standard (`standard_type`), a textual definition of the activity (`definition`), the units in which the standard is measured (`standard_units`), and the normal range for the activity values (`normal_range_min` and `normal_range_max`). This table is crucial for ensuring consistency and accuracy in the interpretation of assay results across various biological and chemical research applications.'
        },
        {
          'name': 'activity_supp',
          'schema': 'public',
          'columns': [
            {
              'name': 'as_id',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'rgid',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Reference Group ID, linking to a specific group of related activities.',
              'LLMComment': 'Identifier for the group of activities this record is associated with.'
            },
            {
              'name': 'smid',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Supplemental Measurement ID, indicating the specific measurement type.',
              'LLMComment': 'Identifier for the specific measurement related to this activity.'
            },
            {
              'name': 'type',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of activity or measurement being recorded.',
              'LLMComment': 'Describes the category or nature of the activity.'
            },
            {
              'name': 'relation',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Describes the relationship of this activity to other activities or measurements.',
              'LLMComment': 'Indicates how this activity is related to others in the dataset.'
            },
            {
              'name': 'value',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Numeric value representing the measurement or result of the activity.',
              'LLMComment': 'The quantitative result of the activity being recorded.'
            },
            {
              'name': 'units',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Units of measurement for the value (e.g., meters, seconds).',
              'LLMComment': 'Specifies the measurement units for the recorded value.'
            },
            {
              'name': 'text_value',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Textual representation of the value, if applicable.',
              'LLMComment': 'A string representation of the value for readability.'
            },
            {
              'name': 'standard_type',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of standard measurement for comparison purposes.',
              'LLMComment': 'Indicates the standard against which the activity is measured.'
            },
            {
              'name': 'standard_relation',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Relationship of the standard measurement to the activity.',
              'LLMComment': 'Describes how the standard measurement relates to the recorded activity.'
            },
            {
              'name': 'standard_value',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Numeric value of the standard measurement for comparison.',
              'LLMComment': 'The quantitative standard value for reference.'
            },
            {
              'name': 'standard_units',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Units of measurement for the standard value.',
              'LLMComment': 'Specifies the measurement units for the standard value.'
            },
            {
              'name': 'standard_text_value',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Textual representation of the standard value, if applicable.',
              'LLMComment': 'A string representation of the standard value for readability.'
            },
            {
              'name': 'comments',
              'table': 'activity_supp',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Additional comments or notes regarding the activity or measurement.',
              'LLMComment': 'Free-text field for any extra information related to the activity.'
            }
          ],
          'rowCount': 1776415,
          'comment': 'This table stores supplementary activity data related to various assays and biological activities, including measurements, types, and standard values.',
          'LLMComment': 'The `activity_supp` table is a crucial component in the biological and chemical research domain, serving as a repository for supplementary data associated with various activities and assays. It includes detailed information such as activity IDs, measurement types, values, and their corresponding units. This data is essential for researchers analyzing biological assays, as it provides context and additional metrics that enhance the understanding of the primary activity data. The table\'s extensive row count indicates its role in capturing a wide array of experimental results, making it a valuable resource for data analysis and interpretation in drug discovery and biotherapeutic development.'
        },
        {
          'name': 'activity_supp_map',
          'schema': 'public',
          'columns': [
            {
              'name': 'actsm_id',
              'table': 'activity_supp_map',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'activity_id',
              'table': 'activity_supp_map',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'smid',
              'table': 'activity_supp_map',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the supporting material associated with the activity.',
              'LLMComment': 'This column represents the unique identifier for the supporting material linked to a specific activity.'
            }
          ],
          'rowCount': 2010125,
          'comment': 'Mapping of activities to their corresponding supplementary materials.',
          'LLMComment': 'The `activity_supp_map` table serves as a crucial link between various activities and their associated supplementary materials in a biological or chemical research context. Each entry in this table connects a specific activity (identified by `activity_id`) to a supplementary material (identified by `smid`), allowing researchers to track and manage the resources that support or enhance the understanding of each activity. This mapping is essential for data integrity and retrieval in complex datasets, particularly in fields like pharmacology, biochemistry, and assay development, where understanding the context and resources related to specific activities is vital for analysis and research.'
        },
        {
          'name': 'assay_class_map',
          'schema': 'public',
          'columns': [
            {
              'name': 'ass_cls_map_id',
              'table': 'assay_class_map',
              'schema': 'public',
              'type': 'bigint',
              'min': 2309321,
              'max': 2526018,
              'isUnique': true
            },
            {
              'name': 'assay_id',
              'table': 'assay_class_map',
              'schema': 'public',
              'type': 'bigint',
              'min': 401,
              'max': 2304852,
              'isUnique': false,
              'comment': 'Identifier for the specific assay being mapped.',
              'LLMComment': 'Unique identifier for the assay associated with this mapping.'
            },
            {
              'name': 'assay_class_id',
              'table': 'assay_class_map',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 580,
              'isUnique': false,
              'comment': 'Identifier for the class to which the assay belongs.',
              'LLMComment': 'Unique identifier representing the classification of the assay.'
            }
          ],
          'rowCount': 216690,
          'comment': 'Mapping of assay classes to assays, linking each assay to its corresponding classification.',
          'LLMComment': 'The `public.assay_class_map` table serves as a crucial link between assays and their respective classifications within a biological or chemical research context. Each entry in this table associates a unique assay (identified by `assay_id`) with a specific assay class (identified by `assay_class_id`), allowing researchers to categorize and analyze assays based on their functional or structural characteristics. This mapping is essential for organizing large datasets in drug discovery and bioassay development, facilitating efficient data retrieval and analysis across related tables in the schema.'
        },
        {
          'name': 'assay_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'assay_class_id',
              'table': 'assay_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 584,
              'isUnique': true,
              'comment': 'Unique identifier for each assay classification.'
            },
            {
              'name': 'l1',
              'table': 'assay_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'MUSCULO-SKELETAL SYSTEM',
                'CARDIOVASCULAR SYSTEM',
                'BLOOD AND BLOOD FORMING ORGANS',
                'SYSTEMIC HORMONAL PREPARATIONS, EXCL. SEX HORMONES AND INSULINS',
                'ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS',
                'SENSORY ORGANS',
                'NERVOUS SYSTEM',
                'DERMATOLOGICALS',
                'ALIMENTARY TRACT AND METABOLISM',
                'GENITO URINARY SYSTEM AND SEX HORMONES',
                'ANTIPARASITIC PRODUCTS, INSECTICIDES AND REPELLENTS',
                'RESPIRATORY SYSTEM'
              ],
              'isUnique': false,
              'comment': 'Primary category of the assay classification, indicating the major biological system involved.',
              'LLMComment': 'This column categorizes the assay into a broad biological system, such as musculoskeletal or cardiovascular.'
            },
            {
              'name': 'l2',
              'table': 'assay_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Melanoma Oncology Models',
                'Anti-Obesity Activity',
                'Acne Models',
                'Calcium Uptake Inhibiting Activity',
                'Intestinal Function',
                'Xenograft Oncology Models',
                'Liver Function',
                'Scleroderma Models',
                'NMRI Methods in Psychoneuropharmacology',
                'Experimental Hypertension',
                'Neoplasm Oncology Models',
                'Anti-Epileptic Activity',
                'Skin Sensitization Testing',
                'Neuroleptic Activity',
                'Test on Salivary Glands',
                'Anti-Arthrotic Activity',
                'Experimental Dermatitis',
                'Induction of Experimental Atherosclerosis',
                'Measurement of Blood Glucose-Lowering and Antidiabetic Activity',
                'Influence on Hair Growth',
                'Pruritus Models',
                'Impaired Renal Function',
                'Pancreatic Function',
                'Endocrine Safety Pharmacology',
                'Anti-Arrhythmic Activity',
                'Metastatic Oncology Models',
                'Influence on Lower Urinary Tract',
                'Uricosuric and Hypo-Uricemic Activity',
                'Cardiovascular Analysis',
                'Genetically Diabetic Animals',
                'Local Anesthetic Activity',
                'Leukemia Oncology Models',
                'Sarcoma Oncology Models',
                'Experimental Obesity',
                'Trypanosomal Activity',
                'Protection against UV light',
                'Diuretic and Saluretic Activity',
                'Anti-Inflammatory Activity',
                'Measurement of Glucose Absorption',
                'General Anesthetics',
                'Intraocular Pressure',
                'Inhibition of Cholesterol Absorption',
                'Malarial Activity',
                'Psoriasis Models',
                'Central Analgesic Activity',
                'Coronary Drugs',
                'Tests for Anxiolytic Activity',
                'Inhibition of Cholesterol Biosynthesis',
                'Thyroid Hormones',
                'Assessment of Renal Function',
                'Ichthyosis Models',
                'Cardiovascular Safety Pharmacology',
                'Monitoring of Diabetic Late Complications',
                'Cardiac Hypertrophy and Insufficiency',
                'Parathyroid Hormone',
                'Peripheral Analgesic Activity',
                'Models of Thrombosis',
                'Gall Bladder Function',
                'Anterior Pituitary Hormones',
                'Learning and Memory',
                'Effects on Different Peptide Hormones',
                'Anti-Pyretic Activity',
                'Genetic Models of Hemostasis and Thrombosis',
                'Analgesic Activity',
                'Ovarian Hormones',
                'Erythropoietic Protoporphyria Models',
                'Pemphigus Models',
                'Lymphoma Oncology Models',
                'Genetically Obese Animals',
                'Neuromuscular Blocking Activity',
                'Hypnotic Activity',
                'Hemostasis Bleeding Models',
                'Effects on Air Ways',
                'Adrenal Steroid Hormones',
                'Biomechanics of Skin',
                'Cutaneous Microcirculation',
                'Effects on Behavior and Muscle Coordination',
                'Transepidermal Water Loss',
                'Effects on Blood Supply and on Arterial and Venous Tonus',
                'Methods for Testing Immunological Factors',
                'Testicular Steroid Hormones',
                'Anti-Depressant Activity',
                'Xeroderma Models',
                'Gastric Function',
                'Carcinoma Oncology Models',
                'Animal Models of Neurological Diseases',
                'Vitiligo Models',
                'Experimental Diabetes Mellitus',
                'Hypothalamic Hormones',
                'Effects on Tracheal Cells and Bronchial Mucus Secretion and Transport',
                'Anti-Parkinsonism Activity',
                'Posterior Pituitary Hormones',
                'Influence on Lipid Metabolism'
              ],
              'isUnique': false,
              'comment': 'Secondary category providing more specific context about the assay\'s focus or application.',
              'LLMComment': 'This column specifies the particular area of research or activity related to the assay, such as oncology or metabolic functions.'
            },
            {
              'name': 'l3',
              'table': 'assay_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Unique label for the assay classification, providing a distinct name for identification.'
            },
            {
              'name': 'class_type',
              'table': 'assay_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'In vivo efficacy'
              ],
              'isUnique': false,
              'comment': 'Type of classification indicating the nature of the assay, specifically its efficacy in vivo.',
              'LLMComment': 'This column describes the classification type, focusing on the effectiveness of the assay in living organisms.'
            },
            {
              'name': 'source',
              'table': 'assay_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'phenotype',
                'Hock_2016',
                'Vogel_2008',
                'Vogel_2013_SafetyPharmacology'
              ],
              'isUnique': false,
              'comment': 'Origin of the assay classification data, indicating the study or dataset it is derived from.',
              'LLMComment': 'This column identifies the source of the classification, which could be a specific phenotype or a published study.'
            }
          ],
          'rowCount': 584,
          'comment': 'This table categorizes various assay types used in biological and chemical research, providing a structured classification system for assays.',
          'LLMComment': 'The `public.assay_classification` table serves as a key reference for categorizing assays in the context of biological and chemical research. It contains 584 entries, each identified by a unique `assay_class_id`. The table includes multiple classification levels (`l1`, `l2`, `l3`) that allow for hierarchical organization of assay types, along with `class_type` to specify the nature of the assay and `source` to indicate the origin of the classification data. This structured classification is essential for researchers to efficiently retrieve and analyze assay-related data across various related tables in the schema, facilitating better understanding and utilization of assay information in drug discovery and development.'
        },
        {
          'name': 'assay_parameters',
          'schema': 'public',
          'columns': [
            {
              'name': 'assay_param_id',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'bigint',
              'min': 3710880,
              'max': 4103145,
              'isUnique': true,
              'comment': 'Unique identifier for the assay parameter.'
            },
            {
              'name': 'assay_id',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'bigint',
              'min': 138,
              'max': 2306542,
              'isUnique': false,
              'comment': 'Identifier for the associated assay.'
            },
            {
              'name': 'type',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The type of the assay parameter, indicating its nature or category.',
              'LLMComment': 'Describes the specific nature or classification of the assay parameter.'
            },
            {
              'name': 'relation',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                '>',
                '='
              ],
              'isUnique': false,
              'comment': 'The relational operator used to compare the parameter value.',
              'LLMComment': 'Indicates how the parameter value relates to a standard or reference value.'
            },
            {
              'name': 'value',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'numeric',
              'min': 0,
              'max': 12000000000,
              'isUnique': false,
              'comment': 'The numeric value of the assay parameter.',
              'LLMComment': 'Represents the measured or calculated value of the parameter.'
            },
            {
              'name': 'units',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'RU',
                'null',
                'Hz',
                'mg',
                'min',
                '%',
                'pM/ mg protein',
                'hr',
                'M-1 s-1',
                'gM',
                'mg/ kg',
                'mg / kg',
                'folds diluted',
                'mgKg-1',
                'mg.Kg-1',
                'C',
                'Week',
                'Day',
                'degree C',
                'minutes',
                'M',
                'Year',
                'mM',
                'mg.m-2',
                'nL',
                'mV',
                'nM',
                'uM',
                'ng.ml-1',
                'ms',
                'mg/kg',
                'day',
                'uL/min',
                'mg per kg',
                'weeks',
                'ug.mL-1',
                'mg kg-1',
                'mgkg',
                'h',
                'hour',
                'kDa',
                's-1',
                'g/mol',
                'µM',
                's',
                'ug.ml-1',
                'ug.eq.Kg-1'
              ],
              'isUnique': false,
              'comment': 'The units of measurement for the parameter value.',
              'LLMComment': 'Specifies the measurement units for the parameter value, such as mg or Hz.'
            },
            {
              'name': 'text_value',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'A textual representation of the parameter value, if applicable.',
              'LLMComment': 'Provides a descriptive text alternative for the numeric value.'
            },
            {
              'name': 'standard_type',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The type of the standard against which the parameter is measured.',
              'LLMComment': 'Indicates the classification of the standard value for comparison.'
            },
            {
              'name': 'standard_relation',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                '>',
                '='
              ],
              'isUnique': false,
              'comment': 'The relational operator used to compare the standard value.',
              'LLMComment': 'Indicates how the standard value relates to the parameter value.'
            },
            {
              'name': 'standard_value',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'numeric',
              'min': 0,
              'max': 12000000000,
              'isUnique': false,
              'comment': 'The numeric standard value for comparison with the parameter value.',
              'LLMComment': 'Represents the reference or standard value for the parameter.'
            },
            {
              'name': 'standard_units',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'RU',
                'null',
                'Hz',
                'mg',
                'min',
                '%',
                'pM/ mg protein',
                'hr',
                'M-1 s-1',
                'gM',
                'mg/ kg',
                'mg / kg',
                'mg.kg-1',
                'folds diluted',
                'mg.Kg-1',
                'C',
                'Week',
                'Day',
                'degree C',
                'Year',
                'mM',
                'mg.m-2',
                'nL',
                'mV',
                'nM',
                'uM',
                'ng.ml-1',
                'ms',
                'uL/min',
                'weeks',
                'ug.mL-1',
                'h',
                'hour',
                'kDa',
                's-1',
                'g/mol',
                'µM',
                's',
                'ug.ml-1',
                'ug.eq.Kg-1'
              ],
              'isUnique': false,
              'comment': 'The units of measurement for the standard value.',
              'LLMComment': 'Specifies the measurement units for the standard value.'
            },
            {
              'name': 'standard_text_value',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'A textual representation of the standard value, if applicable.',
              'LLMComment': 'Provides a descriptive text alternative for the standard numeric value.'
            },
            {
              'name': 'comments',
              'table': 'assay_parameters',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'http://www.bioassayontology.org/bao#BAO_0013016',
                'null',
                'http://www.bioassayontology.org/bao#BAO_0000363',
                'http://www.bioassayontology.org/bao#BAO_0000444',
                'http://www.bioassayontology.org/bao#BAO_0000450',
                'http://www.bioassayontology.org/bao#BAO_0001098',
                'Gut Microbiota Medium (GMM) used for bacterial assays at pH7.',
                'Is the measured interaction considered due to direct binding to target?',
                'http://www.bioassayontology.org/bao#BAO_0010001',
                'http://www.bioassayontology.org/bao#BAO_0000572',
                'http://www.bioassayontology.org/bao#BAO_0000127',
                'http://www.bioassayontology.org/bao#BAO_0000124',
                'http://www.bioassayontology.org/bao#BAO_0010056',
                'http://www.bioassayontology.org/bao#BAO_0000045',
                'http://purl.obolibrary.org/obo/BTO_0000759',
                'Tested drug concentrations in the bacterial assays.',
                'http://www.bioassayontology.org/bao#BAO_0002648',
                'http://www.bioassayontology.org/bao#BAO_0000217',
                'http://www.bioassayontology.org/bao#BAO_0000010',
                'http://www.bioassayontology.org/bao#BAO_0000246',
                'http://www.bioassayontology.org/bao#BAO_0000574',
                'http://purl.obolibrary.org/obo/DOID_934',
                'http://www.bioassayontology.org/bao#BAO_0000535',
                'http://www.bioassayontology.org/bao#BAO_0000219',
                'http://www.bioassayontology.org/bao#BAO_0000030',
                'http://www.bioassayontology.org/bao#BAO_0020008',
                'http://www.bioassayontology.org/bao#BAO_0002764',
                'http://www.bioassayontology.org/bao#BAO_0000701',
                '+/- 0.0000033 (Mca-RPKPVE-Nval-WRK(Dnp)-NH2; ES002 from R&D Systems)'
              ],
              'isUnique': false,
              'comment': 'Additional comments or notes related to the assay parameter.',
              'LLMComment': 'Contains supplementary information or references regarding the assay parameter.'
            }
          ],
          'rowCount': 292750,
          'comment': 'Table containing parameters associated with various assays, including their values, units, and relationships to standards.',
          'LLMComment': 'The `public.assay_parameters` table is a critical component in the domain of bioassays and pharmacology. It stores detailed information about the parameters used in different assays, which are experiments conducted to measure the effects of substances on biological systems. Each entry in this table includes an identifier for the assay parameter, the associated assay ID, the type of parameter, its value, and the units of measurement. Additionally, it captures standard values and relationships, allowing for comparisons and standardization across different assays. This table is essential for researchers and developers in the pharmaceutical and biotechnology fields, as it provides the necessary data to analyze assay results and ensure consistency in experimental methodologies.'
        },
        {
          'name': 'assay_type',
          'schema': 'public',
          'columns': [
            {
              'name': 'assay_type',
              'table': 'assay_type',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'A',
                'B',
                'P',
                'U',
                'T',
                'F'
              ],
              'isUnique': true,
              'comment': 'Type of assay, categorized into various types such as A, B, P, U, T.',
              'LLMComment': 'Represents the classification of the assay, indicating its specific type, with possible values including A, B, P, U, T, etc.'
            },
            {
              'name': 'assay_desc',
              'table': 'assay_type',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Physicochemical',
                'Toxicity',
                'ADME',
                'Binding',
                'Functional',
                'Unassigned'
              ],
              'isUnique': true,
              'comment': 'Description of the assay type, detailing its specific focus such as Physicochemical, Toxicity, ADME, Binding, Functional.',
              'LLMComment': 'Provides a textual description of the assay type, outlining its purpose or area of study, with categories like Physicochemical, Toxicity, ADME, Binding, and Functional.'
            }
          ],
          'rowCount': 6,
          'comment': 'This table defines different types of assays used in biological and chemical research, along with their descriptions.',
          'LLMComment': 'The \'public.assay_type\' table is crucial in the context of biological and chemical research as it categorizes various assay types, which are experimental procedures used to measure the effects of substances on biological systems. Each entry includes a specific assay type and a detailed description, facilitating the understanding and selection of appropriate assays in research and development processes. This table interacts with other related tables, such as \'assays\' and \'assay_parameters\', to provide comprehensive data on experimental methodologies.'
        },
        {
          'name': 'assays',
          'schema': 'public',
          'columns': [
            {
              'name': 'assay_id',
              'table': 'assays',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'doc_id',
              'table': 'assays',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'description',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'semanticType': 'Text',
              'comment': 'Detailed description of the assay, providing context and specifics about the experimental setup or findings.',
              'LLMComment': 'A textual description that elaborates on the assay\'s purpose, methodology, and results.'
            },
            {
              'name': 'assay_type',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'The general classification of the assay, indicating the nature of the test conducted (e.g., biochemical, cellular).',
              'LLMComment': 'Categorizes the assay based on its primary function or methodology.'
            },
            {
              'name': 'assay_test_type',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Specifies the specific type of test performed within the assay (e.g., inhibition, activation).',
              'LLMComment': 'Indicates the precise nature of the test conducted in the assay.'
            },
            {
              'name': 'assay_category',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Broad classification of the assay, often used for grouping similar assays together.',
              'LLMComment': 'Helps in organizing assays into larger categories for easier retrieval and analysis.'
            },
            {
              'name': 'assay_organism',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'The organism from which the assay samples were derived, indicating the biological context of the assay.',
              'LLMComment': 'Identifies the species or model organism used in the assay.'
            },
            {
              'name': 'assay_tax_id',
              'table': 'assays',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Taxonomic identifier for the organism used in the assay, linking to biological classification databases.',
              'LLMComment': 'A unique identifier that corresponds to the organism\'s taxonomic classification.'
            },
            {
              'name': 'assay_strain',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Specific strain of the organism used in the assay, which may have unique genetic or phenotypic characteristics.',
              'LLMComment': 'Details the particular strain of the organism that was tested.'
            },
            {
              'name': 'assay_tissue',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of tissue from which samples were taken for the assay, relevant for understanding biological context.',
              'LLMComment': 'Specifies the tissue type involved in the assay.'
            },
            {
              'name': 'assay_cell_type',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of cells used in the assay, important for interpreting results in a cellular context.',
              'LLMComment': 'Indicates the specific cell type utilized in the assay.'
            },
            {
              'name': 'assay_subcellular_fraction',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Refers to the specific subcellular component from which the assay samples were derived, if applicable.',
              'LLMComment': 'Describes the subcellular location relevant to the assay.'
            },
            {
              'name': 'tid',
              'table': 'assays',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'relationship_type',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Describes the nature of the relationship between the assay and other entities or data points.',
              'LLMComment': 'Defines how this assay relates to other assays or data in the database.'
            },
            {
              'name': 'confidence_score',
              'table': 'assays',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'A numerical score indicating the confidence level in the assay results, often based on data quality or consistency.',
              'LLMComment': 'Quantifies the reliability of the assay findings.'
            },
            {
              'name': 'curated_by',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Identifier for the individual or system that curated the assay data, ensuring data integrity and traceability.',
              'LLMComment': 'Indicates who or what system is responsible for the data curation.'
            },
            {
              'name': 'src_id',
              'table': 'assays',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Identifier for the source of the assay data, linking to external databases or repositories.',
              'LLMComment': 'References the origin of the assay data.'
            },
            {
              'name': 'src_assay_id',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Identifier for the assay in the source database, allowing cross-referencing with external data.',
              'LLMComment': 'Links to the original assay ID in the source database.'
            },
            {
              'name': 'chembl_id',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'semanticType': 'CHEMBL_ID',
              'comment': 'Unique identifier for the assay in the ChEMBL database, facilitating integration with chemical biology data.',
              'LLMComment': 'A reference ID for the assay in the ChEMBL database.'
            },
            {
              'name': 'cell_id',
              'table': 'assays',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the specific cell line or primary cells used in the assay, important for reproducibility.',
              'LLMComment': 'Links to the specific cell type utilized in the assay.'
            },
            {
              'name': 'bao_format',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Format specification for the assay data, often used for standardization in reporting.',
              'LLMComment': 'Indicates the format in which the assay data is structured.'
            },
            {
              'name': 'tissue_id',
              'table': 'assays',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the specific tissue type used in the assay, linking to a tissue database.',
              'LLMComment': 'References the tissue type in a standardized database.'
            },
            {
              'name': 'variant_id',
              'table': 'assays',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for any genetic variants tested in the assay, relevant for genetic studies.',
              'LLMComment': 'Links to specific genetic variants associated with the assay.'
            },
            {
              'name': 'aidx',
              'table': 'assays',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Additional index or identifier for the assay, used for internal tracking or categorization.',
              'LLMComment': 'An auxiliary identifier for internal database management.'
            }
          ],
          'rowCount': 1644390,
          'comment': 'This table contains detailed information about various biological assays, including their types, organisms, and associated metadata.',
          'LLMComment': 'The \'public.assays\' table is a comprehensive repository of biological assay data, crucial for understanding experimental methods in drug discovery and biological research. It includes key identifiers like \'assay_id\' and \'doc_id\', along with descriptive fields such as \'assay_type\', \'assay_category\', and \'assay_organism\'. This table serves as a foundational element for linking assays to their respective activities and outcomes, facilitating insights into the efficacy and mechanisms of compounds in various biological contexts.'
        },
        {
          'name': 'atc_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'who_name',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Name of the World Health Organization (WHO) classification system.',
              'LLMComment': 'The name associated with the WHO classification for the ATC system.'
            },
            {
              'name': 'level1',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'B',
                'V',
                'J',
                'D',
                'N',
                'P',
                'H',
                'G',
                'S',
                'R',
                'A',
                'C',
                'M',
                'L'
              ],
              'isUnique': false,
              'comment': 'First level of the ATC classification, representing broad categories of drugs.',
              'LLMComment': 'Top-level category in the ATC classification, indicating major therapeutic areas.'
            },
            {
              'name': 'level2',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'M03',
                'L04',
                'L01',
                'R06',
                'P01',
                'M04',
                'B01',
                'R02',
                'J04',
                'S03',
                'D10',
                'D09',
                'L03',
                'A03',
                'A02',
                'C04',
                'V03',
                'G01',
                'V04',
                'A16',
                'N02',
                'M02',
                'D02',
                'D11',
                'D07',
                'D03',
                'D05',
                'C01',
                'A05',
                'M01',
                'V08',
                'S02',
                'B03',
                'P03',
                'N07',
                'A10',
                'G04',
                'J06',
                'N03',
                'C05',
                'B02',
                'V06',
                'C09',
                'A04',
                'G03',
                'C08',
                'B06',
                'L02',
                'V10',
                'H04',
                'A07',
                'V09',
                'J07',
                'N01',
                'V01',
                'D01',
                'R07',
                'M05',
                'S01',
                'G02',
                'D04',
                'H01',
                'N05',
                'J05',
                'H02',
                'A14',
                'C03',
                'C07',
                'M09',
                'H05',
                'A11',
                'R03',
                'R05',
                'D08',
                'P02',
                'A06',
                'N06',
                'C10',
                'B05',
                'C02',
                'A01',
                'A09',
                'A08',
                'H03',
                'J01',
                'N04',
                'J02',
                'R01',
                'A12',
                'D06'
              ],
              'isUnique': false,
              'comment': 'Second level of the ATC classification, providing more specific categories within level 1.',
              'LLMComment': 'Sub-category in the ATC classification that further specifies the drug\'s therapeutic use.'
            },
            {
              'name': 'level3',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Third level of the ATC classification, detailing specific pharmacological subcategories.',
              'LLMComment': 'More detailed classification of drugs, indicating specific pharmacological actions.'
            },
            {
              'name': 'level4',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Fourth level of the ATC classification, representing even more specific drug categories.',
              'LLMComment': 'Further refinement of drug classification, indicating specific drug types.'
            },
            {
              'name': 'level5',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Fifth level of the ATC classification, unique identifier for specific drug products.',
              'LLMComment': 'Most specific level in the ATC classification, uniquely identifying individual drug products.'
            },
            {
              'name': 'level1_description',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'BLOOD AND BLOOD FORMING ORGANS',
                'SYSTEMIC HORMONAL PREPARATIONS, EXCL. SEX HORMONES AND INSULINS',
                'ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS',
                'ANTIINFECTIVES FOR SYSTEMIC USE',
                'DERMATOLOGICALS',
                'MUSCULO-SKELETAL SYSTEM',
                'CARDIOVASCULAR SYSTEM',
                'VARIOUS',
                'SENSORY ORGANS',
                'NERVOUS SYSTEM',
                'ALIMENTARY TRACT AND METABOLISM',
                'GENITO URINARY SYSTEM AND SEX HORMONES',
                'ANTIPARASITIC PRODUCTS, INSECTICIDES AND REPELLENTS',
                'RESPIRATORY SYSTEM'
              ],
              'isUnique': false,
              'comment': 'Description of the first level category, providing context for the classification.',
              'LLMComment': 'Descriptive text explaining the major therapeutic area represented by level 1.'
            },
            {
              'name': 'level2_description',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'ANTIINFLAMMATORY AND ANTIRHEUMATIC PRODUCTS',
                'THYROID THERAPY',
                'PANCREATIC HORMONES',
                'DRUGS FOR ACID RELATED DISORDERS',
                'DIAGNOSTIC AGENTS',
                'UROLOGICALS',
                'ANTIEMETICS AND ANTINAUSEANTS',
                'DRUGS FOR TREATMENT OF BONE DISEASES',
                'VACCINES',
                'CORTICOSTEROIDS, DERMATOLOGICAL PREPARATIONS',
                'OTHER HEMATOLOGICAL AGENTS',
                'ANTIVIRALS FOR SYSTEMIC USE',
                'ANTIDIARRHEALS, INTESTINAL ANTIINFLAMMATORY/ANTIINFECTIVE AGENTS',
                'COUGH AND COLD PREPARATIONS',
                'ANESTHETICS',
                'GENERAL NUTRIENTS',
                'SEX HORMONES AND MODULATORS OF THE GENITAL SYSTEM',
                'LIPID MODIFYING AGENTS',
                'CONTRAST MEDIA',
                'PITUITARY AND HYPOTHALAMIC HORMONES AND ANALOGUES',
                'CORTICOSTEROIDS FOR SYSTEMIC USE',
                'DRUGS FOR OBSTRUCTIVE AIRWAY DISEASES',
                'ANTHELMINTICS',
                'AGENTS ACTING ON THE RENIN-ANGIOTENSIN SYSTEM',
                'OTHER RESPIRATORY SYSTEM PRODUCTS',
                'ANTIPSORIATICS',
                'ANTIHEMORRHAGICS',
                'PSYCHOLEPTICS',
                'ANTIBACTERIALS FOR SYSTEMIC USE',
                'THROAT PREPARATIONS',
                'TOPICAL PRODUCTS FOR JOINT AND MUSCULAR PAIN',
                'ECTOPARASITICIDES, INCL. SCABICIDES, INSECTICIDES AND REPELLENTS',
                'ANTIHISTAMINES FOR SYSTEMIC USE',
                'OPHTHALMOLOGICAL AND OTOLOGICAL PREPARATIONS',
                'MUSCLE RELAXANTS',
                'ANTIMYCOTICS FOR SYSTEMIC USE',
                'PSYCHOANALEPTICS',
                'DIGESTIVES, INCL. ENZYMES',
                'ANTIPROTOZOALS',
                'OTHER DERMATOLOGICAL PREPARATIONS',
                'VITAMINS',
                'ALL OTHER THERAPEUTIC PRODUCTS',
                'NASAL PREPARATIONS',
                'IMMUNOSTIMULANTS',
                'GYNECOLOGICAL ANTIINFECTIVES AND ANTISEPTICS',
                'OTHER ALIMENTARY TRACT AND METABOLISM PRODUCTS',
                'PREPARATIONS FOR TREATMENT OF WOUNDS AND ULCERS',
                'IMMUNE SERA AND IMMUNOGLOBULINS',
                'DRUGS FOR FUNCTIONAL GASTROINTESTINAL DISORDERS',
                'IMMUNOSUPPRESSANTS',
                'PERIPHERAL VASODILATORS',
                'OTHER DRUGS FOR DISORDERS OF THE MUSCULO-SKELETAL SYSTEM',
                'ANTISEPTICS AND DISINFECTANTS',
                'ANTIBIOTICS AND CHEMOTHERAPEUTICS FOR DERMATOLOGICAL USE',
                'ANTIGOUT PREPARATIONS',
                'ANTIEPILEPTICS',
                'ANTIANEMIC PREPARATIONS',
                'ANTIMYCOBACTERIALS',
                'ANTIHYPERTENSIVES',
                'ANTIFUNGALS FOR DERMATOLOGICAL USE',
                'OPHTHALMOLOGICALS',
                'CARDIAC THERAPY',
                'MINERAL SUPPLEMENTS',
                'CALCIUM HOMEOSTASIS',
                'ALLERGENS',
                'ANTINEOPLASTIC AGENTS',
                'EMOLLIENTS AND PROTECTIVES',
                'BETA BLOCKING AGENTS',
                'ANTIOBESITY PREPARATIONS, EXCL. DIET PRODUCTS',
                'OTHER NERVOUS SYSTEM DRUGS',
                'DRUGS FOR CONSTIPATION',
                'ANTIPRURITICS, INCL. ANTIHISTAMINES, ANESTHETICS, ETC.',
                'STOMATOLOGICAL PREPARATIONS',
                'ANTI-PARKINSON DRUGS',
                'ANALGESICS',
                'BILE AND LIVER THERAPY',
                'VASOPROTECTIVES',
                'OTHER GYNECOLOGICALS',
                'ANABOLIC AGENTS FOR SYSTEMIC USE',
                'ANTI-ACNE PREPARATIONS',
                'DIAGNOSTIC RADIOPHARMACEUTICALS',
                'ENDOCRINE THERAPY',
                'DIURETICS',
                'THERAPEUTIC RADIOPHARMACEUTICALS',
                'OTOLOGICALS',
                'DRUGS USED IN DIABETES',
                'BLOOD SUBSTITUTES AND PERFUSION SOLUTIONS',
                'MEDICATED DRESSINGS',
                'CALCIUM CHANNEL BLOCKERS',
                'ANTITHROMBOTIC AGENTS'
              ],
              'isUnique': false,
              'comment': 'Description of the second level category, offering details on specific drug types.',
              'LLMComment': 'Descriptive text for level 2, explaining the specific therapeutic uses of the drugs.'
            },
            {
              'name': 'level3_description',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Description of the third level category, detailing pharmacological actions.',
              'LLMComment': 'Text that describes the pharmacological actions associated with level 3 classifications.'
            },
            {
              'name': 'level4_description',
              'table': 'atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Description of the fourth level category, indicating specific drug types.',
              'LLMComment': 'Text that provides details on the specific types of drugs classified at level 4.'
            }
          ],
          'rowCount': 5369,
          'comment': 'This table contains the ATC (Anatomical Therapeutic Chemical) classification system, which categorizes drugs based on their therapeutic use and anatomical and chemical characteristics.',
          'LLMComment': 'The \'public.atc_classification\' table is a comprehensive resource for understanding the ATC classification system, which is essential for categorizing pharmaceuticals. It includes multiple levels of classification (from level 1 to level 5) that detail the anatomical and therapeutic properties of drugs. Each level is accompanied by descriptive fields that provide context and clarity about the classification hierarchy. This table is crucial for researchers, healthcare professionals, and AI systems that analyze drug classifications, therapeutic areas, and pharmacological research.'
        },
        {
          'name': 'binding_sites',
          'schema': 'public',
          'columns': [
            {
              'name': 'site_id',
              'table': 'binding_sites',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 38639,
              'isUnique': true
            },
            {
              'name': 'site_name',
              'table': 'binding_sites',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The name of the binding site, which may describe its function or location.',
              'LLMComment': 'A descriptive label for the binding site, often indicating its biological or chemical significance.'
            },
            {
              'name': 'tid',
              'table': 'binding_sites',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 120294,
              'isUnique': false,
              'comment': 'The identifier for the target associated with the binding site, linking it to a specific entity in the database.',
              'LLMComment': 'Target ID that connects the binding site to its corresponding biological target, such as a protein or molecule.'
            }
          ],
          'rowCount': 4549,
          'comment': 'This table contains information about binding sites, including unique identifiers and names associated with specific targets.',
          'LLMComment': 'The \'public.binding_sites\' table is a crucial component in the pharmacological and biochemical domain, specifically focusing on the interaction points where drugs or compounds bind to biological targets. Each entry in this table is identified by a unique \'site_id\' and is associated with a \'site_name\' that describes the binding site. The \'tid\' column links these binding sites to specific targets, facilitating the understanding of how different compounds interact with biological systems. This table is essential for researchers and developers in drug discovery and development, as it provides foundational data for analyzing binding interactions and optimizing therapeutic efficacy.'
        },
        {
          'name': 'bio_component_sequences',
          'schema': 'public',
          'columns': [
            {
              'name': 'component_id',
              'table': 'bio_component_sequences',
              'schema': 'public',
              'type': 'bigint',
              'min': 20525,
              'max': 22778,
              'isUnique': true,
              'comment': 'Unique identifier for each bio component.'
            },
            {
              'name': 'component_type',
              'table': 'bio_component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Nucleic acid',
                'Nucleic Acid',
                'NUCLEIC ACID',
                'Protein',
                'PROTEIN'
              ],
              'isUnique': false,
              'comment': 'Type of the biological component, indicating whether it is a nucleic acid or protein.',
              'LLMComment': 'Specifies the category of the bio component, which can be a nucleic acid or protein, with variations in case.'
            },
            {
              'name': 'description',
              'table': 'bio_component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'A textual description providing additional details about the bio component.',
              'LLMComment': 'Descriptive text that elaborates on the characteristics or functions of the bio component.'
            },
            {
              'name': 'sequence',
              'table': 'bio_component_sequences',
              'schema': 'public',
              'type': 'text',
              'isUnique': false,
              'comment': 'The actual sequence of nucleotides or amino acids for the bio component.',
              'LLMComment': 'Contains the biological sequence data, which can be a string of nucleotides for nucleic acids or amino acids for proteins.'
            },
            {
              'name': 'sequence_md5sum',
              'table': 'bio_component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'MD5 checksum of the sequence for integrity verification.',
              'LLMComment': 'A hash value used to verify the integrity of the sequence data, ensuring it has not been altered.'
            },
            {
              'name': 'tax_id',
              'table': 'bio_component_sequences',
              'schema': 'public',
              'type': 'bigint',
              'min': 4513,
              'max': 13688,
              'isUnique': false,
              'comment': 'Taxonomic identifier for the organism from which the bio component is derived.',
              'LLMComment': 'A unique identifier that corresponds to the taxonomic classification of the organism, linking to a broader biological taxonomy.'
            },
            {
              'name': 'organism',
              'table': 'bio_component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                'Hordeum vulgare',
                'Mus musculus',
                'Homo sapiens',
                'Novosphingobium capsulatum'
              ],
              'isUnique': false,
              'comment': 'The specific organism from which the bio component is sourced.',
              'LLMComment': 'Indicates the species or strain of the organism, providing context for the biological component\'s origin.'
            }
          ],
          'rowCount': 2254,
          'comment': 'This table stores sequences of biological components along with their associated metadata, such as type, description, and organism information.',
          'LLMComment': 'The `public.bio_component_sequences` table is a crucial part of a biological database that catalogs various biological components, including proteins, nucleic acids, and other biomolecules. Each entry in this table includes a unique identifier (`component_id`), the type of component (`component_type`), a textual description (`description`), the actual sequence of the component (`sequence`), and a checksum for data integrity (`sequence_md5sum`). Additionally, it links to taxonomic information through `tax_id` and provides the name of the organism (`organism`) from which the component is derived. This table is essential for researchers and bioinformaticians who need to analyze biological sequences and their properties in the context of various organisms.'
        },
        {
          'name': 'bioassay_ontology',
          'schema': 'public',
          'columns': [
            {
              'name': 'bao_id',
              'table': 'bioassay_ontology',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true
            },
            {
              'name': 'label',
              'table': 'bioassay_ontology',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true
            }
          ],
          'rowCount': 311,
          'comment': 'This table contains the bioassay ontology, which includes identifiers and labels for various bioassays used in biological and pharmacological research.',
          'LLMComment': 'The `public.bioassay_ontology` table serves as a foundational reference for bioassays in the database, providing a structured vocabulary through `bao_id` and `label` columns. This ontology is crucial for standardizing the classification and description of bioassays, facilitating data integration and analysis across various biological studies and drug development processes. With 311 entries, it supports the understanding of how different bioassays relate to biological activities and therapeutic mechanisms.'
        },
        {
          'name': 'biotherapeutic_components',
          'schema': 'public',
          'columns': [
            {
              'name': 'biocomp_id',
              'table': 'biotherapeutic_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 1240,
              'max': 3701,
              'isUnique': true,
              'comment': 'Unique identifier for each biotherapeutic component.',
              'LLMComment': 'A unique identifier assigned to each biotherapeutic component in the database.'
            },
            {
              'name': 'molregno',
              'table': 'biotherapeutic_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 20668,
              'max': 2832777,
              'isUnique': false,
              'semanticType': 'molregno',
              'comment': 'Molecular registration number associated with the biotherapeutic component.',
              'LLMComment': 'A unique identifier used in molecular databases to reference the specific molecular structure of the component.'
            },
            {
              'name': 'component_id',
              'table': 'biotherapeutic_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 20525,
              'max': 22778,
              'isUnique': false,
              'comment': 'Identifier for the specific component within the biotherapeutic.',
              'LLMComment': 'A unique identifier that distinguishes each component in the biotherapeutic components table.'
            }
          ],
          'rowCount': 2462,
          'comment': 'This table stores information about biotherapeutic components, linking them to their respective identifiers and molecular registration numbers.',
          'LLMComment': 'The \'public.biotherapeutic_components\' table is a crucial part of the biotherapeutics database, containing 2462 entries that represent various biocomponents used in therapeutic applications. Each entry is identified by a unique \'biocomp_id\' and is associated with a \'molregno\' (molecular registration number) and a \'component_id\'. This table serves as a foundational reference for understanding the specific components that make up biotherapeutic products, facilitating research and development in the field of biomedicine. It connects to other tables that provide detailed information on activities, assays, and classifications related to these components, thereby supporting comprehensive analysis and data retrieval in biotherapeutic research.'
        },
        {
          'name': 'biotherapeutics',
          'schema': 'public',
          'columns': [
            {
              'name': 'molregno',
              'table': 'biotherapeutics',
              'schema': 'public',
              'type': 'bigint',
              'min': 197,
              'max': 2832777,
              'isUnique': true,
              'semanticType': 'molregno',
              'comment': 'Unique identifier for each biotherapeutic.',
              'LLMComment': 'A unique numeric identifier assigned to each biotherapeutic in the database.'
            },
            {
              'name': 'description',
              'table': 'biotherapeutics',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Detailed description of the biotherapeutic.',
              'LLMComment': 'A textual description providing information about the biotherapeutic\'s characteristics, uses, or properties.'
            },
            {
              'name': 'helm_notation',
              'table': 'biotherapeutics',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'semanticType': 'Macromolecule',
              'comment': 'HELM notation representing the structure of the macromolecule.',
              'LLMComment': 'A standardized notation used to describe the structure of macromolecules, facilitating their representation and analysis.'
            }
          ],
          'rowCount': 23734,
          'comment': 'Table containing information about biotherapeutics, including their unique identifiers, descriptions, and HELM notation for structural representation.',
          'LLMComment': 'The \'public.biotherapeutics\' table serves as a central repository for biotherapeutic agents, which are biological products used for therapeutic purposes. Each entry is identified by a unique \'molregno\' (molecule registration number), accompanied by a descriptive text that provides insights into the biotherapeutic\'s nature and application. Additionally, the \'helm_notation\' column offers a standardized format for representing the structure of these biotherapeutics, facilitating computational analysis and integration with other biological data. This table is crucial for researchers and developers in the biopharmaceutical domain, as it links to various related tables that provide detailed information on actions, activities, and classifications associated with these therapeutic agents.'
        },
        {
          'name': 'cell_dictionary',
          'schema': 'public',
          'columns': [
            {
              'name': 'cell_id',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 5992,
              'isUnique': true,
              'comment': 'Unique identifier for each cell entry.'
            },
            {
              'name': 'cell_name',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The name of the cell type.'
            },
            {
              'name': 'cell_description',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'A detailed description of the cell type, including characteristics and functions.'
            },
            {
              'name': 'cell_source_tissue',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The specific tissue from which the cell type is derived.'
            },
            {
              'name': 'cell_source_organism',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Xenopus laevis',
                'Cricetulus griseus',
                'Chlorocebus aethiops',
                'Spodoptera frugiperda',
                'Murinae',
                'Drosophila melanogaster',
                'Macaca mulatta',
                'Ovis aries',
                'Trichoplusia ni',
                'Macaca fascicularis',
                'Saguinus oedipus',
                'Mesocricetus auratus',
                'Danio rerio',
                'Bos taurus',
                'null',
                'Mus musculus',
                'Oryctolagus cuniculus',
                'Gallus gallus',
                'Canis lupus familiaris',
                'Manduca',
                'Homo sapiens',
                'Chlorocebus sabaeus',
                'Sus scrofa',
                'Rattus norvegicus',
                'Anura',
                'Equus caballus'
              ],
              'isUnique': false,
              'comment': 'The organism from which the cell type is sourced, with specific categories listed.',
              'LLMComment': 'Examples include Xenopus laevis and Cricetulus griseus.'
            },
            {
              'name': 'cell_source_tax_id',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'bigint',
              'min': 7108,
              'max': 60711,
              'isUnique': false,
              'comment': 'Taxonomic ID corresponding to the source organism of the cell.'
            },
            {
              'name': 'clo_id',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the cell type in the Cell Line Ontology.'
            },
            {
              'name': 'efo_id',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the cell type in the Experimental Factor Ontology.'
            },
            {
              'name': 'cellosaurus_id',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the cell type in the Cellosaurus database.'
            },
            {
              'name': 'cl_lincs_id',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the cell type in the LINCS database.'
            },
            {
              'name': 'chembl_id',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'semanticType': 'CHEMBL_ID',
              'comment': 'Unique identifier for the cell type in the ChEMBL database, related to bioactivity data.',
              'LLMComment': 'Semantic identifier used for linking chemical and biological data.'
            },
            {
              'name': 'cell_ontology_id',
              'table': 'cell_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'CL_0000057',
                'null',
                'CL_0002609',
                'CL_0002586',
                'CL_0000623',
                'CL_0000094',
                'CL_0010012',
                'CL_0000136',
                'CL_0000598',
                'CL_0000359',
                'CL_0000322',
                'CL_0000775',
                'CL_0010004',
                'CL_1000497',
                'CL_0002327',
                'CL_0000584',
                'CL_0000066',
                'CL_2000058',
                'CL_0002131',
                'CL_1001608',
                'CL_0002539',
                'CL_0000232',
                'CL_0000023',
                'CL_0002608',
                'CL_0000893',
                'CL_0000576',
                'CL_0000062',
                'CL_0000148',
                'CL_0002591',
                'CL_0000451',
                'CL_0002548',
                'CL_0000771',
                'CL_0000540',
                'CL_0011103',
                'CL_0000492',
                'CL_0000129',
                'CL_0000092',
                'CL_0000169',
                'CL_0000047',
                'CL_2000004',
                'CL_0002584',
                'CL_0000530',
                'CL_0000746',
                'CL_2000001',
                'CL_0000097',
                'CL_0000233',
                'CL_2000074',
                'CL_0002306',
                'CL_0000312',
                'CL_0002476',
                'CL_0000084',
                'CL_0000182',
                'CL_0000236',
                'CL_0000134',
                'CL_2000008',
                'CL_0000581',
                'CL_0000127',
                'CL_0000056',
                'CL_0002092',
                'CL_0000235',
                'CL_0000542'
              ],
              'isUnique': false,
              'comment': 'Identifier for the cell type in the Cell Ontology, with specific categories listed.',
              'LLMComment': 'Examples include CL_0000057 and CL_0002609.'
            }
          ],
          'rowCount': 2023,
          'comment': 'This table contains a comprehensive dictionary of cell types, including their identifiers, names, descriptions, and sources.',
          'LLMComment': 'The \'public.cell_dictionary\' table serves as a central repository for information about various cell types used in biological and biomedical research. It includes essential identifiers such as \'cell_id\', \'clo_id\', and \'cell_ontology_id\' that link to external databases, facilitating cross-referencing and integration with other biological datasets. The table also provides descriptive fields like \'cell_name\' and \'cell_description\', which help researchers understand the characteristics and origins of each cell type, including the source tissue and organism. This structured information is crucial for studies in cell biology, drug development, and therapeutic research, enabling AI systems to analyze and interpret biological data effectively.'
        },
        {
          'name': 'chembl_id_lookup',
          'schema': 'public',
          'columns': [
            {
              'name': 'chembl_id',
              'table': 'chembl_id_lookup',
              'schema': 'public',
              'type': 'character varying',
              'semanticType': 'CHEMBL_ID'
            },
            {
              'name': 'entity_type',
              'table': 'chembl_id_lookup',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of entity associated with the CHEMBL ID (e.g., compound, target, assay).',
              'LLMComment': 'Describes the category of the entity linked to the CHEMBL ID, providing context for its role in the database.'
            },
            {
              'name': 'entity_id',
              'table': 'chembl_id_lookup',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Unique identifier for the entity within its type.',
              'LLMComment': 'A numeric identifier that uniquely distinguishes each entity within its respective type.'
            },
            {
              'name': 'status',
              'table': 'chembl_id_lookup',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Current status of the entity (e.g., active, inactive).',
              'LLMComment': 'Indicates the operational state of the entity, which can affect its usability in research or applications.'
            },
            {
              'name': 'last_active',
              'table': 'chembl_id_lookup',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Timestamp indicating the last time the entity was active or updated.',
              'LLMComment': 'Represents the last recorded activity or modification date for the entity, useful for tracking changes over time.'
            }
          ],
          'rowCount': 4640251,
          'comment': 'This table provides a mapping of ChEMBL IDs to their corresponding entities, including their types and statuses, along with the last active year of the entity.',
          'LLMComment': 'The `public.chembl_id_lookup` table serves as a crucial reference for identifying and categorizing various entities in the ChEMBL database, which is a comprehensive resource for bioactive drug-like compounds. Each entry in this table links a unique ChEMBL ID to its entity type (such as compound, target, or assay), a unique entity ID, its current status (active, inactive, etc.), and the last year it was active. This information is essential for researchers and AI systems to understand the lifecycle and classification of bioactive compounds and their interactions, facilitating drug discovery and development processes.'
        },
        {
          'name': 'chembl_release',
          'schema': 'public',
          'columns': [
            {
              'name': 'chembl_release_id',
              'table': 'chembl_release',
              'schema': 'public',
              'type': 'integer',
              'min': 1,
              'max': 34,
              'isUnique': true
            },
            {
              'name': 'chembl_release',
              'table': 'chembl_release',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'CHEMBL_10',
                'CHEMBL_29',
                'CHEMBL_3',
                'CHEMBL_14',
                'CHEMBL_11',
                'CHEMBL_1',
                'CHEMBL_24',
                'CHEMBL_19',
                'CHEMBL_16',
                'CHEMBL_5',
                'CHEMBL_2',
                'CHEMBL_15',
                'CHEMBL_13',
                'CHEMBL_26',
                'CHEMBL_22',
                'CHEMBL_7',
                'CHEMBL_28',
                'CHEMBL_32',
                'CHEMBL_33',
                'CHEMBL_23',
                'CHEMBL_17',
                'CHEMBL_21',
                'CHEMBL_18',
                'CHEMBL_20',
                'CHEMBL_34',
                'CHEMBL_27',
                'CHEMBL_12',
                'CHEMBL_6',
                'CHEMBL_8',
                'CHEMBL_30',
                'CHEMBL_9',
                'CHEMBL_4',
                'CHEMBL_31',
                'CHEMBL_25'
              ],
              'isUnique': true,
              'comment': 'Identifier for the specific ChEMBL release version, indicating the dataset\'s version.',
              'LLMComment': 'A unique identifier representing the version of the ChEMBL database release, such as CHEMBL_10 or CHEMBL_29.'
            },
            {
              'name': 'creation_date',
              'table': 'chembl_release',
              'schema': 'public',
              'type': 'timestamp without time zone',
              'isUnique': true,
              'comment': 'The date and time when the ChEMBL release was created.',
              'LLMComment': 'Timestamp indicating when this specific ChEMBL release was generated.'
            }
          ],
          'rowCount': 34,
          'comment': 'This table stores information about different releases of the ChEMBL database, including unique identifiers, release names, and their creation dates.',
          'LLMComment': 'The `public.chembl_release` table is a key component of the ChEMBL database, which is a comprehensive resource for bioactive drug-like small molecules. This table tracks the various versions of the ChEMBL database, providing essential metadata such as the unique identifier for each release, the name of the release, and the date it was created. This information is crucial for researchers and developers who need to reference specific versions of the database for consistency in data analysis and application development.'
        },
        {
          'name': 'component_class',
          'schema': 'public',
          'columns': [
            {
              'name': 'component_id',
              'table': 'component_class',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 20524,
              'isUnique': false
            },
            {
              'name': 'protein_class_id',
              'table': 'component_class',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 1715,
              'isUnique': false,
              'comment': 'Identifier for the protein class associated with the component.',
              'LLMComment': 'This ID links to a specific protein class, categorizing the component based on its protein characteristics.'
            },
            {
              'name': 'comp_class_id',
              'table': 'component_class',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 14064,
              'isUnique': true,
              'comment': 'Unique identifier for the component class.',
              'LLMComment': 'This ID uniquely identifies the class of the component, distinguishing it from other component classes.'
            }
          ],
          'rowCount': 11476,
          'comment': 'This table links components to their respective protein classes, providing a classification system for biological components based on their protein characteristics.',
          'LLMComment': 'The `public.component_class` table serves as a crucial mapping between biological components and their associated protein classes within a larger biological database. Each entry in this table connects a unique component (identified by `component_id`) to a specific protein class (`protein_class_id`) and a classification identifier (`comp_class_id`). This structure is essential for understanding the functional roles of various components in biological processes and aids in the classification and retrieval of data related to proteins, which is vital for research in fields such as biochemistry, pharmacology, and molecular biology.'
        },
        {
          'name': 'component_domains',
          'schema': 'public',
          'columns': [
            {
              'name': 'compd_id',
              'table': 'component_domains',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'domain_id',
              'table': 'component_domains',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the domain associated with the component.',
              'LLMComment': 'Unique identifier representing a specific domain in the context of the component.'
            },
            {
              'name': 'component_id',
              'table': 'component_domains',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the component that this domain is related to.',
              'LLMComment': 'Unique identifier for the component that the domain is part of.'
            },
            {
              'name': 'start_position',
              'table': 'component_domains',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'The starting position of the domain within the component.',
              'LLMComment': 'Indicates where the domain begins in the sequence of the component.'
            },
            {
              'name': 'end_position',
              'table': 'component_domains',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'The ending position of the domain within the component.',
              'LLMComment': 'Indicates where the domain ends in the sequence of the component.'
            }
          ],
          'rowCount': 0,
          'comment': 'This table maps components to their respective domains, indicating the start and end positions of each component within a domain.',
          'LLMComment': 'The `public.component_domains` table serves as a crucial link between biological components and their functional domains. Each entry in this table associates a specific component (identified by `component_id`) with a domain (identified by `domain_id`), detailing the positions (`start_position` and `end_position`) that define the span of the component within the domain. This mapping is essential for understanding the structural and functional relationships in biological data, particularly in the context of biotherapeutics and molecular biology, where the interaction of components and domains can influence drug design and therapeutic efficacy.'
        },
        {
          'name': 'component_go',
          'schema': 'public',
          'columns': [
            {
              'name': 'comp_go_id',
              'table': 'component_go',
              'schema': 'public',
              'type': 'bigint',
              'min': 2547872,
              'max': 2698227,
              'isUnique': true
            },
            {
              'name': 'component_id',
              'table': 'component_go',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 20524,
              'isUnique': false,
              'comment': 'Identifier for the component associated with the GO term.',
              'LLMComment': 'This column links to the unique identifier of a component in the system.'
            },
            {
              'name': 'go_id',
              'table': 'component_go',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Gene Ontology identifier representing a biological concept.',
              'LLMComment': 'This column stores the unique identifier for a Gene Ontology term, which describes a specific biological function, process, or cellular component.'
            }
          ],
          'rowCount': 150356,
          'comment': 'This table links components to Gene Ontology (GO) terms, facilitating the understanding of biological functions associated with various components.',
          'LLMComment': 'The \'public.component_go\' table serves as a crucial mapping between biological components and their corresponding Gene Ontology (GO) identifiers. With 150,356 entries, it connects each component (identified by \'component_id\') to specific GO terms (represented by \'go_id\'), which describe the biological processes, cellular components, and molecular functions associated with those components. This linkage is essential for researchers and AI systems to analyze and interpret the functional roles of various biological entities in the context of genomics and proteomics, enhancing the understanding of their contributions to biological systems.'
        },
        {
          'name': 'component_sequences',
          'schema': 'public',
          'columns': [
            {
              'name': 'component_id',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 20524,
              'isUnique': true,
              'comment': 'Unique identifier for each component.'
            },
            {
              'name': 'component_type',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'DNA',
                'RNA',
                'PROTEIN'
              ],
              'isUnique': false,
              'comment': 'Type of biological macromolecule (DNA, RNA, PROTEIN).',
              'LLMComment': 'Indicates the classification of the component, which can be DNA, RNA, or PROTEIN.'
            },
            {
              'name': 'accession',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Accession number for the component, used for database referencing.',
              'LLMComment': 'A unique identifier assigned to the component in biological databases.'
            },
            {
              'name': 'sequence',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'text',
              'isUnique': false,
              'semanticType': 'Macromolecule',
              'comment': 'The actual sequence of the macromolecule.',
              'LLMComment': 'Represents the linear sequence of nucleotides or amino acids in the macromolecule.'
            },
            {
              'name': 'sequence_md5sum',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'MD5 checksum of the sequence for integrity verification.',
              'LLMComment': 'A hash value used to verify the integrity of the sequence data.'
            },
            {
              'name': 'description',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'A textual description of the component.',
              'LLMComment': 'Provides additional information about the component, such as its function or characteristics.'
            },
            {
              'name': 'tax_id',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'bigint',
              'min': 69,
              'max': 3052230,
              'isUnique': false,
              'comment': 'Taxonomic identifier for the organism from which the component is derived.',
              'LLMComment': 'A unique identifier that corresponds to the organism in taxonomic databases.'
            },
            {
              'name': 'organism',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Name of the organism from which the component originates.',
              'LLMComment': 'The scientific name of the organism associated with the biological component.'
            },
            {
              'name': 'db_source',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Manual',
                'TREMBL',
                'SWISS-PROT'
              ],
              'isUnique': false,
              'comment': 'Source of the database entry (e.g., Manual, TREMBL, SWISS-PROT).',
              'LLMComment': 'Indicates the origin of the data, which can help assess its reliability.'
            },
            {
              'name': 'db_version',
              'table': 'component_sequences',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                '2022_02',
                '2022_05',
                '2024_01',
                '2023_02',
                '2021_04'
              ],
              'isUnique': false,
              'comment': 'Version of the database from which the component data is sourced.',
              'LLMComment': 'Specifies the version of the database that contains the component information.'
            }
          ],
          'rowCount': 11281,
          'comment': 'This table stores sequences of various biological components, including their identifiers, types, and associated metadata.',
          'LLMComment': 'The `public.component_sequences` table is a crucial part of a biological database that catalogs sequences of different biological components, such as proteins or nucleic acids. Each entry includes a unique `component_id`, the type of component (e.g., protein, DNA), its sequence, and metadata like the organism it belongs to and the source database. This table is essential for researchers in bioinformatics and molecular biology, as it provides the foundational data needed for sequence analysis, comparison, and functional studies.'
        },
        {
          'name': 'component_synonyms',
          'schema': 'public',
          'columns': [
            {
              'name': 'compsyn_id',
              'table': 'component_synonyms',
              'schema': 'public',
              'type': 'bigint',
              'min': 860862,
              'max': 1860009,
              'isUnique': true
            },
            {
              'name': 'component_id',
              'table': 'component_synonyms',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 20524,
              'isUnique': false,
              'comment': 'Identifier for the associated component in the database.',
              'LLMComment': 'Unique identifier linking to the specific component.'
            },
            {
              'name': 'component_synonym',
              'table': 'component_synonyms',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Alternative name or term used to refer to the component.',
              'LLMComment': 'Synonym for the component, which may be used in various contexts.'
            },
            {
              'name': 'syn_type',
              'table': 'component_synonyms',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'UNIPROT',
                'GENE_SYMBOL',
                'EC_NUMBER',
                'GENE_SYMBOL_OTHER'
              ],
              'isUnique': false,
              'comment': 'Type of synonym, indicating the source or classification of the synonym.',
              'LLMComment': 'Categorizes the synonym based on its origin, such as UNIPROT or GENE_SYMBOL.'
            }
          ],
          'rowCount': 108004,
          'comment': 'This table stores synonyms for various components used in biological and chemical contexts, allowing for different naming conventions and terminologies to be linked to a single component ID.',
          'LLMComment': 'The `public.component_synonyms` table is crucial for managing the diverse nomenclature associated with biological and chemical components. With 108,004 entries, it provides a mapping between unique component identifiers (`component_id`) and their synonyms (`component_synonym`), which can vary by context or usage (`syn_type`). This is particularly important in fields like drug discovery and bioinformatics, where different databases or research articles may refer to the same component using different names. By standardizing these synonyms, the table facilitates data integration and retrieval across various biological datasets, enhancing the ability to perform comprehensive analyses and ensuring consistency in research and application.'
        },
        {
          'name': 'compound_properties',
          'schema': 'public',
          'columns': [
            {
              'name': 'molregno',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'bigint',
              'semanticType': 'molregno'
            },
            {
              'name': 'mw_freebase',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Molecular weight derived from Freebase data.',
              'LLMComment': 'Represents the molecular weight of the compound as sourced from Freebase, a large collaborative knowledge base.'
            },
            {
              'name': 'alogp',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Calculated logP value based on the compound\'s structure.',
              'LLMComment': 'The calculated logP (partition coefficient) value indicating the hydrophobicity of the compound.'
            },
            {
              'name': 'hba',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Number of hydrogen bond acceptors in the compound.',
              'LLMComment': 'Counts the hydrogen bond acceptor sites in the molecular structure.'
            },
            {
              'name': 'hbd',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Number of hydrogen bond donors in the compound.',
              'LLMComment': 'Counts the hydrogen bond donor sites in the molecular structure.'
            },
            {
              'name': 'psa',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Polar surface area of the compound.',
              'LLMComment': 'The polar surface area, which can influence the compound\'s absorption and permeability.'
            },
            {
              'name': 'rtb',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Number of rotatable bonds in the compound.',
              'LLMComment': 'Indicates the flexibility of the molecule by counting the rotatable bonds.'
            },
            {
              'name': 'ro3_pass',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Indicates if the compound passes the Rule of Three for drug-like properties.',
              'LLMComment': 'A flag indicating whether the compound meets the criteria of the Rule of Three, which is often used in drug discovery.'
            },
            {
              'name': 'num_ro5_violations',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Number of violations of Lipinski\'s Rule of Five.',
              'LLMComment': 'Counts how many of Lipinski\'s Rule of Five criteria the compound does not satisfy.'
            },
            {
              'name': 'cx_most_apka',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Most positive atom contribution to the polar surface area.',
              'LLMComment': 'Represents the most significant contribution of positive atoms to the polar surface area.'
            },
            {
              'name': 'cx_most_bpka',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Most negative atom contribution to the polar surface area.',
              'LLMComment': 'Represents the most significant contribution of negative atoms to the polar surface area.'
            },
            {
              'name': 'cx_logp',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Calculated logP value using the ChemAxon method.',
              'LLMComment': 'The logP value calculated using ChemAxon\'s algorithms, indicating the compound\'s hydrophobicity.'
            },
            {
              'name': 'cx_logd',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Calculated logD value using the ChemAxon method.',
              'LLMComment': 'The logD value calculated using ChemAxon\'s algorithms, representing the distribution of the compound between octanol and water.'
            },
            {
              'name': 'molecular_species',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of molecular species (e.g., ion, neutral).',
              'LLMComment': 'Describes the molecular species classification, such as whether the compound is an ion or neutral.'
            },
            {
              'name': 'full_mwt',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Full molecular weight including isotopes.',
              'LLMComment': 'The complete molecular weight of the compound, accounting for all isotopes present.'
            },
            {
              'name': 'aromatic_rings',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Number of aromatic rings in the compound.',
              'LLMComment': 'Counts the number of aromatic rings, which can affect the compound\'s properties.'
            },
            {
              'name': 'heavy_atoms',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Number of heavy atoms in the compound (non-hydrogen).',
              'LLMComment': 'Counts the heavy atoms (all atoms except hydrogen) in the molecular structure.'
            },
            {
              'name': 'qed_weighted',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Quantitative Estimation of Drug-likeness (weighted).',
              'LLMComment': 'A weighted score representing the drug-likeness of the compound based on various molecular properties.'
            },
            {
              'name': 'mw_monoisotopic',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Monoisotopic molecular weight of the compound.',
              'LLMComment': 'The molecular weight calculated using the most abundant isotopes of each element.'
            },
            {
              'name': 'full_molformula',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Full molecular formula of the compound.',
              'LLMComment': 'The complete molecular formula representing the composition of the compound.'
            },
            {
              'name': 'hba_lipinski',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Number of hydrogen bond acceptors according to Lipinski\'s rules.',
              'LLMComment': 'Counts the hydrogen bond acceptors as per Lipinski\'s criteria for drug-likeness.'
            },
            {
              'name': 'hbd_lipinski',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Number of hydrogen bond donors according to Lipinski\'s rules.',
              'LLMComment': 'Counts the hydrogen bond donors as per Lipinski\'s criteria for drug-likeness.'
            },
            {
              'name': 'num_lipinski_ro5_violations',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Number of Lipinski\'s Rule of Five violations.',
              'LLMComment': 'Counts the total number of violations against Lipinski\'s Rule of Five.'
            },
            {
              'name': 'np_likeness_score',
              'table': 'compound_properties',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Score indicating the likeness to natural products.',
              'LLMComment': 'A score that assesses how closely the compound resembles natural products.'
            }
          ],
          'rowCount': 2412898,
          'comment': 'This table contains detailed properties of chemical compounds, including molecular weight, logP values, and various structural features relevant to drug discovery and development.',
          'LLMComment': 'The \'compound_properties\' table is a comprehensive repository of chemical compound characteristics essential for pharmaceutical research. It includes quantitative measures such as molecular weight (mw_freebase), partition coefficients (alogp, cx_logp), and structural attributes (number of hydrogen bond acceptors (hba), donors (hbd), and aromatic rings). These properties are crucial for assessing the drug-likeness and bioavailability of compounds, aiding in the identification of potential drug candidates. The table also tracks compliance with Lipinski\'s Rule of Five, which is a set of criteria used to evaluate the drug-like properties of compounds. With over 2.4 million entries, this table serves as a foundational resource for researchers in medicinal chemistry and pharmacology.'
        },
        {
          'name': 'compound_records',
          'schema': 'public',
          'columns': [
            {
              'name': 'record_id',
              'table': 'compound_records',
              'schema': 'public',
              'type': 'bigint'
            },
            {
              'name': 'molregno',
              'table': 'compound_records',
              'schema': 'public',
              'type': 'bigint',
              'semanticType': 'molregno',
              'comment': 'Unique identifier for the compound in the database, used for referencing molecular records.',
              'LLMComment': 'Molecular registration number, a unique identifier for compounds.'
            },
            {
              'name': 'doc_id',
              'table': 'compound_records',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the document associated with the compound record, linking to additional information or sources.',
              'LLMComment': 'Document ID that provides context or source for the compound record.'
            },
            {
              'name': 'compound_key',
              'table': 'compound_records',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'A unique key representing the compound, often used for indexing or retrieval purposes.',
              'LLMComment': 'Key used to uniquely identify the compound within the database.'
            },
            {
              'name': 'compound_name',
              'table': 'compound_records',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'The common or systematic name of the compound, used for identification in chemical contexts.',
              'LLMComment': 'Name of the compound, which may be a common or IUPAC name.'
            },
            {
              'name': 'src_id',
              'table': 'compound_records',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Identifier for the source of the compound data, indicating where the information was obtained from.',
              'LLMComment': 'Source ID that indicates the origin of the compound data.'
            },
            {
              'name': 'src_compound_id',
              'table': 'compound_records',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Identifier for the compound as it appears in the source database, allowing cross-referencing.',
              'LLMComment': 'Compound ID from the source database for cross-referencing purposes.'
            },
            {
              'name': 'cidx',
              'table': 'compound_records',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Index or identifier used for internal tracking of the compound within the system.',
              'LLMComment': 'Internal index for tracking the compound in the database.'
            }
          ],
          'rowCount': 3106257,
          'comment': 'This table stores records of chemical compounds, including their identifiers and names, as well as references to related documents and sources.',
          'LLMComment': 'The \'compound_records\' table is a central repository for chemical compound data within the database. It contains over 3 million entries, each uniquely identified by \'record_id\'. Key attributes include \'molregno\' (a unique identifier for the compound), \'compound_key\' (a string representing the compound\'s key), and \'compound_name\' (the name of the compound). This table serves as a foundational element for linking to other related tables, such as \'compound_properties\' and \'compound_structures\', which provide additional details about the compounds\' characteristics and structural information. Understanding this table is crucial for tasks involving chemical data analysis, compound identification, and research in pharmacology and medicinal chemistry.'
        },
        {
          'name': 'compound_structural_alerts',
          'schema': 'public',
          'columns': [
            {
              'name': 'cpd_str_alert_id',
              'table': 'compound_structural_alerts',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Unique identifier for each structural alert in the database.',
              'LLMComment': 'A unique identifier assigned to each compound structural alert.'
            },
            {
              'name': 'molregno',
              'table': 'compound_structural_alerts',
              'schema': 'public',
              'type': 'bigint',
              'semanticType': 'molregno',
              'comment': 'Molecular registration number associated with the compound.',
              'LLMComment': 'A unique identifier for the compound, used for tracking and referencing.'
            },
            {
              'name': 'alert_id',
              'table': 'compound_structural_alerts',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the specific structural alert type.',
              'LLMComment': 'A unique identifier that categorizes the type of structural alert.'
            }
          ],
          'rowCount': 4436020,
          'comment': 'This table stores information about structural alerts associated with chemical compounds, linking each alert to a specific compound and an alert identifier.',
          'LLMComment': 'The `public.compound_structural_alerts` table is a critical component in the domain of cheminformatics and toxicology. It contains over 4.4 million records that map chemical compounds to their respective structural alerts, which are features or patterns in a compound\'s structure that may indicate potential toxicity or undesirable biological activity. Each record includes a unique identifier for the structural alert (`cpd_str_alert_id`), a reference to the compound (`molregno`), and the specific alert identifier (`alert_id`). This table is essential for researchers and developers working on drug discovery and safety assessment, as it helps in identifying compounds that may pose risks based on their structural characteristics.'
        },
        {
          'name': 'compound_structures',
          'schema': 'public',
          'columns': [
            {
              'name': 'molregno',
              'table': 'compound_structures',
              'schema': 'public',
              'type': 'bigint',
              'semanticType': 'molregno',
              'comment': 'Unique identifier for the compound structure.'
            },
            {
              'name': 'molfile',
              'table': 'compound_structures',
              'schema': 'public',
              'type': 'text',
              'semanticType': 'Molecule',
              'comment': 'Text representation of the molecular structure.',
              'LLMComment': 'Contains the molecular structure in a text format, often used for chemical informatics.'
            },
            {
              'name': 'standard_inchi',
              'table': 'compound_structures',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Standard International Chemical Identifier for the compound.',
              'LLMComment': 'A textual identifier that provides a standard way to represent a chemical substance.'
            },
            {
              'name': 'standard_inchi_key',
              'table': 'compound_structures',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Shortened version of the standard InChI for easier referencing.',
              'LLMComment': 'A fixed-length, hashed representation of the InChI, used for quick identification.'
            },
            {
              'name': 'canonical_smiles',
              'table': 'compound_structures',
              'schema': 'public',
              'type': 'character varying',
              'semanticType': 'Molecule',
              'comment': 'Canonical representation of the molecule in SMILES format.',
              'LLMComment': 'A string representation of the molecular structure that encodes the connectivity of atoms.'
            }
          ],
          'rowCount': 2409270,
          'comment': 'This table stores information about chemical compounds, including their unique identifiers and structural representations.',
          'LLMComment': 'The \'compound_structures\' table is a crucial component of a chemical database, containing over 2.4 million entries of chemical compounds. Each entry includes a unique molecular registration number (molregno), a textual representation of the molecular structure (molfile), and standardized formats for chemical identification such as InChI (standard_inchi) and SMILES (canonical_smiles). This table serves as a foundational resource for researchers and AI systems in cheminformatics, enabling the analysis and retrieval of compound data for drug discovery, molecular modeling, and various bioinformatics applications.'
        },
        {
          'name': 'confidence_score_lookup',
          'schema': 'public',
          'columns': [
            {
              'name': 'confidence_score',
              'table': 'confidence_score_lookup',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 9,
              'isUnique': true,
              'comment': 'A numerical score indicating the confidence level of the target assignment, ranging from 0 to 9.',
              'LLMComment': 'Represents the confidence level of the target assignment, with higher values indicating greater confidence.'
            },
            {
              'name': 'description',
              'table': 'confidence_score_lookup',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Default value - Target unknown or has yet to be assigned',
                'Multiple homologous protein targets may be assigned',
                'Homologous single protein target assigned',
                'Target assigned is subcellular fraction',
                'Direct single protein target assigned',
                'Target assigned is molecular non-protein target',
                'Direct protein complex subunits assigned',
                'Multiple direct protein targets may be assigned',
                'Target assigned is non-molecular',
                'Homologous protein complex subunits assigned'
              ],
              'isUnique': true,
              'comment': 'A textual description categorizing the nature of the target assignment, with various possible values indicating the status of target assignment.',
              'LLMComment': 'Describes the specific nature of the target assignment, providing context on whether the target is known, homologous, or assigned to a specific subcellular fraction.'
            },
            {
              'name': 'target_mapping',
              'table': 'confidence_score_lookup',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Protein',
                'Homologous protein',
                'Multiple proteins',
                'Homologous protein complex',
                'Protein complex',
                'Non-molecular',
                'Unassigned',
                'Subcellular fraction',
                'Multiple homologous proteins',
                'Molecular (non-protein)'
              ],
              'isUnique': true,
              'comment': 'A textual categorization of the type of target mapping, indicating whether it is a single protein, homologous proteins, or complexes.',
              'LLMComment': 'Categorizes the type of target mapping, helping to understand the relationship between proteins and their assignments.'
            }
          ],
          'rowCount': 10,
          'comment': 'This table contains lookup values for confidence scores associated with various assessments or predictions in the database.',
          'LLMComment': 'The \'confidence_score_lookup\' table serves as a reference for confidence scores, which are numerical indicators (smallint) reflecting the reliability of certain predictions or assessments within the domain. Each score is accompanied by a descriptive text that explains its significance and a \'target_mapping\' field that likely links to specific entities or assessments in other tables. This table is crucial for interpreting the confidence levels of various actions or activities related to biotherapeutics, drug mechanisms, and other biological assessments, thereby aiding in decision-making processes in research and development.'
        },
        {
          'name': 'curation_lookup',
          'schema': 'public',
          'columns': [
            {
              'name': 'curated_by',
              'table': 'curation_lookup',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Intermediate',
                'Expert',
                'Autocuration'
              ],
              'isUnique': true,
              'comment': 'Indicates the level of expertise of the curator responsible for the entry.',
              'LLMComment': 'Specifies the curator\'s expertise level, which can be Intermediate, Expert, or Autocuration.'
            },
            {
              'name': 'description',
              'table': 'curation_lookup',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Curated against extractor target assignment',
                'Curated against ChEMBL target assignment from original publication',
                'Curated against ChEMBL target assignment from assay description'
              ],
              'isUnique': true,
              'comment': 'Describes the basis for curation, detailing the source of the target assignment.',
              'LLMComment': 'Provides context on how the curation was performed, referencing specific target assignment sources such as extractor targets or ChEMBL.'
            }
          ],
          'rowCount': 3,
          'comment': 'Lookup table for curation details, including who curated the data and a description of the curation.',
          'LLMComment': 'The `public.curation_lookup` table serves as a reference for understanding the curation process within the database. It contains two key columns: `curated_by`, which indicates the individual or entity responsible for the curation, and `description`, which provides additional context or details about the curation. This table is essential for tracking the provenance of data entries and ensuring transparency in the curation process, which is particularly important in domains like biomedicine and pharmacology where data integrity is critical.'
        },
        {
          'name': 'data_validity_lookup',
          'schema': 'public',
          'columns': [
            {
              'name': 'data_validity_comment',
              'table': 'data_validity_lookup',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Potential missing data',
                'Author confirmed error',
                'Manually validated',
                'Outside typical range',
                'Potential author error',
                'Non standard unit for type',
                'Potential transcription error'
              ],
              'isUnique': true,
              'comment': 'Comment indicating the validity status of the data.',
              'LLMComment': 'Describes the reason for data validity status, such as potential missing data or author confirmed error.'
            },
            {
              'name': 'description',
              'table': 'data_validity_lookup',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'semanticType': 'Text',
              'comment': 'A textual description providing additional context about the data validity comment.',
              'LLMComment': 'Unique text that elaborates on the specific data validity status.'
            }
          ],
          'rowCount': 7,
          'comment': 'Lookup table for data validity comments and their descriptions.',
          'LLMComment': 'The `data_validity_lookup` table serves as a reference for understanding the validity of data entries within the database. It contains a list of comments related to data validity, each paired with a detailed description. This is crucial for ensuring data integrity and providing context for users and systems that interact with the database, particularly in domains like biomedicine and pharmacology where data accuracy is paramount.'
        },
        {
          'name': 'defined_daily_dose',
          'schema': 'public',
          'columns': [
            {
              'name': 'atc_code',
              'table': 'defined_daily_dose',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Anatomical Therapeutic Chemical classification code for the drug.',
              'LLMComment': 'A code that classifies the drug based on its anatomical and therapeutic properties.'
            },
            {
              'name': 'ddd_units',
              'table': 'defined_daily_dose',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'mcg',
                'MU',
                'tablet',
                'TU',
                'null',
                'g',
                'ml',
                'U',
                'LSU',
                'mmol',
                'mg'
              ],
              'isUnique': false,
              'comment': 'Units of measurement for the defined daily dose, such as mcg or tablet.',
              'LLMComment': 'Specifies the unit in which the defined daily dose is measured, indicating the form of the drug.'
            },
            {
              'name': 'ddd_admr',
              'table': 'defined_daily_dose',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'ointment',
                'O',
                'Inhal.aerosol',
                'null',
                'Chewing gum',
                'V',
                'lamella',
                'oral aerosol',
                's.c. implant',
                'TD patch',
                'urethral',
                'N',
                'P',
                'Inhal.powder     delivered dose',
                'Inhal.solution',
                'SL',
                'Inhal',
                'intravesical',
                'R',
                'Inh.aerosol',
                'implant',
                'Refers to umeclidinium, delivered dose',
                'TD',
                'Instill.sol.',
                'Inhal.powder'
              ],
              'isUnique': false,
              'comment': 'Administration route or form of the drug, such as ointment or inhalation aerosol.',
              'LLMComment': 'Describes how the drug is administered or its physical form.'
            },
            {
              'name': 'ddd_comment',
              'table': 'defined_daily_dose',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Additional comments or notes regarding the defined daily dose.',
              'LLMComment': 'Any supplementary information or remarks related to the defined daily dose.'
            },
            {
              'name': 'ddd_id',
              'table': 'defined_daily_dose',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 3438,
              'isUnique': true
            },
            {
              'name': 'ddd_value',
              'table': 'defined_daily_dose',
              'schema': 'public',
              'type': 'numeric',
              'min': 0.10000000149011612,
              'max': 800,
              'isUnique': false,
              'comment': 'The numerical value representing the defined daily dose for the drug.',
              'LLMComment': 'The actual dosage amount that is considered a standard daily dose for the drug.'
            }
          ],
          'rowCount': 2672,
          'comment': 'This table stores the defined daily doses (DDD) for various drugs, identified by their ATC codes, along with associated units, administration routes, and additional comments.',
          'LLMComment': 'The \'public.defined_daily_dose\' table is a crucial resource in pharmacology and drug regulation, providing standardized information on the defined daily doses (DDD) for medications. Each entry is linked to an ATC (Anatomical Therapeutic Chemical) code, which categorizes drugs based on their therapeutic use. The table includes details such as the units of measurement for the DDD, the recommended administration route, and any relevant comments that may provide additional context for healthcare professionals. This data is essential for ensuring safe and effective medication dosing in clinical settings.'
        },
        {
          'name': 'docs',
          'schema': 'public',
          'columns': [
            {
              'name': 'doc_id',
              'table': 'docs',
              'schema': 'public',
              'type': 'bigint',
              'min': -1,
              'max': 126266,
              'isUnique': true,
              'comment': 'Document ID'
            },
            {
              'name': 'journal',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Name of the journal where the document is published',
              'LLMComment': 'The journal title associated with the publication.'
            },
            {
              'name': 'year',
              'table': 'docs',
              'schema': 'public',
              'type': 'integer',
              'min': 1974,
              'max': 2023,
              'isUnique': false,
              'comment': 'Year of publication',
              'LLMComment': 'The year in which the document was published.'
            },
            {
              'name': 'volume',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Volume number of the journal',
              'LLMComment': 'The volume number of the journal in which the document appears.'
            },
            {
              'name': 'issue',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Issue number of the journal',
              'LLMComment': 'The specific issue number of the journal for the publication.'
            },
            {
              'name': 'first_page',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'First page number of the document',
              'LLMComment': 'The starting page number of the document in the journal.'
            },
            {
              'name': 'last_page',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Last page number of the document',
              'LLMComment': 'The ending page number of the document in the journal.'
            },
            {
              'name': 'pubmed_id',
              'table': 'docs',
              'schema': 'public',
              'type': 'bigint',
              'min': 1534,
              'max': 38116432,
              'isUnique': false,
              'comment': 'PubMed identifier for the document',
              'LLMComment': 'Unique identifier for the document in the PubMed database.'
            },
            {
              'name': 'doi',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Digital Object Identifier for the document',
              'LLMComment': 'A unique alphanumeric string assigned to the document for identification.'
            },
            {
              'name': 'chembl_id',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'semanticType': 'CHEMBL_ID',
              'comment': 'ChEMBL identifier for the document',
              'LLMComment': 'Unique identifier used in the ChEMBL database for chemical entities.'
            },
            {
              'name': 'title',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Title of the document',
              'LLMComment': 'The title of the publication or document.'
            },
            {
              'name': 'doc_type',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'PATENT',
                'BOOK',
                'PUBLICATION',
                'DATASET'
              ],
              'isUnique': false,
              'comment': 'Type of document (e.g., PATENT, BOOK, PUBLICATION, DATASET)',
              'LLMComment': 'Categorization of the document type.'
            },
            {
              'name': 'authors',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'List of authors of the document',
              'LLMComment': 'Names of the individuals who authored the document.'
            },
            {
              'name': 'abstract',
              'table': 'docs',
              'schema': 'public',
              'type': 'text',
              'isUnique': false,
              'semanticType': 'Text',
              'comment': 'Summary of the document\'s content',
              'LLMComment': 'A brief summary or overview of the document\'s main points.'
            },
            {
              'name': 'patent_id',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the associated patent, if applicable',
              'LLMComment': 'Unique identifier for the patent related to the document.'
            },
            {
              'name': 'ridx',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Reference index for internal use',
              'LLMComment': 'An internal reference index for the document.'
            },
            {
              'name': 'src_id',
              'table': 'docs',
              'schema': 'public',
              'type': 'integer',
              'min': 0,
              'max': 69,
              'isUnique': false,
              'comment': 'Source identifier for the document',
              'LLMComment': 'Identifier indicating the source of the document.'
            },
            {
              'name': 'chembl_release_id',
              'table': 'docs',
              'schema': 'public',
              'type': 'integer',
              'min': 1,
              'max': 34,
              'isUnique': false,
              'comment': 'Release identifier for ChEMBL data',
              'LLMComment': 'Identifier indicating the specific release version of ChEMBL data.'
            },
            {
              'name': 'contact',
              'table': 'docs',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null'
              ],
              'isUnique': false,
              'comment': 'Contact information for the corresponding author or institution',
              'LLMComment': 'Details for reaching the corresponding author or institution.'
            }
          ],
          'rowCount': 89892,
          'comment': 'This table contains metadata for scientific documents, including journal articles, patents, and other publications relevant to the field of chemistry and biomedicine.',
          'LLMComment': 'The \'public.docs\' table serves as a comprehensive repository of scientific literature, encompassing various types of documents such as journal articles and patents. Each entry is identified by a unique \'doc_id\' and includes essential bibliographic details like \'title\', \'authors\', \'abstract\', and publication specifics (e.g., \'journal\', \'year\', \'volume\', \'issue\'). Additionally, it links to external identifiers such as \'pubmed_id\', \'doi\', and \'chembl_id\', facilitating cross-referencing with other databases. This table is crucial for researchers and AI systems seeking to access and analyze scientific knowledge in the domains of chemistry and biomedicine.'
        },
        {
          'name': 'domains',
          'schema': 'public',
          'columns': [
            {
              'name': 'domain_id',
              'table': 'domains',
              'schema': 'public',
              'type': 'bigint',
              'min': 2629,
              'max': 13996,
              'isUnique': true,
              'comment': 'Unique identifier for each domain'
            },
            {
              'name': 'domain_type',
              'table': 'domains',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Pfam-A'
              ],
              'isUnique': false,
              'comment': 'Type of the domain, typically categorized under Pfam-A',
              'LLMComment': 'Specifies the classification of the domain, indicating its family or type within the Pfam database.'
            },
            {
              'name': 'source_domain_id',
              'table': 'domains',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the source domain, if applicable',
              'LLMComment': 'Refers to the original domain ID from which this domain is derived or associated.'
            },
            {
              'name': 'domain_name',
              'table': 'domains',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Name of the domain'
            },
            {
              'name': 'domain_description',
              'table': 'domains',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null'
              ],
              'isUnique': false,
              'comment': 'Description of the domain\'s function or characteristics',
              'LLMComment': 'Provides detailed information about the domain, including its biological significance and role.'
            }
          ],
          'rowCount': 2385,
          'comment': 'This table stores information about various domains related to biological entities, including their types, names, and descriptions.',
          'LLMComment': 'The \'public.domains\' table is a key component in the biological data schema, containing 2385 entries that represent distinct domains associated with biological entities. Each entry includes a unique identifier (domain_id), the type of domain (domain_type), a reference to a source domain (source_domain_id), the name of the domain (domain_name), and a detailed description (domain_description). This table is essential for linking biological data across various related tables, such as assays, components, and mechanisms, facilitating comprehensive analysis and research in bioinformatics and drug discovery.'
        },
        {
          'name': 'drug_indication',
          'schema': 'public',
          'columns': [
            {
              'name': 'drugind_id',
              'table': 'drug_indication',
              'schema': 'public',
              'type': 'bigint',
              'min': 22606,
              'max': 149873,
              'isUnique': true
            },
            {
              'name': 'record_id',
              'table': 'drug_indication',
              'schema': 'public',
              'type': 'bigint',
              'min': 1343055,
              'max': 3956801,
              'isUnique': false
            },
            {
              'name': 'molregno',
              'table': 'drug_indication',
              'schema': 'public',
              'type': 'bigint',
              'min': 97,
              'max': 2832777,
              'isUnique': false,
              'semanticType': 'molregno',
              'comment': 'Unique identifier for the molecule in the database.',
              'LLMComment': 'Molecule registration number used to identify specific compounds.'
            },
            {
              'name': 'max_phase_for_ind',
              'table': 'drug_indication',
              'schema': 'public',
              'type': 'numeric',
              'min': -1,
              'max': 4,
              'isUnique': false,
              'comment': 'The highest clinical trial phase reached for the drug indication, where -1 indicates no trials.',
              'LLMComment': 'Indicates the maximum clinical development phase achieved for the drug indication, ranging from 0 (preclinical) to 4 (marketed).'
            },
            {
              'name': 'mesh_id',
              'table': 'drug_indication',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Unique identifier for the Medical Subject Headings (MeSH) term associated with the drug indication.',
              'LLMComment': 'Identifier linking to the MeSH vocabulary for standardized medical terminology.'
            },
            {
              'name': 'mesh_heading',
              'table': 'drug_indication',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The descriptive term associated with the MeSH ID, representing the drug indication.',
              'LLMComment': 'The MeSH heading that describes the drug indication in standardized medical terms.'
            },
            {
              'name': 'efo_id',
              'table': 'drug_indication',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the Experimental Factor Ontology (EFO) term related to the drug indication.',
              'LLMComment': 'Unique identifier for the EFO term, which provides a structured vocabulary for describing experimental factors.'
            },
            {
              'name': 'efo_term',
              'table': 'drug_indication',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The term associated with the EFO ID, representing the drug indication in the context of experimental factors.',
              'LLMComment': 'Descriptive term from the EFO that relates to the drug indication, aiding in the categorization of experimental data.'
            }
          ],
          'rowCount': 55442,
          'comment': 'This table contains information about drug indications, linking drugs to their therapeutic uses and associated identifiers.',
          'LLMComment': 'The \'public.drug_indication\' table serves as a crucial resource in pharmacology and drug development, cataloging the various indications for which specific drugs are approved or investigated. Each entry is identified by a unique \'drugind_id\' and includes references to drug records (\'record_id\'), molecular registration numbers (\'molregno\'), and the maximum clinical trial phase for the indication (\'max_phase_for_ind\'). Additionally, it incorporates Medical Subject Headings (MeSH) and Experimental Factor Ontology (EFO) identifiers, which facilitate the classification and retrieval of drug-related information in biomedical research. This table is integral for understanding the therapeutic applications of drugs and their regulatory status, making it a key component in drug discovery and clinical research.'
        },
        {
          'name': 'drug_mechanism',
          'schema': 'public',
          'columns': [
            {
              'name': 'mec_id',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'bigint',
              'min': 13,
              'max': 9851,
              'isUnique': true,
              'comment': 'Unique identifier for the drug mechanism entry.'
            },
            {
              'name': 'record_id',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'bigint',
              'min': 1343055,
              'max': 3956801,
              'isUnique': false,
              'comment': 'Identifier for the record associated with this drug mechanism.'
            },
            {
              'name': 'molregno',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'bigint',
              'min': 115,
              'max': 2832777,
              'isUnique': false,
              'semanticType': 'molregno',
              'comment': 'Unique identifier for the molecule in the context of drug mechanisms.'
            },
            {
              'name': 'mechanism_of_action',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Description of how the drug exerts its effects at the molecular level.',
              'LLMComment': 'This field provides insights into the pharmacological action of the drug.'
            },
            {
              'name': 'tid',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 121145,
              'isUnique': false,
              'comment': 'Identifier for the target associated with the drug mechanism.'
            },
            {
              'name': 'site_id',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'bigint',
              'min': 248,
              'max': 38639,
              'isUnique': false,
              'comment': 'Identifier for the specific site related to the drug mechanism.'
            },
            {
              'name': 'action_type',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                'AGONIST',
                'POSITIVE MODULATOR',
                'BLOCKER',
                'NEGATIVE MODULATOR',
                'ACTIVATOR',
                'OTHER',
                'BINDING AGENT',
                'OXIDATIVE ENZYME',
                'ANTISENSE INHIBITOR',
                'EXOGENOUS PROTEIN',
                'POSITIVE ALLOSTERIC MODULATOR',
                'INVERSE AGONIST',
                'EXOGENOUS GENE',
                'PARTIAL AGONIST',
                'CHELATING AGENT',
                'NEGATIVE ALLOSTERIC MODULATOR',
                'STABILISER',
                'REDUCING AGENT',
                'VACCINE ANTIGEN',
                'SEQUESTERING AGENT',
                'DEGRADER',
                'DISRUPTING AGENT',
                'ALLOSTERIC ANTAGONIST',
                'CROSS-LINKING AGENT',
                'INHIBITOR',
                'OPENER',
                'SUBSTRATE',
                'MODULATOR',
                'PROTEOLYTIC ENZYME',
                'HYDROLYTIC ENZYME',
                'RELEASING AGENT',
                'RNAI INHIBITOR',
                'ANTAGONIST'
              ],
              'isUnique': false,
              'comment': 'Type of action the drug has on its target, such as AGONIST or BLOCKER.',
              'LLMComment': 'Categorizes the nature of the drug\'s interaction with its target.'
            },
            {
              'name': 'direct_interaction',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates whether there is a direct interaction with the target (1) or not (0).',
              'LLMComment': 'Binary flag to denote direct engagement with the target.'
            },
            {
              'name': 'molecular_mechanism',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the molecular mechanism is defined (1) or not (0).',
              'LLMComment': 'Binary flag to indicate the presence of a defined molecular mechanism.'
            },
            {
              'name': 'disease_efficacy',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the drug has demonstrated efficacy in treating a disease (1) or not (0).',
              'LLMComment': 'Binary flag to show if the drug is effective against a specific disease.'
            },
            {
              'name': 'mechanism_comment',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Additional comments or notes regarding the mechanism of action.',
              'LLMComment': 'Provides further context or details about the mechanism.'
            },
            {
              'name': 'selectivity_comment',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Comments on the selectivity of the drug for its target.',
              'LLMComment': 'Insights into how selective the drug is for its intended target.'
            },
            {
              'name': 'binding_site_comment',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Comments on the binding site characteristics or specifics.',
              'LLMComment': 'Details about the binding site where the drug interacts.'
            },
            {
              'name': 'variant_id',
              'table': 'drug_mechanism',
              'schema': 'public',
              'type': 'bigint',
              'min': -1,
              'max': 2608,
              'isUnique': false,
              'comment': 'Identifier for any variants related to the drug mechanism.',
              'LLMComment': 'Links to specific variants that may affect the drug\'s action.'
            }
          ],
          'rowCount': 7330,
          'comment': 'This table contains information about various drug mechanisms, detailing how different drugs interact with biological systems to exert their effects.',
          'LLMComment': 'The \'public.drug_mechanism\' table is a comprehensive repository of drug mechanisms, capturing critical data on how drugs operate at the molecular level. It includes fields such as \'mechanism_of_action\' which describes the specific biological processes affected by the drug, and \'action_type\' which categorizes the nature of the drug\'s interaction. This table is essential for understanding the pharmacological profiles of drugs, their efficacy against diseases, and their selectivity towards different biological targets. It serves as a foundational element for drug discovery and development, enabling researchers and AI systems to analyze and predict drug behavior in various therapeutic contexts.'
        },
        {
          'name': 'drug_warning',
          'schema': 'public',
          'columns': [
            {
              'name': 'warning_id',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 3844,
              'isUnique': true,
              'comment': 'Unique identifier for each warning record.'
            },
            {
              'name': 'record_id',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'bigint',
              'min': 1343058,
              'max': 3956801,
              'isUnique': false,
              'comment': 'Identifier for the specific record associated with the warning.'
            },
            {
              'name': 'molregno',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'bigint',
              'min': 146,
              'max': 2832776,
              'isUnique': false,
              'semanticType': 'molregno',
              'comment': 'Unique identifier for the molecule related to the warning, used in chemical databases.'
            },
            {
              'name': 'warning_type',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Black Box Warning',
                'Withdrawn'
              ],
              'isUnique': false,
              'comment': 'Type of warning issued, indicating the severity or nature of the warning.',
              'LLMComment': 'Categorizes the warning as either a Black Box Warning or Withdrawn.'
            },
            {
              'name': 'warning_class',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'immune system toxicity',
                'teratogenicity',
                'vascular toxicity',
                'musculoskeletal toxicity',
                'occular toxicity',
                'hepatotoxicity',
                'misuse',
                'null',
                'carcinogenicity',
                'respiratory toxicity',
                'cardiotoxicity',
                'hematological toxicity',
                'gastrointestinal toxicity',
                'metabolic toxicity',
                'dermatological toxicity',
                'neurotoxicity',
                'nephrotoxicity',
                'psychiatric toxicity',
                'infectious disease'
              ],
              'isUnique': false,
              'comment': 'Classification of the warning based on the type of toxicity or risk it represents.',
              'LLMComment': 'Describes the specific toxicity category such as immune system toxicity or teratogenicity.'
            },
            {
              'name': 'warning_description',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Detailed description of the warning, providing context and specifics.'
            },
            {
              'name': 'warning_country',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'semanticType': 'gis-country',
              'comment': 'Country where the warning is applicable or issued.',
              'LLMComment': 'Indicates the geographical relevance of the warning.'
            },
            {
              'name': 'warning_year',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'integer',
              'min': 1960,
              'max': 2020,
              'isUnique': false,
              'comment': 'Year in which the warning was issued or became relevant.',
              'LLMComment': 'Provides temporal context for the warning.'
            },
            {
              'name': 'efo_term',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Term from the Experimental Factor Ontology related to the warning.',
              'LLMComment': 'Links the warning to a specific term in the EFO.'
            },
            {
              'name': 'efo_id',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier from the Experimental Factor Ontology for the warning.',
              'LLMComment': 'Unique identifier for the EFO term associated with the warning.'
            },
            {
              'name': 'efo_id_for_warning_class',
              'table': 'drug_warning',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'EFO:0009880',
                'EFO:0011057',
                'null',
                'EFO:0011046',
                'EFO:0011056',
                'EFO:0005741',
                'EFO:0020928',
                'EFO:1001482',
                'EFO:0011051',
                'EFO:0011052',
                'EFO:0011055',
                'EFO:0011048',
                'EFO:0011062',
                'EFO:0011059',
                'EFO:0011050',
                'EFO:0011060',
                'EFO:0011049',
                'EFO:0011053',
                'EFO:0011054'
              ],
              'isUnique': false,
              'comment': 'EFO identifier corresponding to the warning class, linking to specific toxicity categories.',
              'LLMComment': 'Categorizes the warning class using EFO identifiers.'
            }
          ],
          'rowCount': 1676,
          'comment': 'This table stores information about drug warnings, including their types, classes, descriptions, and associated metadata such as country and year.',
          'LLMComment': 'The \'public.drug_warning\' table is a critical component of pharmacovigilance, containing 1,676 records of drug warnings. Each entry includes a unique warning ID, a reference to the record it pertains to, and details about the warning such as its type (e.g., side effects, contraindications), class (e.g., serious, moderate), and a description of the warning. Additionally, it captures the geographical context (country) and the year the warning was issued, along with relevant EFO (Experimental Factor Ontology) terms and IDs that link to broader biomedical concepts. This table is essential for understanding the safety profile of drugs and ensuring that healthcare professionals and patients are informed about potential risks associated with drug use.'
        },
        {
          'name': 'formulations',
          'schema': 'public',
          'columns': [
            {
              'name': 'product_id',
              'table': 'formulations',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the product associated with the formulation.',
              'LLMComment': 'A unique identifier representing the product to which this formulation belongs.'
            },
            {
              'name': 'ingredient',
              'table': 'formulations',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Name of the ingredient used in the formulation.',
              'LLMComment': 'The specific component or substance that is part of the formulation.'
            },
            {
              'name': 'strength',
              'table': 'formulations',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Concentration or potency of the ingredient in the formulation.',
              'LLMComment': 'Indicates the strength or dosage of the ingredient present in the formulation.'
            },
            {
              'name': 'record_id',
              'table': 'formulations',
              'schema': 'public',
              'type': 'bigint',
              'min': 1343055,
              'max': 3956801,
              'isUnique': false,
              'comment': 'Unique identifier for the record in the database.',
              'LLMComment': 'A unique numeric identifier for tracking the record within the database.'
            },
            {
              'name': 'molregno',
              'table': 'formulations',
              'schema': 'public',
              'type': 'bigint',
              'min': 115,
              'max': 2832776,
              'isUnique': false,
              'semanticType': 'molregno',
              'comment': 'Unique identifier for the chemical substance in the context of regulatory databases.',
              'LLMComment': 'A unique identifier that corresponds to the chemical substance in regulatory frameworks.'
            },
            {
              'name': 'formulation_id',
              'table': 'formulations',
              'schema': 'public',
              'type': 'bigint',
              'min': 347,
              'max': 201503,
              'isUnique': true,
              'comment': 'Unique identifier for the formulation entry.',
              'LLMComment': 'A unique numeric identifier for the specific formulation record.'
            }
          ],
          'rowCount': 52755,
          'comment': 'This table contains detailed information about various pharmaceutical formulations, including their ingredients and strengths, linked to specific products and formulations.',
          'LLMComment': 'The \'public.formulations\' table is a critical component of the pharmaceutical database, serving as a repository for information on different drug formulations. Each entry in the table represents a unique formulation identified by \'formulation_id\', detailing the \'product_id\' associated with the formulation, the \'ingredient\' used, and its \'strength\'. The \'record_id\' and \'molregno\' provide additional identifiers for tracking and referencing these formulations within the broader context of drug development and regulatory compliance. This table is essential for understanding the composition of pharmaceutical products and is interconnected with various other tables that provide insights into drug mechanisms, indications, and classifications.'
        },
        {
          'name': 'frac_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'frac_class_id',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 470,
              'max': 686,
              'isUnique': true,
              'comment': 'Unique identifier for each fraction classification.'
            },
            {
              'name': 'active_ingredient',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'The primary chemical component responsible for the biological activity.'
            },
            {
              'name': 'level1',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'B',
                'U',
                'E',
                'I',
                'D',
                'A',
                'C',
                'P',
                'H',
                'NC',
                'M',
                'G',
                'F'
              ],
              'isUnique': false,
              'comment': 'Top-level classification category for the active ingredient.',
              'LLMComment': 'Represents the broadest classification of the active ingredient, indicating its primary mode of action.'
            },
            {
              'name': 'level1_description',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'RESPIRATION',
                'MITOSIS AND CELL DIVISION',
                'NOT CLASSIFIED',
                'NUCLEIC ACID SYNTHESIS',
                'SIGNAL TRANSDUCTION',
                'HOST PLANT DEFENCE INDUCTION',
                'MULTI-SITE CONTACT ACTIVITY',
                'STEROL BIOSYNTHESIS IN MEMBRANES',
                'MELANIN SYNTHESIS IN CELL WALL',
                'AMINO ACIDS AND PROTEIN SYNTHESIS',
                'CELL WALL BIOSYNTHESIS',
                'UNKNOWN MODE OF ACTION',
                'LIPID SYNTHESIS AND MEMBRANE INTEGRITY'
              ],
              'isUnique': false,
              'comment': 'Description of the top-level classification category.',
              'LLMComment': 'Provides context for the level1 category, detailing the biological processes affected.'
            },
            {
              'name': 'level2',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'C1',
                'P3',
                'C5',
                'P1',
                'F6',
                'I1',
                'B1',
                'C4',
                'H3',
                'P5',
                'D4',
                'E3',
                'C6',
                'A1',
                'H5',
                'B4',
                'D5',
                'D3',
                'M1',
                'A3',
                'C8',
                'C2',
                'U2',
                'B5',
                'G2',
                'F4',
                'A4',
                'F7',
                'U1',
                'P2',
                'E2',
                'G3',
                'F2',
                'P4',
                'U3',
                'H4',
                'E1',
                'F3',
                'C7',
                'B3',
                'D2',
                'D1',
                'C3',
                'B2',
                'G4',
                'I2',
                'A2',
                'G1'
              ],
              'isUnique': false,
              'comment': 'Second-level classification category for the active ingredient.',
              'LLMComment': 'Refines the classification further, indicating specific mechanisms or targets.'
            },
            {
              'name': 'level2_description',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'ACTIN DISRUPTION (PROPOSED)',
                'DNA TOPOISOMERASE TYPE II (GYRASE)',
                'DELTA14-REDUCTASE AND (DELTA8 TO DELTA7)-ISOMERASE IN STEROL BIOSYNTHESIS (ERG24, ERG2)',
                'SALICYLIC ACID PATHWAY',
                'UNCOUPLERS OF OXIDATIVE PHOSPHORYLATION',
                'BETA-TUBULINE ASSEMBLY IN MITOSIS',
                'null',
                'MAP/HISTIDINE-KINASE IN OSMOTIC SIGNAL TRANSDUCTION (OS-2, HOG1)',
                'DNA/RNA SYNTHESIS (PROPOSED)',
                'SIGNAL TRANSDUCTION (MECHANISM UNKNOWN)',
                'PROTEIN SYNTHESIS',
                'C14-DEMETHYLASE IN STEROL BIOSYNTHESIS (ERG11/CYP51)',
                'CELL DIVISION (PROPOSED)',
                'COMPLEX III: CYTOCHROME BC1 (UBIQUINONE REDUCTASE) AT QX (UNKNOWN) SITE',
                'MULTI-SITE CONTACT ACTIVITY',
                'DELOCALISATION OF SPECTRIN-LIKE PROTEINS',
                'RNA POLYMERASE I',
                'METHIONINE BIOSYNTHESIS (PROPOSED) (CGS GENE)',
                'MAP/HISTIDINE-KINASE IN OSMOTIC SIGNAL TRANSDUCTION (OS-1, DAF1)',
                'PHOSPHOLIPID BIOSYNTHESIS, METHYLTRANSFERASE',
                'DEHYDRATASE IN MELANIN BIOSYNTHESIS',
                'CELL MEMBRANE DISRUPTION (PROPOSED)',
                'COMPLEX II: SUCCINATE-DEHYDROGENASE',
                'CELLULOSE SYNTHASE',
                'MICROBIAL DISRUPTERS OF PATHOGEN CELL MEMBRANES',
                'ATP PRODUCTION (PROPOSED)',
                'ADENOSIN-DEAMINASE',
                'COMPLEX III: CYTOCHROME BC1 (UBIQUINOL OXIDASE) AT QO SITE (CYT B GENE)',
                'LIPID PEROXIDATION (PROPOSED)',
                'COMPLEX III: CYTOCHROME BC1(UBIQUINONE REDUCTASE) AT QI SITE',
                'CHITIN SYNTHASE',
                '3-KETO REDUCTASE, C4-DEMETHYLATION (ERG27)',
                'REDUCTASE IN MELANIN BIOSYNTHESIS',
                'COMPLEX I NADH OXIDO-REDUCTASE',
                'TREHALASE AND INOSITOL-BIOSYNTHESIS',
                'SQUALENE-EPOXIDASE IN STEROL BIOSYNTHESIS (ERG1)',
                'INHIBITORS OF OXIDATIVE PHOSPHORYLATION, ATP SYNTHASE',
                'UNKNOWN',
                'CELL MEMBRANE PERMEABILITY, FATTY ACIDS (PROPOSED)'
              ],
              'isUnique': false,
              'comment': 'Description of the second-level classification category.',
              'LLMComment': 'Elaborates on the level2 category, specifying the biological targets or pathways involved.'
            },
            {
              'name': 'level3',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'M1M3',
                'F428',
                'U133',
                'D541',
                'U2U8',
                'U3U12',
                'G317',
                'M1M1',
                'C139',
                'G25',
                'D425',
                'I116.1',
                'H540',
                'M1M6',
                'G418',
                'C529',
                'U142',
                'C630',
                'M1M4',
                'U127',
                'G13',
                'E212',
                'B322',
                'M1M11',
                'D19',
                'C738',
                'H326',
                'U134',
                'B210',
                'D223',
                'P5P',
                'P4P',
                'F644',
                'U1NC',
                'H419',
                'U1U13',
                'P3P',
                'A431',
                'E32',
                'F26',
                'C421',
                'M1M5',
                'M1M2',
                'D324',
                'U135',
                'U137',
                'B543',
                'P2P',
                'E113',
                'U136',
                'M1M8',
                'A28',
                'M1M9',
                'A14',
                'C311',
                'I216.2',
                'F746',
                'F314',
                'B420',
                'P1P',
                'A332',
                'U1U6',
                'B11',
                'C845',
                'M1M7',
                'U1U14',
                'C27',
                'M1M10'
              ],
              'isUnique': false,
              'comment': 'Third-level classification category for the active ingredient.',
              'LLMComment': 'Further narrows down the classification, often indicating specific compounds or actions.'
            },
            {
              'name': 'level3_description',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'GLUCOPYRANOSYL ANTIBIOTIC',
                'DICARBOXIMIDES',
                'null',
                '(SBI: CLASS III)',
                'BENZOTRIAZINES',
                'BENZENE-SULPHONAMIDES',
                'POLYOXINS',
                'PYRIMIDINAMINES',
                'CARBOXYLIC ACIDS',
                'QUINOXALINES',
                'INORGANIC',
                'HETEROAROMATICS',
                'QXI - FUNGICIDE (QUINONE X INHIBITOR)',
                'AP - FUNGICIDES (ANILINO-PYRIMIDINES)',
                'ARYL-PHENYL-KETONE',
                'AZA-NAPHTHALENES',
                'THIOCARBAMATE',
                'CHLORONITRILES (PHTHALONITRILES)',
                'HYDROXY-(2-AMINO-)PYRIMIDINES',
                'PLANT EXTRACT',
                'MBC - FUNGICIDES (METHYL BENZIMIDAZOLE CARBAMATES)',
                'N-PHENYL CARBAMATES',
                'SULFAMIDES',
                'PHTHALAMIC ACIDS',
                'SDHI (SUCCINATE DEHYDROGENASE INHIBITORS)',
                'DMI-FUNGICIDES (DEMETHYLATION INHIBITORS) (SBI: CLASS I)',
                'THIOPHENE-CARBOXAMIDES',
                'PP-FUNGICIDES (PHENYLPYRROLES)',
                'BENZAMIDES',
                'HEXOPYRANOSYL ANTIBIOTIC',
                'DITHIOCARBAMATES AND RELATIVES',
                'QII - FUNGICIDES (QUINONE INSIDE INHIBITORS)',
                'THIAZOLIDINE',
                'PYRIDAZINONES',
                'THIAZOLE CARBOXAMIDE',
                'MBI-R (MELANIN BIOSYNTHESIS INHIBITORS - REDUCTASE)',
                'BENZO-THIADIAZOLE BTH',
                'PA - FUNGICIDES (PHENYLAMIDES)',
                'MBI-D (MELANIN BIOSYNTHESIS INHIBITORS - DEHYDRATASE)',
                'CYANOACETAMIDE-OXIME',
                'NATURAL COMPOUND',
                'PHOSPHORO-THIOLATES',
                'MALEIMIDE',
                'PYRIMIDINONE-HYDRAZONES',
                'TRIAZINES',
                'PHOSPHONATES',
                'PHENYL-ACETAMIDE',
                'DIVERSE',
                'CARBAMATES',
                '(SBI CLASS IV)',
                'BENZISOTHIAZOLE',
                'CAA-FUNGICIDES (CARBOXYLIC ACID AMIDES)',
                'GUANIDINES',
                'MICROBIAL (BACILLUS SP.)',
                'ORGANO TIN COMPOUNDS',
                'THIADIAZOLE-CARBOXAMIDE',
                'PHTHALIMIDES',
                'QUINONES (ANTHRAQUINONES)',
                'AH-FUNGICIDES (AROMATIC HYDROCARBONS) (CHLOROPHENYLS, NITROANILINES)',
                'TETRACYCLINE ANTIBIOTIC',
                'ENOPYRANURONIC ACID ANTIBIOTIC',
                'AMINES ("MORPHOLINES") (SBI: CLASS II)',
                'PHENYLUREAS',
                'DITHIOLANES',
                'PYRAZOLE-MET1',
                'QOI-FUNGICIDES (QUINONE OUTSIDE INHIBITORS)'
              ],
              'isUnique': false,
              'comment': 'Description of the third-level classification category.',
              'LLMComment': 'Details the specific nature of the compounds or actions represented in level3.'
            },
            {
              'name': 'level4',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Fourth-level classification category for the active ingredient.',
              'LLMComment': 'Additional classification that may provide more granularity in categorization.'
            },
            {
              'name': 'level4_description',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Description of the fourth-level classification category.',
              'LLMComment': 'Offers insights into the specifics of the level4 classification.'
            },
            {
              'name': 'level5',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Fifth-level classification category for the active ingredient.',
              'LLMComment': 'The most detailed classification, often indicating specific subcategories.'
            },
            {
              'name': 'frac_code',
              'table': 'frac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'M1',
                'M5',
                '13',
                '28',
                'P',
                '39',
                '19',
                '23',
                '20',
                '14',
                '25',
                '34',
                '6',
                '18',
                '2',
                '38',
                '32',
                'U14',
                'M3',
                '37',
                '9',
                '4',
                'M6',
                'M4',
                '33',
                'M9',
                '16.1',
                '30',
                '45',
                '16.2',
                '46',
                '36',
                '27',
                '11',
                '31',
                'U8',
                '43',
                '42',
                '12',
                'M8',
                '10',
                '7',
                '5',
                'U13',
                '40',
                'M2',
                'U12',
                '44',
                'U6',
                'M7',
                'M10',
                '29',
                'NC',
                '1',
                '21',
                '24',
                '17',
                '35',
                '41',
                '22',
                '3',
                '26',
                'M11',
                '8'
              ],
              'isUnique': false,
              'comment': 'Code representing the fraction classification, often used for quick reference.',
              'LLMComment': 'A shorthand notation for the classification, useful for identification and categorization.'
            }
          ],
          'rowCount': 217,
          'comment': 'This table contains classifications for various fractions used in pharmacological contexts, detailing their active ingredients and hierarchical classification levels.',
          'LLMComment': 'The `public.frac_classification` table serves as a comprehensive reference for categorizing different fractions in the pharmaceutical domain. It includes a unique identifier for each fraction (`frac_class_id`), the active ingredient involved, and a multi-level classification system (from level 1 to level 5) that provides detailed descriptions at each level. This hierarchical structure allows for nuanced understanding and organization of fractions, which is crucial for drug development, regulatory compliance, and research purposes. The table is interconnected with various other tables in the schema, facilitating a rich dataset for analysis and insights into drug mechanisms, indications, and classifications.'
        },
        {
          'name': 'go_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'go_id',
              'table': 'go_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true
            },
            {
              'name': 'parent_go_id',
              'table': 'go_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                'GO:0004872',
                'GO:0044281',
                'GO:0008233',
                'GO:0009636',
                'GO:0004104',
                'GO:0005694',
                'GO:0005102',
                'GO:0016485',
                'GO:0048856',
                'GO:0009058',
                'GO:0005575',
                'GO:0006810',
                'GO:0005198',
                'GO:0030312',
                'GO:0016772',
                'GO:0042802',
                'GO:0000166',
                'GO:0042562',
                'GO:0004620',
                'GO:0005125',
                'GO:0019838',
                'GO:0016740',
                'GO:0051400',
                'GO:0008150',
                'GO:0003824',
                'GO:0016788',
                'GO:0003674',
                'GO:0042277',
                'GO:0019904',
                'GO:0043226',
                'GO:0031638',
                'GO:0016020',
                'GO:0030246',
                'GO:0005515',
                'GO:0005520',
                'GO:0005773',
                'GO:0005085',
                'GO:0008146',
                'GO:0065003',
                'GO:0032403',
                'GO:0008152',
                'GO:0061024',
                'GO:0009308',
                'GO:0017111',
                'GO:0005215',
                'GO:0050790',
                'GO:0016787',
                'GO:0019835',
                'GO:0034641',
                'GO:0003723',
                'GO:0016298',
                'GO:0006508',
                'GO:0003676',
                'GO:0003823',
                'GO:0022607',
                'GO:0009056',
                'GO:0005179',
                'GO:0004623',
                'GO:0040011'
              ],
              'isUnique': false,
              'comment': 'Identifier for the parent Gene Ontology term, indicating hierarchical relationships.',
              'LLMComment': 'This column represents the parent term in the Gene Ontology hierarchy, linking to broader categories.'
            },
            {
              'name': 'pref_name',
              'table': 'go_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Preferred name for the Gene Ontology term, used for display purposes.',
              'LLMComment': 'The commonly used name for the Gene Ontology term, which may be more user-friendly than the ID.'
            },
            {
              'name': 'class_level',
              'table': 'go_classification',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 7,
              'isUnique': false,
              'comment': 'Numerical representation of the classification level within the Gene Ontology hierarchy, where lower numbers indicate higher-level categories.',
              'LLMComment': 'This small integer indicates the depth of the term in the Gene Ontology hierarchy, with 0 being the most general.'
            },
            {
              'name': 'aspect',
              'table': 'go_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'C',
                'P',
                'F'
              ],
              'isUnique': false,
              'comment': 'Indicates the aspect of the biological process represented by the term: Cellular Component (C), Biological Process (P), or Molecular Function (F).',
              'LLMComment': 'This character categorizes the Gene Ontology term into one of three aspects: C for Cellular Component, P for Biological Process, and F for Molecular Function.'
            },
            {
              'name': 'path',
              'table': 'go_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Unique path representation of the Gene Ontology term\'s location in the hierarchy, often used for traversal or visualization.',
              'LLMComment': 'This unique string describes the path to the term within the Gene Ontology structure, useful for understanding its context.'
            }
          ],
          'rowCount': 309,
          'comment': 'This table contains information about Gene Ontology (GO) classifications, including unique identifiers, parent-child relationships, preferred names, classification levels, aspects, and paths.',
          'LLMComment': 'The `public.go_classification` table is a crucial component of the biological ontology domain, specifically focusing on Gene Ontology (GO) classifications. It provides a structured representation of biological concepts, where each entry is identified by a unique `go_id` and may have a hierarchical relationship with other GO terms through `parent_go_id`. The `pref_name` column offers a human-readable name for each GO term, while `class_level` indicates the depth of the term in the ontology hierarchy. The `aspect` column categorizes the GO terms into three main aspects: biological process, molecular function, and cellular component, and the `path` column outlines the lineage of the term within the ontology. This table is essential for researchers and AI systems that analyze biological data, enabling them to understand and utilize the relationships and classifications of genes and proteins effectively.'
        },
        {
          'name': 'hrac_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'hrac_class_id',
              'table': 'hrac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 273,
              'isUnique': true,
              'comment': 'Unique identifier for each HRAC classification.'
            },
            {
              'name': 'active_ingredient',
              'table': 'hrac_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The chemical compound that acts as the active substance in the classification.'
            },
            {
              'name': 'level1',
              'table': 'hrac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'B',
                'Z',
                'K1',
                'C3',
                'C2',
                'I',
                'D',
                'K2',
                'N',
                'P',
                'H',
                'G',
                'C1',
                'E',
                'F2',
                'A',
                'F1',
                'O',
                'M',
                'F3',
                'L',
                'K3'
              ],
              'isUnique': false,
              'comment': 'The primary classification level indicating the mode of action of the active ingredient.',
              'LLMComment': 'Represents the broadest category of herbicide action, such as synthetic auxins or inhibitors.'
            },
            {
              'name': 'level1_description',
              'table': 'hrac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'ACTION LIKE INDOLE ACETIC ACID (SYNTHETIC AUXINS)',
                'MICROTUBULE ASSEMBLY INHIBITION',
                'BLEACHING: INHIBITION OF CAROTENOID BIOSYNTHESIS AT THE PHYTOENE DESATURASE STEP (PDS)',
                'INHIBITION OF LIPID SYNTHESIS - NOT ACCASE INHIBITION',
                'INHIBITION OF GLUTAMINE SYNTHETASE',
                'BLEACHING: INHIBITION OF 4-HYDROXYPHENYL-PYRUVATE-DIOXYGENASE (4-HPPD)',
                'UNCOUPLING (MEMBRANE DISRUPTION)',
                'INHIBITION OF ACETYL COA CARBOXYLASE (ACCASE)',
                'INHIBITION OF AUXIN TRANSPORT',
                'INHIBITION OF PROTOPORPHYRINOGEN OXIDASE (PPO)',
                'INHIBITION OF ACETOLACTATE SYNTHASE ALS (ACETOHYDROXYACID SYNTHASE AHAS)',
                'INHIBITION OF DHP (DIHYDROPTEROATE) SYNTHASE',
                'PHOTOSYSTEM-I-ELECTRON DIVERSION',
                'INHIBITION OF CELL WALL (CELLULOSE) SYNTHESIS',
                'UNKNOWN',
                'INHIBITION OF MITOSIS / MICROTUBULE ORGANISATION',
                'INHIBITION OF PHOTOSYNTHESIS AT PHOTOSYSTEM II',
                'BLEACHING: INHIBITION OF CAROTENOID BIOSYNTHESIS (UNKNOWN TARGET)',
                'INHIBITION OF VLCFAS (INHIBITION OF CELL DIVISION)',
                'INHIBITION OF EPSP SYNTHASE'
              ],
              'isUnique': false,
              'comment': 'Detailed description of the action associated with level1 classification.',
              'LLMComment': 'Describes the specific biochemical processes affected by the active ingredient, such as inhibition of carotenoid biosynthesis.'
            },
            {
              'name': 'level2',
              'table': 'hrac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'K11',
                'C16',
                'Z2',
                'K12',
                'C33',
                'C11',
                'O2',
                'C32',
                'E4',
                'H1',
                'E8',
                'K36',
                'E3',
                'F24',
                'A1',
                'M1',
                'E9',
                'E6',
                'E5',
                'C14',
                'C13',
                'N2',
                'F21',
                'L4',
                'K15',
                'F23',
                'C21',
                'O5',
                'L2',
                'B3',
                'C15',
                'F12',
                'B2',
                'Z4',
                'K14',
                'G1',
                'A2',
                'N1',
                'O4',
                'K31',
                'C31',
                'P1',
                'K33',
                'I1',
                'B1',
                'O1',
                'N4',
                'K34',
                'F32',
                'O3',
                'B4',
                'K21',
                'L1',
                'A3',
                'K13',
                'F22',
                'B5',
                'K35',
                'F34',
                'F13',
                'K32',
                'C12',
                'E2',
                'N3',
                'C22',
                'Z1',
                'F33',
                'F11',
                'Z3',
                'F31',
                'E1',
                'D1',
                'L3',
                'E7'
              ],
              'isUnique': false,
              'comment': 'The secondary classification level providing more specific categorization.',
              'LLMComment': 'Refines the classification further, indicating a more specific mode of action.'
            },
            {
              'name': 'level2_description',
              'table': 'hrac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'ARYLAMINOPROPIONIC ACID',
                'OXYACETAMIDE',
                'null',
                'PHENYL-CARBAMATE',
                'TRIAZOLINONE',
                'TRIAZOLOCARBOXAMIDE',
                'BENZOIC ACID',
                'PYRAZOLE',
                'OTHER',
                'DIPHENYLETHER',
                'TRIAZOLOPYRIMIDINE',
                'PHENOXY-CARBOXYLIC-ACID',
                'PYRIDINECARBOXAMIDE',
                'GLYCINE',
                'IMIDAZOLINONE',
                'TRIAZOLE',
                'ISOXAZOLIDINONE',
                'ISOXAZOLE',
                'OXADIAZOLE',
                'THIADIAZOLE',
                'PHTHALAMATE SEMICARBAZONE',
                'TRIAZINE',
                'THIOCARBAMATE',
                'SULFONYLUREA',
                'BENZOFURAN',
                'TETRAZOLINONE',
                'PHENYLPYRAZOLINE \'DEN\'',
                'TRIAZINONE',
                'PYRIDAZINONE',
                'URACIL',
                'NITRILE',
                'CHLOROACETAMIDE',
                'OXAZOLIDINEDIONE',
                'PYRIDINE',
                'SULFONYLAMINOCARBONYL-TRIAZOLINONE',
                'BENZAMIDE',
                'DINITROANILINE',
                'N-PHENYLPHTHALIMIDE',
                'BENZOTHIADIAZINONE',
                'UREA',
                'DINITROPHENOL',
                'PHENYLPYRAZOLE',
                'AMIDE',
                'PYRIDINE CARBOXYLIC ACID',
                'PHOSPHINIC ACID',
                'BIPYRIDYLIUM',
                'CARBAMATE',
                'CHLORO-CARBONIC-ACID',
                'CYCLOHEXANEDIONE \'DIMS\'',
                'QUINOLINE CARBOXYLIC ACID',
                'ORGANOARSENICAL',
                'PYRIMIDINYL(THIO)BENZOATE',
                'PYRAZOLIUM',
                'ARYLOXYPHENOXY-PROPIONATE \'FOPS\'',
                'ACETAMIDE',
                'PHOSPHORODITHIOATE',
                'TRIKETONE',
                'PHOSPHOROAMIDATE',
                'PHENYL-PYRIDAZINE',
                'PYRIMIDINDIONE'
              ],
              'isUnique': false,
              'comment': 'Detailed description of the action associated with level2 classification.',
              'LLMComment': 'Provides insight into the specific chemical nature or mechanism of the active ingredient.'
            },
            {
              'name': 'level3',
              'table': 'hrac_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'The tertiary classification level, providing the most specific categorization.',
              'LLMComment': 'Indicates a unique classification that may relate to specific chemical structures or actions.'
            },
            {
              'name': 'hrac_code',
              'table': 'hrac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'B',
                'Z',
                'K1',
                'C3',
                'C2',
                'I',
                'D',
                'K2',
                'N',
                'P',
                'H',
                'G',
                'C1',
                'E',
                'F2',
                'A',
                'F1',
                'O',
                'M',
                'F3',
                'L',
                'K3'
              ],
              'isUnique': false,
              'comment': 'A short code representing the HRAC classification.',
              'LLMComment': 'A concise alphanumeric code that summarizes the classification for easy reference.'
            }
          ],
          'rowCount': 273,
          'comment': 'This table contains classifications of herbicides based on their active ingredients and their corresponding HRAC (Herbicide Resistance Action Committee) codes.',
          'LLMComment': 'The `public.hrac_classification` table is essential for understanding the classification of herbicides in agricultural science. It categorizes herbicides by their active ingredients and assigns them HRAC codes, which are crucial for identifying their mechanisms of action and resistance management. The table includes hierarchical levels of classification (level1, level2, level3) along with descriptions, providing a structured way to analyze herbicide properties and their applications in pest management.'
        },
        {
          'name': 'indication_refs',
          'schema': 'public',
          'columns': [
            {
              'name': 'indref_id',
              'table': 'indication_refs',
              'schema': 'public',
              'type': 'bigint',
              'min': 593841,
              'max': 682330,
              'isUnique': true,
              'comment': 'Unique identifier for each indication reference.'
            },
            {
              'name': 'drugind_id',
              'table': 'indication_refs',
              'schema': 'public',
              'type': 'bigint',
              'min': 22606,
              'max': 149873,
              'isUnique': false,
              'comment': 'Identifier linking to the specific drug indication.'
            },
            {
              'name': 'ref_type',
              'table': 'indication_refs',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'USAN',
                'ClinicalTrials',
                'EMA',
                'INN',
                'DailyMed',
                'ATC',
                'FDA'
              ],
              'isUnique': false,
              'comment': 'Type of reference, indicating the source or category of the information.',
              'LLMComment': 'Categorizes the reference source, such as USAN, ClinicalTrials, EMA, INN, DailyMed, etc.'
            },
            {
              'name': 'ref_id',
              'table': 'indication_refs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the specific reference within its source.',
              'LLMComment': 'Unique ID associated with the reference in its respective database or source.'
            },
            {
              'name': 'ref_url',
              'table': 'indication_refs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'semanticType': 'URL',
              'comment': 'URL linking to the reference source for more information.',
              'LLMComment': 'Web address providing access to the detailed reference information.'
            }
          ],
          'rowCount': 88432,
          'comment': 'This table stores references related to drug indications, including their types and associated URLs.',
          'LLMComment': 'The \'indication_refs\' table is a crucial component of the pharmacological database, containing 88,432 entries that link drug indications to various reference types. Each entry includes a unique identifier for the reference (\'indref_id\'), the associated drug indication (\'drugind_id\'), the type of reference (\'ref_type\'), a unique reference identifier (\'ref_id\'), and a URL for further information (\'ref_url\'). This structure allows for the organization and retrieval of external resources and documentation related to drug indications, enhancing the understanding of drug uses and their supporting literature.'
        },
        {
          'name': 'irac_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'irac_class_id',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 251,
              'isUnique': true,
              'comment': 'Unique identifier for each IRAC classification entry.'
            },
            {
              'name': 'active_ingredient',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'The primary chemical compound that exhibits biological activity.'
            },
            {
              'name': 'level1',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'B',
                'U',
                'E',
                'D',
                'A',
                'C',
                'M'
              ],
              'isUnique': false,
              'comment': 'The first level of classification indicating broad categories of action.',
              'LLMComment': 'Represents the highest level of classification, categorizing active ingredients into major functional groups.'
            },
            {
              'name': 'level1_description',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'ENERGY METABOLISM',
                'GROWTH REGULATION',
                'UNKNOWN',
                'LIPID SYNTHESIS, GROWTH REGULATION',
                'NERVE ACTION',
                'MISCELLANEOUS',
                'NERVE AND MUSCLE ACTION'
              ],
              'isUnique': false,
              'comment': 'Description of the first level category, providing context for the classification.',
              'LLMComment': 'Describes the specific biological function or mechanism associated with the level1 category.'
            },
            {
              'name': 'level2',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'C16',
                'C17',
                'B28',
                'D25',
                'A14',
                'C7',
                'C15',
                'A1',
                'M8',
                'UUN',
                'E23',
                'C18',
                'B6',
                'A3',
                'A9',
                'U11',
                'D13',
                'A19',
                'A22',
                'D21',
                'C10',
                'D20',
                'A4',
                'D24',
                'D12',
                'A5',
                'A2'
              ],
              'isUnique': false,
              'comment': 'The second level of classification providing more specific categories.',
              'LLMComment': 'Refines the classification further into more specific action categories.'
            },
            {
              'name': 'level2_description',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'CHLORIDE CHANNEL ACTIVATORS',
                'RYANODINE RECEPTOR MODULATORS',
                'VOLTAGE-DEPENDENT SODIUM CHANNEL BLOCKERS',
                'ECDYSONE RECEPTOR AGONISTS',
                'JUVENILE HORMONE MIMICS',
                'NICOTINIC ACETYLCHOLINE RECEPTOR (NACHR) AGONISTS',
                'OCTOPAMINE RECEPTOR AGONISTS',
                'MITOCHONDRIAL COMPLEX IV ELECTRON TRANSPORT INHIBITORS',
                'MITOCHONDRIAL COMPLEX III ELECTRON TRANSPORT INHIBITORS',
                'UNCOUPLERS OF OXIDATIVE PHOSPHORYLATION VIA DISRUPTION OF THE PROTON GRADIENT',
                'MOULTING DISRUPTOR, DIPTERAN',
                'MITE GROWTH INHIBITORS',
                'GABA-GATED CHLORIDE CHANNEL ANTAGONISTS',
                'COMPOUNDS OF UNKNOWN OR UNCERTAIN MOA',
                'INHIBITORS OF MITOCHONDRIAL ATP SYNTHASE',
                'NICOTINIC ACETYLCHOLINE RECEPTOR (NACHR) ALLOSTERIC ACTIVATORS',
                'MICROBIAL DISRUPTORS OF INSECT MIDGUT MEMBRANES',
                'INHIBITORS OF CHITIN BIOSYNTHESIS, TYPE 1',
                'INHIBITORS OF CHITIN BIOSYNTHESIS, TYPE 0',
                'MITOCHONDRIAL COMPLEX II ELECTRON TRANSPORT INHIBITORS',
                'MITOCHONDRIAL COMPLEX I ELECTRON TRANSPORT INHIBITORS',
                'NICOTINIC ACETYLCHOLINE RECEPTOR (NACHR) CHANNEL BLOCKERS',
                'MISCELLANEOUS NONSPECIFIC (MULTI-SITE) INHIBITORS',
                'SODIUM CHANNEL MODULATORS',
                'SELECTIVE HOMOPTERAN FEEDING BLOCKERS',
                'INHIBITORS OF ACETYL COA CARBOXYLASE',
                'ACETYLCHOLINESTERASE (ACHE) INHIBITORS'
              ],
              'isUnique': false,
              'comment': 'Description of the second level category, detailing the specific action type.',
              'LLMComment': 'Provides detailed information about the specific mechanisms of action for the level2 category.'
            },
            {
              'name': 'level3',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'A22B',
                'A11B',
                'UUNUNG',
                'D2424A',
                'A44C',
                'D2020C',
                'C1818A',
                'UUNUNH',
                'M88A',
                'C1515A',
                'UUNUNB',
                'B2828A',
                'E2323A',
                'M88E',
                'A2222A',
                'C1010B',
                'C77B',
                'A44B',
                'C1616A',
                'D2525A',
                'D1212B',
                'A99B',
                'A44A',
                'D2424B',
                'D1212A',
                'A55A',
                'D2020A',
                'C77A',
                'B66A',
                'A33A',
                'C1010A',
                'UUNUNC',
                'U1111A',
                'M88B',
                'D1212D',
                'A11A',
                'M88C',
                'UUNUNA',
                'A1919A',
                'A1414A',
                'C1717A',
                'D2020B',
                'D1212C',
                'A33B',
                'D2121A',
                'M88D',
                'A2222B',
                'A22A',
                'U1111B',
                'D2121B',
                'D1313A',
                'A99C',
                'C77C',
                'UUNUNF',
                'UUNUNE',
                'UUNUND'
              ],
              'isUnique': false,
              'comment': 'The third level of classification, offering even more specificity.',
              'LLMComment': 'Indicates a more granular classification of the active ingredient\'s action.'
            },
            {
              'name': 'level3_description',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'CHLORFENAPYR, DNOC, SULFLURAMID',
                'SULFURYL FLUORIDE',
                'BACILLUS SPHAERICUS',
                'PYRIFLUQUINAZON',
                'CLOFENTEZINE, HEXYTHIAZOX, DIFLOVIDAZIN',
                'ACEQUINOCYL',
                'ETOXAZOLE',
                'FENOXYCARB',
                'PROPARGITE',
                'PYMETROZINE',
                'ALKYL HALIDES',
                'FLONICAMID',
                'BENZOXIMATE',
                'DDT, METHOXYCHLOR',
                'NICOTINE',
                'NEONICOTINOIDS',
                'DICOFOL',
                'AZADIRACHTIN',
                'TETRONIC AND TETRAMIC ACID DERIVATIVES',
                'BUPROFEZIN',
                'CYCLODIENE ORGANOCHLORINES',
                'BETA-KETONITRILE DERIVATIVES',
                'NEREISTOXIN ANALOGUES',
                'TARTAR EMETIC',
                'TETRADIFON',
                'JUVENILE HORMONE ANALOGUES',
                'SULFOXAFLOR',
                'CRYOLITE',
                'DIAMIDES',
                'FLUACRYPYRIM',
                'ROTENONE',
                'BENZOYLUREAS',
                'METAFLUMIZONE',
                'HYDRAMETHYLNON',
                'PYRETHROIDS, PYRETHRINS',
                'METI ACARICIDES AND INSECTICIDES',
                'ORGANOTIN MITICIDES',
                'PHOSPHINE',
                'DIAFENTHIURON',
                'CHLOROPICRIN',
                'BIFENAZATE',
                'PYRIPROXYFEN',
                'CYROMAZINE',
                'PHENYLPYRAZOLES (FIPROLES)',
                'PYRIDALYL',
                'CHINOMETHIONAT',
                'AVERMECTINS, MILBEMYCINS',
                'SPINOSYNS',
                'AMITRAZ',
                'ORGANOPHOSPHATES',
                'CARBAMATES',
                'CYANIDE',
                'BACILLUS THURINGIENSIS AND THE INSECTICIDAL PROTEINS THEY PRODUCE',
                'DIACYLHYDRAZINES',
                'BORAX',
                'INDOXACARB'
              ],
              'isUnique': false,
              'comment': 'Description of the third level category, identifying specific compounds or actions.',
              'LLMComment': 'Lists specific active ingredients or compounds that fall under the level3 classification.'
            },
            {
              'name': 'level4',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'The fourth level of classification, unique to each entry.',
              'LLMComment': 'Provides the most detailed classification level, often used for regulatory or identification purposes.'
            },
            {
              'name': 'irac_code',
              'table': 'irac_classification',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                '2A',
                '9C',
                '21A',
                '12B',
                '11A',
                '1B',
                '3B',
                '8C',
                '16',
                '8B',
                '20B',
                '4B',
                '7C',
                '10B',
                '12C',
                '4C',
                '24B',
                '10A',
                '13',
                '3A',
                '28',
                '15',
                '5',
                '20A',
                '12A',
                '8D',
                '19',
                '1A',
                '23',
                '25',
                '14',
                '8E',
                '21B',
                '9B',
                '11B',
                '6',
                '22B',
                '7B',
                '7A',
                '18',
                '2B',
                '22A',
                'UN',
                '4A',
                '12D',
                '8A',
                '17',
                '24A',
                '20C'
              ],
              'isUnique': false,
              'comment': 'A standardized code representing the classification of the active ingredient.',
              'LLMComment': 'Used for regulatory purposes to identify the classification of pesticides and other chemicals.'
            }
          ],
          'rowCount': 251,
          'comment': 'This table contains classifications for insecticides based on the IRAC (Insecticide Resistance Action Committee) system, detailing various levels of classification and their descriptions.',
          'LLMComment': 'The `public.irac_classification` table is essential for understanding the classification of insecticides according to the IRAC system. It includes fields such as `irac_class_id`, which uniquely identifies each classification, and `active_ingredient`, which specifies the chemical component of the insecticide. The table is structured into multiple levels (level1 to level4), each with corresponding descriptions, allowing for a hierarchical understanding of insecticide categories. This classification is crucial for researchers and practitioners in agriculture and pest management to identify and manage insecticide resistance effectively.'
        },
        {
          'name': 'ligand_eff',
          'schema': 'public',
          'columns': [
            {
              'name': 'activity_id',
              'table': 'ligand_eff',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Unique identifier for the activity associated with the ligand.'
            },
            {
              'name': 'bei',
              'table': 'ligand_eff',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Binding energy index, representing the energy required for the ligand to bind to its target.',
              'LLMComment': 'A numeric value indicating the binding energy, which is crucial for understanding ligand-target interactions.'
            },
            {
              'name': 'sei',
              'table': 'ligand_eff',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Selectivity energy index, indicating the selectivity of the ligand for its target over other potential targets.',
              'LLMComment': 'A numeric measure of how selectively the ligand binds to its intended target compared to others.'
            },
            {
              'name': 'le',
              'table': 'ligand_eff',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Ligand efficiency, a metric that relates the binding affinity to the molecular weight of the ligand.',
              'LLMComment': 'A calculated value that helps assess the effectiveness of a ligand in relation to its size.'
            },
            {
              'name': 'lle',
              'table': 'ligand_eff',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Ligand lipophilicity efficiency, indicating the balance between lipophilicity and binding affinity.',
              'LLMComment': 'A measure that evaluates the lipophilicity of the ligand in relation to its binding efficiency.'
            }
          ],
          'rowCount': 1558068,
          'comment': 'This table stores ligand efficiency metrics for various activities, including activity IDs and numerical values representing different efficiency measures.',
          'LLMComment': 'The `ligand_eff` table is a crucial component in pharmacological databases, specifically designed to capture and analyze the efficiency of ligands in biological assays. It contains over 1.5 million records, each linked to a unique `activity_id`, which corresponds to specific biological activities. The columns `bei`, `sei`, `le`, and `lle` represent various ligand efficiency indices, providing insights into the effectiveness of compounds in drug discovery and development. This table is essential for researchers and AI models focused on understanding ligand-receptor interactions and optimizing drug candidates.'
        },
        {
          'name': 'mechanism_refs',
          'schema': 'public',
          'columns': [
            {
              'name': 'mecref_id',
              'table': 'mechanism_refs',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 16558,
              'isUnique': true,
              'comment': 'Unique identifier for each mechanism reference.'
            },
            {
              'name': 'mec_id',
              'table': 'mechanism_refs',
              'schema': 'public',
              'type': 'bigint',
              'min': 13,
              'max': 9851,
              'isUnique': false,
              'comment': 'Identifier linking to a specific mechanism.'
            },
            {
              'name': 'ref_type',
              'table': 'mechanism_refs',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'PubMed',
                'Expert',
                'PMC',
                'KEGG',
                'Wikipedia',
                'InterPro',
                'PubChem',
                'UniProt',
                'DOI',
                'FDA',
                'ISBN',
                'BNF',
                'Patent',
                'ClinicalTrials',
                'EMA',
                'Other',
                'DailyMed',
                'HMA',
                'IUPHAR',
                'PMDA'
              ],
              'isUnique': false,
              'comment': 'Type of reference, indicating the source of the information.',
              'LLMComment': 'Categories include PubMed, Expert, PMC, KEGG, Wikipedia, etc.'
            },
            {
              'name': 'ref_id',
              'table': 'mechanism_refs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the specific reference within its source.'
            },
            {
              'name': 'ref_url',
              'table': 'mechanism_refs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'semanticType': 'URL',
              'comment': 'URL linking to the reference source.',
              'LLMComment': 'Semantic type indicating a web address for accessing the reference.'
            }
          ],
          'rowCount': 13343,
          'comment': 'This table stores references related to various mechanisms in the database, linking them to specific types and identifiers.',
          'LLMComment': 'The \'mechanism_refs\' table is a crucial component of the database that catalogs references associated with different biological mechanisms. Each entry includes a unique identifier for the reference (\'mecref_id\'), a link to the mechanism it pertains to (\'mec_id\'), the type of reference (\'ref_type\'), a specific identifier for the reference (\'ref_id\'), and a URL for further information (\'ref_url\'). This structure allows for the organization and retrieval of diverse references, such as publications or online resources, that provide additional context or data about the mechanisms of action for various drugs or biological processes, thereby supporting research and development in pharmacology and biomedicine.'
        },
        {
          'name': 'metabolism',
          'schema': 'public',
          'columns': [
            {
              'name': 'met_id',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'bigint',
              'min': 119,
              'max': 2810,
              'isUnique': true,
              'comment': 'Unique identifier for each metabolism record.'
            },
            {
              'name': 'drug_record_id',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'bigint',
              'min': 2468083,
              'max': 3460335,
              'isUnique': false,
              'comment': 'Identifier for the drug associated with the metabolism record.'
            },
            {
              'name': 'substrate_record_id',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'bigint',
              'min': 2468083,
              'max': 3460343,
              'isUnique': false,
              'comment': 'Identifier for the substrate involved in the metabolic process.'
            },
            {
              'name': 'metabolite_record_id',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'bigint',
              'min': 2468090,
              'max': 3460343,
              'isUnique': false,
              'comment': 'Identifier for the metabolite produced or consumed in the metabolic process.'
            },
            {
              'name': 'pathway_id',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 5,
              'isUnique': false,
              'comment': 'Identifier for the metabolic pathway to which this record belongs.'
            },
            {
              'name': 'pathway_key',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Fig. 2.1.1.1, p.35',
                'null',
                'fig 4_2',
                'Fig 6D',
                'Fig. 8',
                'fig 1',
                'Fig. 2, p.19',
                'fig 9C',
                'Scheme 3',
                'fig 8_E_human',
                'fig 8_S_human',
                'Fig 5',
                'Fig. 2.1.4.1, p.46',
                'Fig 9',
                'Figure 3.6 and scheme on p. 85',
                'fig 8_S_rat',
                'Fig 6',
                'Scheme 1 and Table 4',
                'fig 10',
                'fig 9B',
                'Section 2.6.4.1, p. 33 and Figure 2.2.6.1',
                'Fig. 9, p164/165',
                'fig 2',
                'fig 3S_Mur',
                'Figure 2',
                'Fig. 7',
                'Fig 2',
                'Figure 8.3.7',
                'fig 4_3',
                'Fig 6C',
                'Fig. 1',
                'Fig2',
                'Figure 3.6',
                'fig 9',
                'Scheme 1',
                'fig 4_5',
                'Fig. 4',
                'Figure 4',
                'fig 3',
                'Figure 5 and section 2.2.5.6',
                'fig 5',
                'Section 2.6.4.1, p. 33',
                'Scheme 2',
                'fig 6',
                'SCHEME 1',
                'fig 9D',
                'Fig 1A',
                'fig 4_4',
                'Fig 1',
                'fig 8',
                'Fig5',
                'Fig 1B',
                'Figure 1',
                'fig 4',
                'Fig 4',
                'Scheme 4',
                'fig 8_E_rat',
                'Fig S1',
                'fig 3S_Pel',
                'Fig 3',
                'fig 11',
                'fig 9A',
                'Fig3',
                'fig 4_1',
                'Fig 6A',
                'Fig 6B',
                'fig 7',
                'Fig 7'
              ],
              'isUnique': false,
              'comment': 'Key representing the specific metabolic pathway, often linked to figures in literature.',
              'LLMComment': 'This key helps in identifying the metabolic pathway visually represented in scientific figures.'
            },
            {
              'name': 'enzyme_name',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Name of the enzyme that catalyzes the metabolic reaction.'
            },
            {
              'name': 'enzyme_tid',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'bigint',
              'min': 30,
              'max': 109924,
              'isUnique': false,
              'comment': 'Identifier for the enzyme, possibly linked to a taxonomy or database.'
            },
            {
              'name': 'met_conversion',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Description of the conversion process occurring in the metabolism.',
              'LLMComment': 'This field outlines the specific transformation of substrates into metabolites.'
            },
            {
              'name': 'organism',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                'Mus musculus',
                'Oryctolagus cuniculus',
                'Canis lupus familiaris',
                'Homo sapiens',
                'Sus scrofa',
                'Rattus norvegicus',
                'Macaca fascicularis',
                'Callithrix jacchus'
              ],
              'isUnique': false,
              'comment': 'The organism from which the metabolic data is derived.',
              'LLMComment': 'This indicates the species relevant to the metabolic process, which can affect the interpretation of the data.'
            },
            {
              'name': 'tax_id',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'bigint',
              'min': 9483,
              'max': 10116,
              'isUnique': false,
              'comment': 'Taxonomic identifier for the organism associated with the metabolism record.'
            },
            {
              'name': 'met_comment',
              'table': 'metabolism',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Additional comments or notes regarding the metabolism record.',
              'LLMComment': 'This field can contain supplementary information that provides context or clarifications about the metabolic process.'
            }
          ],
          'rowCount': 2147,
          'comment': 'This table stores information about the metabolism of various substances, including drugs and their metabolites, along with associated pathways and enzymes.',
          'LLMComment': 'The \'public.metabolism\' table is a crucial component in pharmacology and biochemistry, detailing the metabolic processes of drugs and their interactions within biological systems. It includes identifiers for drugs, substrates, and metabolites, as well as information on the enzymes involved in these metabolic pathways. This data is essential for understanding drug metabolism, efficacy, and safety, making it valuable for researchers and developers in the pharmaceutical industry.'
        },
        {
          'name': 'metabolism_refs',
          'schema': 'public',
          'columns': [
            {
              'name': 'metref_id',
              'table': 'metabolism_refs',
              'schema': 'public',
              'type': 'bigint',
              'min': 39,
              'max': 13313,
              'isUnique': true,
              'comment': 'Unique identifier for each metabolism reference.'
            },
            {
              'name': 'met_id',
              'table': 'metabolism_refs',
              'schema': 'public',
              'type': 'bigint',
              'min': 119,
              'max': 2810,
              'isUnique': false,
              'comment': 'Identifier linking to a specific metabolism entry.'
            },
            {
              'name': 'ref_type',
              'table': 'metabolism_refs',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'DailyMed',
                'OTHER',
                'DOI',
                'FDA',
                'DAILYMED',
                'ISBN',
                'PMID'
              ],
              'isUnique': false,
              'comment': 'Type of reference, indicating the source or nature of the reference (e.g., DailyMed, DOI).',
              'LLMComment': 'Categorizes the source of the metabolism reference, helping to identify the context of the information.'
            },
            {
              'name': 'ref_id',
              'table': 'metabolism_refs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the specific reference, which may vary by source type.',
              'LLMComment': 'A unique identifier associated with the reference, used to retrieve or cite the source.'
            },
            {
              'name': 'ref_url',
              'table': 'metabolism_refs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'semanticType': 'URL',
              'comment': 'URL linking to the reference source for further information.',
              'LLMComment': 'A web link that directs to the online resource or document related to the metabolism reference.'
            }
          ],
          'rowCount': 3296,
          'comment': 'This table stores references related to metabolic processes, linking various metabolic entities to their corresponding references, such as articles or databases.',
          'LLMComment': 'The `public.metabolism_refs` table is a crucial component of the metabolic database schema, containing 3,296 entries that connect metabolic entities (identified by `met_id`) to their references. Each entry includes a unique identifier (`metref_id`), the type of reference (`ref_type`), a reference identifier (`ref_id`), and a URL (`ref_url`) for accessing the reference. This table facilitates the integration of external information sources, enhancing the understanding of metabolic processes by providing context and supporting literature for researchers and AI systems analyzing metabolic data.'
        },
        {
          'name': 'molecule_atc_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'mol_atc_id',
              'table': 'molecule_atc_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 91486,
              'max': 95888,
              'isUnique': true,
              'comment': 'Unique identifier for the ATC classification of the molecule.',
              'LLMComment': 'A unique identifier assigned to each molecule\'s ATC classification, ranging from 91486 to 95888.'
            },
            {
              'name': 'level5',
              'table': 'molecule_atc_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Fifth level of the ATC classification, representing a specific therapeutic subgroup.',
              'LLMComment': 'This column indicates the fifth level of the Anatomical Therapeutic Chemical (ATC) classification, which specifies a particular therapeutic subgroup.'
            },
            {
              'name': 'molregno',
              'table': 'molecule_atc_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 97,
              'max': 2832777,
              'isUnique': false,
              'semanticType': 'molregno',
              'comment': 'Unique identifier for the molecule in the database.',
              'LLMComment': 'A unique identifier for the molecule, known as molregno, with a range from 97 to 2832777.'
            }
          ],
          'rowCount': 4403,
          'comment': 'This table stores the classification of molecules based on the Anatomical Therapeutic Chemical (ATC) classification system, linking each molecule to its corresponding ATC code at level 5.',
          'LLMComment': 'The `molecule_atc_classification` table is crucial for understanding how various molecules are categorized within the ATC classification system, which is widely used for drug classification. Each entry in this table associates a unique molecule identifier (`molregno`) with its specific ATC classification code at the fifth level (`level5`). This classification helps in identifying the therapeutic use of the molecules, facilitating drug discovery and development processes by providing insights into the pharmacological properties and potential applications of the compounds.'
        },
        {
          'name': 'molecule_dictionary',
          'schema': 'public',
          'columns': [
            {
              'name': 'molregno',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'bigint',
              'semanticType': 'molregno'
            },
            {
              'name': 'pref_name',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Preferred name of the molecule as recognized in the database.',
              'LLMComment': 'The commonly used name for the molecule.'
            },
            {
              'name': 'chembl_id',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'semanticType': 'CHEMBL_ID',
              'comment': 'Unique identifier for the molecule in the ChEMBL database.',
              'LLMComment': 'Identifier used to reference the molecule in the ChEMBL database.'
            },
            {
              'name': 'max_phase',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'numeric',
              'comment': 'Highest clinical development phase reached by the molecule.',
              'LLMComment': 'Indicates the maximum phase of clinical trials the molecule has undergone.'
            },
            {
              'name': 'therapeutic_flag',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule is approved for therapeutic use.',
              'LLMComment': 'Flag indicating whether the molecule is considered a therapeutic agent.'
            },
            {
              'name': 'dosed_ingredient',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule is a dosed ingredient in formulations.',
              'LLMComment': 'Flag to show if the molecule is used as an active ingredient in dosages.'
            },
            {
              'name': 'structure_type',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Type of chemical structure (e.g., small molecule, biologic).',
              'LLMComment': 'Describes the structural classification of the molecule.'
            },
            {
              'name': 'chebi_par_id',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Identifier for the molecule in the ChEBI database.',
              'LLMComment': 'Unique identifier for the molecule in the Chemical Entities of Biological Interest database.'
            },
            {
              'name': 'molecule_type',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Classification of the molecule (e.g., small molecule, peptide).',
              'LLMComment': 'Type of molecule based on its chemical structure and properties.'
            },
            {
              'name': 'first_approval',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Year the molecule was first approved for use.',
              'LLMComment': 'The year when the molecule received its first regulatory approval.'
            },
            {
              'name': 'oral',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule can be administered orally.',
              'LLMComment': 'Flag indicating the oral bioavailability of the molecule.'
            },
            {
              'name': 'parenteral',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule can be administered parenterally.',
              'LLMComment': 'Flag indicating the suitability of the molecule for parenteral administration.'
            },
            {
              'name': 'topical',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule can be applied topically.',
              'LLMComment': 'Flag indicating the molecule\'s use in topical applications.'
            },
            {
              'name': 'black_box_warning',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule has a black box warning from regulatory agencies.',
              'LLMComment': 'Flag indicating the presence of serious safety warnings associated with the molecule.'
            },
            {
              'name': 'natural_product',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule is derived from natural sources.',
              'LLMComment': 'Flag indicating whether the molecule is a natural product.'
            },
            {
              'name': 'first_in_class',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule is the first of its kind in its therapeutic class.',
              'LLMComment': 'Flag indicating if the molecule is a pioneer in its therapeutic category.'
            },
            {
              'name': 'chirality',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule has chiral centers.',
              'LLMComment': 'Flag indicating the presence of chirality in the molecule.'
            },
            {
              'name': 'prodrug',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule is a prodrug.',
              'LLMComment': 'Flag indicating whether the molecule is administered in an inactive form that is metabolized into an active form.'
            },
            {
              'name': 'inorganic_flag',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule is inorganic.',
              'LLMComment': 'Flag indicating whether the molecule is classified as inorganic.'
            },
            {
              'name': 'usan_year',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'integer',
              'comment': 'Year the molecule received its USAN (United States Adopted Name).',
              'LLMComment': 'The year when the molecule was assigned a USAN.'
            },
            {
              'name': 'availability_type',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Type of availability for the molecule (e.g., prescription, over-the-counter).',
              'LLMComment': 'Indicates the availability status of the molecule in the market.'
            },
            {
              'name': 'usan_stem',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Stem used in the USAN for the molecule.',
              'LLMComment': 'The stem that forms the basis of the molecule\'s USAN.'
            },
            {
              'name': 'polymer_flag',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule is a polymer.',
              'LLMComment': 'Flag indicating whether the molecule is classified as a polymer.'
            },
            {
              'name': 'usan_substem',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Substem used in the USAN for the molecule.',
              'LLMComment': 'The sub-stem that provides additional classification in the USAN.'
            },
            {
              'name': 'usan_stem_definition',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Definition of the USAN stem used for the molecule.',
              'LLMComment': 'Describes the meaning or significance of the USAN stem.'
            },
            {
              'name': 'indication_class',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'comment': 'Class of indications for which the molecule is used.',
              'LLMComment': 'Categorizes the therapeutic indications associated with the molecule.'
            },
            {
              'name': 'withdrawn_flag',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule has been withdrawn from the market.',
              'LLMComment': 'Flag indicating whether the molecule is no longer available for use.'
            },
            {
              'name': 'chemical_probe',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule is classified as a chemical probe.',
              'LLMComment': 'Flag indicating whether the molecule is used as a chemical probe in research.'
            },
            {
              'name': 'orphan',
              'table': 'molecule_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'comment': 'Indicates if the molecule is designated for orphan drug status.',
              'LLMComment': 'Flag indicating whether the molecule is recognized as an orphan drug.'
            }
          ],
          'rowCount': 2431025,
          'comment': 'This table contains a comprehensive dictionary of molecules, including their identifiers, names, types, approval statuses, and various pharmacological properties.',
          'LLMComment': 'The `molecule_dictionary` table serves as a central repository for information about various chemical entities, including drugs and biologics. It includes key identifiers like `chembl_id` and `molregno`, along with attributes such as `max_phase` indicating the highest clinical trial phase reached, `therapeutic_flag` for categorizing therapeutic uses, and `first_approval` to track the year of initial market approval. This table is crucial for researchers and developers in the pharmaceutical domain, as it provides essential data for drug discovery, development, and regulatory compliance. The extensive row count suggests a rich dataset that can be leveraged for various analyses, including drug efficacy, safety profiles, and market trends.'
        },
        {
          'name': 'molecule_frac_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'mol_frac_id',
              'table': 'molecule_frac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 147,
              'isUnique': true,
              'comment': 'Unique identifier for each molecular fraction entry.',
              'LLMComment': 'Unique identifier for each molecular fraction entry.'
            },
            {
              'name': 'frac_class_id',
              'table': 'molecule_frac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 470,
              'max': 685,
              'isUnique': false,
              'comment': 'Identifier for the classification of the molecular fraction, linking to a classification system.',
              'LLMComment': 'Identifier for the classification of the molecular fraction, linking to a classification system.'
            },
            {
              'name': 'molregno',
              'table': 'molecule_frac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 6738,
              'max': 1678348,
              'isUnique': true,
              'semanticType': 'molregno',
              'comment': 'Unique identifier for the molecule, used to reference specific molecules in the database.',
              'LLMComment': 'Unique identifier for the molecule, used to reference specific molecules in the database.'
            }
          ],
          'rowCount': 145,
          'comment': 'This table stores the classification of molecular fractions, linking each fraction to its corresponding molecular registration number and classification ID.',
          'LLMComment': 'The \'public.molecule_frac_classification\' table is a key component in the domain of molecular biology and pharmacology. It contains 145 records that associate molecular fractions with their respective classifications and registration numbers. Each entry is identified by a unique \'mol_frac_id\' and links to a \'frac_class_id\' that categorizes the fraction, as well as a \'molregno\' that references the specific molecular entity. This table is crucial for understanding how different molecular fractions are classified and can be used in various analyses related to drug development, compound properties, and biological activities.'
        },
        {
          'name': 'molecule_hierarchy',
          'schema': 'public',
          'columns': [
            {
              'name': 'molregno',
              'table': 'molecule_hierarchy',
              'schema': 'public',
              'type': 'bigint',
              'semanticType': 'molregno',
              'comment': 'Molecule registration number, uniquely identifies a molecule.'
            },
            {
              'name': 'parent_molregno',
              'table': 'molecule_hierarchy',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Registration number of the parent molecule in the hierarchy.',
              'LLMComment': 'This column indicates the hierarchical relationship by referencing the parent molecule\'s registration number.'
            },
            {
              'name': 'active_molregno',
              'table': 'molecule_hierarchy',
              'schema': 'public',
              'type': 'bigint',
              'comment': 'Registration number of the active form of the molecule, if applicable.',
              'LLMComment': 'This column is used to identify the active variant of the molecule, which may differ from its base form.'
            }
          ],
          'rowCount': 2342837,
          'comment': 'This table represents the hierarchy of molecules, detailing their relationships through unique identifiers.',
          'LLMComment': 'The \'molecule_hierarchy\' table is a crucial component of the molecular database, containing over 2.3 million entries that define the hierarchical relationships between different molecules. Each entry includes a unique \'molregno\' (molecule registration number), a \'parent_molregno\' indicating the parent molecule in the hierarchy, and an \'active_molregno\' which may represent an active form or variant of the molecule. This structure allows for the organization and retrieval of molecular data, facilitating the understanding of molecular relationships and classifications within the broader context of drug discovery and bioinformatics.'
        },
        {
          'name': 'molecule_hrac_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'mol_hrac_id',
              'table': 'molecule_hrac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 205,
              'isUnique': true,
              'comment': 'Unique identifier for each HRAC classification entry.'
            },
            {
              'name': 'hrac_class_id',
              'table': 'molecule_hrac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 272,
              'isUnique': false,
              'comment': 'Identifier for the HRAC classification category to which the molecule belongs.'
            },
            {
              'name': 'molregno',
              'table': 'molecule_hrac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 4361,
              'max': 1768763,
              'isUnique': false,
              'semanticType': 'molregno',
              'comment': 'Unique registration number assigned to the molecule, used for tracking and identification.'
            }
          ],
          'rowCount': 205,
          'comment': 'This table stores the classification of molecules based on the HRAC (Human Receptor Activity Classification) system, linking each molecule to its corresponding HRAC class.',
          'LLMComment': 'The `public.molecule_hrac_classification` table is essential for categorizing molecules according to the HRAC system, which is used to classify human receptor activities. It contains three key columns: `mol_hrac_id`, which uniquely identifies the classification entry; `hrac_class_id`, which links to the specific HRAC class; and `molregno`, which references the molecule\'s registration number. This classification is crucial for understanding the pharmacological properties and potential therapeutic uses of various molecules in drug development and research.'
        },
        {
          'name': 'molecule_irac_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'mol_irac_id',
              'table': 'molecule_irac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 186,
              'isUnique': true,
              'comment': 'Unique identifier for each molecule IRAC classification entry.'
            },
            {
              'name': 'irac_class_id',
              'table': 'molecule_irac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 251,
              'isUnique': false,
              'comment': 'Identifier for the specific IRAC classification category assigned to the molecule.'
            },
            {
              'name': 'molregno',
              'table': 'molecule_irac_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 115,
              'max': 1765082,
              'isUnique': true,
              'semanticType': 'molregno',
              'comment': 'Unique registration number for the molecule, used for tracking and identification.'
            }
          ],
          'rowCount': 184,
          'comment': 'This table links molecules to their respective IRAC (Insecticide Resistance Action Committee) classification IDs, providing a way to categorize molecules based on their mode of action against pests.',
          'LLMComment': 'The `public.molecule_irac_classification` table serves as a crucial link between chemical entities (molecules) and their corresponding IRAC classification IDs, which are essential for understanding the mechanisms of action of insecticides. With 184 entries, this table allows researchers and developers in the agrochemical domain to identify and categorize molecules based on their resistance management strategies. The `mol_irac_id` serves as a unique identifier for the classification, while `irac_class_id` connects to the specific IRAC category, and `molregno` refers to the unique registration number of the molecule. This classification is vital for developing effective pest control strategies and ensuring sustainable agricultural practices.'
        },
        {
          'name': 'molecule_synonyms',
          'schema': 'public',
          'columns': [
            {
              'name': 'molregno',
              'table': 'molecule_synonyms',
              'schema': 'public',
              'type': 'bigint',
              'min': 23,
              'max': 2832784,
              'isUnique': false,
              'semanticType': 'molregno',
              'comment': 'Molecule registration number, uniquely identifies a molecule.'
            },
            {
              'name': 'syn_type',
              'table': 'molecule_synonyms',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'INN',
                'USP',
                'OTHER',
                'SYSTEMATIC',
                'E_NUMBER',
                'TRADE_NAME',
                'FDA',
                'USAN',
                'BNF',
                'JAN',
                'ATC',
                'DCF',
                'RESEARCH_CODE',
                'NATIONAL_FORMULARY',
                'MERCK_INDEX',
                'BAN'
              ],
              'isUnique': false,
              'comment': 'Type of synonym, indicating the classification such as INN (International Nonproprietary Name), USP (United States Pharmacopeia), etc.',
              'LLMComment': 'Describes the category of the synonym, which helps in understanding its usage and context.'
            },
            {
              'name': 'molsyn_id',
              'table': 'molecule_synonyms',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 3124464,
              'isUnique': true,
              'comment': 'Unique identifier for the synonym entry.'
            },
            {
              'name': 'res_stem_id',
              'table': 'molecule_synonyms',
              'schema': 'public',
              'type': 'bigint',
              'min': 3,
              'max': 671,
              'isUnique': false,
              'comment': 'Identifier for the root stem of the synonym, used for categorizing related synonyms.',
              'LLMComment': 'Links to the root structure of the molecule, aiding in the classification of synonyms.'
            },
            {
              'name': 'synonyms',
              'table': 'molecule_synonyms',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The actual synonym name for the molecule.'
            }
          ],
          'rowCount': 211047,
          'comment': 'This table stores synonyms for various molecules, allowing for different naming conventions and identifiers used in chemical and biological research.',
          'LLMComment': 'The \'molecule_synonyms\' table is a crucial component of the database that facilitates the identification and retrieval of molecules by their various synonyms. With over 211,000 entries, it includes columns for the unique molecule registration number (molregno), the type of synonym (syn_type), a unique identifier for the synonym (molsyn_id), a reference to the research stem (res_stem_id), and the actual synonym text (synonyms). This table supports the integration of diverse nomenclature used in scientific literature, ensuring that researchers can find and link related molecular data across different studies and databases.'
        },
        {
          'name': 'organism_class',
          'schema': 'public',
          'columns': [
            {
              'name': 'oc_id',
              'table': 'organism_class',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 4545,
              'isUnique': true
            },
            {
              'name': 'tax_id',
              'table': 'organism_class',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 3052763,
              'isUnique': true
            },
            {
              'name': 'l1',
              'table': 'organism_class',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Bacteria',
                'Viruses',
                'Fungi',
                'Unclassified',
                'Archaea',
                'Eukaryotes'
              ],
              'isUnique': false,
              'comment': 'Primary classification level indicating the broad category of the organism, such as Bacteria or Viruses.',
              'LLMComment': 'This column represents the highest taxonomic rank for the organism, categorizing it into major groups like Bacteria, Viruses, and Archaea.'
            },
            {
              'name': 'l2',
              'table': 'organism_class',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                'Kinetoplastida',
                'Rhizaria',
                'Echinodermata',
                'Gram-Positive',
                'Amphibia',
                'Gram-positive',
                'Arthropoda',
                'Chondrichthyes',
                'Other',
                'Fornicata',
                'Mollusca',
                'Microsporidia',
                'Platyhelminthes',
                'ssRNA',
                'Teleostei',
                'Mammalia',
                'retro-transcribing',
                'dsDNA',
                'Apicomplexa',
                'Lepidosauria',
                'Gram-Negative',
                'dsRNA',
                'Nematoda',
                'Annelida ',
                'Parabasalia',
                'Ascomycota',
                'Alveolata',
                'Aves',
                'Viridiplantae',
                'Amoebozoa',
                'Basidiomycota',
                'Mucorales',
                'Oomycetes'
              ],
              'isUnique': false,
              'comment': 'Secondary classification level providing more specific categories within the primary classification.',
              'LLMComment': 'This column offers a finer classification of the organism, detailing subcategories like Kinetoplastida or Echinodermata.'
            },
            {
              'name': 'l3',
              'table': 'organism_class',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Tertiary classification level that may include additional details or specific classifications not covered in l1 or l2.',
              'LLMComment': 'This column can contain further distinctions or classifications of the organism, providing additional context beyond the first two levels.'
            }
          ],
          'rowCount': 4190,
          'comment': 'This table contains classifications of organisms based on their taxonomic identifiers and hierarchical levels.',
          'LLMComment': 'The \'public.organism_class\' table serves as a key reference for categorizing various organisms within a biological or ecological context. It includes a unique identifier for each organism (oc_id), a taxonomic identifier (tax_id) that links to broader taxonomic databases, and three levels of classification (l1, l2, l3) that provide a hierarchical structure for understanding the relationships between different organisms. This table is essential for researchers and AI models working in fields such as bioinformatics, ecology, and pharmacology, as it helps in organizing and retrieving organism-related data across various biological studies and applications.'
        },
        {
          'name': 'patent_use_codes',
          'schema': 'public',
          'columns': [
            {
              'name': 'patent_use_code',
              'table': 'patent_use_codes',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Unique code representing a specific patent use.',
              'LLMComment': 'A unique identifier for each type of patent use, allowing for easy reference and categorization.'
            },
            {
              'name': 'definition',
              'table': 'patent_use_codes',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Description of the patent use associated with the code.',
              'LLMComment': 'A detailed explanation of what the patent use code signifies, providing context and understanding of its application.'
            }
          ],
          'rowCount': 3737,
          'comment': 'This table contains codes used to classify the various uses of patents in the context of biotherapeutics and related research.',
          'LLMComment': 'The `public.patent_use_codes` table serves as a reference for categorizing the different applications of patents within the biotherapeutics domain. Each entry consists of a unique `patent_use_code` and its corresponding `definition`, which provides clarity on how each code is utilized in the context of drug development, research activities, and regulatory compliance. This classification is essential for understanding the landscape of intellectual property in biomedicine, aiding researchers and developers in navigating patent-related information effectively.'
        },
        {
          'name': 'predicted_binding_domains',
          'schema': 'public',
          'columns': [
            {
              'name': 'predbind_id',
              'table': 'predicted_binding_domains',
              'schema': 'public',
              'type': 'bigint',
              'min': 7699241,
              'max': 8521966,
              'isUnique': true
            },
            {
              'name': 'activity_id',
              'table': 'predicted_binding_domains',
              'schema': 'public',
              'type': 'bigint',
              'min': 31931,
              'max': 22994980,
              'isUnique': false
            },
            {
              'name': 'site_id',
              'table': 'predicted_binding_domains',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 38445,
              'isUnique': false
            },
            {
              'name': 'prediction_method',
              'table': 'predicted_binding_domains',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Multi domain',
                'Single domain'
              ],
              'isUnique': false,
              'comment': 'Indicates whether the binding domain prediction is based on multiple domains or a single domain.',
              'LLMComment': 'Categorizes the prediction approach as either \'Multi domain\' or \'Single domain\'.'
            },
            {
              'name': 'confidence',
              'table': 'predicted_binding_domains',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'medium',
                'high'
              ],
              'isUnique': false,
              'comment': 'Represents the confidence level of the prediction, indicating the reliability of the binding domain prediction.',
              'LLMComment': 'Classifies the confidence in the prediction as either \'medium\' or \'high\'.'
            }
          ],
          'rowCount': 822666,
          'comment': 'This table stores predicted binding domains for various biological entities, including their associated confidence levels and prediction methods.',
          'LLMComment': 'The `public.predicted_binding_domains` table is a crucial component in the domain of bioinformatics and drug discovery. It contains 822,666 entries that represent predicted binding domains, which are essential for understanding how different biological molecules interact with each other. Each entry includes a unique identifier (`predbind_id`), links to specific activities (`activity_id`), and binding sites (`site_id`). The `prediction_method` column indicates the algorithm or approach used to make the prediction, while the `confidence` column provides a measure of certainty regarding the prediction\'s accuracy. This table is integral for researchers analyzing molecular interactions, developing therapeutics, and studying biological pathways.'
        },
        {
          'name': 'product_patents',
          'schema': 'public',
          'columns': [
            {
              'name': 'prod_pat_id',
              'table': 'product_patents',
              'schema': 'public',
              'type': 'bigint',
              'min': 257,
              'max': 33399,
              'isUnique': true,
              'comment': 'Unique identifier for each product patent record.'
            },
            {
              'name': 'product_id',
              'table': 'product_patents',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the product associated with the patent.'
            },
            {
              'name': 'patent_no',
              'table': 'product_patents',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The number assigned to the patent.'
            },
            {
              'name': 'patent_expire_date',
              'table': 'product_patents',
              'schema': 'public',
              'type': 'timestamp without time zone',
              'isUnique': false,
              'comment': 'The date when the patent protection expires.',
              'LLMComment': 'Indicates the expiration date of the patent, after which the patent holder may lose exclusive rights.'
            },
            {
              'name': 'drug_substance_flag',
              'table': 'product_patents',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the patent is related to a drug substance (1) or not (0).',
              'LLMComment': 'A binary flag indicating whether the patent pertains to the active pharmaceutical ingredient.'
            },
            {
              'name': 'drug_product_flag',
              'table': 'product_patents',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the patent is related to a drug product (1) or not (0).',
              'LLMComment': 'A binary flag indicating whether the patent pertains to the final drug formulation.'
            },
            {
              'name': 'patent_use_code',
              'table': 'product_patents',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Code representing the specific use of the patent.',
              'LLMComment': 'A code that categorizes the use or application of the patent in the pharmaceutical context.'
            },
            {
              'name': 'delist_flag',
              'table': 'product_patents',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the patent has been delisted (1) or not (0).',
              'LLMComment': 'A binary flag indicating whether the patent has been removed from active status.'
            },
            {
              'name': 'submission_date',
              'table': 'product_patents',
              'schema': 'public',
              'type': 'timestamp without time zone',
              'isUnique': false,
              'comment': 'The date when the patent application was submitted.',
              'LLMComment': 'Records the date the patent application was filed, important for tracking patent timelines.'
            }
          ],
          'rowCount': 19031,
          'comment': 'This table stores information about patents associated with various products, including their expiration dates and flags indicating their status as drug substances or products.',
          'LLMComment': 'The `public.product_patents` table is a crucial component in the pharmaceutical and biotechnology domain, as it links products to their respective patents. Each entry includes a unique product patent ID, the associated product ID, the patent number, and the expiration date of the patent. Additionally, it contains flags that indicate whether the patent pertains to a drug substance or a drug product, as well as a patent use code that categorizes the patent\'s application. This table is essential for tracking the intellectual property status of pharmaceutical products, which is vital for regulatory compliance, market exclusivity, and strategic planning in drug development.'
        },
        {
          'name': 'products',
          'schema': 'public',
          'columns': [
            {
              'name': 'dosage_form',
              'table': 'products',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Form in which the product is administered (e.g., tablet, liquid).',
              'LLMComment': 'Describes the physical form of the medication, indicating how it is delivered to the patient.'
            },
            {
              'name': 'route',
              'table': 'products',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'URETHRAL',
                'ENTERAL',
                'INTRAVENOUS, INTRAOCULAR, INTRAMUSCULAR, SUBCUTANEOUS',
                'INTRAVENOUS, SUBCUTANEOUS',
                'RECTAL',
                'ORAL',
                'URETERAL',
                'INHALATION, INTRAVENOUS',
                'PERFUSION, BILIARY',
                'RECTAL, SUBLINGUAL',
                'INHALATION, ORAL',
                'IM-IV',
                'INTRAMUSCULAR, SUBCUTANEOUS',
                'INHALATION',
                'INTRA-ANAL',
                'INTRAMUSCULAR, INTRAVENOUS, SUBCUTANEOUS',
                'INTRA-ARTICULAR, INTRAMUSCULAR, INTRAVITREAL',
                'BUCCAL, SUBLINGUAL',
                'SPINAL',
                'VAGINAL',
                'IRRIGATION, URETHRAL',
                'DENTAL',
                'INJECTION, ORAL',
                'SUBLINGUAL',
                'IV (INFUSION), SUBCUTANEOUS',
                'PERFUSION, CARDIAC',
                'N/A',
                'ORAL-21',
                'TOPICAL',
                'IMPLANTATION',
                'TOPICAL, VAGINAL',
                'INHALATION, INJECTION',
                'INTRAMUSCULAR',
                'ENDOCERVICAL',
                'INTRAVENOUS',
                'OPHTHALMIC',
                'INTRACRANIAL',
                'TRANSDERMAL',
                'INTRAOCULAR',
                'INTRAMUSCULAR, INTRAVENOUS',
                'IRRIGATION',
                'INTRAPLEURAL',
                'INTRAVENOUS, INTERSTITIAL',
                'EPIDURAL',
                'NASAL',
                'OTIC',
                'PERIODONTAL',
                'SUBCUTANEOUS',
                'INTRADERMAL',
                'OPHTHALMIC, OTIC',
                'PYELOCALYCEAL',
                'INTRAVENOUS, INTRAMUSCULAR, SUBCUTANEOUS, INTRAOSSEOUS, ENDOTRACHEAL',
                'INJECTION, ORAL, RECTAL',
                'INTRA-ARTICULAR',
                'INTRACAVITARY, INTRAVENOUS, INTRAVESICAL',
                'IONTOPHORESIS',
                'PERIARTICULAR',
                'BUCCAL',
                'INTRAVESICAL',
                'INTRALYMPHATIC, INTRAUTERINE',
                'ORAL-28',
                'INTRAPERITONEAL',
                'INTRAVENOUS, INTRAVENOUS',
                'INTRAVENOUS, INTRAVESICULAR, OPHTHALMIC',
                'INJECTION',
                'INTRATRACHEAL',
                'INTRAMUSCULAR, ORAL',
                'ORAL, RECTAL',
                'INTRA-ARTERIAL',
                'INFILTRATION',
                'TRANSMUCOSAL',
                'INTRAUTERINE',
                'IONTOPHORESIS, TRANSDERMAL',
                'FOR RX COMPOUNDING',
                'IV (INFUSION)',
                'IONTOPHORESIS, TOPICAL',
                'ORAL-20',
                'INTRATHECAL',
                'INTRAVITREAL',
                'INTRAVENOUS, ORAL'
              ],
              'isUnique': false,
              'comment': 'Method of administration for the product, such as URETHRAL or INTRAVENOUS.',
              'LLMComment': 'Specifies the pathway through which the medication is delivered into the body, crucial for understanding its use.'
            },
            {
              'name': 'trade_name',
              'table': 'products',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Commercial name under which the product is marketed.',
              'LLMComment': 'The brand name given to the product, which may differ from its generic name.'
            },
            {
              'name': 'approval_date',
              'table': 'products',
              'schema': 'public',
              'type': 'timestamp without time zone',
              'isUnique': false,
              'comment': 'Date when the product received regulatory approval.',
              'LLMComment': 'Indicates when the product was officially approved for use by regulatory authorities.'
            },
            {
              'name': 'ad_type',
              'table': 'products',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'RX',
                'DISCN',
                'OTC'
              ],
              'isUnique': false,
              'comment': 'Type of advertisement classification for the product, such as RX or OTC.',
              'LLMComment': 'Categorizes the product based on its marketing and prescription status.'
            },
            {
              'name': 'oral',
              'table': 'products',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the product can be taken orally (1 for yes, 0 for no).',
              'LLMComment': 'Binary indicator showing whether the product is suitable for oral administration.'
            },
            {
              'name': 'topical',
              'table': 'products',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the product is applied to the skin (1 for yes, 0 for no).',
              'LLMComment': 'Binary indicator showing whether the product is intended for topical use.'
            },
            {
              'name': 'parenteral',
              'table': 'products',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the product is administered by injection (1 for yes, 0 for no).',
              'LLMComment': 'Binary indicator showing whether the product is suitable for parenteral administration.'
            },
            {
              'name': 'black_box_warning',
              'table': 'products',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the product has a black box warning (1 for yes, 0 for no).',
              'LLMComment': 'Binary indicator that signifies if the product carries a serious safety warning.'
            },
            {
              'name': 'applicant_full_name',
              'table': 'products',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Full name of the entity that applied for product approval.',
              'LLMComment': 'Identifies the organization or individual responsible for submitting the product for approval.'
            },
            {
              'name': 'innovator_company',
              'table': 'products',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates if the product is from an innovator company (1 for yes, 0 for no).',
              'LLMComment': 'Binary indicator showing whether the product is developed by the original innovator.'
            },
            {
              'name': 'product_id',
              'table': 'products',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true
            },
            {
              'name': 'nda_type',
              'table': 'products',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'A',
                'null',
                'N'
              ],
              'isUnique': false,
              'comment': 'Type of New Drug Application, such as A or N.',
              'LLMComment': 'Categorizes the application type for the product\'s approval process.'
            }
          ],
          'rowCount': 45008,
          'comment': 'This table contains detailed information about pharmaceutical products, including their dosage forms, routes of administration, trade names, approval dates, and various classifications related to their use and safety.',
          'LLMComment': 'The \'public.products\' table serves as a comprehensive repository for pharmaceutical product data, crucial for understanding the landscape of approved medications. It includes key attributes such as dosage forms (e.g., oral, topical, parenteral), trade names, and approval dates, which are essential for regulatory compliance and market analysis. Additionally, the table captures safety information like black box warnings and categorizes products by their application types, aiding in the assessment of therapeutic options and potential risks associated with each product. This data is vital for researchers, healthcare professionals, and regulatory bodies in the pharmaceutical domain.'
        },
        {
          'name': 'protein_class_synonyms',
          'schema': 'public',
          'columns': [
            {
              'name': 'protclasssyn_id',
              'table': 'protein_class_synonyms',
              'schema': 'public',
              'type': 'bigint',
              'min': 33666,
              'max': 44911,
              'isUnique': true,
              'comment': 'Unique identifier for each protein class synonym entry.'
            },
            {
              'name': 'protein_class_id',
              'table': 'protein_class_synonyms',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 1715,
              'isUnique': false,
              'comment': 'Identifier linking to the corresponding protein class.'
            },
            {
              'name': 'protein_class_synonym',
              'table': 'protein_class_synonyms',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Alternative name or term used to refer to the protein class.',
              'LLMComment': 'This column contains synonyms that may be used interchangeably with the primary protein class name.'
            },
            {
              'name': 'syn_type',
              'table': 'protein_class_synonyms',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'MESH_XREF',
                'UMLS',
                'CW_XREF',
                'CONCEPT_WIKI'
              ],
              'isUnique': false,
              'comment': 'Type of synonym, indicating the source or classification of the synonym.',
              'LLMComment': 'This column categorizes the synonym based on its origin, such as MESH_XREF for Medical Subject Headings or UMLS for Unified Medical Language System.'
            }
          ],
          'rowCount': 7539,
          'comment': 'This table stores synonyms for various protein classes, allowing for alternative naming and categorization of proteins in biological research.',
          'LLMComment': 'The `public.protein_class_synonyms` table is essential for managing the diverse nomenclature associated with protein classes in biological databases. It contains 7,539 entries that link unique identifiers (`protclasssyn_id`) to specific protein classes (`protein_class_id`) and their synonyms (`protein_class_synonym`). The `syn_type` column categorizes the type of synonym, which can be crucial for researchers and AI systems to understand the relationships and variations in protein classification. This table supports the integration and retrieval of protein-related data across various biological and biochemical applications, enhancing the accuracy and comprehensiveness of protein classification in research.'
        },
        {
          'name': 'protein_classification',
          'schema': 'public',
          'columns': [
            {
              'name': 'protein_class_id',
              'table': 'protein_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 0,
              'max': 1715,
              'isUnique': true,
              'comment': 'Unique identifier for each protein class.'
            },
            {
              'name': 'parent_id',
              'table': 'protein_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 0,
              'max': 1713,
              'isUnique': false,
              'comment': 'Identifier for the parent protein class, indicating hierarchical relationships.'
            },
            {
              'name': 'pref_name',
              'table': 'protein_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Preferred name of the protein class, used for display purposes.'
            },
            {
              'name': 'short_name',
              'table': 'protein_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Abbreviated name for the protein class, often used in concise contexts.'
            },
            {
              'name': 'protein_class_desc',
              'table': 'protein_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Detailed description of the protein class, providing context and information.'
            },
            {
              'name': 'definition',
              'table': 'protein_classification',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'semanticType': 'Text',
              'comment': 'A semantic text providing a formal definition of the protein class.'
            },
            {
              'name': 'class_level',
              'table': 'protein_classification',
              'schema': 'public',
              'type': 'bigint',
              'min': 0,
              'max': 6,
              'isUnique': false,
              'comment': 'Numerical representation of the classification level, indicating the hierarchy of the protein class.'
            }
          ],
          'rowCount': 905,
          'comment': 'This table contains classifications of proteins, including their identifiers, names, descriptions, and hierarchical relationships.',
          'LLMComment': 'The `public.protein_classification` table serves as a comprehensive repository for categorizing proteins based on various attributes. It includes unique identifiers for each protein class, hierarchical relationships through parent-child links, and descriptive fields that provide insights into the protein\'s function and classification. This table is crucial for bioinformatics and molecular biology research, as it helps in organizing and retrieving information about protein families, facilitating studies on protein functions, interactions, and their roles in biological processes.'
        },
        {
          'name': 'relationship_type',
          'schema': 'public',
          'columns': [
            {
              'name': 'relationship_type',
              'table': 'relationship_type',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'N',
                'S',
                'H',
                'U',
                'M',
                'D'
              ],
              'isUnique': true,
              'comment': 'Type of relationship, represented by a single character code.',
              'LLMComment': 'A unique character code representing the type of relationship, such as N for \'Nuclear\', S for \'Subcellular\', etc.'
            },
            {
              'name': 'relationship_desc',
              'table': 'relationship_type',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Direct protein target assigned',
                'Molecular target other than protein assigned',
                'Subcellular target assigned',
                'Default value - Target has yet to be curated',
                'Homologous protein target assigned',
                'Non-molecular target assigned'
              ],
              'isUnique': true,
              'comment': 'Description of the relationship type, providing more context.',
              'LLMComment': 'A detailed description of the relationship type, explaining what each code signifies, such as \'Direct protein target assigned\' or \'Homologous protein target assigned\'.'
            }
          ],
          'rowCount': 6,
          'comment': 'This table defines different types of relationships that can exist between entities in the database, along with descriptions for each type.',
          'LLMComment': 'The `public.relationship_type` table serves as a reference for categorizing various relationships within the database. It contains two columns: `relationship_type`, which specifies the type of relationship (e.g., \'parent\', \'child\', \'sibling\'), and `relationship_desc`, which provides a detailed description of what that relationship entails. This table is crucial for understanding how different entities interact with one another, facilitating data integrity and clarity in relational mappings across the schema.'
        },
        {
          'name': 'research_companies',
          'schema': 'public',
          'columns': [
            {
              'name': 'co_stem_id',
              'table': 'research_companies',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 812,
              'isUnique': true,
              'comment': 'Unique identifier for the research company.'
            },
            {
              'name': 'res_stem_id',
              'table': 'research_companies',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 671,
              'isUnique': false,
              'comment': 'Identifier linking to the research study or project associated with the company.'
            },
            {
              'name': 'company',
              'table': 'research_companies',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Name of the research company.'
            },
            {
              'name': 'country',
              'table': 'research_companies',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Switzerland',
                'New Zealand',
                'Italy',
                'Hungary',
                'China',
                'Korea',
                'Sweden',
                'Norway',
                'USA',
                'Netherlands',
                'Austria',
                'UK',
                'Australia',
                'Ireland',
                'Germany',
                'Singapore',
                'Canada',
                'Finland',
                'Spain',
                'Slovenia',
                'India',
                'Belgium',
                'France',
                'italy',
                'Israel',
                'Iceland',
                'Japan',
                'Denmark'
              ],
              'isUnique': false,
              'semanticType': 'gis-country',
              'comment': 'Country where the research company is based, categorized by GIS standards.',
              'LLMComment': 'Represents the geographical location of the company, useful for regional analysis.'
            },
            {
              'name': 'previous_company',
              'table': 'research_companies',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Xenova',
                'Lederle',
                'Athenagen',
                'Ciba-Geigy',
                'Schering-Plough',
                'null',
                'Medimmune',
                'Parke Davis',
                'Ciba Geigy',
                'Astra',
                'Schering',
                'Igeneon',
                'Sterix',
                'AbGenix',
                'Imperial Chemical Industries',
                'DOV Pharmaceuticals',
                'Saegis',
                'Daiichi',
                'Kudos',
                'Chiron',
                'Altus',
                'Rinat',
                'SmithKlineBeecham',
                'Sanofi Research',
                'Fujisawa',
                'Hoescht Marrion Roussell',
                'Solvay',
                'Pharmacia and Upjohn',
                'Agouron',
                'Searle',
                'Rousell-Uclaf',
                'Raven',
                'Aventis',
                'Sirtris',
                'Corgentech',
                'Rhone-Poulenc-Rohrer',
                'Tibotec',
                'Millennium',
                'Rhone-Poulenc',
                'Borroughs Wellcome',
                'IDEC',
                'Ayerst Laboratories',
                'Organon',
                'Piramed',
                'Structural Genomix',
                'Tularik',
                'Sugen',
                'Byk Gulden Lomberg Chemische Fabrik',
                'William H. orerR',
                'Winthrop',
                'Texas Biotechnology Corp',
                'ICOS',
                'McNeil',
                'Cita NeuroPharmaceuticals',
                'Upjohn',
                'Sandoz',
                'Allen & Hanburys',
                'Atugen',
                'SmithKlineFrench',
                'Agensys',
                'Yamanouchi',
                'Ecopia',
                'CAT',
                'Sankyo',
                'Mead Johnson',
                'Warner-Lambert',
                'Oxford GlycoSciences',
                'Chugai',
                'Wyeth',
                'Parke-Davis',
                'May & Baker',
                'Pharmacia',
                'Fisons',
                'Hoescht',
                'Chinoin Pharm.',
                'British Biotech'
              ],
              'isUnique': false,
              'comment': 'Name of the previous company the research company was associated with, if applicable.',
              'LLMComment': 'Indicates the prior affiliations of the company, which may provide insights into its history and expertise.'
            }
          ],
          'rowCount': 812,
          'comment': 'This table contains information about research companies involved in various scientific and pharmaceutical activities, including their identifiers, names, and country of operation.',
          'LLMComment': 'The \'public.research_companies\' table serves as a key resource for understanding the landscape of research entities in the pharmaceutical and biotechnology sectors. It includes unique identifiers for both companies and their associated research stems, along with the names of the companies and their countries of operation. This information is crucial for linking research activities to specific organizations, facilitating analysis of global research trends, collaborations, and the impact of these companies on drug development and innovation.'
        },
        {
          'name': 'research_stem',
          'schema': 'public',
          'columns': [
            {
              'name': 'res_stem_id',
              'table': 'research_stem',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 671,
              'isUnique': true,
              'comment': 'Unique identifier for each research stem entry.'
            },
            {
              'name': 'research_stem',
              'table': 'research_stem',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Descriptive name of the research stem, representing a specific area of research.'
            }
          ],
          'rowCount': 670,
          'comment': 'This table stores information about various research stems used in scientific studies, identified by a unique ID.',
          'LLMComment': 'The \'public.research_stem\' table is a key component in the database that catalogs different research stems, which are essential for categorizing and organizing research activities in the scientific domain. Each entry in this table is identified by a unique \'res_stem_id\' and includes a descriptive \'research_stem\' field that provides context for the type of research being conducted. This table is likely used in conjunction with other tables related to activities, assays, and drug mechanisms, facilitating a comprehensive understanding of research trends and methodologies in the life sciences.'
        },
        {
          'name': 'site_components',
          'schema': 'public',
          'columns': [
            {
              'name': 'sitecomp_id',
              'table': 'site_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 51624,
              'isUnique': true
            },
            {
              'name': 'site_id',
              'table': 'site_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 3,
              'max': 38631,
              'isUnique': false
            },
            {
              'name': 'component_id',
              'table': 'site_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 5,
              'max': 19723,
              'isUnique': false
            },
            {
              'name': 'domain_id',
              'table': 'site_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 2629,
              'max': 4420,
              'isUnique': false,
              'comment': 'Identifier for the domain associated with the site component, linking to domain-specific data.',
              'LLMComment': 'This column represents the unique identifier for the domain that the site component belongs to, facilitating domain-related queries.'
            },
            {
              'name': 'site_residues',
              'table': 'site_components',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null'
              ],
              'isUnique': false,
              'comment': 'A string representing specific residues or attributes related to the site component, potentially used for categorization or analysis.',
              'LLMComment': 'This column contains a variable-length character string that may include specific residues or attributes pertinent to the site component, useful for detailed analysis.'
            }
          ],
          'rowCount': 2537,
          'comment': 'This table stores the relationships between various site components and their associated sites, components, and domains, along with any relevant site residues.',
          'LLMComment': 'The `public.site_components` table is a crucial part of the database schema that links specific components to their respective sites and domains within a biological or chemical context. Each entry in this table represents a unique association identified by `sitecomp_id`, which connects a `site_id` (indicating the site where the component is used), a `component_id` (referring to the specific component), and a `domain_id` (which may represent a biological domain or functional area). The `site_residues` column provides additional information about the residues present at the site, which can be important for understanding the functional characteristics of the components in various biological assays or therapeutic applications. This table is essential for researchers and developers working with biological data, as it helps in mapping out the interactions and functionalities of different components within a given site.'
        },
        {
          'name': 'source',
          'schema': 'public',
          'columns': [
            {
              'name': 'src_id',
              'table': 'source',
              'schema': 'public',
              'type': 'integer',
              'min': 0,
              'max': 69,
              'isUnique': true
            },
            {
              'name': 'src_description',
              'table': 'source',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Millipore Kinase Screening (DEPRECATED - MERGED WITH SRC_ID = 1',
                'PDBe Ligands',
                'BindingDB Database',
                'Curated Drug Metabolism Pathways',
                'TP-search Transporter Database',
                'Deposited Supplementary Bioactivity Data',
                'HESi',
                'GSK Kinetoplastid Screening',
                'AstraZeneca Deposited Data',
                'Zimmerman Lab Biotransformation data Dec 2023',
                'Winzeler Lab Plasmodium Screening Data',
                'MMV Malaria HGL',
                'FDA Approved Drug Products with Therapeutic Equivalence Evaluations (Orange Book)',
                'International Nonproprietary Names',
                'EUbOPEN Chemogenomic Library',
                'St Jude Malaria Screening',
                'Open TG-GATEs',
                'Open Source Malaria Screening',
                'Drugs for Neglected Diseases Initiative (DNDi)',
                'Cardiff Schistosomiasis Dataset 2023',
                'St Jude Leishmania Screening',
                'University of Dundee - T. cruzi data',
                'Patent Bioactivity Data',
                'IMI-CARE SARS-CoV-2 Data',
                'HeCaToS Compounds',
                'Curated Drug Pharmacokinetic Data',
                'PubChem BioAssays',
                'Literature data from EUbOPEN Chemogenomic Library',
                'External Project Compounds',
                'Active Ingredient of a Prodrug',
                'SARS-CoV-2 Screening Data 2020-21',
                'British National Formulary',
                'Gene Expression Atlas Compounds',
                'EU-OPENSCREEN dataset',
                'GSK Tuberculosis Screening',
                'GSK Published Kinase Inhibitor Set',
                'Undefined',
                'Fraunhofer HDAC6',
                'FDA New Molecular Entity and New Therapeutic Biological Product Approvals (New FDA Drugs)',
                'K4DD Project',
                'Scientific Literature',
                'Kuster lab chemical proteomics drug profiling',
                'CO-ADD antimicrobial screening data',
                'Donated Chemical Probes - SGC Frankfurt',
                'GSK Malaria Screening',
                'RESOLUTE - Research Empowerment\non Solute Carriers',
                'USP Dictionary of USAN and International Drug Names',
                'FDA Approval Packages',
                'WHO-TDR Malaria Screening',
                'Clinical Candidates',
                'Published Kinase Inhibitor Set 2',
                'Withdrawn Drugs',
                'Gates Library compound collection',
                'WHO Anatomical Therapeutic Chemical Classification',
                'DrugMatrix',
                'Novartis Malaria Screening',
                'Salvensis and LSHTM Schistosomiasis screening data',
                'European Medicines Agency',
                'Harvard Malaria Screening',
                'MMV Malaria Box',
                'Sanger Institute Genomics of Drug Sensitivity in Cancer',
                'MMV Pathogen Box',
                'Guide to Receptors and Channels (DEPRECATED)'
              ],
              'isUnique': true,
              'comment': 'A detailed description of the source, including its purpose and any relevant notes about its status or categorization.',
              'LLMComment': 'This column provides a comprehensive overview of the source\'s function and historical context, including any deprecated statuses or mergers with other sources.'
            },
            {
              'name': 'src_short_name',
              'table': 'source',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'HESI',
                'USP/USAN',
                'EXT. PROJECT CPDS',
                'TG_GATES',
                'WITHDRAWN',
                'BINDINGDB',
                'GSK_PKIS',
                'PATENT',
                'GATES_LIBRARY',
                'CSD23',
                'K4DD',
                'BNF',
                'LITERATURE',
                'EU-OPENSCREEN',
                'EMA',
                'LIT_EUBOPEN_CGL',
                'RESOLUTE',
                'SUPPLEMENTARY',
                'MMV_PBOX',
                'ASTRAZENECA',
                'MMV_MALARIA_HGL',
                'HDAC6',
                'FDA_ORANGE_BOOK',
                'UNDEFINED',
                'COADD',
                'GSK_TB',
                'PDBE',
                'EUBOPEN_CGL',
                'DUNDEE_T_CRUZI',
                'FDA_NEW_DRUGS',
                'SARS_COV_2',
                'ATC',
                'HECATOS',
                'HARVARD',
                'DONATED_PROBES',
                'NOVARTIS',
                'INN',
                'MILLIPORE',
                'FDA_APPROVAL',
                'MMV_MBOX',
                'DNDI',
                'ST_JUDE',
                'DRUG_PK',
                'SANGER',
                'SALVENSIS_LSHTM',
                'CANDIDATES',
                'GSK_TCAKS',
                'ZIMM_BT_12_23',
                'DRUGMATRIX',
                'ST_JUDE_LEISH',
                'PKIS2',
                'PRODRUG_ACTIVE',
                'CARE',
                'WINZ_PLASMO',
                'OSM',
                'WHO_TDR',
                'GRAC',
                'METABOLISM',
                'GSK_TCMDC',
                'PUBCHEM_BIOASSAY',
                'ATLAS',
                'TUM_PROTEOMIC_KUSTER',
                'TP_TRANSPORTER'
              ],
              'isUnique': true,
              'comment': 'A concise identifier for the source, often used for quick reference or categorization.',
              'LLMComment': 'This column contains abbreviated names for the sources, which are useful for quick identification and reference in various contexts.'
            }
          ],
          'rowCount': 63,
          'comment': 'This table contains information about various sources used in the database, including their unique identifiers, descriptions, and short names.',
          'LLMComment': 'The \'public.source\' table serves as a reference for different sources of data within the database. Each entry is identified by a unique \'src_id\' and includes a detailed \'src_description\' and a concise \'src_short_name\'. This table is crucial for understanding the origin of data entries across various related tables, such as assays, compounds, and biotherapeutics, thereby providing context for data provenance and facilitating data integration and analysis.'
        },
        {
          'name': 'structural_alert_sets',
          'schema': 'public',
          'columns': [
            {
              'name': 'alert_set_id',
              'table': 'structural_alert_sets',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 6,
              'isUnique': true
            },
            {
              'name': 'set_name',
              'table': 'structural_alert_sets',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'BMS',
                'Glaxo',
                'Dundee',
                'MLSMR',
                'PAINS'
              ],
              'isUnique': true,
              'comment': 'Name of the alert set, categorized into specific groups.',
              'LLMComment': 'Descriptive name for the alert set, indicating its category such as BMS, Glaxo, etc.'
            },
            {
              'name': 'priority',
              'table': 'structural_alert_sets',
              'schema': 'public',
              'type': 'smallint',
              'min': 3,
              'max': 8,
              'isUnique': true,
              'comment': 'Numerical value indicating the urgency or importance of the alert set, with a defined range.',
              'LLMComment': 'Priority level assigned to the alert set, indicating its significance on a scale from 3 to 8.'
            }
          ],
          'rowCount': 5,
          'comment': 'This table stores sets of structural alerts used in the context of chemical and biological research, including their identifiers, names, and associated priorities.',
          'LLMComment': 'The \'structural_alert_sets\' table is a key component in the domain of chemical and biological research, specifically focusing on the categorization of structural alerts. Each entry in this table represents a unique set of alerts that signal potential issues or noteworthy features in molecular structures. The \'alert_set_id\' serves as a unique identifier for each alert set, while \'set_name\' provides a descriptive label for easy reference. The \'priority\' column indicates the importance or urgency of the alert set, which can be crucial for researchers and developers in assessing the safety and efficacy of compounds. This table interacts with various other tables related to compounds, assays, and biological activities, making it integral to the overall data architecture in drug discovery and development.'
        },
        {
          'name': 'structural_alerts',
          'schema': 'public',
          'columns': [
            {
              'name': 'alert_id',
              'table': 'structural_alerts',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 1103,
              'isUnique': true
            },
            {
              'name': 'alert_set_id',
              'table': 'structural_alerts',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 6,
              'isUnique': false,
              'comment': 'Identifier for the set to which this alert belongs, indicating grouping or categorization of alerts.',
              'LLMComment': 'This column represents the unique identifier for a collection of alerts, allowing for organization and management of related alerts.'
            },
            {
              'name': 'alert_name',
              'table': 'structural_alerts',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Descriptive name of the alert, providing context about the nature of the alert.',
              'LLMComment': 'This field contains the name of the alert, which helps users understand the specific issue or condition being flagged.'
            },
            {
              'name': 'smarts',
              'table': 'structural_alerts',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'A string representing the SMARTS notation for the chemical structure associated with the alert.',
              'LLMComment': 'This column holds the SMARTS (SMiles ARbitrary Target Specification) representation, which is a way to describe molecular structures in a format that can be processed by cheminformatics software.'
            }
          ],
          'rowCount': 936,
          'comment': 'This table stores information about structural alerts related to compounds, including unique identifiers and descriptions of the alerts.',
          'LLMComment': 'The \'structural_alerts\' table is a key component in the domain of cheminformatics and drug discovery. It contains 936 entries that detail various structural alerts, which are specific molecular features that may indicate potential safety or efficacy issues in drug candidates. Each alert is identified by a unique \'alert_id\' and is associated with an \'alert_set_id\' that groups related alerts. The \'alert_name\' provides a descriptive label for the alert, while the \'smarts\' column contains the SMARTS representation, a notation used to describe molecular structures. This table is crucial for researchers and developers in the pharmaceutical industry as it helps in identifying and mitigating risks associated with chemical compounds during the drug development process.'
        },
        {
          'name': 'target_components',
          'schema': 'public',
          'columns': [
            {
              'name': 'tid',
              'table': 'target_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 121145,
              'isUnique': false
            },
            {
              'name': 'component_id',
              'table': 'target_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 20524,
              'isUnique': false
            },
            {
              'name': 'targcomp_id',
              'table': 'target_components',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 23949,
              'isUnique': true,
              'comment': 'Unique identifier for each target component in the database.',
              'LLMComment': 'This is a unique identifier assigned to each target component, ensuring that each entry can be distinctly referenced.'
            },
            {
              'name': 'homologue',
              'table': 'target_components',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 2,
              'isUnique': false,
              'comment': 'Indicates the homology status of the component, where 0 = no homologue, 1 = partial homologue, 2 = full homologue.',
              'LLMComment': 'This column represents the homology classification of the component, providing insights into its evolutionary relationships with other components.'
            }
          ],
          'rowCount': 14394,
          'comment': 'This table stores information about target components associated with various biological entities, including their identifiers and homologous relationships.',
          'LLMComment': 'The \'public.target_components\' table is a crucial part of a biological database that catalogs target components, which are essential for understanding interactions in drug discovery and development. It contains 14,394 rows, each representing a unique target component identified by \'tid\' (target ID), \'component_id\' (specific component identifier), \'targcomp_id\' (target component identifier), and \'homologue\' (indicating related components). This table is interconnected with various other tables in the schema, such as \'assays\', \'biotherapeutics\', and \'drug_indication\', providing a comprehensive view of how these components relate to biological activities and therapeutic applications.'
        },
        {
          'name': 'target_dictionary',
          'schema': 'public',
          'columns': [
            {
              'name': 'tid',
              'table': 'target_dictionary',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 121168,
              'isUnique': true,
              'comment': 'Unique identifier for each target in the dictionary.'
            },
            {
              'name': 'target_type',
              'table': 'target_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'PROTEIN NUCLEIC-ACID COMPLEX',
                'SELECTIVITY GROUP',
                'LIPID',
                'ADMET',
                'NO TARGET',
                'ORGANISM',
                'OLIGOSACCHARIDE',
                'NUCLEIC-ACID',
                'CELL-LINE',
                'SMALL MOLECULE',
                'TISSUE',
                'UNCHECKED',
                'PROTEIN FAMILY',
                'NON-MOLECULAR',
                'CHIMERIC PROTEIN',
                'PROTEIN COMPLEX',
                'PHENOTYPE',
                'SINGLE PROTEIN',
                'UNKNOWN',
                'PROTEIN-PROTEIN INTERACTION',
                'SUBCELLULAR',
                'MACROMOLECULE',
                'METAL',
                'PROTEIN COMPLEX GROUP'
              ],
              'isUnique': false,
              'comment': 'Categorizes the type of biological target, such as protein or nucleic acid.',
              'LLMComment': 'Indicates the biological classification of the target, which can influence drug design and research.'
            },
            {
              'name': 'pref_name',
              'table': 'target_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Preferred name for the target, used for identification.',
              'LLMComment': 'The commonly accepted name for the target, which may differ from other identifiers.'
            },
            {
              'name': 'tax_id',
              'table': 'target_dictionary',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 3052763,
              'isUnique': false,
              'comment': 'Taxonomic identifier for the organism associated with the target.',
              'LLMComment': 'Links the target to a specific organism in biological classification systems.'
            },
            {
              'name': 'organism',
              'table': 'target_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Name of the organism from which the target is derived.',
              'LLMComment': 'Specifies the biological source of the target, important for understanding its function.'
            },
            {
              'name': 'chembl_id',
              'table': 'target_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'semanticType': 'CHEMBL_ID',
              'comment': 'Unique identifier for the target in the ChEMBL database.',
              'LLMComment': 'Facilitates cross-referencing with the ChEMBL database for drug discovery data.'
            },
            {
              'name': 'species_group_flag',
              'table': 'target_dictionary',
              'schema': 'public',
              'type': 'smallint',
              'min': 0,
              'max': 1,
              'isUnique': false,
              'comment': 'Indicates whether the target belongs to a specific species group (0 or 1).',
              'LLMComment': 'Used to categorize targets based on their species affiliation, which can affect drug interactions.'
            }
          ],
          'rowCount': 15598,
          'comment': 'This table contains a dictionary of biological targets, including their identifiers, types, and associated organisms.',
          'LLMComment': 'The `public.target_dictionary` table serves as a comprehensive reference for biological targets in the database. It includes essential information such as the target\'s unique identifier (`tid`), its type (`target_type`), preferred name (`pref_name`), and the associated organism\'s taxonomy ID (`tax_id`). This table is crucial for linking various biological activities and assays to specific targets, facilitating research in pharmacology and drug discovery. The presence of a `chembl_id` indicates integration with the ChEMBL database, which is vital for accessing chemical and biological data related to drug development. The `species_group_flag` helps categorize targets based on their biological classification, enhancing the understanding of target diversity across different organisms.'
        },
        {
          'name': 'target_relations',
          'schema': 'public',
          'columns': [
            {
              'name': 'tid',
              'table': 'target_relations',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 121140,
              'isUnique': false
            },
            {
              'name': 'relationship',
              'table': 'target_relations',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'EQUIVALENT TO',
                'OVERLAPS WITH',
                'SUBSET OF',
                'SUPERSET OF'
              ],
              'isUnique': false,
              'comment': 'Describes the type of relationship between the target and the related target, with specific categories indicating the nature of the connection.',
              'LLMComment': 'This column specifies the nature of the relationship, such as whether one target is equivalent to, overlaps with, is a subset of, or is a superset of another target.'
            },
            {
              'name': 'related_tid',
              'table': 'target_relations',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 121140,
              'isUnique': false,
              'comment': 'References the target ID that is related to the primary target identified by \'tid\'.',
              'LLMComment': 'This column holds the ID of the target that has a specified relationship with the primary target.'
            },
            {
              'name': 'targrel_id',
              'table': 'target_relations',
              'schema': 'public',
              'type': 'bigint',
              'min': 489398,
              'max': 544717,
              'isUnique': true,
              'comment': 'A unique identifier for each relationship entry in the table, ensuring that each relationship can be distinctly referenced.',
              'LLMComment': 'This is a unique identifier for the relationship record, allowing for precise identification of each relationship in the dataset.'
            }
          ],
          'rowCount': 49763,
          'comment': 'This table stores relationships between different targets in the database, allowing for the mapping of related targets based on various types of relationships.',
          'LLMComment': 'The \'public.target_relations\' table is a crucial component in the biological and pharmacological domain, as it captures the relationships between various biological targets (identified by \'tid\') and their related targets (identified by \'related_tid\'). Each relationship is characterized by a type (stored in the \'relationship\' column), which can denote different kinds of interactions or associations. This table facilitates the understanding of how different targets are interconnected, which is essential for drug discovery, understanding biological pathways, and analyzing the effects of compounds on various targets. The \'targrel_id\' serves as a unique identifier for each relationship entry, ensuring data integrity and ease of reference.'
        },
        {
          'name': 'target_type',
          'schema': 'public',
          'columns': [
            {
              'name': 'target_type',
              'table': 'target_type',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'METAL',
                'PROTEIN COMPLEX GROUP',
                'UNDEFINED',
                'PROTEIN NUCLEIC-ACID COMPLEX',
                'SELECTIVITY GROUP',
                'LIPID',
                'ADMET',
                'ORGANISM',
                'NO TARGET',
                'PROTEIN',
                'OLIGOSACCHARIDE',
                'NUCLEIC-ACID',
                'CELL-LINE',
                'MOLECULAR',
                'SMALL MOLECULE',
                'UNCHECKED',
                'TISSUE',
                'PROTEIN FAMILY',
                'NON-MOLECULAR',
                'PROTEIN COMPLEX',
                'CHIMERIC PROTEIN',
                'SINGLE PROTEIN',
                'PHENOTYPE',
                'UNKNOWN',
                'PROTEIN-PROTEIN INTERACTION',
                'SUBCELLULAR',
                'MACROMOLECULE'
              ],
              'isUnique': true,
              'comment': 'Type of target, indicating the category it belongs to.',
              'LLMComment': 'Defines the classification of the target, such as METAL or PROTEIN COMPLEX GROUP.'
            },
            {
              'name': 'target_desc',
              'table': 'target_type',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Description of the target type, providing additional context or details.',
              'LLMComment': 'A unique textual description that elaborates on the target type.'
            },
            {
              'name': 'parent_type',
              'table': 'target_type',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                'UNDEFINED',
                'PROTEIN',
                'MOLECULAR',
                'NON-MOLECULAR'
              ],
              'isUnique': false,
              'comment': 'Higher-level category that the target type falls under, indicating its broader classification.',
              'LLMComment': 'Specifies the parent category of the target, such as PROTEIN or MOLECULAR.'
            }
          ],
          'rowCount': 27,
          'comment': 'This table defines various types of biological targets used in drug discovery and development, including their descriptions and hierarchical relationships.',
          'LLMComment': 'The \'public.target_type\' table is essential for categorizing biological targets in pharmacology and drug development. It contains 27 entries that specify different target types, their descriptions, and any parent-child relationships among them. This hierarchical structure aids in understanding how various targets relate to one another, which is crucial for researchers and AI systems analyzing drug interactions, mechanisms of action, and therapeutic strategies. By organizing target types, this table supports the broader schema of drug discovery, linking to other tables that detail actions, activities, and compounds associated with these targets.'
        },
        {
          'name': 'tissue_dictionary',
          'schema': 'public',
          'columns': [
            {
              'name': 'tissue_id',
              'table': 'tissue_dictionary',
              'schema': 'public',
              'type': 'bigint',
              'min': 2,
              'max': 10000192,
              'isUnique': true,
              'comment': 'Unique identifier for each tissue entry.'
            },
            {
              'name': 'uberon_id',
              'table': 'tissue_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier from the Uberon ontology, used for anatomical classification.',
              'LLMComment': 'This ID links to the Uberon ontology, which provides a standardized vocabulary for anatomical structures.'
            },
            {
              'name': 'pref_name',
              'table': 'tissue_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Preferred name for the tissue, used for display purposes.',
              'LLMComment': 'The primary name used to refer to the tissue in various contexts.'
            },
            {
              'name': 'efo_id',
              'table': 'tissue_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier from the Experimental Factor Ontology, used for categorizing experimental factors.',
              'LLMComment': 'This ID is used to classify tissues in relation to experimental factors in biomedical research.'
            },
            {
              'name': 'chembl_id',
              'table': 'tissue_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'semanticType': 'CHEMBL_ID',
              'comment': 'Identifier from the ChEMBL database, linking to chemical entities related to the tissue.',
              'LLMComment': 'This ID connects to the ChEMBL database, which contains information on bioactive drug-like small molecules.'
            },
            {
              'name': 'bto_id',
              'table': 'tissue_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier from the BRENDA Tissue Ontology, used for tissue classification.',
              'LLMComment': 'This ID is part of the BRENDA Tissue Ontology, which provides a systematic classification of tissues.'
            },
            {
              'name': 'caloha_id',
              'table': 'tissue_dictionary',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier from the CALOHA database, used for linking to histological data.',
              'LLMComment': 'This ID links to the CALOHA database, which focuses on histological and anatomical data.'
            }
          ],
          'rowCount': 782,
          'comment': 'This table contains a dictionary of tissues, providing unique identifiers and names for various tissues used in biological and medical research.',
          'LLMComment': 'The `public.tissue_dictionary` table serves as a comprehensive reference for tissues in the context of biomedical research. It includes key identifiers such as `tissue_id`, `uberon_id`, and various other IDs (e.g., `efo_id`, `chembl_id`, `bto_id`, `caloha_id`) that link tissues to different ontologies and databases. This table is crucial for ensuring consistent terminology and facilitating data integration across various biological datasets, making it easier for researchers to access and analyze tissue-related information in studies involving drug discovery, disease mechanisms, and therapeutic development.'
        },
        {
          'name': 'usan_stems',
          'schema': 'public',
          'columns': [
            {
              'name': 'usan_stem_id',
              'table': 'usan_stems',
              'schema': 'public',
              'type': 'bigint',
              'min': 3341,
              'max': 4134,
              'isUnique': true
            },
            {
              'name': 'stem',
              'table': 'usan_stems',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'The actual stem of the drug name, representing a specific chemical structure or pharmacological action.',
              'LLMComment': 'This column contains the stem of the drug name, which is crucial for identifying the chemical structure or pharmacological action associated with the drug.'
            },
            {
              'name': 'subgroup',
              'table': 'usan_stems',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'A classification that further categorizes the stem into a more specific group within the broader category.',
              'LLMComment': 'This column indicates the subgroup classification of the stem, providing additional context for its pharmacological properties.'
            },
            {
              'name': 'annotation',
              'table': 'usan_stems',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Additional notes or comments regarding the stem, which may include historical or usage information.',
              'LLMComment': 'This column contains annotations that provide extra information about the stem, such as its usage, history, or other relevant details.'
            },
            {
              'name': 'stem_class',
              'table': 'usan_stems',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'infix',
                'prefix',
                'combined prefix and suffix',
                'suffix'
              ],
              'isUnique': false,
              'comment': 'The type of stem based on its morphological characteristics, indicating how it is used in drug naming conventions.',
              'LLMComment': 'This column categorizes the stem into types such as infix, prefix, combined prefix and suffix, or suffix, which are important for understanding its role in drug nomenclature.'
            },
            {
              'name': 'major_class',
              'table': 'usan_stems',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'null',
                'GPCR',
                'protease',
                'NR',
                'PDE',
                'ion channel',
                'kinase'
              ],
              'isUnique': false,
              'comment': 'The primary classification of the stem based on its biological target or mechanism of action, such as GPCR or protease.',
              'LLMComment': 'This column specifies the major class of the stem, indicating its biological target or mechanism of action, which is essential for pharmacological categorization.'
            }
          ],
          'rowCount': 698,
          'comment': 'Table containing information about various drug stems used in pharmacological research and development.',
          'LLMComment': 'The \'public.usan_stems\' table is a key component in pharmacological databases, specifically focusing on the classification and categorization of drug stems. Each entry includes a unique identifier (\'usan_stem_id\'), the stem itself, its subgroup, and various annotations that provide context about its classification (\'stem_class\' and \'major_class\'). This table is essential for researchers and AI models to understand the relationships and classifications of drug components, aiding in drug discovery and development processes.'
        },
        {
          'name': 'variant_sequences',
          'schema': 'public',
          'columns': [
            {
              'name': 'variant_id',
              'table': 'variant_sequences',
              'schema': 'public',
              'type': 'bigint',
              'min': -1,
              'max': 2702,
              'isUnique': true,
              'comment': 'Unique identifier for each variant sequence.'
            },
            {
              'name': 'mutation',
              'table': 'variant_sequences',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Describes the specific mutation present in the variant sequence.',
              'LLMComment': 'A textual representation of the mutation, indicating changes in the nucleotide or amino acid sequence.'
            },
            {
              'name': 'accession',
              'table': 'variant_sequences',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Accession number associated with the variant sequence, used for database referencing.',
              'LLMComment': 'A unique identifier assigned to the sequence in biological databases.'
            },
            {
              'name': 'version',
              'table': 'variant_sequences',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 6,
              'isUnique': false,
              'comment': 'Indicates the version of the variant sequence, reflecting updates or changes.',
              'LLMComment': 'A numeric value representing the version of the sequence, useful for tracking revisions.'
            },
            {
              'name': 'isoform',
              'table': 'variant_sequences',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 4,
              'isUnique': false,
              'comment': 'Represents the specific isoform of the protein or gene associated with the variant.',
              'LLMComment': 'A numeric identifier for different isoforms, indicating variations in protein structure.'
            },
            {
              'name': 'sequence',
              'table': 'variant_sequences',
              'schema': 'public',
              'type': 'text',
              'isUnique': false,
              'semanticType': 'Macromolecule',
              'comment': 'The actual nucleotide or amino acid sequence of the variant.',
              'LLMComment': 'A text representation of the biological macromolecule sequence, crucial for analysis.'
            },
            {
              'name': 'organism',
              'table': 'variant_sequences',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'Rhizobium leguminosarum',
                'Mycobacterium leprae',
                'Influenza A virus',
                'Pyrococcus horikoshii',
                'Candida glabrata',
                'Streptococcus pneumoniae TIGR4',
                'Bacillus subtilis',
                'Influenza B virus (strain B/Lee/1940)',
                'null',
                'Human immunodeficiency virus type 1 group M subtype B (isolate HXB2) (HIV-1)',
                'Escherichia coli',
                'Human respiratory syncytial virus',
                'Lactococcus lactis subsp. cremoris (strain MG1363)',
                'Zaire ebolavirus',
                'Pneumocystis carinii',
                'Thermus thermophilus',
                'Human herpesvirus 8',
                'Chlamydomonas reinhardtii',
                'Staphylococcus aureus subsp. aureus MW2',
                'Influenza A virus (strain A/Wilson-Smith/1933 H1N1) (Influenza A virus (strain A/WS/1933 H1N1))',
                'Vibrio cholerae',
                'Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv)',
                'Mycobacterium tuberculosis H37Rv',
                'Human immunodeficiency virus type 1 group M subtype B (isolate HXB2)(HIV-1)',
                'Human T-lymphotropic virus 1',
                'Pseudomonas aeruginosa',
                'Magnetospirillum gryphiswaldense',
                'Cryptosporidium hominis',
                'Arabidopsis thaliana',
                'Escherichia coli K-12',
                'Gallus gallus',
                'Magnaporthe oryzae',
                'Salinispora tropica CNB-440',
                'Rhodococcus sp.',
                'Bacteroides thetaiotaomicron',
                'Candida tropicalis',
                'Escherichia coli O127:H6 str. E2348/69',
                'Acinetobacter baumannii',
                'Hepatitis C virus subtype 1b',
                'Plasmodium falciparum',
                'Leishmania donovani',
                'Plasmodium vivax',
                'Streptomyces antibioticus',
                'Streptococcus pneumoniae R6',
                'Synechocystis sp. PCC 6803',
                'Anopheles gambiae (African malaria mosquito)',
                'Bos taurus',
                'Helicobacter pylori',
                'Candida albicans',
                'Influenza A virus (A/udorn/1972(H3N2))',
                'Plasmodium falciparum (isolate 3D7)',
                'Homo sapiens',
                'Paramecium bursaria Chlorella virus 1',
                'Enterobacteria phage T4',
                'Human immunodeficiency virus',
                'Hepacivirus C',
                'Mycobacterium tuberculosis (strain CDC 1551 / Oshkosh)',
                'Escherichia coli SMS-3-5',
                'Escherichia coli (strain K12)',
                'Bombyx mori',
                'Human enterovirus 71',
                'Rous sarcoma virus',
                'Mus musculus',
                'Hepatitis C virus',
                'Mycobacterium tuberculosis',
                'Aspergillus fumigatus',
                'Saccharomyces cerevisiae',
                'Francisella tularensis subsp. tularensis SCHU S4',
                'Rattus norvegicus',
                'Bacillus anthracis',
                'Staphylococcus aureus',
                'Leishmania major',
                'Issatchenkia orientalis',
                'Zaire ebolavirus (strain Mayinga-76) (ZEBOV) (Zaire Ebola virus)',
                'Influenza A virus (strain A/Udorn/1972 H3N2)',
                'Mycobacterium smegmatis',
                'Human immunodeficiency virus 1'
              ],
              'isUnique': false,
              'comment': 'The organism from which the variant sequence is derived.',
              'LLMComment': 'Categorical data indicating the species, important for understanding the biological context.'
            },
            {
              'name': 'tax_id',
              'table': 'variant_sequences',
              'schema': 'public',
              'type': 'bigint',
              'min': 287,
              'max': 1111708,
              'isUnique': false,
              'comment': 'Taxonomic identifier for the organism, linking to taxonomic databases.',
              'LLMComment': 'A numeric identifier that provides information about the organism\'s classification in biological taxonomy.'
            }
          ],
          'rowCount': 2473,
          'comment': 'This table stores information about various genetic variants, including their sequences and associated metadata.',
          'LLMComment': 'The `variant_sequences` table is a crucial component in the domain of genomics and bioinformatics. It contains detailed records of genetic variants, identified by a unique `variant_id`, along with their specific mutations, sequences, and the organisms they are associated with. This table serves as a reference for researchers and developers working on genetic analysis, variant interpretation, and understanding the biological implications of mutations across different organisms. The inclusion of fields like `accession`, `version`, and `tax_id` provides essential context for each variant, linking them to broader biological databases and ensuring accurate data retrieval and analysis.'
        },
        {
          'name': 'version',
          'schema': 'public',
          'columns': [
            {
              'name': 'name',
              'table': 'version',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'ChEMBL_34',
                'Bioassay Ontology 2.0',
                'RDKit 2022.09.4',
                'Gene Ontology 2024-02-22',
                'InChI v1.06',
                'EFO v3.56.0',
                'ChEMBL_Structure_Pipeline 1.2.0',
                'Swiss-Prot 2024_01',
                'MeSH 2024'
              ],
              'isUnique': true
            },
            {
              'name': 'creation_date',
              'table': 'version',
              'schema': 'public',
              'type': 'timestamp without time zone',
              'isUnique': false,
              'comment': 'The date and time when the version was created.',
              'LLMComment': 'Timestamp indicating when this version entry was generated.'
            },
            {
              'name': 'comments',
              'table': 'version',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': true,
              'comment': 'Additional remarks or notes regarding the version.',
              'LLMComment': 'Unique comments that provide context or details about the version.'
            }
          ],
          'rowCount': 9,
          'comment': 'This table stores version information for various entities in the database, including their names, creation dates, and any associated comments.',
          'LLMComment': 'The \'public.version\' table is crucial for tracking the different versions of entities within the database. It contains three key columns: \'name\', which identifies the version; \'creation_date\', indicating when the version was created; and \'comments\', which can provide additional context or notes about the version. This table helps maintain a historical record of changes and updates across various biological and chemical data entities, ensuring that users can reference specific versions and understand the evolution of the data.'
        },
        {
          'name': 'warning_refs',
          'schema': 'public',
          'columns': [
            {
              'name': 'warnref_id',
              'table': 'warning_refs',
              'schema': 'public',
              'type': 'bigint',
              'min': 11,
              'max': 8376,
              'isUnique': true
            },
            {
              'name': 'warning_id',
              'table': 'warning_refs',
              'schema': 'public',
              'type': 'bigint',
              'min': 1,
              'max': 3784,
              'isUnique': false
            },
            {
              'name': 'ref_type',
              'table': 'warning_refs',
              'schema': 'public',
              'type': 'character varying',
              'categoryValues': [
                'PubMed',
                'USGPO',
                'HIS',
                'HPRA',
                'TGA',
                'DOI',
                'FDA',
                'NICE',
                'Other',
                'EMA',
                'DailyMed',
                'MEDSAFE',
                'Health Canada',
                'MHRA',
                'WHO'
              ],
              'isUnique': false,
              'comment': 'Type of reference, indicating the source of the warning information.',
              'LLMComment': 'Categorical identifier for the source of the reference, such as PubMed or USGPO.'
            },
            {
              'name': 'ref_id',
              'table': 'warning_refs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'comment': 'Identifier for the specific reference within the source type.',
              'LLMComment': 'Unique identifier associated with the reference, used to link to the source.'
            },
            {
              'name': 'ref_url',
              'table': 'warning_refs',
              'schema': 'public',
              'type': 'character varying',
              'isUnique': false,
              'semanticType': 'URL',
              'comment': 'URL linking to the reference source for further information.',
              'LLMComment': 'Web address that directs to the online resource or document related to the warning.'
            }
          ],
          'rowCount': 3900,
          'comment': 'This table stores references related to various warnings associated with drugs or compounds, including their types and URLs for more information.',
          'LLMComment': 'The `public.warning_refs` table is a critical component in the pharmacovigilance domain, serving as a repository for references to warnings linked to drugs or compounds. Each entry includes a unique identifier for the warning reference (`warnref_id`), the associated warning (`warning_id`), the type of reference (`ref_type`), a specific identifier for the reference (`ref_id`), and a URL (`ref_url`) that provides additional context or information about the warning. This structure allows for efficient tracking and retrieval of safety information, which is essential for regulatory compliance and ensuring patient safety.'
        }
      ],
      'comment': 'The public schema contains a comprehensive set of tables related to biological and chemical data, including assays, compounds, drug mechanisms, and various classifications. It serves as a central repository for research and development in biotherapeutics and pharmacology.',
      'LLMComment': 'This schema is designed to support queries related to biological assays, drug classifications, and compound properties, facilitating research in pharmacology and biotherapeutics.'
    }
  ],
  'relations': [
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'src_id'
      ],
      'toSchema': 'public',
      'toTable': 'source',
      'toColumns': [
        'src_id'
      ],
      'comment': 'Each activity is associated with a source, indicating where the activity data originates from.',
      'LLMComment': 'Each activity belongs to one source; one source can have many activities associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'data_validity_comment'
      ],
      'toSchema': 'public',
      'toTable': 'data_validity_lookup',
      'toColumns': [
        'data_validity_comment'
      ],
      'comment': 'The activities table references the data validity lookup table through the data_validity_comment field.',
      'LLMComment': 'Each activity has a data validity comment that corresponds to a description in the data validity lookup table; one data validity comment can be associated with many activities.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The activities table references the molecule_dictionary table through the molregno field, indicating a relationship between activities and molecules.',
      'LLMComment': 'Each activity in the activities table is associated with one molecule from the molecule_dictionary table; one molecule can have many activities associated with it.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'record_id'
      ],
      'toSchema': 'public',
      'toTable': 'compound_records',
      'toColumns': [
        'record_id'
      ],
      'comment': 'Each activity is associated with one compound record through the record_id.',
      'LLMComment': 'Each activity belongs to one compound record; one compound record can have many activities associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'doc_id'
      ],
      'toSchema': 'public',
      'toTable': 'docs',
      'toColumns': [
        'doc_id'
      ],
      'comment': 'Each activity is associated with one document; one document can have many activities.',
      'LLMComment': 'Each activity belongs to one document; one document can have many activities associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'assay_id'
      ],
      'toSchema': 'public',
      'toTable': 'assays',
      'toColumns': [
        'assay_id'
      ],
      'comment': 'Each activity is associated with a specific assay, indicating a one-to-many relationship between assays and activities.',
      'LLMComment': 'Each activity belongs to one assay; one assay can have many activities associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'bao_endpoint'
      ],
      'toSchema': 'public',
      'toTable': 'bioassay_ontology',
      'toColumns': [
        'bao_id'
      ],
      'comment': 'The activities table references the bioassay ontology table through the bao_endpoint and bao_id fields.',
      'LLMComment': 'Each activity in the activities table is associated with a specific bioassay ontology entry, indicating the type of assay or endpoint being measured; one bioassay ontology entry can be referenced by many activities.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'action_type'
      ],
      'toSchema': 'public',
      'toTable': 'action_type',
      'toColumns': [
        'action_type'
      ],
      'comment': 'The activities table references action types through the action_type field.',
      'LLMComment': 'Each activity has a specific action type that categorizes it; one action type can be associated with many activities.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activity_properties',
      'fromColumns': [
        'activity_id'
      ],
      'toSchema': 'public',
      'toTable': 'activities',
      'toColumns': [
        'activity_id'
      ],
      'comment': 'Each activity_property is associated with one activity; one activity can have multiple activity_properties.',
      'LLMComment': 'Each activity_property belongs to one activity; one activity can have many activity_properties.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activity_supp',
      'fromColumns': [
        'smid'
      ],
      'toSchema': 'public',
      'toTable': 'activity_smid',
      'toColumns': [
        'smid'
      ],
      'comment': 'The activity_supp table references the activity_smid table through the smid foreign key.',
      'LLMComment': 'Each activity_supp entry is associated with one activity_smid; one activity_smid can be referenced by many activity_supp entries.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activity_supp_map',
      'fromColumns': [
        'activity_id'
      ],
      'toSchema': 'public',
      'toTable': 'activities',
      'toColumns': [
        'activity_id'
      ],
      'comment': 'The activity_supp_map table maps activities to supplementary data.',
      'LLMComment': 'Each activity_supp_map entry links one activity to supplementary data; one activity can have many supplementary mappings.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activity_supp_map',
      'fromColumns': [
        'smid'
      ],
      'toSchema': 'public',
      'toTable': 'activity_smid',
      'toColumns': [
        'smid'
      ],
      'comment': 'The activity_supp_map table maps activities to supplementary information using a shared smid.',
      'LLMComment': 'Each activity_supp_map entry links to one activity_smid entry; one activity_smid can be associated with many activity_supp_map entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assay_class_map',
      'fromColumns': [
        'assay_id'
      ],
      'toSchema': 'public',
      'toTable': 'assays',
      'toColumns': [
        'assay_id'
      ],
      'comment': 'The assay_class_map table maps assays to their respective classes.',
      'LLMComment': 'Each assay_class_map entry links one assay to a specific assay class; one assay can belong to multiple assay classes through different entries in the assay_class_map.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assay_class_map',
      'fromColumns': [
        'assay_class_id'
      ],
      'toSchema': 'public',
      'toTable': 'assay_classification',
      'toColumns': [
        'assay_class_id'
      ],
      'comment': 'The assay_class_map table links assays to their classifications through the assay_class_id.',
      'LLMComment': 'Each assay_class_map entry associates one assay with one assay classification; one assay classification can be associated with many assay_class_map entries.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assay_parameters',
      'fromColumns': [
        'assay_id'
      ],
      'toSchema': 'public',
      'toTable': 'assays',
      'toColumns': [
        'assay_id'
      ],
      'comment': 'Each assay_parameter is associated with one assay; one assay can have multiple assay_parameters.',
      'LLMComment': 'Each assay_parameter belongs to one assay; one assay has many assay_parameters.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'cell_id'
      ],
      'toSchema': 'public',
      'toTable': 'cell_dictionary',
      'toColumns': [
        'cell_id'
      ],
      'comment': 'Each assay references a specific cell from the cell dictionary.',
      'LLMComment': 'Each assay belongs to one cell; one cell can be associated with many assays.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'bao_format'
      ],
      'toSchema': 'public',
      'toTable': 'bioassay_ontology',
      'toColumns': [
        'bao_id'
      ],
      'comment': 'The assays table references the bioassay ontology table through the bao_format and bao_id fields.',
      'LLMComment': 'Each assay has a specific format defined in the bioassay ontology; one bioassay ontology can be associated with many assays.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'tissue_id'
      ],
      'toSchema': 'public',
      'toTable': 'tissue_dictionary',
      'toColumns': [
        'tissue_id'
      ],
      'comment': 'Each assay references a tissue type from the tissue dictionary.',
      'LLMComment': 'Each assay belongs to one tissue type; one tissue type can be associated with many assays.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'variant_id'
      ],
      'toSchema': 'public',
      'toTable': 'variant_sequences',
      'toColumns': [
        'variant_id'
      ],
      'comment': 'Each assay references a variant sequence through the variant_id.',
      'LLMComment': 'Each assay belongs to one variant sequence; one variant sequence can have many assays associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'src_id'
      ],
      'toSchema': 'public',
      'toTable': 'source',
      'toColumns': [
        'src_id'
      ],
      'comment': 'Each assay belongs to one source; one source can have many assays.',
      'LLMComment': 'Each assay is associated with a single source, indicating that the source provides the context or origin for the assay. Conversely, a single source can be linked to multiple assays, reflecting that a source may contribute to various assays over time.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'tid'
      ],
      'toSchema': 'public',
      'toTable': 'target_dictionary',
      'toColumns': [
        'tid'
      ],
      'comment': 'Each assay belongs to one target; one target can have many assays.',
      'LLMComment': 'Each assay is associated with a single target, and a single target can be linked to multiple assays, indicating a one-to-many relationship between targets and assays.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'relationship_type'
      ],
      'toSchema': 'public',
      'toTable': 'relationship_type',
      'toColumns': [
        'relationship_type'
      ],
      'comment': 'The assays table references the relationship_type table through the relationship_type field.',
      'LLMComment': 'Each assay belongs to one relationship type; one relationship type can be associated with many assays.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'confidence_score'
      ],
      'toSchema': 'public',
      'toTable': 'confidence_score_lookup',
      'toColumns': [
        'confidence_score'
      ],
      'comment': 'The assays table references the confidence_score_lookup table through the confidence_score field.',
      'LLMComment': 'Each assay has a confidence score that is defined in the confidence_score_lookup table; one confidence score can be associated with many assays.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'curated_by'
      ],
      'toSchema': 'public',
      'toTable': 'curation_lookup',
      'toColumns': [
        'curated_by'
      ],
      'comment': 'The curated_by field in the assays table references the curated_by field in the curation_lookup table.',
      'LLMComment': 'Each assay is curated by a specific curator, and each curator can curate multiple assays.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'chembl_id'
      ],
      'toSchema': 'public',
      'toTable': 'chembl_id_lookup',
      'toColumns': [
        'chembl_id'
      ],
      'comment': 'The assays table references the chembl_id_lookup table through the chembl_id field.',
      'LLMComment': 'Each assay has a unique chembl_id that corresponds to an entry in the chembl_id_lookup table; one chembl_id_lookup can be associated with many assays.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'doc_id'
      ],
      'toSchema': 'public',
      'toTable': 'docs',
      'toColumns': [
        'doc_id'
      ],
      'comment': 'Each assay belongs to one document; one document can have many assays.',
      'LLMComment': 'Each assay belongs to one document; one document can have many assays.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'assays',
      'fromColumns': [
        'assay_type'
      ],
      'toSchema': 'public',
      'toTable': 'assay_type',
      'toColumns': [
        'assay_type'
      ],
      'comment': 'Each assay is associated with one assay type; one assay type can be associated with many assays.',
      'LLMComment': 'Each assay belongs to one assay type; one assay type can have many assays associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'binding_sites',
      'fromColumns': [
        'tid'
      ],
      'toSchema': 'public',
      'toTable': 'target_dictionary',
      'toColumns': [
        'tid'
      ],
      'comment': 'Each binding_site is associated with one target; one target can have many binding_sites.',
      'LLMComment': 'Each binding_site belongs to one target; one target can have many binding_sites.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'biotherapeutic_components',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'biotherapeutics',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The biotherapeutic_components table references the biotherapeutics table through the molregno field, indicating a relationship between biotherapeutic components and their corresponding biotherapeutics.',
      'LLMComment': 'Each biotherapeutic_component belongs to one biotherapeutic; one biotherapeutic can have many biotherapeutic_components.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'biotherapeutic_components',
      'fromColumns': [
        'component_id'
      ],
      'toSchema': 'public',
      'toTable': 'bio_component_sequences',
      'toColumns': [
        'component_id'
      ],
      'comment': 'The biotherapeutic_components table references the bio_component_sequences table through the component_id field.',
      'LLMComment': 'Each biotherapeutic_component is associated with one bio_component_sequence; one bio_component_sequence can be associated with many biotherapeutic_components.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'biotherapeutics',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The biotherapeutics table references the molecule_dictionary table through the molregno field.',
      'LLMComment': 'Each biotherapeutic is associated with one molecule in the molecule dictionary; one molecule can be referenced by many biotherapeutics.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'cell_dictionary',
      'fromColumns': [
        'chembl_id'
      ],
      'toSchema': 'public',
      'toTable': 'chembl_id_lookup',
      'toColumns': [
        'chembl_id'
      ],
      'comment': 'The cell_dictionary table references chembl_id_lookup through the chembl_id field.',
      'LLMComment': 'Each entry in the cell_dictionary is associated with one entry in the chembl_id_lookup; one chembl_id_lookup can be referenced by many entries in the cell_dictionary.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'component_class',
      'fromColumns': [
        'protein_class_id'
      ],
      'toSchema': 'public',
      'toTable': 'protein_classification',
      'toColumns': [
        'protein_class_id'
      ],
      'comment': 'The component_class table references the protein_classification table through the protein_class_id foreign key.',
      'LLMComment': 'Each component_class is associated with one protein_classification; one protein_classification can have many component_classes.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'component_class',
      'fromColumns': [
        'component_id'
      ],
      'toSchema': 'public',
      'toTable': 'component_sequences',
      'toColumns': [
        'component_id'
      ],
      'comment': 'The component_class table references the component_sequences table through the component_id foreign key.',
      'LLMComment': 'Each component_class entry is associated with one component_sequence; one component_sequence can belong to many component_classes.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'component_domains',
      'fromColumns': [
        'component_id'
      ],
      'toSchema': 'public',
      'toTable': 'component_sequences',
      'toColumns': [
        'component_id'
      ],
      'comment': 'The component_domains table references the component_sequences table through the component_id foreign key.',
      'LLMComment': 'Each component_domain is associated with one component_sequence; one component_sequence can have many component_domains.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'component_domains',
      'fromColumns': [
        'domain_id'
      ],
      'toSchema': 'public',
      'toTable': 'domains',
      'toColumns': [
        'domain_id'
      ],
      'comment': 'The component_domains table references the domains table through the domain_id foreign key.',
      'LLMComment': 'Each component_domain is associated with one domain; one domain can have many component_domains associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'component_go',
      'fromColumns': [
        'component_id'
      ],
      'toSchema': 'public',
      'toTable': 'component_sequences',
      'toColumns': [
        'component_id'
      ],
      'comment': 'The component_go table references the component_sequences table using component_id as a foreign key.',
      'LLMComment': 'Each component_go entry is associated with one component_sequence; one component_sequence can have many component_go entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'component_go',
      'fromColumns': [
        'go_id'
      ],
      'toSchema': 'public',
      'toTable': 'go_classification',
      'toColumns': [
        'go_id'
      ],
      'comment': 'The component_go table references the go_classification table using go_id as a foreign key.',
      'LLMComment': 'Each component_go entry is associated with one go_classification entry; one go_classification can be associated with many component_go entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'component_synonyms',
      'fromColumns': [
        'component_id'
      ],
      'toSchema': 'public',
      'toTable': 'component_sequences',
      'toColumns': [
        'component_id'
      ],
      'comment': 'The component_synonyms table references the component_sequences table through the component_id foreign key.',
      'LLMComment': 'Each component_synonym is associated with one component_sequence; one component_sequence can have many component_synonyms.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'compound_properties',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The compound_properties table references the molecule_dictionary table through the molregno field.',
      'LLMComment': 'Each compound_property entry corresponds to one molecule in the molecule_dictionary; one molecule can have many compound_properties associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'compound_records',
      'fromColumns': [
        'src_id'
      ],
      'toSchema': 'public',
      'toTable': 'source',
      'toColumns': [
        'src_id'
      ],
      'comment': 'The compound_records table references the source table through the src_id foreign key.',
      'LLMComment': 'Each compound_record is associated with one source; one source can have many compound_records.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'compound_records',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Each compound_record is associated with one molecule in the molecule_dictionary.',
      'LLMComment': 'Each compound_record belongs to one molecule; one molecule can have many compound_records associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'compound_records',
      'fromColumns': [
        'doc_id'
      ],
      'toSchema': 'public',
      'toTable': 'docs',
      'toColumns': [
        'doc_id'
      ],
      'comment': 'Each compound_record is associated with one document; one document can have many compound_records.',
      'LLMComment': 'Each compound_record belongs to one doc; one doc can have many compound_records.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'compound_structural_alerts',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The compound_structural_alerts table references the molecule_dictionary table through the molregno field.',
      'LLMComment': 'Each compound_structural_alert belongs to one molecule; one molecule can have many compound_structural_alerts.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'compound_structural_alerts',
      'fromColumns': [
        'alert_id'
      ],
      'toSchema': 'public',
      'toTable': 'structural_alerts',
      'toColumns': [
        'alert_id'
      ],
      'comment': 'The compound_structural_alerts table references the structural_alerts table through the alert_id foreign key.',
      'LLMComment': 'Each compound_structural_alert belongs to one structural_alert; one structural_alert can have many compound_structural_alerts associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'compound_structures',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The compound_structures table references the molecule_dictionary table through the molregno field.',
      'LLMComment': 'Each compound_structure corresponds to one molecule in the molecule_dictionary; one molecule can have many compound_structures.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'defined_daily_dose',
      'fromColumns': [
        'atc_code'
      ],
      'toSchema': 'public',
      'toTable': 'atc_classification',
      'toColumns': [
        'level5'
      ],
      'comment': 'The defined_daily_dose table references the atc_classification table through the atc_code field.',
      'LLMComment': 'Each defined_daily_dose entry is associated with one atc_classification entry; one atc_classification can have many defined_daily_dose entries.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'docs',
      'fromColumns': [
        'chembl_id'
      ],
      'toSchema': 'public',
      'toTable': 'chembl_id_lookup',
      'toColumns': [
        'chembl_id'
      ],
      'comment': 'The docs table references chembl_id_lookup through the chembl_id field.',
      'LLMComment': 'Each document in the docs table is associated with one entry in the chembl_id_lookup table, and each chembl_id_lookup entry can be referenced by multiple documents.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'docs',
      'fromColumns': [
        'src_id'
      ],
      'toSchema': 'public',
      'toTable': 'source',
      'toColumns': [
        'src_id'
      ],
      'comment': 'The docs table references the source table through the src_id foreign key.',
      'LLMComment': 'Each document in the docs table is associated with one source from the source table; one source can be referenced by many documents.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'docs',
      'fromColumns': [
        'chembl_release_id'
      ],
      'toSchema': 'public',
      'toTable': 'chembl_release',
      'toColumns': [
        'chembl_release_id'
      ],
      'comment': 'Each document is associated with a specific ChEMBL release.',
      'LLMComment': 'Each doc belongs to one chembl_release; one chembl_release can have many docs.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_indication',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The drug_indication table references the molecule_dictionary table through the molregno field, indicating a relationship between drug indications and molecules.',
      'LLMComment': 'Each drug indication is associated with one molecule; one molecule can have many drug indications.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_indication',
      'fromColumns': [
        'record_id'
      ],
      'toSchema': 'public',
      'toTable': 'compound_records',
      'toColumns': [
        'record_id'
      ],
      'comment': 'The drug_indication table references the compound_records table through the record_id field.',
      'LLMComment': 'Each drug indication is associated with one compound record; one compound record can have many drug indications.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_mechanism',
      'fromColumns': [
        'variant_id'
      ],
      'toSchema': 'public',
      'toTable': 'variant_sequences',
      'toColumns': [
        'variant_id'
      ],
      'comment': 'The drug_mechanism table references the variant_sequences table through the variant_id field.',
      'LLMComment': 'Each drug_mechanism entry is associated with one variant_sequence; one variant_sequence can be linked to many drug_mechanism entries.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_mechanism',
      'fromColumns': [
        'action_type'
      ],
      'toSchema': 'public',
      'toTable': 'action_type',
      'toColumns': [
        'action_type'
      ],
      'comment': 'The drug_mechanism table references action_type to define the type of action associated with a drug mechanism.',
      'LLMComment': 'Each drug_mechanism has one action_type; one action_type can be associated with many drug_mechanisms.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_mechanism',
      'fromColumns': [
        'site_id'
      ],
      'toSchema': 'public',
      'toTable': 'binding_sites',
      'toColumns': [
        'site_id'
      ],
      'comment': 'The drug_mechanism table references binding_sites through site_id, indicating a relationship between drug mechanisms and their corresponding binding sites.',
      'LLMComment': 'Each drug_mechanism entry is associated with one binding_site; one binding_site can be associated with many drug_mechanism entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_mechanism',
      'fromColumns': [
        'tid'
      ],
      'toSchema': 'public',
      'toTable': 'target_dictionary',
      'toColumns': [
        'tid'
      ],
      'comment': 'The drug_mechanism table references the target_dictionary table through the tid field, indicating a relationship between drug mechanisms and targets.',
      'LLMComment': 'Each drug_mechanism entry is associated with one target from the target_dictionary; one target can have many drug_mechanism entries associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_mechanism',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The drug_mechanism table references the molecule_dictionary table through the molregno field, indicating a relationship between drug mechanisms and their corresponding molecules.',
      'LLMComment': 'Each drug_mechanism entry is associated with one molecule from the molecule_dictionary; one molecule can have many drug_mechanism entries.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_mechanism',
      'fromColumns': [
        'record_id'
      ],
      'toSchema': 'public',
      'toTable': 'compound_records',
      'toColumns': [
        'record_id'
      ],
      'comment': 'The drug_mechanism table references the compound_records table through the record_id field.',
      'LLMComment': 'Each drug_mechanism entry is associated with one compound_record; one compound_record can have many drug_mechanism entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_warning',
      'fromColumns': [
        'record_id'
      ],
      'toSchema': 'public',
      'toTable': 'compound_records',
      'toColumns': [
        'record_id'
      ],
      'comment': 'The drug_warning table references compound_records through record_id, indicating a relationship where each drug warning is associated with a specific compound record.',
      'LLMComment': 'Each drug_warning belongs to one compound_record; one compound_record can have many drug_warnings.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'formulations',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Each formulation is associated with one molecule; one molecule can have many formulations.',
      'LLMComment': 'Each formulation belongs to one molecule, and one molecule can have multiple formulations associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'formulations',
      'fromColumns': [
        'record_id'
      ],
      'toSchema': 'public',
      'toTable': 'compound_records',
      'toColumns': [
        'record_id'
      ],
      'comment': 'Each formulation is associated with one compound record; one compound record can have many formulations.',
      'LLMComment': 'Each formulation belongs to one compound record; one compound record can have many formulations.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'formulations',
      'fromColumns': [
        'product_id'
      ],
      'toSchema': 'public',
      'toTable': 'products',
      'toColumns': [
        'product_id'
      ],
      'comment': 'Each formulation is associated with one product; one product can have many formulations.',
      'LLMComment': 'Each formulation belongs to one product; one product can have many formulations associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'indication_refs',
      'fromColumns': [
        'drugind_id'
      ],
      'toSchema': 'public',
      'toTable': 'drug_indication',
      'toColumns': [
        'drugind_id'
      ],
      'comment': 'Each indication_ref is associated with one drug_indication; one drug_indication can have many indication_refs.',
      'LLMComment': 'Each indication_ref belongs to one drug_indication; one drug_indication can have many indication_refs.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'ligand_eff',
      'fromColumns': [
        'activity_id'
      ],
      'toSchema': 'public',
      'toTable': 'activities',
      'toColumns': [
        'activity_id'
      ],
      'comment': 'Each ligand_eff record is associated with one activity; one activity can have many ligand_eff records.',
      'LLMComment': 'Each ligand_eff belongs to one activity; one activity has many ligand_eff records associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'mechanism_refs',
      'fromColumns': [
        'mec_id'
      ],
      'toSchema': 'public',
      'toTable': 'drug_mechanism',
      'toColumns': [
        'mec_id'
      ],
      'comment': 'The mechanism_refs table references the drug_mechanism table through the mec_id foreign key.',
      'LLMComment': 'Each mechanism_ref belongs to one drug_mechanism; one drug_mechanism can have many mechanism_refs.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'metabolism',
      'fromColumns': [
        'enzyme_tid'
      ],
      'toSchema': 'public',
      'toTable': 'target_dictionary',
      'toColumns': [
        'tid'
      ],
      'comment': 'The metabolism table references the target_dictionary table through the enzyme_tid foreign key.',
      'LLMComment': 'Each metabolism entry is associated with one target from the target_dictionary; one target can be associated with many metabolism entries.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'metabolism',
      'fromColumns': [
        'drug_record_id'
      ],
      'toSchema': 'public',
      'toTable': 'compound_records',
      'toColumns': [
        'record_id'
      ],
      'comment': 'The metabolism table references the compound_records table through drug_record_id and record_id.',
      'LLMComment': 'Each metabolism entry is associated with one compound record; one compound record can be referenced by many metabolism entries.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'metabolism',
      'fromColumns': [
        'substrate_record_id'
      ],
      'toSchema': 'public',
      'toTable': 'compound_records',
      'toColumns': [
        'record_id'
      ],
      'comment': 'The metabolism table references compound_records through substrate_record_id, indicating a relationship where each metabolism entry is associated with a specific compound.',
      'LLMComment': 'Each metabolism entry is associated with one compound as a substrate; one compound can be involved in many metabolism entries.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'metabolism',
      'fromColumns': [
        'metabolite_record_id'
      ],
      'toSchema': 'public',
      'toTable': 'compound_records',
      'toColumns': [
        'record_id'
      ],
      'comment': 'Each metabolism record can reference one compound record as a metabolite.',
      'LLMComment': 'Each metabolism entry is associated with one compound record that represents the metabolite involved in the metabolism process.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'metabolism_refs',
      'fromColumns': [
        'met_id'
      ],
      'toSchema': 'public',
      'toTable': 'metabolism',
      'toColumns': [
        'met_id'
      ],
      'comment': 'The metabolism_refs table references the metabolism table using the met_id foreign key.',
      'LLMComment': 'Each metabolism_ref belongs to one metabolism; one metabolism can have many metabolism_refs.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_atc_classification',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The molecule_atc_classification table references the molecule_dictionary table through the molregno field.',
      'LLMComment': 'Each molecule_atc_classification entry corresponds to one molecule in the molecule_dictionary; one molecule can have multiple ATC classifications.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_atc_classification',
      'fromColumns': [
        'level5'
      ],
      'toSchema': 'public',
      'toTable': 'atc_classification',
      'toColumns': [
        'level5'
      ],
      'comment': 'The molecule_atc_classification table references the atc_classification table through the level5 field.',
      'LLMComment': 'Each molecule_atc_classification entry corresponds to one atc_classification entry based on the level5 classification; one atc_classification can have many molecule_atc_classification entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_dictionary',
      'fromColumns': [
        'chembl_id'
      ],
      'toSchema': 'public',
      'toTable': 'chembl_id_lookup',
      'toColumns': [
        'chembl_id'
      ],
      'comment': 'The chembl_id in molecule_dictionary is a foreign key that references the chembl_id in chembl_id_lookup.',
      'LLMComment': 'Each molecule in the molecule_dictionary is associated with a unique chembl_id, which can be looked up in the chembl_id_lookup table for additional information about that chembl_id.',
      'cardinality': 'one-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_frac_classification',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The molecule_frac_classification table references the molecule_dictionary table through the molregno field.',
      'LLMComment': 'Each molecule_frac_classification entry corresponds to one molecule in the molecule_dictionary; one molecule can have many classifications in molecule_frac_classification.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_frac_classification',
      'fromColumns': [
        'frac_class_id'
      ],
      'toSchema': 'public',
      'toTable': 'frac_classification',
      'toColumns': [
        'frac_class_id'
      ],
      'comment': 'The molecule_frac_classification table references the frac_classification table through the frac_class_id foreign key.',
      'LLMComment': 'Each molecule_frac_classification entry is associated with one frac_classification; one frac_classification can have many molecule_frac_classification entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_hierarchy',
      'fromColumns': [
        'parent_molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The molecule_hierarchy table references the molecule_dictionary table to establish a parent-child relationship between molecules.',
      'LLMComment': 'Each molecule_hierarchy entry references a parent molecule from the molecule_dictionary, indicating that one molecule can have multiple child molecules, while each child molecule has one parent.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_hierarchy',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The molecule_hierarchy table references the molecule_dictionary table through the molregno field, indicating a hierarchical relationship between molecules.',
      'LLMComment': 'Each molecule_hierarchy entry corresponds to a molecule in the molecule_dictionary; one molecule can have multiple hierarchical relationships.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_hierarchy',
      'fromColumns': [
        'active_molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The molecule_hierarchy table references the molecule_dictionary table through the molregno field.',
      'LLMComment': 'Each active molecule in the molecule_hierarchy belongs to one molecule in the molecule_dictionary; one molecule in the molecule_dictionary can have many active molecules in the molecule_hierarchy.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_hrac_classification',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The molecule_hrac_classification table references the molecule_dictionary table through the molregno field.',
      'LLMComment': 'Each molecule_hrac_classification entry corresponds to one molecule in the molecule_dictionary; one molecule can have many hrac classifications.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_hrac_classification',
      'fromColumns': [
        'hrac_class_id'
      ],
      'toSchema': 'public',
      'toTable': 'hrac_classification',
      'toColumns': [
        'hrac_class_id'
      ],
      'comment': 'The molecule_hrac_classification table references the hrac_classification table through hrac_class_id.',
      'LLMComment': 'Each molecule_hrac_classification entry is associated with one hrac_classification entry; one hrac_classification can be associated with many molecule_hrac_classification entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_irac_classification',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The molecule_irac_classification table references the molecule_dictionary table using molregno as a foreign key.',
      'LLMComment': 'Each molecule_irac_classification entry corresponds to one molecule in the molecule_dictionary; one molecule can have many IRAC classifications.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_irac_classification',
      'fromColumns': [
        'irac_class_id'
      ],
      'toSchema': 'public',
      'toTable': 'irac_classification',
      'toColumns': [
        'irac_class_id'
      ],
      'comment': 'The molecule_irac_classification table references the irac_classification table through the irac_class_id foreign key.',
      'LLMComment': 'Each molecule_irac_classification entry is associated with one irac_classification; one irac_classification can be associated with many molecule_irac_classification entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_synonyms',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'The molecule_synonyms table references the molecule_dictionary table using molregno as a foreign key.',
      'LLMComment': 'Each molecule_synonym belongs to one molecule_dictionary; one molecule_dictionary can have many molecule_synonyms.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_synonyms',
      'fromColumns': [
        'res_stem_id'
      ],
      'toSchema': 'public',
      'toTable': 'research_stem',
      'toColumns': [
        'res_stem_id'
      ],
      'comment': 'The molecule_synonyms table references the research_stem table through the res_stem_id foreign key.',
      'LLMComment': 'Each molecule_synonym is associated with one research_stem; one research_stem can have many molecule_synonyms.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'predicted_binding_domains',
      'fromColumns': [
        'activity_id'
      ],
      'toSchema': 'public',
      'toTable': 'activities',
      'toColumns': [
        'activity_id'
      ],
      'comment': 'The predicted_binding_domains table references the activities table through the activity_id foreign key.',
      'LLMComment': 'Each predicted_binding_domain entry is associated with one activity; one activity can have many predicted_binding_domains.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'predicted_binding_domains',
      'fromColumns': [
        'site_id'
      ],
      'toSchema': 'public',
      'toTable': 'binding_sites',
      'toColumns': [
        'site_id'
      ],
      'comment': 'The predicted_binding_domains table references binding_sites through site_id.',
      'LLMComment': 'Each predicted_binding_domain is associated with one binding_site; one binding_site can have many predicted_binding_domains.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'product_patents',
      'fromColumns': [
        'patent_use_code'
      ],
      'toSchema': 'public',
      'toTable': 'patent_use_codes',
      'toColumns': [
        'patent_use_code'
      ],
      'comment': 'The product_patents table references the patent_use_codes table through the patent_use_code field.',
      'LLMComment': 'Each product_patent is associated with one patent_use_code; one patent_use_code can be used by many product_patents.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'product_patents',
      'fromColumns': [
        'product_id'
      ],
      'toSchema': 'public',
      'toTable': 'products',
      'toColumns': [
        'product_id'
      ],
      'comment': 'The product_patents table references the products table through product_id, indicating a relationship between patents and products.',
      'LLMComment': 'Each product_patent is associated with one product; one product can have many product_patents.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'protein_class_synonyms',
      'fromColumns': [
        'protein_class_id'
      ],
      'toSchema': 'public',
      'toTable': 'protein_classification',
      'toColumns': [
        'protein_class_id'
      ],
      'comment': 'The protein_class_synonyms table references the protein_classification table through the protein_class_id foreign key.',
      'LLMComment': 'Each protein_class_synonym belongs to one protein_classification; one protein_classification can have many protein_class_synonyms.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'research_companies',
      'fromColumns': [
        'res_stem_id'
      ],
      'toSchema': 'public',
      'toTable': 'research_stem',
      'toColumns': [
        'res_stem_id'
      ],
      'comment': 'The research_companies table references the research_stem table through the res_stem_id foreign key.',
      'LLMComment': 'Each research company is associated with one research stem; one research stem can have many research companies associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'site_components',
      'fromColumns': [
        'domain_id'
      ],
      'toSchema': 'public',
      'toTable': 'domains',
      'toColumns': [
        'domain_id'
      ],
      'comment': 'The site_components table references the domains table through the domain_id foreign key.',
      'LLMComment': 'Each site_component belongs to one domain; one domain can have many site_components.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'site_components',
      'fromColumns': [
        'component_id'
      ],
      'toSchema': 'public',
      'toTable': 'component_sequences',
      'toColumns': [
        'component_id'
      ],
      'comment': 'The site_components table references the component_sequences table through the component_id foreign key.',
      'LLMComment': 'Each site_component is associated with one component_sequence; one component_sequence can be associated with many site_components.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'site_components',
      'fromColumns': [
        'site_id'
      ],
      'toSchema': 'public',
      'toTable': 'binding_sites',
      'toColumns': [
        'site_id'
      ],
      'comment': 'The site_components table references the binding_sites table through the site_id foreign key.',
      'LLMComment': 'Each site_component belongs to one binding_site; one binding_site can have many site_components.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'structural_alerts',
      'fromColumns': [
        'alert_set_id'
      ],
      'toSchema': 'public',
      'toTable': 'structural_alert_sets',
      'toColumns': [
        'alert_set_id'
      ],
      'comment': 'Each structural_alert belongs to one structural_alert_set; one structural_alert_set can have many structural_alerts.',
      'LLMComment': 'Each structural_alert belongs to one structural_alert_set; one structural_alert_set can have many structural_alerts.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'target_components',
      'fromColumns': [
        'tid'
      ],
      'toSchema': 'public',
      'toTable': 'target_dictionary',
      'toColumns': [
        'tid'
      ],
      'comment': 'The target_components table references the target_dictionary table through the tid foreign key.',
      'LLMComment': 'Each target_component belongs to one target_dictionary; one target_dictionary can have many target_components.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'target_components',
      'fromColumns': [
        'component_id'
      ],
      'toSchema': 'public',
      'toTable': 'component_sequences',
      'toColumns': [
        'component_id'
      ],
      'comment': 'The target_components table references the component_sequences table through the component_id field.',
      'LLMComment': 'Each target_component is associated with one component_sequence; one component_sequence can be associated with many target_components.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'target_dictionary',
      'fromColumns': [
        'target_type'
      ],
      'toSchema': 'public',
      'toTable': 'target_type',
      'toColumns': [
        'target_type'
      ],
      'comment': 'The target_dictionary table references the target_type table through the target_type field.',
      'LLMComment': 'Each entry in the target_dictionary is associated with one target_type; one target_type can be associated with many entries in the target_dictionary.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'target_dictionary',
      'fromColumns': [
        'chembl_id'
      ],
      'toSchema': 'public',
      'toTable': 'chembl_id_lookup',
      'toColumns': [
        'chembl_id'
      ],
      'comment': 'The target_dictionary table references chembl_id_lookup through the chembl_id field.',
      'LLMComment': 'Each target in the target_dictionary is associated with one entry in the chembl_id_lookup; one chembl_id_lookup can be referenced by many targets in the target_dictionary.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'target_relations',
      'fromColumns': [
        'tid'
      ],
      'toSchema': 'public',
      'toTable': 'target_dictionary',
      'toColumns': [
        'tid'
      ],
      'comment': 'The target_relations table references the target_dictionary table using the tid column.',
      'LLMComment': 'Each target_relation is associated with one target from the target_dictionary; one target can have many target_relations.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'target_relations',
      'fromColumns': [
        'related_tid'
      ],
      'toSchema': 'public',
      'toTable': 'target_dictionary',
      'toColumns': [
        'tid'
      ],
      'comment': 'The target_relations table references the target_dictionary table through the related_tid foreign key.',
      'LLMComment': 'Each target_relation references one target in the target_dictionary; one target can have many target_relations.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'tissue_dictionary',
      'fromColumns': [
        'chembl_id'
      ],
      'toSchema': 'public',
      'toTable': 'chembl_id_lookup',
      'toColumns': [
        'chembl_id'
      ],
      'comment': 'The tissue_dictionary table references chembl_id_lookup through the chembl_id field.',
      'LLMComment': 'Each tissue in the tissue_dictionary is associated with a unique chembl_id in the chembl_id_lookup, indicating that each tissue can be linked to its corresponding entry in the chembl_id_lookup.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'warning_refs',
      'fromColumns': [
        'warning_id'
      ],
      'toSchema': 'public',
      'toTable': 'drug_warning',
      'toColumns': [
        'warning_id'
      ],
      'comment': 'Each warning_ref is associated with one drug_warning; one drug_warning can have many warning_refs.',
      'LLMComment': 'Each warning_ref belongs to one drug_warning; one drug_warning can have many warning_refs.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'biotherapeutics',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Both tables reference the same molecule identifier (molregno), indicating a potential relationship between activities and biotherapeutics.',
      'LLMComment': 'Activities may involve biotherapeutics identified by the same molregno.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'public',
      'fromTable': 'activities',
      'fromColumns': [
        'assay_id'
      ],
      'toSchema': 'public',
      'toTable': 'assay_parameters',
      'toColumns': [
        'assay_id'
      ],
      'comment': 'Both tables reference assay_id, indicating that activities can be linked to their corresponding assay parameters.',
      'LLMComment': 'Activities can be linked to assay parameters through the assay_id.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_indication',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Both tables reference the same molecule identifier (molregno), indicating a potential relationship between drug indications and molecules.',
      'LLMComment': 'Drug indications can be linked to molecules through the molregno.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_mechanism',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Both tables reference the same molecule identifier (molregno), indicating a potential relationship between drug mechanisms and molecules.',
      'LLMComment': 'Drug mechanisms can be linked to molecules through the molregno.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'public',
      'fromTable': 'drug_warning',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Both tables reference the same molecule identifier (molregno), indicating a potential relationship between drug warnings and molecules.',
      'LLMComment': 'Drug warnings can be linked to molecules through the molregno.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'public',
      'fromTable': 'formulations',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Both tables reference the same molecule identifier (molregno), indicating a potential relationship between formulations and molecules.',
      'LLMComment': 'Formulations can be linked to molecules through the molregno.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_hierarchy',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Both tables reference the same molecule identifier (molregno), indicating a potential relationship between molecule hierarchy and molecules.',
      'LLMComment': 'Molecule hierarchy can be linked to molecules through the molregno.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_frac_classification',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Both tables reference the same molecule identifier (molregno), indicating a potential relationship between molecule fraction classification and molecules.',
      'LLMComment': 'Molecule fraction classification can be linked to molecules through the molregno.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_hrac_classification',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Both tables reference the same molecule identifier (molregno), indicating a potential relationship between molecule HRAC classification and molecules.',
      'LLMComment': 'Molecule HRAC classification can be linked to molecules through the molregno.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'public',
      'fromTable': 'molecule_irac_classification',
      'fromColumns': [
        'molregno'
      ],
      'toSchema': 'public',
      'toTable': 'molecule_dictionary',
      'toColumns': [
        'molregno'
      ],
      'comment': 'Both tables reference the same molecule identifier (molregno), indicating a potential relationship between molecule IRAC classification and molecules.',
      'LLMComment': 'Molecule IRAC classification can be linked to molecules through the molregno.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': false
    }
  ],
  'comment': 'Chembl is a comprehensive database of bioactive drug-like small molecules, providing detailed information on their biological activities, chemical properties, and pharmacological data. It serves as a valuable resource for researchers in drug discovery and development.',
  'LLMComment': 'The Chembl database consists of a single schema named \'public\' and contains 79 tables that store diverse data related to bioactive compounds, including their chemical structures, biological activities, and target interactions. With 102 defined relations, the database facilitates complex queries and data retrieval, enabling users to explore relationships between compounds, targets, and their associated activities effectively.'
};
