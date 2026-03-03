/* eslint-disable max-len */
import {DBConnectionMeta} from '../db-index-tools';

export const biologicsIndex: DBConnectionMeta = {
  'name': 'biologics',
  'schemas': [
    {
      'name': 'biologics',
      'tables': [
        {
          'name': 'adc',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'adc',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 600,
              'isUnique': true
            },
            {
              'name': 'identifier',
              'table': 'adc',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'semanticType': 'DG_BIOLOGICS_ADC_ID',
              'comment': 'Unique identifier for the ADC, following the DG_BIOLOGICS_ADC_ID semantic.',
              'LLMComment': 'A unique text identifier for the ADC, adhering to the DG_BIOLOGICS_ADC_ID format.'
            },
            {
              'name': 'name',
              'table': 'adc',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'comment': 'The name of the ADC, which is unique within the table.',
              'LLMComment': 'A unique text name assigned to the ADC.'
            },
            {
              'name': 'antibody_id',
              'table': 'adc',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 200,
              'isUnique': false,
              'comment': 'Foreign key referencing the unique ID of the associated antibody.',
              'LLMComment': 'An integer representing the ID of the antibody linked to this ADC.'
            },
            {
              'name': 'linker_id',
              'table': 'adc',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 10,
              'isUnique': false,
              'comment': 'Foreign key referencing the unique ID of the linker used in the ADC.',
              'LLMComment': 'An integer representing the ID of the linker component in the ADC.'
            },
            {
              'name': 'drug_id',
              'table': 'adc',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 200,
              'isUnique': false,
              'comment': 'Foreign key referencing the unique ID of the drug associated with the ADC.',
              'LLMComment': 'An integer representing the ID of the drug that the ADC delivers.'
            },
            {
              'name': 'glyph',
              'table': 'adc',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': false,
              'semanticType': 'rawPng',
              'comment': 'Raw PNG representation of the ADC\'s graphical representation.',
              'LLMComment': 'A text field containing the raw PNG data for the ADC\'s glyph.'
            }
          ],
          'rowCount': 600,
          'comment': 'Table containing information about antibody-drug conjugates (ADCs), including their identifiers, names, and associated components like antibodies, linkers, and drugs.',
          'LLMComment': 'The biologics.adc table is a crucial component in the domain of biopharmaceuticals, specifically focusing on antibody-drug conjugates (ADCs). This table stores detailed records of ADCs, which are targeted cancer therapies that combine an antibody with a cytotoxic drug. Each entry includes unique identifiers, the names of the ADCs, and foreign keys linking to other essential components such as antibodies, linkers, and drugs. This structure allows for comprehensive tracking and analysis of ADCs in relation to their performance in various assays, their chemical properties, and their interactions with target organisms, thereby supporting research and development in targeted therapies.'
        },
        {
          'name': 'assay_results',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'assay_results',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 3542,
              'isUnique': true
            },
            {
              'name': 'assay_id',
              'table': 'assay_results',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 9,
              'isUnique': false,
              'comment': 'Identifier for the specific assay conducted, linking to assay metadata.',
              'LLMComment': 'Unique identifier for the assay, which provides context for the results.'
            },
            {
              'name': 'target_organism_id',
              'table': 'assay_results',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 3,
              'isUnique': false,
              'comment': 'Identifier for the organism that is the target of the assay, linking to organism details.',
              'LLMComment': 'References the target organism for the assay, essential for biological context.'
            },
            {
              'name': 'result_value',
              'table': 'assay_results',
              'schema': 'biologics',
              'type': 'real',
              'min': 0.009999999776482582,
              'max': 4978.5,
              'isUnique': false,
              'comment': 'The quantitative result of the assay, indicating the measured effect or concentration.',
              'LLMComment': 'Numerical value representing the outcome of the assay, crucial for analysis.'
            },
            {
              'name': 'units',
              'table': 'assay_results',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'h',
                'nM',
                'µg·h/mL',
                'ratio',
                'min',
                'µg/mL',
                'RFU'
              ],
              'isUnique': false,
              'comment': 'The measurement units for the result_value, indicating the scale or type of measurement.',
              'LLMComment': 'Specifies the units of measurement for the result_value, important for interpretation.'
            },
            {
              'name': 'measured_at',
              'table': 'assay_results',
              'schema': 'biologics',
              'type': 'timestamp with time zone',
              'isUnique': false,
              'comment': 'Timestamp indicating when the assay result was recorded, including time zone information.',
              'LLMComment': 'Date and time when the measurement was taken, important for temporal analysis.'
            },
            {
              'name': 'adc_id',
              'table': 'assay_results',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 600,
              'isUnique': false,
              'comment': 'Identifier for the antibody-drug conjugate associated with the assay result.',
              'LLMComment': 'Links the result to a specific antibody-drug conjugate, relevant for drug development.'
            },
            {
              'name': 'peptide_id',
              'table': 'assay_results',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 91,
              'isUnique': false,
              'comment': 'Identifier for the peptide used in the assay, linking to peptide details.',
              'LLMComment': 'References the specific peptide involved in the assay, important for understanding the assay context.'
            }
          ],
          'rowCount': 3542,
          'comment': 'This table stores the results of various biological assays, including the measured values, units, and associated metadata.',
          'LLMComment': 'The biologics.assay_results table is a crucial component in the domain of biological research and drug development. It contains 3542 records of assay results, where each entry represents the outcome of a specific assay linked to a target organism. The table includes key columns such as \'result_value\' (the quantitative outcome of the assay), \'units\' (the measurement units for the result), and timestamps indicating when the measurements were taken. Additionally, it references other tables like \'adc\' and \'peptides\', which provide context about the assays and the biological entities involved. This structure allows researchers to analyze the efficacy and characteristics of biologics in relation to various target organisms.'
        },
        {
          'name': 'assay_types',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'assay_types',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 9,
              'isUnique': true
            },
            {
              'name': 'name',
              'table': 'assay_types',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'AUC',
                'Average Drug-to-Antibody Ratio (DAR)',
                'Binding affinity (KD / KDapp / EC50 / IC50)',
                'Tmax',
                'Internalization half-life (t½)',
                'Cell binding constant (KD / EC50)',
                'Caspase activity',
                'Cmax',
                'IC50'
              ],
              'isUnique': true,
              'comment': 'Type of assay used in biologics evaluation, indicating the specific measurement or characteristic being assessed.',
              'LLMComment': 'Describes the assay type, such as AUC or Binding affinity, which are critical for understanding drug behavior.'
            },
            {
              'name': 'description',
              'table': 'assay_types',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'Time at which Cmax is observed',
                'Flow cytometry binding curve-derived apparent affinity',
                'Affinity from SPR, BLI, or ELISA dose-response experiments',
                'Cytotoxic potency (half-maximal inhibitory concentration)',
                'Apoptosis induction via caspase activation readout',
                'Average number of drug molecules per antibody (MS / HIC / UV/Vis)',
                'Maximum observed concentration in PK profile',
                'Area under the concentration–time curve (systemic exposure)',
                'Time to reach 50% internalization of the construct'
              ],
              'isUnique': true,
              'comment': 'Detailed explanation of what the assay measures, providing context for its application in biologics research.',
              'LLMComment': 'Offers insights into the significance of the assay type, such as the meaning of Tmax or Internalization half-life.'
            }
          ],
          'rowCount': 9,
          'comment': 'This table stores different types of assays used in the evaluation of biologics, including their unique identifiers, names, and descriptions.',
          'LLMComment': 'The \'biologics.assay_types\' table is crucial for categorizing and defining various assay methodologies employed in the development and testing of biologic drugs. Each entry includes an \'id\' for unique identification, a \'name\' that specifies the type of assay, and a \'description\' providing detailed information about the assay\'s purpose and application. This table interacts with other tables in the schema, such as \'assay_results\' and \'drugs\', to facilitate comprehensive data analysis and reporting in the biologics domain.'
        },
        {
          'name': 'drugs',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'drugs',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 200,
              'isUnique': true
            },
            {
              'name': 'identifier',
              'table': 'drugs',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'semanticType': 'DG_BIOLOGICS_DRUG_ID',
              'comment': 'Unique identifier for the drug in the biologics database, following the DG_BIOLOGICS_DRUG_ID format.',
              'LLMComment': 'A unique textual identifier assigned to each drug, used for referencing in the biologics domain.'
            },
            {
              'name': 'name',
              'table': 'drugs',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'comment': 'The official name of the drug, which is unique within the database.',
              'LLMComment': 'The distinct name by which the drug is known, ensuring no duplicates exist.'
            },
            {
              'name': 'smiles',
              'table': 'drugs',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'semanticType': 'Molecule',
              'comment': 'A unique representation of the drug\'s molecular structure using the SMILES notation.',
              'LLMComment': 'A string that encodes the molecular structure of the drug in SMILES format, allowing for easy interpretation by cheminformatics tools.'
            }
          ],
          'rowCount': 200,
          'comment': 'This table contains information about biologic drugs, including their unique identifiers, names, and chemical structure representations in SMILES format.',
          'LLMComment': 'The \'biologics.drugs\' table serves as a central repository for biologic drug data within the biologics schema. It includes essential attributes such as \'id\' for unique identification, \'identifier\' for external references, \'name\' for the drug\'s common name, and \'smiles\' for the chemical structure representation. This table is crucial for linking to other related tables like \'adc\' (antibody-drug conjugates), \'assay_results\', and \'target_organisms\', facilitating comprehensive analysis and research in the field of biologics.'
        },
        {
          'name': 'expression_batches',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'expression_batches',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 600,
              'isUnique': true
            },
            {
              'name': 'identifier',
              'table': 'expression_batches',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'semanticType': 'DG_BIOLOGICS_EXPRESSION_ID',
              'comment': 'Unique identifier for the expression batch, following the DG_BIOLOGICS_EXPRESSION_ID format.',
              'LLMComment': 'A unique text identifier for the expression batch, adhering to the DG_BIOLOGICS_EXPRESSION_ID semantic.'
            },
            {
              'name': 'sequence_id',
              'table': 'expression_batches',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 200,
              'isUnique': false,
              'comment': 'Identifier for the sequence associated with the expression batch, ranging from 1 to 200.',
              'LLMComment': 'An integer representing the sequence ID linked to this expression batch, within a defined range.'
            },
            {
              'name': 'expression_system',
              'table': 'expression_batches',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'CHO',
                'E.coli',
                'HEK293'
              ],
              'isUnique': false,
              'comment': 'Type of expression system used for the batch, such as CHO, E.coli, or HEK293.',
              'LLMComment': 'Categorical text indicating the expression system utilized, with options including CHO, E.coli, and HEK293.'
            },
            {
              'name': 'yield_mg',
              'table': 'expression_batches',
              'schema': 'biologics',
              'type': 'real',
              'min': 0.5099999904632568,
              'max': 149.50999450683594,
              'isUnique': false,
              'comment': 'Amount of product yielded from the expression batch, measured in milligrams.',
              'LLMComment': 'A real number representing the yield of the expression batch in milligrams, within a specified range.'
            },
            {
              'name': 'name',
              'table': 'expression_batches',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'comment': 'Unique name assigned to the expression batch.',
              'LLMComment': 'A unique text name for the expression batch, ensuring distinct identification.'
            },
            {
              'name': 'notes',
              'table': 'expression_batches',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'Auto-generated expression batch'
              ],
              'isUnique': false,
              'comment': 'Additional notes regarding the expression batch, typically auto-generated.',
              'LLMComment': 'Text field for supplementary information about the expression batch, often auto-generated.'
            },
            {
              'name': 'created_at',
              'table': 'expression_batches',
              'schema': 'biologics',
              'type': 'timestamp with time zone',
              'isUnique': false,
              'comment': 'Timestamp indicating when the expression batch was created, including time zone information.',
              'LLMComment': 'A timestamp with time zone that records the creation date and time of the expression batch.'
            }
          ],
          'rowCount': 600,
          'comment': 'This table stores information about batches of biologics expression, including their identifiers, yield, and associated metadata.',
          'LLMComment': 'The biologics.expression_batches table is crucial for tracking the production of biologics through various expression batches. Each entry represents a unique batch of biologics, identified by an ID and a specific identifier, and includes details such as the sequence used, the expression system employed, the yield in milligrams, and any relevant notes. The created_at timestamp helps in managing the lifecycle of these batches. This table interacts with other tables in the schema, such as sequences and purification_batches, to provide a comprehensive view of the biologics production process.'
        },
        {
          'name': 'linkers',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'linkers',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 10,
              'isUnique': true
            },
            {
              'name': 'identifier',
              'table': 'linkers',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'GROKLINKER-000001',
                'GROKLINKER-000002',
                'GROKLINKER-000003',
                'GROKLINKER-000004',
                'GROKLINKER-000005',
                'GROKLINKER-000006',
                'GROKLINKER-000007',
                'GROKLINKER-000008',
                'GROKLINKER-000009',
                'GROKLINKER-000010'
              ],
              'isUnique': true,
              'semanticType': 'DG_BIOLOGICS_LINKER_ID',
              'comment': 'Unique identifier for the linker, following a specific naming convention.',
              'LLMComment': 'A unique text identifier for each linker, formatted as DG_BIOLOGICS_LINKER_ID, such as GROKLINKER-000001.'
            },
            {
              'name': 'linker_type',
              'table': 'linkers',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'PROTEIN',
                'SMALL'
              ],
              'isUnique': false,
              'comment': 'Type of linker, indicating its classification as either a protein or small molecule.',
              'LLMComment': 'Categorizes the linker as either a PROTEIN or SMALL molecule.'
            },
            {
              'name': 'linker_molecule_smiles',
              'table': 'linkers',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'null',
                'CCOCCOCCOCCOCCOCCOCCOCCO',
                'CCNCCNCCNCCNCCNCCNCCN',
                'CNC(=O)CCNC(=O)CCNC(=O)CCNC(=O)C',
                'COCOCOCOCOCOCOCO',
                'NCC(=O)NCC(=O)NCC(=O)NCC=O'
              ],
              'isUnique': false,
              'semanticType': 'Molecule',
              'comment': 'SMILES representation of the linker molecule, used for chemical structure encoding.',
              'LLMComment': 'A text representation of the chemical structure of the linker in SMILES format, such as CCOCCOCCOCCO.'
            },
            {
              'name': 'linker_sequence',
              'table': 'linkers',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'null',
                'GPGPG',
                'GGSGGS',
                'GSGSG',
                'GGGGS',
                'GGGSGGGS'
              ],
              'isUnique': false,
              'comment': 'A sequence representation of the linker, often used in biological contexts.',
              'LLMComment': 'A text sequence that represents the linker, which may include repeating units like GPGPG or GGSGGS.'
            }
          ],
          'rowCount': 10,
          'comment': 'This table stores information about linkers used in biologics, including their identifiers, types, and molecular structures.',
          'LLMComment': 'The \'biologics.linkers\' table is a crucial component in the domain of biologics research, specifically focusing on the linkers that connect various components in antibody-drug conjugates (ADCs). Each entry includes an \'id\' for unique identification, \'identifier\' for referencing the linker, \'linker_type\' to categorize the type of linker, \'linker_molecule_smiles\' which provides the SMILES representation of the linker molecule for computational analysis, and \'linker_sequence\' that details the specific sequence of the linker. This table interacts with other tables such as \'adc\' for linking to specific drug conjugates, \'assay_results\' for evaluating linker performance, and \'peptides\' for understanding how linkers interact with peptide sequences.'
        },
        {
          'name': 'peptides',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'peptides',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 91,
              'isUnique': true
            },
            {
              'name': 'identifier',
              'table': 'peptides',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'GROKPEP-000008',
                'GROKPEP-000041',
                'GROKPEP-000036',
                'GROKPEP-000079',
                'GROKPEP-000031',
                'GROKPEP-000033',
                'GROKPEP-000021',
                'GROKPEP-000071',
                'GROKPEP-000003',
                'GROKPEP-000034',
                'GROKPEP-000016',
                'GROKPEP-000006',
                'GROKPEP-000004',
                'GROKPEP-000020',
                'GROKPEP-000039',
                'GROKPEP-000010',
                'GROKPEP-000053',
                'GROKPEP-000048',
                'GROKPEP-000042',
                'GROKPEP-000035',
                'GROKPEP-000024',
                'GROKPEP-000089',
                'GROKPEP-000076',
                'GROKPEP-000017',
                'GROKPEP-000032',
                'GROKPEP-000026',
                'GROKPEP-000018',
                'GROKPEP-000007',
                'GROKPEP-000067',
                'GROKPEP-000052',
                'GROKPEP-000046',
                'GROKPEP-000014',
                'GROKPEP-000029',
                'GROKPEP-000019',
                'GROKPEP-000005',
                'GROKPEP-000038',
                'GROKPEP-000087',
                'GROKPEP-000086',
                'GROKPEP-000065',
                'GROKPEP-000059',
                'GROKPEP-000040',
                'GROKPEP-000074',
                'GROKPEP-000001',
                'GROKPEP-000002',
                'GROKPEP-000043',
                'GROKPEP-000080',
                'GROKPEP-000060',
                'GROKPEP-000064',
                'GROKPEP-000091',
                'GROKPEP-000013',
                'GROKPEP-000068',
                'GROKPEP-000057',
                'GROKPEP-000044',
                'GROKPEP-000055',
                'GROKPEP-000061',
                'GROKPEP-000023',
                'GROKPEP-000062',
                'GROKPEP-000077',
                'GROKPEP-000054',
                'GROKPEP-000075',
                'GROKPEP-000069',
                'GROKPEP-000049',
                'GROKPEP-000088',
                'GROKPEP-000073',
                'GROKPEP-000056',
                'GROKPEP-000066',
                'GROKPEP-000058',
                'GROKPEP-000063',
                'GROKPEP-000011',
                'GROKPEP-000083',
                'GROKPEP-000082',
                'GROKPEP-000012',
                'GROKPEP-000078',
                'GROKPEP-000051',
                'GROKPEP-000025',
                'GROKPEP-000085',
                'GROKPEP-000045',
                'GROKPEP-000037',
                'GROKPEP-000028',
                'GROKPEP-000070',
                'GROKPEP-000050',
                'GROKPEP-000090',
                'GROKPEP-000084',
                'GROKPEP-000081',
                'GROKPEP-000022',
                'GROKPEP-000030',
                'GROKPEP-000015',
                'GROKPEP-000027',
                'GROKPEP-000009',
                'GROKPEP-000072',
                'GROKPEP-000047'
              ],
              'isUnique': true,
              'semanticType': 'DG_BIOLOGICS_PEPTIDE_ID',
              'comment': 'Unique identifier for each peptide, following the DG_BIOLOGICS_PEPTIDE_ID semantic format.',
              'LLMComment': 'A unique text identifier for the peptide, formatted as DG_BIOLOGICS_PEPTIDE_ID, used for tracking and referencing specific peptides.'
            },
            {
              'name': 'name',
              'table': 'peptides',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'Peptide_82',
                'Peptide_61',
                'Peptide_57',
                'Peptide_34',
                'Peptide_74',
                'Peptide_25',
                'Peptide_72',
                'Peptide_11',
                'Peptide_73',
                'Peptide_29',
                'Peptide_12',
                'Peptide_66',
                'Peptide_3',
                'Peptide_49',
                'Peptide_20',
                'Peptide_70',
                'Peptide_53',
                'Peptide_65',
                'Peptide_10',
                'Peptide_30',
                'Peptide_27',
                'Peptide_16',
                'Peptide_54',
                'Peptide_71',
                'Peptide_26',
                'Peptide_7',
                'Peptide_35',
                'Peptide_15',
                'Peptide_83',
                'Peptide_38',
                'Peptide_32',
                'Peptide_67',
                'Peptide_21',
                'Peptide_55',
                'Peptide_51',
                'Peptide_18',
                'Peptide_75',
                'Peptide_37',
                'Peptide_31',
                'Peptide_78',
                'Peptide_69',
                'Peptide_24',
                'Peptide_87',
                'Peptide_91',
                'Peptide_90',
                'Peptide_52',
                'Peptide_62',
                'Peptide_56',
                'Peptide_44',
                'Peptide_13',
                'Peptide_4',
                'Peptide_88',
                'Peptide_48',
                'Peptide_5',
                'Peptide_77',
                'Peptide_63',
                'Peptide_28',
                'Peptide_59',
                'Peptide_45',
                'Peptide_60',
                'Peptide_8',
                'Peptide_39',
                'Peptide_6',
                'Peptide_46',
                'Peptide_23',
                'Peptide_81',
                'Peptide_33',
                'Peptide_85',
                'Peptide_19',
                'Peptide_80',
                'Peptide_36',
                'Peptide_68',
                'Peptide_58',
                'Peptide_41',
                'Peptide_22',
                'Peptide_1',
                'Peptide_42',
                'Peptide_40',
                'Peptide_9',
                'Peptide_47',
                'Peptide_76',
                'Peptide_84',
                'Peptide_14',
                'Peptide_86',
                'Peptide_64',
                'Peptide_50',
                'Peptide_43',
                'Peptide_79',
                'Peptide_2',
                'Peptide_89',
                'Peptide_17'
              ],
              'isUnique': true,
              'comment': 'Unique name assigned to each peptide, used for identification purposes.',
              'LLMComment': 'A unique text name for the peptide, which helps in identifying and categorizing the peptide within the database.'
            },
            {
              'name': 'helm',
              'table': 'peptides',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'semanticType': 'Macromolecule',
              'comment': 'HELM notation representing the macromolecular structure of the peptide.',
              'LLMComment': 'A unique text representation of the peptide\'s structure in HEuristics for Lifelike Macromolecules (HELM) format.'
            },
            {
              'name': 'molecule_structure',
              'table': 'peptides',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'null'
              ],
              'isUnique': false,
              'comment': 'Describes the structural composition of the peptide, specific categories may apply.',
              'LLMComment': 'Textual representation of the molecular structure of the peptide, detailing its chemical composition and arrangement.'
            }
          ],
          'rowCount': 91,
          'comment': 'This table stores information about peptides used in biologics research, including their identifiers, names, and structural representations.',
          'LLMComment': 'The \'biologics.peptides\' table is a key component in the biologics research domain, specifically focusing on peptides, which are short chains of amino acids that play crucial roles in various biological functions. This table contains 91 entries, each representing a unique peptide characterized by an integer ID, a textual identifier, a name, a HELM notation (a standard for representing complex biomolecules), and a molecular structure description. The data in this table is essential for researchers working with biologics, as it provides foundational information that can be linked to other tables such as \'adc\' (antibody-drug conjugates), \'assay_results\', and \'target_organisms\', facilitating comprehensive analysis and development of peptide-based therapeutics.'
        },
        {
          'name': 'purification_batches',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'purification_batches',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 600,
              'isUnique': true
            },
            {
              'name': 'identifier',
              'table': 'purification_batches',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'semanticType': 'DG_BIOLOGICS_PURIFICATION_ID',
              'comment': 'Unique identifier for the purification batch, following the DG_BIOLOGICS_PURIFICATION_ID semantic.',
              'LLMComment': 'A unique text identifier for each purification batch, adhering to the DG_BIOLOGICS standard.'
            },
            {
              'name': 'sequence_id',
              'table': 'purification_batches',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 200,
              'isUnique': false,
              'comment': 'Sequential identifier for tracking the order of purification batches, ranging from 1 to 200.',
              'LLMComment': 'An integer representing the order of the purification batch in a sequence, useful for tracking purposes.'
            },
            {
              'name': 'name',
              'table': 'purification_batches',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': false,
              'comment': 'Descriptive name of the purification batch, providing context or identification.',
              'LLMComment': 'A text field that gives a descriptive name to the purification batch for easier identification.'
            },
            {
              'name': 'notes',
              'table': 'purification_batches',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'Auto-generated purification batch'
              ],
              'isUnique': false,
              'comment': 'Additional information or comments related to the purification batch, typically auto-generated.',
              'LLMComment': 'A text field for storing supplementary notes about the purification batch, often auto-generated.'
            },
            {
              'name': 'created_at',
              'table': 'purification_batches',
              'schema': 'biologics',
              'type': 'timestamp with time zone',
              'isUnique': false,
              'comment': 'Timestamp indicating when the purification batch record was created, including time zone information.',
              'LLMComment': 'A timestamp that records the creation date and time of the purification batch entry, with time zone details.'
            }
          ],
          'rowCount': 600,
          'comment': 'This table stores information about purification batches of biologics, including unique identifiers, associated sequences, and metadata such as creation timestamps and notes.',
          'LLMComment': 'The \'biologics.purification_batches\' table is crucial in the context of biologics development, specifically for tracking the purification processes of various biologic products. Each entry represents a distinct purification batch, identified by a unique ID and associated with a specific sequence (via \'sequence_id\'). The \'name\' and \'notes\' fields provide additional context about each batch, while the \'created_at\' timestamp helps in tracking the timeline of purification activities. This table interacts with other tables in the schema, such as \'adc\' and \'drugs\', to provide a comprehensive view of the biologics lifecycle from purification to final product.'
        },
        {
          'name': 'sequences',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'sequences',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 200,
              'isUnique': true
            },
            {
              'name': 'identifier',
              'table': 'sequences',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'semanticType': 'DG_BIOLOGICS_SEQUENCE_ID',
              'comment': 'Unique identifier for the biologics sequence, following the DG_BIOLOGICS_SEQUENCE_ID semantic.',
              'LLMComment': 'A unique text identifier for each biologics sequence, adhering to the DG_BIOLOGICS_SEQUENCE_ID standard.'
            },
            {
              'name': 'name',
              'table': 'sequences',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'comment': 'Unique name assigned to the biologics sequence.',
              'LLMComment': 'A distinct name for the biologics sequence, ensuring no duplicates.'
            },
            {
              'name': 'sequence_type',
              'table': 'sequences',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'PROTEIN'
              ],
              'isUnique': false,
              'comment': 'Type of sequence, categorized as PROTEIN.',
              'LLMComment': 'Indicates the type of sequence, specifically categorized under PROTEIN.'
            },
            {
              'name': 'sequence',
              'table': 'sequences',
              'schema': 'biologics',
              'type': 'text',
              'isUnique': true,
              'semanticType': 'Macromolecule',
              'comment': 'The actual macromolecule sequence, unique to each entry.',
              'LLMComment': 'The unique text representation of the macromolecule sequence for the biologics.'
            }
          ],
          'rowCount': 200,
          'comment': 'This table stores information about biological sequences related to various biologics, including their identifiers, names, types, and the actual sequence data.',
          'LLMComment': 'The \'biologics.sequences\' table is a crucial component of the biologics database, containing detailed records of biological sequences used in the development and analysis of biologics. Each entry includes a unique identifier, the name of the sequence, the type of sequence (such as DNA, RNA, or protein), and the actual sequence data. This table serves as a foundational reference for other related tables in the schema, such as \'adc\' for antibody-drug conjugates, \'assay_results\' for experimental outcomes, and \'drugs\' for pharmacological information, thereby facilitating comprehensive research and development in the field of biologics.'
        },
        {
          'name': 'target_organisms',
          'schema': 'biologics',
          'columns': [
            {
              'name': 'id',
              'table': 'target_organisms',
              'schema': 'biologics',
              'type': 'integer',
              'min': 1,
              'max': 3,
              'isUnique': true
            },
            {
              'name': 'identifier',
              'table': 'target_organisms',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'GROKORG-000001',
                'GROKORG-000002',
                'GROKORG-000003'
              ],
              'isUnique': true,
              'semanticType': 'DG_BIOLOGICS_ORGANISM_ID',
              'comment': 'Unique identifier for the organism in the biologics database, following the DG_BIOLOGICS_ORGANISM_ID format.',
              'LLMComment': 'A unique text identifier for each target organism, formatted as DG_BIOLOGICS_ORGANISM_ID.'
            },
            {
              'name': 'name',
              'table': 'target_organisms',
              'schema': 'biologics',
              'type': 'text',
              'categoryValues': [
                'CHO',
                'S. cerevisiae',
                'E. coli'
              ],
              'isUnique': true,
              'comment': 'Common name of the target organism, representing the species used in biologics research.',
              'LLMComment': 'The common name of the organism, indicating the species utilized in biologics applications.'
            }
          ],
          'rowCount': 3,
          'comment': 'This table stores information about target organisms relevant to biologics research, including their unique identifiers and names.',
          'LLMComment': 'The \'biologics.target_organisms\' table is crucial for managing data related to the organisms that are the focus of biologics studies. It contains three columns: \'id\' serves as a unique integer identifier for each organism, \'identifier\' provides a textual reference for the organism, and \'name\' gives the common name of the organism. This table is interconnected with other tables in the schema, such as \'assay_results\' and \'drugs\', which likely reference these target organisms in the context of drug development and biological assays.'
        }
      ],
      'comment': 'The biologics schema is designed to manage and store data related to biologic drugs, including their components, assay results, and related processes.',
      'LLMComment': 'This schema encompasses various tables that track biologic drug development, including details on antibody-drug conjugates (adc), assay results, types of assays, drug information, expression batches, linkers, peptides, purification processes, sequences, and target organisms.'
    }
  ],
  'relations': [
    {
      'fromSchema': 'biologics',
      'fromTable': 'adc',
      'fromColumns': [
        'linker_id'
      ],
      'toSchema': 'biologics',
      'toTable': 'linkers',
      'toColumns': [
        'id'
      ],
      'comment': 'Each ADC (Antibody-Drug Conjugate) is associated with one linker; one linker can be used in multiple ADCs.',
      'LLMComment': 'Each adc record in the biologics.adc table is linked to one linker record in the biologics.linkers table through the linker_id foreign key; conversely, one linker can be associated with multiple adc records.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'adc',
      'fromColumns': [
        'drug_id'
      ],
      'toSchema': 'biologics',
      'toTable': 'drugs',
      'toColumns': [
        'id'
      ],
      'comment': 'Each ADC (Antibody-Drug Conjugate) is associated with one drug; one drug can be associated with multiple ADCs.',
      'LLMComment': 'Each adc entry in the biologics.adc table is linked to one drug in the biologics.drugs table through the drug_id foreign key; conversely, one drug can be linked to multiple adc entries.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'adc',
      'fromColumns': [
        'antibody_id'
      ],
      'toSchema': 'biologics',
      'toTable': 'sequences',
      'toColumns': [
        'id'
      ],
      'comment': 'Each ADC (antibody-drug conjugate) is associated with one sequence, identified by antibody_id.',
      'LLMComment': 'Each ADC (antibody-drug conjugate) is associated with one sequence, identified by antibody_id, and one sequence can be linked to many ADCs.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'assay_results',
      'fromColumns': [
        'adc_id'
      ],
      'toSchema': 'biologics',
      'toTable': 'adc',
      'toColumns': [
        'id'
      ],
      'comment': 'Each assay_result is associated with one ADC; one ADC can have many assay_results.',
      'LLMComment': 'Each assay_result belongs to one ADC; one ADC can have many assay_results, indicating a relationship where multiple results can be linked to a single ADC entry.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'assay_results',
      'fromColumns': [
        'peptide_id'
      ],
      'toSchema': 'biologics',
      'toTable': 'peptides',
      'toColumns': [
        'id'
      ],
      'comment': 'Each assay_result is associated with one peptide; one peptide can have many assay_results.',
      'LLMComment': 'Each assay_result belongs to one peptide; one peptide can have many assay_results.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'assay_results',
      'fromColumns': [
        'target_organism_id'
      ],
      'toSchema': 'biologics',
      'toTable': 'target_organisms',
      'toColumns': [
        'id'
      ],
      'comment': 'Each assay_result is associated with one target_organism; one target_organism can have many assay_results.',
      'LLMComment': 'Each assay_result belongs to one target organism; one target organism can have many assay results associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'assay_results',
      'fromColumns': [
        'assay_id'
      ],
      'toSchema': 'biologics',
      'toTable': 'assay_types',
      'toColumns': [
        'id'
      ],
      'comment': 'Each assay_result is linked to one assay_type; one assay_type can have many assay_results.',
      'LLMComment': 'Each assay_result belongs to one assay_type; one assay_type can have many assay_results.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'expression_batches',
      'fromColumns': [
        'sequence_id'
      ],
      'toSchema': 'biologics',
      'toTable': 'sequences',
      'toColumns': [
        'id'
      ],
      'comment': 'Each expression_batch is associated with one sequence; one sequence can be linked to many expression_batches.',
      'LLMComment': 'Each expression_batch belongs to one sequence; one sequence can have many expression_batches associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'purification_batches',
      'fromColumns': [
        'sequence_id'
      ],
      'toSchema': 'biologics',
      'toTable': 'sequences',
      'toColumns': [
        'id'
      ],
      'comment': 'Each purification_batch is associated with one sequence; one sequence can be associated with many purification_batches.',
      'LLMComment': 'Each purification_batch belongs to one sequence; one sequence can have many purification_batches associated with it.',
      'cardinality': 'one-to-many',
      'IsPrimaryPath': true
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'adc',
      'fromColumns': [
        'name'
      ],
      'toSchema': 'biologics',
      'toTable': 'drugs',
      'toColumns': [
        'name'
      ],
      'comment': 'The name of ADCs may correspond to the name of drugs, indicating a potential relationship.',
      'LLMComment': 'ADC names might be similar or related to drug names, suggesting a logical join for queries involving drug information.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'adc',
      'fromColumns': [
        'identifier'
      ],
      'toSchema': 'biologics',
      'toTable': 'drugs',
      'toColumns': [
        'identifier'
      ],
      'comment': 'The identifier of ADCs may correspond to the identifier of drugs, indicating a potential relationship.',
      'LLMComment': 'Identifiers for ADCs and drugs may follow a similar naming convention, suggesting a logical join for queries involving drug information.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': false
    },
    {
      'fromSchema': 'biologics',
      'fromTable': 'assay_results',
      'fromColumns': [
        'result_value'
      ],
      'toSchema': 'biologics',
      'toTable': 'assay_types',
      'toColumns': [
        'name'
      ],
      'comment': 'Assay results may be categorized by assay types, indicating a potential relationship based on result values and assay types.',
      'LLMComment': 'Result values in assay results could be linked to specific assay types, suggesting a logical join for queries involving assay performance.',
      'cardinality': 'many-to-one',
      'IsPrimaryPath': false
    }
  ],
  'comment': 'The biologics database is designed to store and manage data related to biologic products, including their development, manufacturing, and regulatory information. It serves as a central repository for researchers and manufacturers in the biotechnology and pharmaceutical industries.',
  'LLMComment': 'The biologics database consists of a single schema named \'biologics\' which contains 10 tables interconnected by 9 relationships. This structure allows for efficient data retrieval and management of biologic product information, including details on product types, manufacturing processes, clinical trials, and regulatory compliance. The relationships between tables facilitate complex queries and data analysis, making it a valuable resource for stakeholders in the biologics field.'
};
