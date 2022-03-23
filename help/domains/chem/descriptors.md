<!-- TITLE: Molecular descriptors -->
<!-- SUBTITLE: -->

# Molecular descriptors

The molecular descriptor is the final result of a logic and mathematical procedure which transforms chemical information
encoded within a symbolic representation of a molecule into a useful number or the result of some standardized
experiment.

![descriptors](descriptors.gif)

## Available descriptors modules

### EState.EState_VSA

Hybrid EState-VSA descriptors (like the MOE VSA descriptors)

| Descriptor   | Description                                   |
|--------------|-----------------------------------------------|
| EState_VSA1  | EState VSA Descriptor 1 (-inf < x < -0.39)    |
| EState_VSA10 | EState VSA Descriptor 10 ( 9.17 <= x < 15.00) |
| EState_VSA11 | EState VSA Descriptor 11 ( 15.00 <= x < inf)  |
| EState_VSA2  | EState VSA Descriptor 2 ( -0.39 <= x < 0.29)  |
| EState_VSA3  | EState VSA Descriptor 3 ( 0.29 <= x < 0.72)   |
| EState_VSA4  | EState VSA Descriptor 4 ( 0.72 <= x < 1.17)   |
| EState_VSA5  | EState VSA Descriptor 5 ( 1.17 <= x < 1.54)   |
| EState_VSA6  | EState VSA Descriptor 6 ( 1.54 <= x < 1.81)   |
| EState_VSA7  | EState VSA Descriptor 7 ( 1.81 <= x < 2.05)   |
| EState_VSA8  | EState VSA Descriptor 8 ( 2.05 <= x < 4.69)   |
| EState_VSA9  | EState VSA Descriptor 9 ( 4.69 <= x < 9.17)   |
| VSA_EState1  | VSA EState Descriptor 1 (-inf < x < 4.78)     |
| VSA_EState10 | VSA EState Descriptor 10 ( 11.00 <= x < inf)  |
| VSA_EState2  | VSA EState Descriptor 2 ( 4.78 <= x < 5.00)   |
| VSA_EState3  | VSA EState Descriptor 3 ( 5.00 <= x < 5.41)   |
| VSA_EState4  | VSA EState Descriptor 4 ( 5.41 <= x < 5.74)   |
| VSA_EState5  | VSA EState Descriptor 5 ( 5.74 <= x < 6.00)   |
| VSA_EState6  | VSA EState Descriptor 6 ( 6.00 <= x < 6.07)   |
| VSA_EState7  | VSA EState Descriptor 7 ( 6.07 <= x < 6.45)   |
| VSA_EState8  | VSA EState Descriptor 8 ( 6.45 <= x < 7.00)   |
| VSA_EState9  | VSA EState Descriptor 9 ( 7.00 <= x < 11.00)  |

### QED

QED stands for quantitative estimation of drug-likeness

| Descriptor | Description                           |
|------------|---------------------------------------|
| qed        | Weighted sum of ADS mapped properties |

### GraphDescriptors

Topological/topochemical descriptors

| Descriptor    | Description                                                                                                                                       |
|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| BalabanJ      | Balaban's J value                                                                                                                                 |
| BertzCT       | A topological index meant to quantify "complexity"                                                                                                |
| Chi0          | From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)                                                                        |
| Chi0n         | Similar to Hall Kier Chi0v, but uses nVal instead of valence. This makes a big difference after we get out of the first row                       |
| Chi0v         | From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)                                                                        |
| Chi1          | From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)                                                                       |
| Chi1n         | Similar to Hall Kier Chi1v, but uses nVal instead of valence                                                                                      |
| Chi1v         | From equations (5),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)                                                                       |
| Chi2n         | Similar to Hall Kier Chi2v, but uses nVal instead of valence. This makes a big difference after we get out of the first row                       |
| Chi2v         | From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991)                                                                       |
| Chi3n         | Similar to Hall Kier Chi3v, but uses nVal instead of valence. This makes a big difference after we get out of the first row                       |
| Chi3v         | From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991)                                                                       |
| Chi4n         | Similar to Hall Kier Chi4v, but uses nVal instead of valence. This makes a big difference after we get out of the first row                       |
| Chi4v         | From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991).                                                                      |
| HallKierAlpha | The Hall-Kier alpha value for a molecule                                                                                                          |
| Ipc           | The information content of the coefficients of the characteristic polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule |
| Kappa1        | Hall-Kier Kappa1 value                                                                                                                            |
| Kappa2        | Hall-Kier Kappa2 value                                                                                                                            |
| Kappa3        | Hall-Kier Kappa3 value                                                                                                                            |

### Descriptors

General descriptors

| Descriptor          | Description                                                                      |
|---------------------|----------------------------------------------------------------------------------|
| MinAbsPartialCharge | Minimal absolute partial charge                                                  |
| NumRadicalElectrons | The number of radical electrons the molecule has (says nothing about spin state) |
| FpDensityMorgan2    | Morgan fingerprint, radius 2                                                     |
| FpDensityMorgan3    | Morgan fingerprint, radius 3                                                     |
| FpDensityMorgan1    | Morgan fingerprint, radius 1                                                     |
| HeavyAtomMolWt      | The average molecular weight of the molecule ignoring hydrogens                  |
| MaxAbsPartialCharge | Maximum absolute partial charge                                                  |
| MinPartialCharge    | Minimal partial charge                                                           |
| ExactMolWt          | The exact molecular weight of the molecule                                       |
| MolWt               | The average molecular weight of the molecule                                     |
| NumValenceElectrons | The number of valence electrons the molecule has                                 |
| MaxPartialCharge    | Maximum partial charge                                                           |

### Crippen

Atom-based calculation of LogP and MR using Crippen's approach

| Descriptor | Description                |
|------------|----------------------------|
| MolLogP    | Wildman-Crippen LogP value |
| MolMR      | Wildman-Crippen MR value   |

### MolSurf

MOE-like approximate molecular surface area descriptors

| Descriptor  | Description                                                |
|-------------|------------------------------------------------------------|
| LabuteASA   | Labute's Approximate Surface Area (ASA from MOE)           |
| PEOE_VSA1   | MOE Charge VSA Descriptor 1 (-inf < x < -0.30)             |
| PEOE_VSA10  | MOE Charge VSA Descriptor 10 ( 0.10 <= x < 0.15)           |
| PEOE_VSA11  | MOE Charge VSA Descriptor 11 ( 0.15 <= x < 0.20)           |
| PEOE_VSA12  | MOE Charge VSA Descriptor 12 ( 0.20 <= x < 0.25)           |
| PEOE_VSA13  | MOE Charge VSA Descriptor 13 ( 0.25 <= x < 0.30)           |
| PEOE_VSA14  | MOE Charge VSA Descriptor 14 ( 0.30 <= x < inf)            |
| PEOE_VSA2   | MOE Charge VSA Descriptor 2 (-0.30 <= x < -0.25)           |
| PEOE_VSA3   | MOE Charge VSA Descriptor 3 (-0.25 <= x < -0.20)           |
| PEOE_VSA4   | MOE Charge VSA Descriptor 4 (-0.20 <= x < -0.15)           |
| PEOE_VSA5   | MOE Charge VSA Descriptor 5 (-0.15 <= x < -0.10)           |
| PEOE_VSA6   | MOE Charge VSA Descriptor 6 (-0.10 <= x < -0.05)           |
| PEOE_VSA7   | MOE Charge VSA Descriptor 7 (-0.05 <= x < 0.00)            |
| PEOE_VSA8   | MOE Charge VSA Descriptor 8 ( 0.00 <= x < 0.05)            |
| PEOE_VSA9   | MOE Charge VSA Descriptor 9 ( 0.05 <= x < 0.10)            |
| SMR_VSA1    | MOE MR VSA Descriptor 1 (-inf < x < 1.29)                  |
| SMR_VSA10   | MOE MR VSA Descriptor 10 ( 4.00 <= x < inf)                |
| SMR_VSA2    | MOE MR VSA Descriptor 2 ( 1.29 <= x < 1.82)                |
| SMR_VSA3    | MOE MR VSA Descriptor 3 ( 1.82 <= x < 2.24)                |
| SMR_VSA4    | MOE MR VSA Descriptor 4 ( 2.24 <= x < 2.45)                |
| SMR_VSA5    | MOE MR VSA Descriptor 5 ( 2.45 <= x < 2.75)                |
| SMR_VSA6    | MOE MR VSA Descriptor 6 ( 2.75 <= x < 3.05)                |
| SMR_VSA7    | MOE MR VSA Descriptor 7 ( 3.05 <= x < 3.63)                |
| SMR_VSA8    | MOE MR VSA Descriptor 8 ( 3.63 <= x < 3.80)                |
| SMR_VSA9    | MOE MR VSA Descriptor 9 ( 3.80 <= x < 4.00)                |
| SlogP_VSA1  | MOE logP VSA Descriptor 1 (-inf < x < -0.40)               |
| SlogP_VSA10 | MOE logP VSA Descriptor 10 ( 0.40 <= x < 0.50)             |
| SlogP_VSA11 | MOE logP VSA Descriptor 11 ( 0.50 <= x < 0.60)             |
| SlogP_VSA12 | MOE logP VSA Descriptor 12 ( 0.60 <= x < inf)              |
| SlogP_VSA2  | MOE logP VSA Descriptor 2 (-0.40 <= x < -0.20)             |
| SlogP_VSA3  | MOE logP VSA Descriptor 3 (-0.20 <= x < 0.00)              |
| SlogP_VSA4  | MOE logP VSA Descriptor 4 ( 0.00 <= x < 0.10)              |
| SlogP_VSA5  | MOE logP VSA Descriptor 5 ( 0.10 <= x < 0.15)              |
| SlogP_VSA6  | MOE logP VSA Descriptor 6 ( 0.15 <= x < 0.20)              |
| SlogP_VSA7  | MOE logP VSA Descriptor 7 ( 0.20 <= x < 0.25)              |
| SlogP_VSA8  | MOE logP VSA Descriptor 8 ( 0.25 <= x < 0.30)              |
| SlogP_VSA9  | MOE logP VSA Descriptor 9 ( 0.30 <= x < 0.40)              |
| TPSA        | The polar surface area of a molecule based upon fragments. |

### Descriptors3D

Descriptors derived from a molecule's 3D structure

| Descriptor          | Description                                   |
|---------------------|-----------------------------------------------|
| PMI1                | First (smallest) principal moment of inertia  |
| PMI2                | Second principal moment of inertia            |
| PMI3                | Third (largest) principal moment of inertia   |
| NPR1                | Normalized principal moments ratio 1 (=I1/I3) |
| NPR2                | Normalized principal moments ratio 2 (=I2/I3) |
| RadiusOfGyration    | Radius of gyration                            |
| InertialShapeFactor | Inertial shape factor                         |
| Eccentricity        | Molecular eccentricity                        |
| Asphericity         | Molecular asphericity                         |
| SpherocityIndex     | Molecular spherocity index                    |

### EState.EState

Basic EState descriptors

| Descriptor        | Description                   |
|-------------------|-------------------------------|
| MaxAbsEStateIndex | Maximum absolute EState index |
| MaxEStateIndex    | Maximum EState index          |
| MinEStateIndex    | Minimum EState index          |
| MinAbsEStateIndex | Minimum absolute EState index |

### Fragments

Bunch of fragment descriptors from a file

| Descriptor             | Description                                                                       |
|------------------------|-----------------------------------------------------------------------------------|
| fr_Al_COO              | Number of aliphatic carboxylic acids                                              |
| fr_Al_OH               | Number of aliphatic hydroxyl groups                                               |
| fr_Al_OH_noTert        | Number of aliphatic hydroxyl groups excluding tert-OH                             |
| fr_ArN                 | Number of N functional groups attached to aromatics                               |
| fr_Ar_COO              | Number of Aromatic carboxylic acide                                               |
| fr_Ar_N                | Number of aromatic nitrogens                                                      |
| fr_Ar_NH               | Number of aromatic amines                                                         |
| fr_Ar_OH               | Number of aromatic hydroxyl groups                                                |
| fr_COO                 | Number of carboxylic acids                                                        |
| fr_COO2                | Number of carboxylic acids                                                        |
| fr_C_O                 | Number of carbonyl O                                                              |
| fr_C_O_noCOO           | Number of carbonyl O, excluding COOH                                              |
| fr_C_S                 | Number of thiocarbonyl                                                            |
| fr_HOCCN               | Number of C(OH)CCN-Ctert-alkyl or C(OH)CCNcyclic                                  |
| fr_Imine               | Number of Imines                                                                  |
| fr_NH0                 | Number of Tertiary amines                                                         |
| fr_NH1                 | Number of Secondary amines                                                        |
| fr_NH2                 | Number of Primary amines                                                          |
| fr_N_O                 | Number of hydroxylamine groups                                                    |
| fr_Ndealkylation1      | Number of XCCNR groups                                                            |
| fr_Ndealkylation2      | Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)      |
| fr_Nhpyrrole           | Number of H-pyrrole nitrogens                                                     |
| fr_SH                  | Number of thiol groups                                                            |
| fr_aldehyde            | Number of aldehydes                                                               |
| fr_alkyl_carbamate     | Number of alkyl carbamates (subject to hydrolysis)                                |
| fr_alkyl_halide        | Number of alkyl halides                                                           |
| fr_allylic_oxid        | Number of allylic oxidation sites excluding steroid dienone                       |
| fr_amide               | Number of amides                                                                  |
| fr_amidine             | Number of amidine groups                                                          |
| fr_aniline             | Number of anilines                                                                |
| fr_aryl_methyl         | Number of aryl methyl sites for hydroxylation                                     |
| fr_azide               | Number of azide groups                                                            |
| fr_azo                 | Number of azo groups                                                              |
| fr_barbitur            | Number of barbiturate groups                                                      |
| fr_benzene             | Number of benzene rings                                                           |
| fr_benzodiazepine      | Number of benzodiazepines with no additional fused rings                          |
| fr_bicyclic            | Bicyclic                                                                          |
| fr_diazo               | Number of diazo groups                                                            |
| fr_dihydropyridine     | Number of dihydropyridines                                                        |
| fr_epoxide             | Number of epoxide rings                                                           |
| fr_ester               | Number of esters                                                                  |
| fr_ether               | Number of ether oxygens (including phenoxy)                                       |
| fr_furan               | Number of furan rings                                                             |
| fr_guanido             | Number of guanidine groups                                                        |
| fr_halogen             | Number of halogens                                                                |
| fr_hdrzine             | Number of hydrazine groups                                                        |
| fr_hdrzone             | Number of hydrazone groups                                                        |
| fr_imidazole           | Number of imidazole rings                                                         |
| fr_imide               | Number of imide groups                                                            |
| fr_isocyan             | Number of isocyanates                                                             |
| fr_isothiocyan         | Number of isothiocyanates                                                         |
| fr_ketone              | Number of ketones                                                                 |
| fr_ketone_Topliss      | Number of ketones excluding diaryl, a,b-unsat                                     |
| fr_lactam              | Number of beta lactams                                                            |
| fr_lactone             | Number of cyclic esters (lactones)                                                |
| fr_methoxy             | Number of methoxy groups -OCH3                                                    |
| fr_morpholine          | Number of morpholine rings                                                        |
| fr_nitrile             | Number of nitriles                                                                |
| fr_nitro               | Number of nitro groups                                                            |
| fr_nitro_arom          | Number of nitro benzene ring substituent                                          |
| fr_nitro_arom_nonortho | Number of non-ortho nitro benzene ring substituents                               |
| fr_nitroso             | Number of nitroso groups, excluding NO2                                           |
| fr_oxazole             | Number of oxazole rings                                                           |
| fr_oxime               | Number of oxime groups                                                            |
| fr_para_hydroxylation  | Number of para-hydroxylation sites                                                |
| fr_phenol              | Number of phenols                                                                 |
| fr_phenol_noOrthoHbond | Number of phenolic OH excluding ortho intramolecular Hbond substituents           |
| fr_phos_acid           | Number of phosphoric acid groups                                                  |
| fr_phos_ester          | Number of phosphoric ester groups                                                 |
| fr_piperdine           | Number of piperdine rings                                                         |
| fr_piperzine           | Number of piperzine rings                                                         |
| fr_priamide            | Number of primary amides                                                          |
| fr_prisulfonamd        | Number of primary sulfonamides                                                    |
| fr_pyridine            | Number of pyridine rings                                                          |
| fr_quatN               | Number of quarternary nitrogens                                                   |
| fr_sulfide             | Number of thioether                                                               |
| fr_sulfonamd           | Number of sulfonamides                                                            |
| fr_sulfone             | Number of sulfone groups                                                          |
| fr_term_acetylene      | Number of terminal acetylenes                                                     |
| fr_tetrazole           | Number of tetrazole rings                                                         |
| fr_thiazole            | Number of tetrazole rings                                                         |
| fr_thiocyan            | Number of thiocyanates                                                            |
| fr_thiophene           | Number of thiophene rings                                                         |
| fr_unbrch_alkane       | Number of unbranched alkanes of at least 4 members (excludes halogenated alkanes) |
| fr_urea                | Number of urea groups                                                             |

### Lipinski

Lipinski parameters for molecules

| Descriptor               | Description                                                                                     |
|--------------------------|-------------------------------------------------------------------------------------------------|
| FractionCSP3             | The fraction of C atoms that are SP3 hybridized                                                 |
| HeavyAtomCount           | The number of heavy atoms a molecule                                                            |
| NHOHCount                | The number of NHs or OHs                                                                        |
| NOCount                  | The number of Nitrogens and Oxygens                                                             |
| NumAliphaticCarbocycles  | The number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule  |
| NumAliphaticHeterocycles | The number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule |
| NumAliphaticRings        | The number of aliphatic (containing at least one non-aromatic bond) rings for a molecule        |
| NumAromaticCarbocycles   | The number of aromatic carbocycles for a molecule                                               |
| NumAromaticHeterocycles  | The number of aromatic heterocycles for a molecule                                              |
| NumAromaticRings         | The number of aromatic rings for a molecule                                                     |
| NumHAcceptors            | The number of Hydrogen Bond Acceptors                                                           |
| NumHDonors               | The number of Hydrogen Bond Donors                                                              |
| NumHeteroatoms           | The number of Heteroatoms                                                                       |
| NumRotatableBonds        | The number of Rotatable Bonds                                                                   |
| NumSaturatedCarbocycles  | The number of saturated carbocycles for a molecule                                              |
| NumSaturatedHeterocycles | The number of saturated heterocycles for a molecule                                             |
| NumSaturatedRings        | The number of saturated rings for a molecule                                                    |
| RingCount                | Ring count                                                                                      |

See also:

* [Molecular descriptor](https://en.wikipedia.org/wiki/Molecular_descriptor)
* [RDKit](http://rdkit.org)
