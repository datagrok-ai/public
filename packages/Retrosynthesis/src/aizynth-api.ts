
export interface ReactionMetadata {
    template_hash: string;
    classification: string;
    library_occurence: number;
    policy_probability: number;
    policy_probability_rank: number;
    policy_name: string;
    template_code: number;
    template: string;
    mapped_reaction_smiles: string;
}

export interface Scores {
    'state score': number;
    'number of reactions': number;
    'number of pre-cursors': number;
    'number of pre-cursors in stock': number;
    'average template occurrence': number;
}

export interface TreeNode {
    type: string;
    hide: boolean;
    smiles: string;
    is_chemical?: boolean;
    is_reaction?: boolean;
    in_stock?: boolean;
    metadata?: ReactionMetadata;
    children?: TreeNode[];
}

export interface TreeMetadata {
    created_at_iteration: number;
    is_solved: boolean;
}

export interface Tree {
    type: string;
    hide: boolean;
    smiles: string;
    is_chemical: boolean;
    in_stock: boolean;
    children?: TreeNode[];
    scores: Scores;
    metadata: TreeMetadata;
}

export interface DataEntry {
    index: number;
    target: string;
    search_time: number;
    first_solution_time: number;
    first_solution_iteration: number;
    number_of_nodes: number;
    max_transforms: number;
    max_children: number;
    number_of_routes: number;
    number_of_solved_routes: number;
    top_score: number;
    is_solved: boolean;
    number_of_steps: number;
    number_of_precursors: number;
    number_of_precursors_in_stock: number;
    precursors_in_stock: string;
    precursors_not_in_stock: string;
    precursors_availability: string;
    policy_used_counts: PolicyUsedCounts;
    profiling: Profiling;
    stock_info: StockInfo;
    trees: Tree[];
}

export interface SchemaField {
    name: string;
    type: string;
}

export interface Schema {
    fields: SchemaField[];
    primaryKey: string[];
    pandas_version: string;
}

export interface PolicyUsedCounts {
    [key: string]: number;
}

export interface Profiling {
    expansion_calls: number;
    reactants_generations: number;
    iterations: number;
}

export interface StockInfo {
    [key: string]: string[];
}

export interface ReactionData {
    schema: Schema;
    data: DataEntry[];
}

