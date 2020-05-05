
class chem {
    static similaritySearch(column, molecule, metric = METRIC_TANIMOTO, limit = 10, minScore = 0.7) {
        return new Promise((resolve, reject) => grok_Chem_SimilaritySearch(column.d, molecule, metric,
            limit, minScore, (t) => resolve(new DataFrame(t))));
    }

    static diversitySearch(column, metric = METRIC_TANIMOTO, limit = 10) {
        return new Promise((resolve, reject) => grok_Chem_DiversitySearch(column.d, metric, limit, (mols) => resolve(mols)));
    }

    static substructureSearch(column, pattern, isSmarts = true) {
        return new Promise((resolve, reject) => grok_Chem_SubstructureSearch(column.d, pattern, isSmarts, (bs) => resolve(new BitSet(bs))));
    }

    static rGroup(table, column, core) {
        return new Promise((resolve, reject) => grok_Chem_RGroup(table.d, column, core, () => resolve(table)));
    }

    static mcs(column) {
        return new Promise((resolve, reject) => grok_Chem_MCS(column.d, (mcs) => resolve(mcs)));
    }

    static descriptors(table, column, descriptors) {
        return new Promise((resolve, reject) => grok_Chem_Descriptors(table.d, column, descriptors, () => resolve(table)));
    }
}

