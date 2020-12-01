// Through RDKit we operate on all acceptable strings representing

function _morganFP(molString, fp_length = 128, fp_radius = 2) {
    let mol = Module.get_mol(molString);
    let mfp = mol.get_morgan_fp(fp_radius, fp_length);
    mol.delete();
    return mfp;
}

function moleculesToFingerprints(molStringsColumn, settings) {
    const len = molStringsColumn.length;
    const fpLength = settings.hasOwnProperty('fpLength') ? settings.fpLength : 128;
    const fpRadius = settings.hasOwnProperty('fpRadius') ? settings.fpRadius : 2;
    let fingerprints = molStringsColumn.toList().map((molString) =>
      DG.BitSet.fromString(_morganFP(molString, fpLength, fpRadius)));
    return DG.Column.fromList('object', 'fingerprints', fingerprints);
}

function _foldFingerprint(bitsetFp, newLength) {
    let result = DG.BitSet.create(newLength);
    for (let idx in bitsetFp.getSelectedIndexes())
        result.set(idx % newLength, true, false);
    return result;
}

function fingerprintSimilarity(bitsetFp1, bitsetFp2) {
    const len1 = bitsetFp1.length;
    const len2 = bitsetFp2.length;
    if (len1 < len2)
        bitsetFp2 = _foldFingerprint(bitsetFp2, len1);
    else if (len2 < len1)
        bitsetFp1 = _foldFingerprint(bitsetFp1, len2);
    return bitsetFp1.similarityTo(bitsetFp2, 'tanimoto'); // tanimotoSimilarity(fp1, fp2);
}

// molString can be any string type RDKit supports: smiles, MolBlock
// This also applies to cells content of molColumn
function chemSimilarityScoring(molStringsColumn, molString, settings) {
    
    const fingerprint = moleculesToFingerprints(DG.Column.fromStrings('molecules', [molString]), settings).get(0);
    const fingerprintCol = moleculesToFingerprints(molStringsColumn, settings);
    const len = molStringsColumn.length;
    
    let distances = DG.Column.fromType(DG.TYPE.FLOAT, 'distances', len);
    for (let row = 0; row < len; ++row) {
        const fp = fingerprintCol.get(row);
        distances.set(row, fp === null ? 100.0 : fingerprintSimilarity(fingerprint, fp));
    }
    
    if (settings.hasOwnProperty('sorted') && settings.sorted === true) {
        
        const limit = Math.min((settings.hasOwnProperty('limit') ? settings.limit : len), len);
        const minScore = settings.hasOwnProperty('minScore') ? settings.minScore : 0.0;
        let sortedIndices = Array.from(Array(len).keys()).sort((i1, i2) => {
            const a1 = distances.get(i1);
            const a2 = distances.get(i2);
            if (a2 < a1) return -1;
            if (a2 > a1) return +1;
            return 0; // a2.compareTo(a1)
          });
        let sortedMolStrings = DG.Column.fromType(DG.TYPE.STRING, 'molecule', limit);
        sortedMolStrings.semType = DG.SEMTYPE.MOLECULE;
        let sortedScores = DG.Column.fromType(DG.TYPE.FLOAT, 'score', limit);
    
        for (let n = 0; n < limit; n++) {
            const idx = sortedIndices[n];
            const score = distances.get(idx);
            if (score < minScore) {
                sortedMolStrings.removeAt(n, limit - n);
                sortedScores.removeAt(n, limit - n);
                break;
            }
            sortedMolStrings.set(n, molStringsColumn.get(idx));
            sortedScores.set(n, score);
        }
    
        return DG.DataFrame.fromColumns([sortedMolStrings, sortedScores]);
        
    } else {
        
        return distances;
        
    }
}

function chemSubstructureSearchGraph(molStringsColumn, molString) {
    
    const len = molStringsColumn.length;
    let result = DG.BitSet.create(len);
    let subMol = Module.get_mol(molString);
    for (let i = 0; i < len; ++i) {
        let mol = Module.get_mol(molStringsColumn.get(i));
        let match = mol.get_substruct_match(subMol);
        if (match !== "{}")
            result.set(i, true, false);
        mol.delete();
    }
    subMol.delete();
    return result;
}

function chemSubstructureSearchLibrary(molStringsColumn, molString) {
    
    let df = molStringsColumn.dataFrame;
    let col = molStringsColumn;
    let library = null;
    if (df == null) {
        // no caching mode
        library = new Module.SubstructLibrary();
    } else {
        // caching mode
        df.onValuesChanged.subscribe((_) => {
            if (col != null && col.temp.rdkitLirary != null) {
                col.temp.rdkitLirary.delete();
                col.temp.rdkitLirary = null;
            }
        });
        if (col.temp.rdkitLirary == null) {
            col.temp.rdkitLirary = new Module.SubstructLibrary();
            for (let i = 0; i < col.length; ++i) {
                const smiles = molStringsColumn.get(i);
                col.temp.rdkitLirary.add_trusted_smiles(smiles);
            }
        }
        library = col.temp.rdkitLirary;
    }
  
    var query = Module.get_qmol(molString);
    const matches = JSON.parse(library.get_matches(query));
    query.delete();
    let result = DG.BitSet.create(col.length);
    for (let m of matches) {
        result.set(m, true, false);
    }
    
    if (df == null) {
        library.delete();
    }
    
    return result;
}