class ProteomicsPackageDetectors extends DG.Package {
  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectProteinId(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;
    const name = col.name.toLowerCase();
    if (name === 'protein ids' || name === 'protein id' || name === 'majority protein ids' ||
      name === 'primary protein id' || name === 'leading protein' || name === 'leading razor protein' ||
      name === 'uniprot' || name === 'accession') {
      const sample = col.toList().slice(0, 20).filter((v) => v != null);
      const uniprotPattern = /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$/;
      if (sample.some((v) => uniprotPattern.test(v.split(';')[0]))) {
        col.semType = 'Proteomics-ProteinId';
        return col.semType;
      }
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectGeneSymbol(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;
    const name = col.name.toLowerCase();
    if (name === 'gene names' || name === 'gene name' || name === 'gene symbol' || name === 'gene') {
      col.semType = 'Proteomics-GeneSymbol';
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectSubcellularLocation(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;
    const name = col.name.toLowerCase();
    if (name === 'subcellular location' || name === 'subcellular location [cc]' ||
      name === 'subcellular' || name.includes('subcellular location') || name.includes('subcellular')) {
      col.semType = 'Proteomics-SubcellularLocation';
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectLog2FC(col) {
    if (col.type !== DG.TYPE.FLOAT && col.type !== DG.TYPE.INT)
      return null;
    const name = col.name.toLowerCase();
    if (name.includes('log2') && (name.includes('fc') || name.includes('fold'))) {
      col.semType = 'Proteomics-Log2FC';
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectPValue(col) {
    if (col.type !== DG.TYPE.FLOAT && col.type !== DG.TYPE.INT)
      return null;
    const name = col.name.toLowerCase();
    if ((name.includes('p-value') || name.includes('pvalue') || name.includes('p.value') ||
      name.includes('adj.p') || name.includes('fdr') || name.includes('q-value') || name.includes('qvalue')) &&
      col.min >= 0 && col.max <= 1) {
      col.semType = 'Proteomics-PValue';
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectIntensity(col) {
    if (col.type !== DG.TYPE.FLOAT && col.type !== DG.TYPE.INT)
      return null;
    const name = col.name.toLowerCase();
    const isRawIntensity = name.startsWith('intensity') || name.startsWith('lfq intensity') ||
      name.startsWith('reporter intensity') || name.startsWith('ibaq') ||
      name.endsWith(' maxlfq intensity') || name.endsWith(' razor intensity') ||
      name.endsWith(' intensity');
    const isLog2Intensity = (name.startsWith('log2(') || name.startsWith('log2 ')) &&
      (name.includes('intensity') || name.includes('ibaq'));
    if (isRawIntensity || isLog2Intensity) {
      col.semType = 'Proteomics-Intensity';
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectDisplayName(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;
    if (col.name.toLowerCase() === 'display name') {
      col.semType = 'Proteomics-DisplayName';
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectSourceId(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;
    if (col.name.toLowerCase() === 'source id') {
      col.semType = 'Proteomics-SourceId';
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectNumeratorMean(col) {
    if (col.type !== DG.TYPE.FLOAT && col.type !== DG.TYPE.INT)
      return null;
    if (col.name.toLowerCase() === 'numerator mean') {
      col.semType = 'Proteomics-NumeratorMean';
      return col.semType;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectDenominatorMean(col) {
    if (col.type !== DG.TYPE.FLOAT && col.type !== DG.TYPE.INT)
      return null;
    if (col.name.toLowerCase() === 'denominator mean') {
      col.semType = 'Proteomics-DenominatorMean';
      return col.semType;
    }
    return null;
  }
}
