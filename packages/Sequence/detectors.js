class SequencePackageDetectors extends DG.Package {

  //tags: autostart
  autostartTest() {
    console.log('sequence autostarted.');
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    const rex = /^[ATGCNM-]{30,}$/;
    if (DG.Detector.sampleCategories(col, (s) => rex.test(s)))
      return 'nucleotides';
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectPdb(col) {
    if (col.name === "pdb") {
      col.semType = 'pdb_id';
      return 'pdb_id';
    }

    return null;
  }

  //input: string s
  //output: string ss
  dupXX(s) {
    return s + '_' + s + ' (sequence)';
  }


  //input: string s
  //output: string ss
  dup(s) {
    return s + '_' + s + ' (sequence)';
  }

  //input: column col { semType: nucleotides }
  sequenceNucToRobot(nuc) {
    nuc.dataFrame.columns
      .addNewString('robot')
      .init((i) => 'converted ' + nuc.get(i));
  }

  //input: string seq { semType: sequence }
  //input: string seqUnits { semType: sequenceUnits }
  //output: string seq { semType: sequence }
  // sequenceTranslator(nuc) {
  //   nuc.dataFrame.columns
  //     .addNewString('robot')
  //     .init((i) => 'converted ' + nuc.get(i));
  // }
}