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
    if (col.name === "pdb")
      return 'pdb_id';
  }

  //input: column col { semType: nucleotides }
  sequenceNucToRobot(nuc) {
    nuc.dataFrame.columns
      .addNewString('robot')
      .init((i) => 'converted ' + nuc.get(i));
  }
}