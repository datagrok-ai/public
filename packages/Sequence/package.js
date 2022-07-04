class SequencePackage extends DG.Package {

  //description: Returns a string containing the Amino Acid sequence represented by the nucleotide sequence
  //tags: converter
  //input: string nucleotides {semType: nucleotides}
  //output: string aminoacids {semType: aminoacids}
  nucleotidesToAminoacids(nucleotides) {
    let seq = new Nt.Seq();
    seq.read(nucleotides);
    return seq.translate();
  }

  toFasta(sequences, labels) {
    let result = '';
    for (let i = 0; i < sequences.length; i++) {
      let label = labels === null ? `seq${i}` : labels[i];
      result += `>${label}\n${sequences[i]}\n`;
    }
    return result;
  }

  //description: returns complementary sequence
  //input: string nucleotides {semType: nucleotides}
  //output: string result {semType: nucleotides}
  complement(nucleotides) {
    let seq = new Nt.Seq();
    seq.read(nucleotides);
    return seq.complement().sequence();
  }

  //description: returns a sequence widget
  //input: string nucleotides {semType: nucleotides}
  //tags: panel
  //output: widget result
  sequenceWidget(nucleotides) {
    let e = document.createElement('DIV');
    e.id = 'sq77';
    let aa = this.nucleotidesToAminoacids(nucleotides);

    setTimeout(function () {
      new FeatureViewer(aa,
        '#sq77',
        {
          showAxis: true,
          showSequence: true,
          brushActive: true, //zoom
          toolbar: true, //current zoom & mouse position
          bubbleHelp: true,
          zoomMax: 50 //define the maximum range of the zoom
        });
    }, 1000);

    return new DG.Widget(e);
  }

  //name: sequenceCellRenderer
  //tags: cellRenderer
  //meta-cell-renderer-sem-type: nucleotides
  //output: grid_cell_renderer result
  sequenceCellRenderer() {
    return new SequenceCellRenderer();
  }

  //name: Sequence
  //description: creates a sequence viewer
  //tags: viewer
  //output: viewer result
  sequenceViewer() {
    return new SequenceViewer();
  }

  //name: WebLogoViewer
  //tags: viewer,panel
  //output: viewer result
  webLogoViewer() {
    return new WebLogoViewer();
  }

  // TODO: condition detector to be connected with the file-handler (currently impossible)
  //
  //input: file file
  //output: bool
  //
  isFastaFile(file) {
    let extensions = ["fasta", "aln"];
    return extensions.includes(file.extension);
  }

  //input: string content
  //output: list tables
  //tags: file-handler
  //meta.ext: aln
  fastaAlnFileHandler(fastaString) {
    let fr = new FastaRead();
    fr.read(fastaString);
    if (fr.molCnt <= 0) {
      return;
    }

    let t = DG.DataFrame.create(fr.molCnt);
    t.columns.add(DG.Column.fromStrings('accession', fr.acc_));
    t.columns.add(DG.Column.fromStrings('description', fr.descr_));

    let seqCol = DG.Column.fromStrings('sequence', fr.seq_);
    seqCol.setTag(".mono", "1"); // add tags to mark-up sequence to display monospace
    t.columns.add(seqCol);

    let frc = fr.free_columns_;
    for (let [key, value] of frc) {
      t.columns.add(DG.Column.fromStrings(key, value));
    }
    return [t];
  }

  //input: string content
  //output: list tables
  //tags: file-handler
  //meta.ext: fasta
  fastaFileHandler(fastaString) {
    return this.fastaAlnFileHandler(fastaString);
  }

}

let _seq = new SequencePackage();

class SequenceViewer extends DG.JsViewer {
  onTableAttached() {
    let seqCol = this.dataFrame.columns.toList().find((c) => c.semType === 'nucleotides');
    let fasta = _seq.toFasta(Array.from(seqCol.values()), null);
    let seqs = msa.io.fasta.parse(fasta);

    this.msa = msa({
      el: this.root,
      seqs: seqs
    });
    this.msa.render();

    this.dataFrame.selection.onChanged(() => this.render());
    this.dataFrame.filter.onChanged(() => this.render());
    this.render();
  }

  // reflect changes made to filter/selection
  render() {
  }
}
