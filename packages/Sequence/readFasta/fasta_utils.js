/// Simple FASTA reader class
///
class FastaRead {
  constructor() {
    this.acc_ = []; this.descr_ = []; this.seq_ = [];
    this.re_ = new RegExp(/\[.*?\]/g);
    this.free_columns_ = new Map();
    this.molCnt = 0;
  }

  reset()      { this.molCnt = 0; }
  resetAcc()   { this.acc_ = [];   }
  resetDescr() { this.descr_ = []; }
  resetSeq()   { this.seq_ = [];   }

  add(acc, descr, seq) {
    if (!acc || acc == '') // accession must be defined
      return;
    this.acc_.push(acc);
    this.seq_.push(seq == '' ? null : seq);
    if (descr != null && descr != '') {
      this.descr_.push(descr);
      let keyValArr = this.parseDescr(descr);
      if (keyValArr != null && keyValArr.length > 0) {
        for (let i = 0; i < keyValArr.length; i++) {
          const colName = keyValArr[i].key;
          const featValue = keyValArr[i].val;
          if (this.free_columns_.has(colName)) {
            this.free_columns_.get(colName).push(featValue);
          } else { // create new column
            let valueArr = [];
            for (let j = 0; j < this.molCnt; j++) {
              valueArr.push(null);
            }
            valueArr.push(featValue);
            this.free_columns_.set(colName, valueArr);
          }
        } // for i
        // make sure all free columns has the same element count
        let rowCnt = this.molCnt+1;
        for (let [k, colArr] of this.free_columns_) {
          if (colArr.length < rowCnt)
            colArr.push(null);
        }

      }
    }
    else { // no description (add NULL to descr and all derived columns)
      this.descr.push(null);
      for (let [key, valueArr] of this.free_columns_) {
        valueArr.push(null);
      }
    }
    this.molCnt++;
  }

  /// parse a feature as [key=value], return key-value object
  parseFeature(featStr) {
    featStr = featStr.trim();
    let keyStr = new String();
    let valStr = new String();
    let keyStartIdx = 0; let keyEndIdx = 0; let valStartIdx = 0; let valEndIdx = 0;
    let matchEq = false;
    for (let i = 0; i < featStr.length; i++) {
      let ch = featStr[i];
      switch (ch) {
        case '[':
          keyStartIdx = i;
          break;
        case ']':
          valEndIdx = i;
          break;
        case '=':
          keyEndIdx = i; valStartIdx = i;
          matchEq = true;
          break;
      } // switch
      if (ch == ']')
        break;
    } // for i

    if (!matchEq)
      return null;
    keyStr = featStr.substring(keyStartIdx+1, keyEndIdx);
    valStr = featStr.substring(valStartIdx+1, valEndIdx);

    return { key: keyStr.trim(), val: valStr.trim() } ;
  }

  parseDescr(desrStr) {
    let keyValArr = [];
    let featIter = desrStr.matchAll(this.re_);
    for (const feat of featIter) {
      let featStr = feat[0];
      const featPair = this.parseFeature(featStr);
      if (featPair != null) {
        keyValArr.push(featPair);
      }
    } // for feat
    return keyValArr;
  }


  /// return number of molecules parsed
  readLen() { return this.acc_.length; }

  read(fasta) {
    if (!fasta)
      return false;

    // split the input on newlines... TODO: read it line by line
    let lines = fasta.split('\n');

    let acc = new String();
    let defline = new String(); // FASTA defline
    let seqArr = []; // accumulator for sequence strings belonging to one mol

    for (var i = 0; i < lines.length; i++) {
      let str = lines[i];
      str = str.replaceAll('\r', '').trim();
      if (str == '')
        continue;

      if (str[0] == '>') { // defline
        let seq = seqArr.join('');
        this.add(acc, defline, seq);
        acc = ''; defline = ''; seqArr = []; seq = '';

        const indSpace = str.indexOf(' ');
        if (indSpace > 0) { //found
          acc = str.substring(0, indSpace);
          defline = str.substring(indSpace+1);
        } else {
          acc = str;
          defline = '';
        }
        acc = acc.substring(1); // remove '>'
        seq = '';
      } else { // new sequence line
        seqArr.push(str); // add to the collection of sequence strings
      }
    } // for i

    // add the last accumulated FASTA element
    //
    let seq = seqArr.join('');
    this.add(acc, defline, seq);
  }
}
