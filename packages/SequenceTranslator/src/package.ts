/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

//name: Sequence Translator
//tags: app
export function sequenceTranslator(): void {
  let inputField = ui.textInput('', "");
  let outputValues = {type: 'Input', classic: '', bioSpring: '', axolabs: '', gcrs: '', complement: ''}
  draw(inputField, outputValues);
  inputField.onChanged(() => {
    let pi = DG.TaskBarProgressIndicator.create('Converting sequences...');
    outputValues = convertSequence(inputField.value);
    draw(inputField, outputValues);
    pi.close();
  });
}

function copyToClipboard(text: string) {
  let dummy = document.createElement("textarea");
  document.body.appendChild(dummy);
  dummy.value = text;
  dummy.select();
  document.execCommand("copy");
  document.body.removeChild(dummy);
}

function draw(inputSeq: DG.InputBase, outputValues: {type: string; classic: string; bioSpring: string; complement: string; gcrs: string; axolabs: string;}) {
  let view = grok.shell.newView('Sequence Translator: ' + inputSeq.value, []);
  view.append(
    ui.divV([
      ui.divText('Paste sequence in text field below', 'grok-datajob-publish-alert'),
      ui.divH([
        ui.div([
          ui.h1("Input description: " + outputValues.type),
          ui.inputs([
            inputSeq,
          ], {})
        ], 'sequenceInput'),
        ui.div([
          ui.divH([ui.h1('Classic Code'), ui.iconFA('copy', () => copyToClipboard(outputValues.classic))]),
          ui.divText(outputValues.classic),
          ui.divH([ui.h1('BioSpring Code'), ui.iconFA('copy', () => copyToClipboard(outputValues.bioSpring))]),
          ui.divText(outputValues.bioSpring),
          ui.divH([ui.h1('Axolabs Code'), ui.iconFA('copy', () => copyToClipboard(outputValues.axolabs))]),
          ui.divText(outputValues.axolabs),
          ui.divH([ui.h1('GCRS Code'), ui.iconFA('copy', () => copyToClipboard(outputValues.gcrs))]),
          ui.divText(outputValues.gcrs),
          ui.divH([ui.h1('Complement'), ui.iconFA('copy', () => copyToClipboard(outputValues.complement))]),
          ui.divText(outputValues.complement)
        ],'sequenceOutput')
      ], 'sequence')
    ])
  );
  // @ts-ignore
  $('.sequence').css('padding','20px 0');
  // @ts-ignore
  $('.sequence').children().css('padding','0 10px');
  // @ts-ignore
  $('.sequenceOutput').attr('style','word-break:break-word');
  // @ts-ignore
  $('.sequenceInput .input-base').css('margin','15px 0');
  // @ts-ignore
  $('.sequenceInput textarea').css('resize','none');
  // @ts-ignore
  $('.sequenceInput textarea').css('min-height','100px');
  // @ts-ignore
  $('.sequenceInput textarea').css('width','100%');
  // @ts-ignore
  $('.sequenceInput select').css('width','100%');
  // @ts-ignore
  $('.sequenceInput textarea').attr('placeholder','Type here');
  // @ts-ignore
  $('.sequenceInput select').attr('placeholder','');
}

export function isClassicCode(sequence: string): boolean {return /^[ATGCU]{10,}$/.test(sequence);}

export function isAsoGapmerBioSpringCode(sequence: string): boolean {return /^[*56789ATGC]{30,}$/.test(sequence);}

export function isAsoGapmerGCRSCode(sequence: string): boolean {return /^(?=.*moe)(?=.*5mC)(?=.*ps){30,}/.test(sequence);}

export function isSiRnaBioSpringCode(sequence: string): boolean {return /^[*1-8]{30,}$/.test(sequence);}

export function isSiRnaAxolabsCode(sequence: string): boolean {return /^[fsACGUacgu]{30,}$/.test(sequence);}

export function isSiRnaGCRSCode(sequence: string): boolean {return /^[fmpsACGU]{30,}$/.test(sequence);}

function convertSequence(seq: string) {
  if (isClassicCode(seq))
    return {
      type: "ASO Gapmers / Classic Code",
      classic: seq,
      bioSpring: asoGapmersClassicToBioSpring(seq),
      axolabs: "No code accordance provided",
      gcrs: asoGapmersClassicToGCRS(seq),
      complement: asoGapmersClassicComplement(seq)
  };
  if (isAsoGapmerBioSpringCode(seq))
    return {
      type: "ASO Gapmers / BioSpring Code",
      classic: asoGapmersBioSpringToClassic(seq),
      bioSpring: seq,
      axolabs: "No code accordance provided",
      gcrs: asoGapmersBioSpringToGCRS(seq),
      complement: asoGapmersBioSpringComplement(seq)
    };
  if (isAsoGapmerGCRSCode(seq))
    return {
      type: "ASO Gapmers / GCRS Code",
      classic: asoGapmersGCRSToClassic(seq),
      bioSpring: asoGapmersGCRSToBioSpring(seq),
      axolabs: "No code accordance provided",
      gcrs: seq,
      complement: asoGapmersGCRSComplement(seq)
    };
  if (isSiRnaBioSpringCode(seq))
    return {
      type: "siRNA / bioSpring Code",
      classic: "coming soon",
      bioSpring: seq,
      axolabs: "coming soon",
      gcrs: "coming soon",
      complement: "coming soon"
    };
  if (isSiRnaAxolabsCode(seq))
    return {
      type: "siRNA / Axolabs Code",
      classic: siRnaAxolabsToClassic(seq),
      bioSpring: siRnaAxolabsToBioSpring(seq),
      axolabs: seq,
      gcrs: siRnaAxolabsToGCRS(seq),
      complement: "coming soon"
    };
  if (isSiRnaGCRSCode(seq))
    return {
      type: "siRNA / GCRS Code",
      classic: "coming soon",
      bioSpring: "coming soon",
      axolabs: siRnaGCRSToAxolabs(seq),
      gcrs: seq,
      complement: "coming soon"
    };
  return {
    type: "Type of input sequence is undefined",
    classic: "Type of input sequence is undefined",
    bioSpring: "Type of input sequence is undefined",
    axolabs: "Type of input sequence is undefined",
    gcrs: "Type of input sequence is undefined",
    complement: "Type of input sequence is undefined"
  };
}

//name: asoGapmersClassicToBioSpring
//input: string nucleotides {semType: nucleotides}
//output: string result {semType: BioSpring / Gapmers}
export function asoGapmersClassicToBioSpring(nucleotides: string) {
  let count: number = -1;
  const objForEdges: {[index: string]: string} = {"T": "5*", "A": "6*", "C": "7*", "G": "8*"};
  const objForCenter: {[index: string]: string} = {"C": "9*", "A": "A*", "T": "T*", "G": "G*"};
  return nucleotides.replace(/[ATCG]/g, function (x: string) {
    count++;
    if (count < 5) return objForEdges[x];
    if (count < 15) return objForCenter[x];
    return objForEdges[x];
  }).slice(0, 2 * count + 1);
}

//name: asoGapmersClassicToGCRS
//input: string nucleotides {semType: nucleotides}
//output: string result {semType: GCRS / Gapmers}
export function asoGapmersClassicToGCRS(nucleotides: string) {
  let count: number = -1;
  const objForEdges: {[index: string]: string} = {"T": "moeUnps", "A": "moeAnps", "C": "moe5mCnps", "G": "moeGnps"};
  const objForCenter: {[index: string]: string} = {"C": "Cps", "A": "Aps", "T": "Tps", "G": "Gps"};
  return nucleotides.replace(/[ATCG]/g, function (x: string) {
    count++;
    if (count < 5) return objForEdges[x];
    if (count < 15) return objForCenter[x];
    return objForEdges[x];
  }).slice(0, -3);
}

//name: asoGapmersBioSpringToClassic
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: nucleotides}
export function asoGapmersBioSpringToClassic(nucleotides: string) {
  const obj: {[index: string]: string} = {"*": "", "5": "T", "6": "A", "7": "C", "8": "G", "9": "C"};
  return nucleotides.replace(/[*56789]/g, function (x: string) {return obj[x];});
}

//name: asoGapmersBioSpringToGCRS
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: GCRS / Gapmers}
export function asoGapmersBioSpringToGCRS(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "5*": "moeUnps", "6*": "moeAnps", "7*": "moe5mCnps", "8*": "moeGnps", "9*": "5mCps", "A*": "Aps", "T*": "Tps",
    "G*": "Gps", "C*": "Cps", "5": "moeU", "6": "moeA", "7": "moe5mC", "8": "moeG"
  };
  return nucleotides.replace(/(5\*|6\*|7\*|8\*|9\*|A\*|T\*|G\*|C\*|5|6|7|8])/g, function (x: string) {return obj[x];});
}

//name: asoGapmersGCRSToBioSpring
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: BioSpring / Gapmers}
export function asoGapmersGCRSToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "moeT": "5", "moeA": "6", "moe5mC": "7", "moeG": "8", "moeU": "5", "5mC": "9", "nps": "*", "ps": "*", "U": "T"
  };
  return nucleotides.replace(/(moeT|moeA|moe5mC|moeG|moeU|5mC|nps|ps|U)/g, function (x: string) {return obj[x];});
}

//name: asoGapmersGCRSToClassic
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: nucleotides}
export function asoGapmersGCRSToClassic(nucleotides: string) {
  const obj: {[index: string]: string} = {"moe": "", "5m": "", "n": "", "ps": "", "U": "T"};
  return nucleotides.replace(/(moe|5m|n|ps|U)/g, function (x: string) {return obj[x];});
}

//name: siRnaGCRSToAxolabs
//input: string nucleotides {semType: GCRS / siRNA}
//output: string result {semType: Axolabs / siRNA}
export function siRnaGCRSToAxolabs(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "fU": "Uf", "fA": "Af", "fC": "Cf", "fG": "Gf", "mU": "u", "mA": "a", "mC": "c", "mG": "g", "ps": "s"
  };
  return nucleotides.replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x: string) {return obj[x];});
}

//name: siRnaAxolabsToGCRS
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: GCRS / siRNA}
export function siRnaAxolabsToGCRS(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "Uf": "fU", "Af": "fA", "Cf": "fC", "Gf": "fG", "u": "mU", "a": "mA", "c": "mC", "g": "mG", "s": "ps"
  };
  return nucleotides.replace(/(Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: siRnaAxolabsToBioSpring
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: BioSpring / siRNA}
export function siRnaAxolabsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "Uf": "1", "Af": "2", "Cf": "3", "Gf": "4", "u": "5", "a": "6", "c": "7", "g": "8", "s": "*"
  };
  return nucleotides.replace(/(Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: siRnaAxolabsToClassic
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: nucleotides}
export function siRnaAxolabsToClassic(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "Uf": "U", "Af": "A", "Cf": "C", "Gf": "G", "u": "U", "a": "A", "c": "C", "g": "G", "s": ""
  };
  return nucleotides.replace(/(Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: asoGapmersClassicComplement
//input: string nucleotides {semType: nucleotides}
//output: string result {semType: nucleotides}
export function asoGapmersClassicComplement(nucleotides: string) {
  const obj: {[index: string]: string} = {"A": "T", "T": "A", "G": "C", "C": "G"};
  return nucleotides.replace(/[ATGC]/g, function (x: string) {return obj[x];});
}

//name: asoGapmersBioSpringComplement
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: BioSpring / Gapmers}
export function asoGapmersBioSpringComplement(nucleotides: string) {
  const obj: {[index: string]: string} = {"A": "T", "T": "A", "G": "C", "C": "G", "5": "6", "6": "5", "7": "8", "8": "9", "9": "8"};
  return nucleotides.replace(/[56789ATGC]/g, function (x: string) {return obj[x];});
}

//name: asoGapmersGCRSComplement
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: GCRS / Gapmers}
export function asoGapmersGCRSComplement(nucleotides: string) {
  const obj: {[index: string]: string} = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"};
  return nucleotides.replace(/[ACGTU]/g, function (x: string) {return obj[x];});
}