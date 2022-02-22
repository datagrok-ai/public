/** HELM associated sdf libraries with monomer processing*/
export class MonomerLibrary {
  private static _name = 'monomerLibrary';

  private monomerFields: string[] = [
    'molecule', 'MonomerType', 'MonomerNaturalAnalogCode', 'MonomerName', 'MonomerCode', 'MonomerCaps', 'BranchMonomer',
  ];

  private library: {
    [name: string]: {
      mol: string,
      type: string,
      analogueCode: string,
      linkages: { [link: string]: { atomNumber: number, type: string } }
    }
  } = {};

  private monomers: string[] = [];

  constructor(sdf: string) {
    const sdfReader = new SDFReader();
    const data = sdfReader.getColls(sdf);
    this.monomerFields.forEach((f) => {
      if (!(f in data))
        throw new Error(`Monomer library was not compiled: ${f} field is absent in provided file`);

      if (data[f].length != data.molecule.length)
        throw new Error(`Monomer library was not compiled: ${f} field is not presented for each monomer`);
    });

    for (let i = 0; i < data.molecule.length; i++) {
      const linkData = this.getLinkData(data.molecule[i], data.MonomerCaps[i], data.MonomerName[i]);
      const entry = {
        mol: data.molecule[i],
        type: 'Peptide',
        code: data.MonomerCode[i],
        analogueCode: data.MonomerNaturalAnalogCode[i],
        linkages: linkData,
      };

      const name = data.MonomerCode[i] !== '.' ? data.MonomerCode[i] : data.MonomerName[i];
      this.library[name] = entry;
      this.monomers.push(name);
    }
  }

  /** getting full monomer information from monomer library*/
  public getMonomerEntry(name: string) {
    if (!this.monomers.includes(name))
      throw new Error(`Monomer library do not contain ${name} monomer`);


    return this.library[name];
  }

  /** getting mol as string for monomer*/
  public getMonomerMol(name: string) {
    if (!this.monomers.includes(name))
      throw new Error(`Monomer library do not contain ${name} monomer`);


    const entry = this.library[name];
    let monomerMol = entry.mol.replace(/M  RGP  .+\n/, '');

    //order matters
    const links = Object.keys(entry.linkages);
    for (let i = 0; i < links.length; i++)
      monomerMol = monomerMol.replace('R#', entry.linkages[links[i]].type + ' ');


    return monomerMol;
  }

  /** getting the list of the minomers available in library*/
  get monomerNames() {
    return this.monomers;
  }

  static get id(): string {
    return MonomerLibrary._name;
  }

  private getLinkData(mol: string, caps: string, name: string) {
    const rawData = mol.match(/M  RGP  .+/);
    if (rawData === null)
      throw new Error(`Monomer library was not compiled: ${name} entry has no RGP`);


    const types: { [code: string]: string } = {};
    caps.split('\n')?.forEach((e) => {
      types[e.match(/\d+/)![0]] = e.match(/(?<=\])\w+/)![0];
    });

    const data = rawData![0].replace('M  RGP  ', '').split(/\s+/);
    const res: { [link: string]: { atomNumber: number, type: string } } = {};
    for (let i = 0; i < parseInt(data[0]); i++) {
      const code = parseInt(data[2 * i + 2]);
      let type = '';
      switch (code) {
      case 1:
        type = 'N-terminal';
        break;
      case 2:
        type = 'C-terminal';
        break;
      case 3:
        type = 'branch';
        break;
      default:
        break;
      }
      res[type] = {atomNumber: parseInt(data[2 * i + 1]), type: types[code]};
    }

    return res;
  }
}

//TODO: merge with Chem version
class SDFReader {
  dataColls: { [_: string]: any };

  constructor() {
    this.dataColls = {'molecule': []};
  }

  getColls(content: string) {
    this.read(content);
    return this.dataColls;
  }

  read(content: string) {
    content = content.replaceAll('\r', ''); //equalize old and new sdf standards
    let startIndex = content.indexOf('$$$$', 0);
    this.parse(content, 0, startIndex, (name: string, val: any) => { // TODO: type
      this.dataColls[name] = [];
      this.dataColls[name].push(val);
    });
    startIndex += 5;
    while (startIndex > -1 && startIndex < content.length)
      startIndex = this.readNext(content, startIndex);
  }

  readNext(content: string, startIndex: number) {
    const nextStartIndex = content.indexOf('$$$$', startIndex);
    if (nextStartIndex === -1)
      return -1;
    else {
      this.parse(content, startIndex, nextStartIndex,
        (name: string, val: number) => this.dataColls[name].push(val));
    }

    if (nextStartIndex > -1)
      return nextStartIndex + 5;


    return nextStartIndex;
  }

  parse(content: string, start: number, end: number, handler: any) {
    const molEnd = +content.indexOf('M  END\n', start) + 7;
    let localEnd = start;
    this.dataColls['molecule'].push(content.substr(start, molEnd - start));

    start = molEnd;
    while (localEnd < end) {
      start = content.indexOf('> <', localEnd);
      if (start === -1)
        return;


      start += 3;
      localEnd = content.indexOf('>\n', start);
      if (localEnd === -1)
        return;


      const propertyName = content.substring(start, localEnd);
      start = localEnd + 2;

      localEnd = content.indexOf('\n', start);
      if (localEnd === -1)
        localEnd = end;
      else if (content[localEnd + 1] != '\n')
        localEnd = content.indexOf('\n', ++localEnd);
      ;

      handler(propertyName, content.substring(start, localEnd));
      localEnd += 2;
    }
  }
}
