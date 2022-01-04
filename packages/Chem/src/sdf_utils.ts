class SDFReader {
  dataColls: { [_: string]: any };

  constructor() {
    this.dataColls = {'molecule': []};
  }

  get_colls(content: string) {
    this.read(content);
    return this.dataColls;
  }

  read(content: string) {
    let startIndex = content.indexOf('$$$$', 0);
    this.parse(content, 0, startIndex, (name: string, val: any) => { // TODO: type
      this.dataColls[name] = [];
      this.dataColls[name].push(val);
    });
    startIndex += 5;
    while (startIndex > -1 && startIndex < content.length) {
      startIndex = this.readNext(content, startIndex);
    }
  }

  readNext(content: string, startIndex: number) {
    const nextStartIndex = content.indexOf('$$$$', startIndex);
    if (nextStartIndex === -1) {
      return -1;
    } else {
      this.parse(content, startIndex, nextStartIndex,
        (name: string, val: number) => this.dataColls[name].push(val));
    }

    if (nextStartIndex > -1) {
      return nextStartIndex + 5;
    }
    return nextStartIndex;
  }

  parse(content: string, start: number, end: number, handler: any) {
    const molEnd = +content.indexOf('M  END\n', start) + 7;
    let localEnd = start;
    this.dataColls['molecule'].push(content.substr(start, molEnd - start));

    start = molEnd;
    while (localEnd < end) {
      start = content.indexOf('> <', localEnd);
      if (start === -1) {
        return;
      }
      start += 3;
      localEnd = content.indexOf('>\n', start);
      if (localEnd === -1) {
        return;
      }
      const propertyName = content.substr(start, localEnd - start);
      start = localEnd + 2;

      localEnd = content.indexOf('\n', start);
      if (localEnd === -1) {
        localEnd = end;
      }
      handler(propertyName, content.substr(start, localEnd - start));
      localEnd += 2;
    }
  }
}
