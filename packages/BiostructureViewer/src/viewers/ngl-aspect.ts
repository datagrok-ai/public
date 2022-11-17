import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PdbEntry} from '../pdb-entry.js';

// export let _package = new DG.Package();

const getCircularReplacer = () => {
  const seen = new WeakSet();
  return (key, value) => {
    if (typeof value === 'object' && value !== null) {
      if (seen.has(value)) {
        return;
      }
      seen.add(value);
    }
    return value;
  };
};

export class NglAspect {
  stage: any;
  entry: string;
  colorScheme: any;
  nglHost: HTMLElement;
  repChoice: DG.InputBase;
  selection: any;
  ligandSelection: { [key: string]: any };
  ligandRepresentations: { [key: string]: any } = {};

  public async init(view: DG.TableView, entry: string, colorScheme: {}, nglHost: HTMLElement,
    repChoice: DG.InputBase, twinSelections: {}, ligandSelection: { [key: string]: boolean }) {
    this.colorScheme = colorScheme;
    const col_background = this.colorScheme['col_background'];
    const col_helix = this.colorScheme[''];

    nglHost.style.backgroundColor = col_background;
    view.box = true;
    this.entry = entry;
    this.repChoice = repChoice;
    this.selection = twinSelections;
    //@ts-ignore
    this.stage = new NGL.Stage(nglHost);
    // let originalRender = this.stage.viewer.renderer.render;
    // this.stage.viewer.renderer.render = function(scene: any, camera: any) {
    //     //@ts-ignore
    //     let va = this.stage;
    //     console.log(va.viewer.renderer.domElement.toDataURL('image/png'));
    //   //@ts-ignore
    //   originalRender.bind(this.stage.viewer.renderer)(scene, camera);
    //   //@ts-ignore
    //   va = this.stage;
    //   console.log(va.viewer.renderer.domElement.toDataURL('image/png'));
    // }
    // .bind(this);

    const originalRender = this.stage.viewer.render;
    this.stage.viewer.render = function() {
      //@ts-ignore
      originalRender.bind(this.stage.viewer)();
      console.log('ABABABABABABABABABABAB');
      console.debug('JSON stage: ' + JSON.stringify(this.stage, getCircularReplacer()));
    }
      .bind(this);

    // let listener = () => {
    //   let va = this.stage;
    //   console.log(va.viewer.renderer.domElement.toDataURL('image/png'));
    // };
    // this.stage.viewer.signals.rendered.add(listener);

    await this.render(true, ligandSelection);
    await new Promise((r) => setTimeout(r, 4000));
    await this.nglResize(nglHost);
    let va = this.stage;
    console.log(va.viewer.renderer.domElement.toDataURL('image/png'));
  }

  public async render(reload: boolean, ligandSelection: { [key: string]: boolean }) {
    let switchObj = this.selection;
    this.ligandSelection = ligandSelection;

    let col_chain = this.colorScheme['col_chain'];
    let col_helix = this.colorScheme['col_helix'];
    let col_highlight = this.colorScheme['col_highlight'];

    //highlights in NGL
    let scheme_buffer: string[][] = [];

    Object.keys(switchObj).forEach((keyChain) => {
      Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {

        let index = parseInt(keyFtStart) + 1;
        if (switchObj[keyChain][keyFtStart]['state'] === true) {
          scheme_buffer.push([col_highlight, `${index} and :${keyChain}`]);
        }
      });
    });

    let schemeId;

    // this.entry.entities[0].chains.forEach(chain => {
    //   let str_buffer = '';
    //   Object.keys(chain.tracks).forEach(track => {
    //     if (track === "HELIX_P") {
    //       let index1 = chain.tracks[track][0];
    //       let index2 = chain.tracks[track].slice(-1);
    //       str_buffer = str_buffer + ` or ${index1}-${index2} and :${chain.id}`;
    //     }
    //   });

    //   scheme_buffer.push([col_helix, str_buffer]);
    //   scheme_buffer.push([col_chain, `* and :${chain.id}`]);


    // });

    //TODO: remove - demo purpose only
    scheme_buffer.push(['#0069a7', `* and :A`]);
    scheme_buffer.push(['#f1532b', `* and :B`]);
    scheme_buffer.push(['green', `* and :R`]);
    scheme_buffer.push(['green', `* and :M`]);
    //@ts-ignore
    schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);


    if (reload) {
      this.stage.removeAllComponents();
      await this.loadPdb(this.entry, this.repChoice, {color: schemeId});

      // //recovering ball and stick residual at cartoon view
      // if (this.repChoice.value === "cartoon") {
      //   Object.keys(switchObj).forEach((keyChain) => {
      //     Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {

      //       let index = keyFtStart;

      //       let selection = `${index} and :${keyChain} and (not backbone or .CA or (PRO and .N))`;
      //       let addedRepresentation = this.stage.compList[0].addRepresentation("ball+stick", { sele: selection });
      //       switchObj[keyChain][keyFtStart]['representation'] = addedRepresentation;

      //       if (switchObj[keyChain][keyFtStart]['state'] !== true) {
      //         addedRepresentation.setVisibility(false);
      //       }
      //     });
      //   });
      // }

    } else {
      this.stage.compList[0].addRepresentation(this.repChoice.value, {color: schemeId});

      if (this.repChoice.value === 'cartoon') {
        Object.keys(switchObj).forEach((keyChain) => {
          Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {

            let index = parseInt(keyFtStart) + 1;

            if (typeof switchObj[keyChain][keyFtStart]['representation'] === 'undefined') {
              let selection = `${index} and :${keyChain} and (not backbone or .CA or (PRO and .N))`;
              let addedRepresentation = this.stage.compList[0].addRepresentation('ball+stick', {sele: selection});
              switchObj[keyChain][keyFtStart]['representation'] = addedRepresentation;
            }

            if (switchObj[keyChain][keyFtStart]['state'] !== true) {
              switchObj[keyChain][keyFtStart]['representation'].setVisibility(false);
            } else {
              switchObj[keyChain][keyFtStart]['representation'].setVisibility(true);
            }
          });
        });
      }

      Object.keys(this.ligandSelection).forEach(ligand => {
        if (this.ligandSelection[ligand][0] === true) {
          if (!Object.keys(this.ligandRepresentations).includes(ligand)) {
            this.ligandRepresentations[ligand] =
              this.stage.compList[0].addRepresentation('licorice', {sele: `${this.ligandSelection[ligand][1]} and :${ligand}`});
          }
          this.ligandRepresentations[ligand].setVisibility(true);
        } else if (Object.keys(this.ligandRepresentations).includes(ligand)) {
          this.ligandRepresentations[ligand].setVisibility(false);
        }
      });
      // let addedRepresentation = this.stage.compList[0].addRepresentation("licorice", { sele: '400 and :R' });
      // let addedRepresentation1 = this.stage.compList[0].addRepresentation("licorice", { sele: '401 and :M' });

    }
  }

  // load the 3D model
  private async loadPdb(pdbStr: string, repChoice: DG.InputBase, schemeObj: any) {
    var stringBlob = new Blob([pdbStr], {type: 'text/plain'});
    let a = await this.stage.loadFile(stringBlob, {ext: 'pdb'});
    // .then((o: any) => {
    // //await this.stage.loadFile(pdbStr).then(function (o: any) {
    //   let vv = this.stage;
    //   console.log(vv.viewer.renderer.domElement.toDataURL('image/png'));
    //   o.addRepresentation(repChoice.value, schemeObj);
    //   vv = this.stage;
    //   console.log(vv.viewer.renderer.domElement.toDataURL('image/png'));
    //   o.autoView();
    //   vv = this.stage;
    //   console.log(vv.viewer.renderer.domElement.toDataURL('image/png'));
    // });

    await this.stage.compList[0].addRepresentation(repChoice.value, schemeObj);
    await this.stage.compList[0].autoView();

    let va = this.stage;
    console.log(va.viewer.renderer.domElement.toDataURL('image/png'));
  }

  // viewer resize
  private async _resize(host: HTMLElement) {
    let canvas = host.querySelector('canvas');
    canvas!.width = Math.floor(host.clientWidth * window.devicePixelRatio);
    canvas!.height = Math.floor(host.clientWidth * window.devicePixelRatio);
    await this.stage.handleResize();
    console.log(this.stage.viewer.renderer.domElement.toDataURL('image/png'));
  }

  private async nglResize(host: HTMLElement) {
    ui.onSizeChanged(host).subscribe((_) => {
      this._resize(host);
    });
    await this._resize(host);
  }
}