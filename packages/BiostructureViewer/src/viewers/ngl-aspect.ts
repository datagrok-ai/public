import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as ngl from 'NGL';

const getCircularReplacer = (): (this: any, key: string, value: any) => any => {
  const seen = new WeakSet();
  return (_key: string, value: any) => {
    if (typeof value === 'object' && value !== null) {
      if (seen.has(value)) {
        //
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
    const colBackground = this.colorScheme['col_background'];
    // const colHelix = this.colorScheme[''];

    nglHost.style.backgroundColor = colBackground;
    view.box = true;
    this.entry = entry;
    this.repChoice = repChoice;
    this.selection = twinSelections;

    this.stage = new ngl.Stage(nglHost);
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
    this.stage.viewer.render = () => {
      //@ts-ignore
      originalRender.bind(this.stage.viewer)();
      console.debug('ABABABABABABABABABABAB');
      //Getting JSON causes circular reference exception
      //console.debug('JSON stage: ' + JSON.stringify(this.stage, getCircularReplacer()));
    };

    // let listener = () => {
    //   let va = this.stage;
    //   console.log(va.viewer.renderer.domElement.toDataURL('image/png'));
    // };
    // this.stage.viewer.signals.rendered.add(listener);

    await this.render(true, ligandSelection);
    await new Promise((r) => setTimeout(r, 4000));
    await this.nglResize(nglHost);
    // let va = this.stage;
    // console.log(va.viewer.renderer.domElement.toDataURL('image/png'));
  }

  public async render(reload: boolean, ligandSelection: { [key: string]: boolean }) {
    const switchObj = this.selection;
    this.ligandSelection = ligandSelection;

    // const colChain = this.colorScheme['col_chain'];
    // const colHelix = this.colorScheme['col_helix'];
    const colHighlight = this.colorScheme['col_highlight'];

    //highlights in NGL
    const schemeBuffer: string[][] = [];

    Object.keys(switchObj).forEach((keyChain) => {
      Object.keys(switchObj[keyChain]).forEach((keyFtStart) => {
        const index = parseInt(keyFtStart) + 1;
        if (switchObj[keyChain][keyFtStart]['state'] === true)
          schemeBuffer.push([colHighlight, `${index} and :${keyChain}`]);
      });
    });

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
    schemeBuffer.push(['#0069a7', `* and :A`]);
    schemeBuffer.push(['#f1532b', `* and :B`]);
    schemeBuffer.push(['green', `* and :R`]);
    schemeBuffer.push(['green', `* and :M`]);
    //@ts-ignore
    const schemeId = ngl.ColormakerRegistry.addSelectionScheme(scheme_buffer);

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
            const index = parseInt(keyFtStart) + 1;

            if (typeof switchObj[keyChain][keyFtStart]['representation'] === 'undefined') {
              const selection = `${index} and :${keyChain} and (not backbone or .CA or (PRO and .N))`;
              const addedRepresentation = this.stage.compList[0].addRepresentation('ball+stick', {sele: selection});
              switchObj[keyChain][keyFtStart]['representation'] = addedRepresentation;
            }

            if (switchObj[keyChain][keyFtStart]['state'] !== true)
              switchObj[keyChain][keyFtStart]['representation'].setVisibility(false);
            else
              switchObj[keyChain][keyFtStart]['representation'].setVisibility(true);
          });
        });
      }

      Object.keys(this.ligandSelection).forEach((ligand) => {
        if (this.ligandSelection[ligand][0] === true) {
          if (!Object.keys(this.ligandRepresentations).includes(ligand)) {
            this.ligandRepresentations[ligand] =
              this.stage.compList[0].addRepresentation('licorice',
                {sele: `${this.ligandSelection[ligand][1]} and :${ligand}`});
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
    const stringBlob = new Blob([pdbStr], {type: 'text/plain'});
    await this.stage.loadFile(stringBlob, {ext: 'pdb'});
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

    const va = this.stage;
    console.log(va.viewer.renderer.domElement.toDataURL('image/png'));
  }

  // viewer resize
  private async _resize(host: HTMLElement) {
    const canvas = host.querySelector('canvas');
    canvas!.width = Math.floor(host.clientWidth * window.devicePixelRatio);
    canvas!.height = Math.floor(host.clientWidth * window.devicePixelRatio);
    await this.stage.handleResize();
    // console.log(this.stage.viewer.renderer.domElement.toDataURL('image/png'));
  }

  private async nglResize(host: HTMLElement) {
    ui.onSizeChanged(host).subscribe((_) => {
      this._resize(host);
    });
    await this._resize(host);
  }
}
