class NglViewerPackage extends DG.Package {


  //tags: fileViewer, fileViewer-mol, fileViewer-mol2, fileViewer-cif, fileViewer-mcif, fileViewer-mmcif, fileViewer-gro, fileViewer-pdb, fileViewer-ent, fileViewer-pqr, fileViewer-mmtf, fileViewer-mtl, fileViewer-sdf, fileViewer-sd
  //input: file file
  //output: view v
  nglStructureViewer(file) {
    return this.nglViewer(file);
  }

  //tags: fileViewer, fileViewer-ply, fileViewer-obj
  //input: file file
  //output: view v
  nglSurfaceViewer(file) {
    return this.nglViewer(file);
  }

  //tags: fileViewer, fileViewer-prmtop, fileViewer-parm7, fileViewer-psf, fileViewer-top
  //input: file file
  //output: view v
  nglTopologyViewer(file) {
    return this.nglViewer(file);
  }

  //tags: fileViewer, fileViewer-dsn6, fileViewer-brix, fileViewer-cube, fileViewer-cub, fileViewer-dx, fileViewer-dxbin, fileViewer-xplor, fileViewer-cns, fileViewer-mrc, fileViewer-map, fileViewer-ccp4
  //input: file file
  //output: view v
  nglDensityViewer(file) {
    return this.nglViewer(file);
  }

  //name: PDB Viewer
  //tags: panel
  //input: string pdbId {semType: pdb_id}
  //output: widget w
  nglPdbPanelWidget(pdbId) {
    var host = ui.div([], 'd4-ngl-viewer');
    var stage = new NGL.Stage(host);
    NglViewerPackage.handleResize(host, stage);
    stage.loadFile(`rcsb://${pdbId}`, { defaultRepresentation: true });
    return DG.Widget.fromRoot(host);
  }

  static handleResize(host, stage) {
    let canvas = host.querySelector('canvas');
    function resize() {
      canvas.width = Math.floor(canvas.clientWidth * window.devicePixelRatio);
      canvas.height = Math.floor(canvas.clientHeight * window.devicePixelRatio);
      stage.handleResize();
    }

    ui.onSizeChanged(host).subscribe((_) => resize());
    resize();
  }


  /** @returns {DG.View} */
  nglViewer(file) {
    let view = DG.View.create();
    var host = ui.div([], 'd4-ngl-viewer');
    var stage = new NGL.Stage(host);
    NglViewerPackage.handleResize(host, stage);

    function loadBytes(bytes) {
      var blob = new Blob([bytes], {type: 'application/octet-binary'});
      stage.loadFile(blob, { defaultRepresentation: true, ext: file.extension} );
    }

    file
      .readAsBytes()
      .then(loadBytes);

    view.append(host);
    return view;
  }
}
