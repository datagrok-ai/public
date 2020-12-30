class NglViewerPackage extends DG.Package {

  //tags: fileViewer, fileViewer-sdf
  //input: file file
  //output: view v
  sdfViewer(file) {
    return this.nglViewer(file);
  }

  //tags: fileViewer, fileViewer-cif
  //input: file file
  //output: view v
  cifViewer(file) {
    return this.nglViewer(file);
  }

  //tags: fileViewer, fileViewer-mol, fileViewer-mol2, fileViewer-mmcif
  //input: file file
  //output: view v
  molViewer(file) {
    return this.nglViewer(file);
  }

  /** @returns {DG.View} */
  nglViewer(file) {
    let view = DG.View.create();
    var host = ui.div([], 'd4-ngl-viewer');
    var stage = new NGL.Stage(host);
    let canvas = host.querySelector('canvas');

    function init() {
      function resize() {
        canvas.width = Math.floor(canvas.clientWidth * window.devicePixelRatio);
        canvas.height = Math.floor(canvas.clientHeight * window.devicePixelRatio);
        stage.handleResize();
      }

      ui.onSizeChanged(host).subscribe((_) => resize());
      resize();

      function loadBytes(bytes) {
        var blob = new Blob([bytes], {type: 'application/octet-binary'});
        stage.loadFile(blob, { ext: 'cif'} );
      }

      // function loadString(bytes) {
      //   var blob = new Blob([bytes], {type: 'text/plain'});
      //   stage.loadFile(blob, { defaultRepresentation: true, ext: file.extension} );
      // }

      file
        .readAsBytes()
        .then(loadBytes);
    }

    setTimeout(init, 200);
    view.append(host);
    return view;
  }
}
