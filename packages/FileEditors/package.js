class FileEditorsPackage extends DG.Package {

  //tags: fileViewer, fileViewer-pdf
  //input: file file
  //output: view v
  viewPdf(file) {
    let view = DG.View.create();
    let canvas = ui.canvas();
    view.append(canvas);

    file.readAsBytes()
      .then((pdf) => {
        pdfjsLib.getDocument( { data: pdf })
          .then((pdf) => pdf.getPage(1))
          .then((page) => {
            let scale = 1.0;
            let viewport = page.getViewport(canvas.width / page.getViewport(scale).width);
            let context = canvas.getContext('2d');

            canvas.style.width =  viewport.width + 'px';
            canvas.style.height = viewport.height + 'px';
            canvas.height = viewport.height;
            canvas.width = viewport.width;

            page.render({
              canvasContext: context,
              viewport: viewport
            });
          });
    });

    return view;
  }

  /** @returns {DG.View} */
  nglViewer(file) {
    let view = DG.View.create();
    var host = ui.div([], 'd4-ngl-viewer');
    var stage = new NGL.Stage(host);
    let canvas = host.querySelector('canvas');

    function init() {
      canvas.width = 300;
      canvas.height = 300;
      canvas.style.height = '300px';
      canvas.style.width = '300px';
      stage.handleResize();

      function loadString(bytes) {
        var blob = new Blob([bytes], {type: 'text/plain'});
        stage.loadFile(blob, { ext: file.extension} );
      }

      file
        .readAsString()
        .then(loadString);
    }

    setTimeout(init, 200);
    view.append(host);
    return view;
  }

  //tags: fileViewer, fileViewer-sdf
  //input: file file
  //output: view v
  viewSdf(file) {
    return this.nglViewer(file);
  }

  //tags: fileViewer, fileViewer-cif
  //input: file file
  //output: view v
  viewCif(file) {
    let view = DG.View.create();

    //view.append(ui.bigButton('foo', () => grok.shell.info('Foo!')));
    var host = ui.div([], 'd4-ngl-viewer');
    var stage = new NGL.Stage(host);
    let canvas = host.querySelector('canvas');

    function init() {
      canvas.width = 300;
      canvas.height = 300;
      canvas.style.height = '300px';
      canvas.style.width = '300px';
      stage.handleResize();

      function loadBytes(bytes) {
        var blob = new Blob([bytes], {type: 'application/octet-binary'});
        stage.loadFile(blob, { ext: 'cif'} );
      }

      function loadString(bytes) {
        var blob = new Blob([bytes], {type: 'text/plain'});
        stage.loadFile(blob, { ext: 'cif'} );
      }

      //stage.loadFile("http://files.rcsb.org/download/5IOS.cif", {defaultRepresentation: true});
      file
        .readAsBytes()
        .then(loadBytes);

      //stage.loadFile(url, {defaultRepresentation: true});
    }

    setTimeout(init, 200);

    //stage.loadFile("rcsb://1crn", {defaultRepresentation: true});
    //stage.loadFile(url, {defaultRepresentation: true});

    view.append(host);
    return view;
  }
}
