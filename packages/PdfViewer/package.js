class PdfViewerPackage extends DG.Package {

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
}
