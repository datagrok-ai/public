class FileEditorsPackage extends GrokPackage {

    //tags: file-viewer, file-viewer-pdf
    //input: string url {semType: url}
    //output: view v
    viewPdf(url) {
        let view = View.create();
        let canvas = ui.canvas();
        canvas.style.width = '100%';
        canvas.style.height = '100%';
        view.append(canvas);

        pdfjsLib.getDocument(url)
            .then(function(pdf) {
                return pdf.getPage(1);
            })
            .then(function(page) {
                // Set scale (zoom) level
                var scale = 1.5;

                // Get viewport (dimensions)
                var viewport = page.getViewport(scale);

                // Fetch canvas' 2d context
                var context = canvas.getContext('2d');

                // Set dimensions to Canvas
                canvas.height = viewport.height;
                canvas.width = viewport.width;

                // Prepare object needed by render method
                var renderContext = {
                    canvasContext: context,
                    viewport: viewport
                };

                // Render PDF page
                page.render(renderContext);
            });


        return view;
    }

    //tags: file-viewer, file-viewer-pdf
    //input: string url {semType: url}
    testPdf(url) {
        return url;
    }
}
