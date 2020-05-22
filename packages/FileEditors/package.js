class FileEditorsPackage extends DG.Package {

    //tags: fileViewer, fileViewer-pdf
    //input: string url {semType: url}
    //output: view v
    viewPdf(url) {
        let view = grok.View.create();
        let canvas = ui.canvas();
        canvas.style.width = '100%';
        canvas.style.height = '100%';
        view.append(canvas);

        pdfjsLib.getDocument(url)
            .then((pdf) => pdf.getPage(1))
            .then((page) => {
                let scale = 1.5;
                let viewport = page.getViewport(scale);
                let context = canvas.getContext('2d');
                canvas.height = viewport.height;
                canvas.width = viewport.width;
                page.render({
                    canvasContext: context,
                    viewport: viewport
                });
            });

        return view;
    }

    //tags: fileViewer, fileViewer-pdf
    //input: string url {semType: url}
    testPdf(url) {
        return url;
    }
}
