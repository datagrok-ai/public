class NlpPackageDetectors extends DG.Package {
    // `isTextFile` appears in JavaScript panel conditions
    // `supportedExt` is used in Python scripts

    //input: file file
    //output: bool result
    isTextFile(file) {
        let extensions = ["doc", "docx", "odt", "pdf", "rtf", "txt"];
        return extensions.includes(file.extension);
    }

    //input: string filename
    //output: bool result
    supportedExt(filename) {
        let extensions = ["doc", "docx", "odt", "pdf", "rtf", "txt"];
        return extensions.some(ext => filename.endsWith(ext));
    }

    //input: string filename
    //output: bool result
    supportedExtExtract(filename) {
        let extensions = ["csv", "doc", "docx", "eml", "epub", "gif", "htm", "html",
                        "jpeg", "jpg", "json", "log", "mp3", "msg", "odt", "ogg",
                        "pdf", "png", "pptx", "ps", "psv", "rtf", "tff", "tif",
                        "tiff", "tsv", "txt", "wav", "xls", "xlsx"];
        return extensions.some(ext => filename.endsWith(ext));
    }
}
