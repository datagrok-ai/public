//tags: Files, read
//See more examples in dapi/files

(async () => {
    // Using System:DemoFiles connection
    let filePath = 'System:DemoFiles/cars.csv';

    // Read file as bytes
    let bytes = await grok.dapi.files.readAsBytes(filePath);
    grok.shell.info(`Read file cars.csv as bytes: length ${bytes.length}`);

    // Read file as string. File content should be UTF-8 decodable or error will be thrown
    let text = await grok.dapi.files.readAsText(filePath);
    grok.shell.info(`Read file cars.csv as string: length ${text.length}`);

    // Read csv file and convert it to DG.DataFrame automatically. File content should be valid csv or error will be thrown
    let dataFrame = await grok.dapi.files.readCsv(filePath);
    grok.shell.info(`Read file cars.csv as DataFrame: rows ${dataFrame.rowCount}`);
})();
