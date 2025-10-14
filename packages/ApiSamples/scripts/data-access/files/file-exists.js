//tags: Files, exists
//See more examples in dapi/files

(async () => {
    // Using System:DemoFiles connection
    let filePath = 'System:DemoFiles/cars.csv';

    let exists = await grok.dapi.files.exists(filePath);
    grok.shell.info(`File ${filePath} ${exists ? 'exists' : 'doesn\'t exist'}`);
})();
