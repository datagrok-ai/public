import {after, before, category, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('Dapi: files', () => {
    let filePrefix = 'Demo:TestJobs:Files:DemoFiles/';
    let testTextFileName = 'js-api-testTextFile'
    let testTextFilePath = filePrefix + testTextFileName;

    before(async () => {
        await grok.dapi.files.writeAsText(testTextFilePath, 'testString');
    });

    test('Dapi: files - exists', async () => {
        if (!await grok.dapi.files.exists(testTextFilePath))
            throw "File doesn't exist";
    });

    test('Dapi: files - write/read text', async () => {
        let filePath = filePrefix + 'Dapi: files - write, read text.txt';
        let fileText = 'testString';

        try {
            await grok.dapi.files.writeAsText(filePath, fileText);
            if (fileText !== await grok.dapi.files.readAsText(filePath))
                throw "Сontent is wrong";
        }
        finally {
            await grok.dapi.files.delete(filePath);
        }
    });

    test('Dapi: files - write/read blob', async () => {
        let filePath = filePrefix + 'Dapi: files - write, read blob.txt';
        let content = [0, 1, 2, 3];

        try {
            await grok.dapi.files.write(filePath, content);
            if (content.toString() !== (await grok.dapi.files.readAsBytes(filePath)).toString())
                throw "Сontent is wrong";
        }
        finally {
            await grok.dapi.files.delete(filePath);
        }
    });

    test('Dapi: files - search', async () => {
        if((await grok.dapi.files.list(filePrefix, false, testTextFileName)).length !== 1)
            throw "Can't find the file";
    });

    // test('Dapi: files - move', async () => {
    //     let fileName = 'Dapi: files - move.txt';
    //     let filePath = filePrefix + fileName;
    //     let newFilePrefix = 'texts';
    //
    //     await grok.dapi.files.writeAsText(filePath, 'testString');
    //     await grok.dapi.files.move([filePath], newFilePrefix);
    //
    //     if(await grok.dapi.files.exists(filePath))
    //         throw "Old file exists";
    //
    //     if(!await grok.dapi.files.exists(filePrefix + newFilePrefix + '/' + fileName))
    //         throw "New file doesn't exist";
    //
    //     await grok.dapi.files.delete(filePath);
    // });

    test('Dapi: files - delete', async () => {
        let filePath = filePrefix + 'Dapi: files - delete.txt';

        try {
            await grok.dapi.files.writeAsText(filePath, 'testString');
            await grok.dapi.files.delete(filePath);

            if (await grok.dapi.files.exists(filePath))
                throw "File exists";
        }
        finally {
            if (await grok.dapi.files.exists(filePath))
                await grok.dapi.files.delete(filePath);
        }
    });

    after(async () => {
        await grok.dapi.files.delete(testTextFilePath);
    });

});
