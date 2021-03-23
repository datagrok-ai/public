//name: SaveLoadRowsViewer
//language: javascript


class SaveLoadRows extends DG.JsViewer {

    constructor() {
        super();
        this.uniqueId = null;
        this.fileToSave = ui.stringInput('FileName', 'filename');
        this.savedFilesList = null;
        this.files = null;
        // put here connection to folder, where saveRows json files will be stored
        this.connection = 'Demo:TestJobs:Files:DemoFiles/';
        this.init();

    }

    async saveSelectedRows() {

        let indexes = this.dataFrame
            .groupBy([`${this.uniqueId.value}`])
            .whereRowMask(this.dataFrame.selection)
            .aggregate();

        let data = `${this.dataFrame.toString()}\n${this.uniqueId.stringValue}\n${indexes.col(0).toList()}`;
        await grok.dapi.files.writeAsText(`${this.connection}${this.fileToSave.value}.txt`, data);
    }



    async loadSelectedRows() {

        let res = await grok.dapi.files.readAsText(`${this.connection}${this.savedFilesList.value}`);
        res = res.split("\n");
        let uniqueColumnName = res[1];
        let values = JSON.parse("[" + res[2] + "]");
        values = values.map((e) => parseInt(e));
        this.dataFrame.rows.select((row) => values.includes(row[`${uniqueColumnName}`]));

    }


    saveDialog = () => {
        ui.dialog('Save rows to file')
            .add(this.uniqueId).add(this.fileToSave)
            .onOK(() => this.saveSelectedRows()).show();
    };

    loadDialog = async () => {
        this.files = await grok.dapi.files.list(this.connection, false, '');
        this.files = this.files.map((e) => e.path);
        this.savedFilesList = await ui.choiceInput('Saved Rows', ' ', this.files)

        ui.dialog('Load rows from file')
            .add(this.savedFilesList)
            .onOK(() => this.loadSelectedRows()).show();
    };

    async init() {

        this.files = await grok.dapi.files.list(this.connection, false, '');
        this.files = this.files.map((e) => e.path);
        this.savedFilesList = await ui.choiceInput('Saved Rows', ' ', this.files)
        // please define here you primary key column
        this.uniqueId = ui.columnInput('Unique id column', this.dataFrame, this.dataFrame.col('subj'))



        // initial UI
        //Buttons yo have to integrate yo your UI

        let saveRowsButton = ui.button('Save rows');
        saveRowsButton.addEventListener("click", this.saveDialog);

        let loadRowsButton = ui.button('Load rows')
        loadRowsButton.addEventListener("click", this.loadDialog);


        this.root.appendChild(ui.div(ui.inputs([saveRowsButton, loadRowsButton])));

    }

}


// uncomment for testing
//grok.shell.registerViewer('SaveLoadRows', 'JavaScript-based viewer', () => new SaveLoadRows());

//demog = grok.data.demo.demog();
//view = grok.shell.addTableView(demog);
//hist = view.addViewer('SaveLoadRows');