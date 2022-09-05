import {AnalysisLoader} from "../service/nx/analysis/AnalysisLoader";

export class MetadataViewer extends DG.JsViewer {

    onTableAttached()
    {

           /**
         *
         * Updates the content of the viewer based on currently selected cell in the active dataframe.
         * @param viewer a reference to this viewer.
         * @param eCanvasCpd a reference to the DOM canvas element where the structure should be displayed.
         * @param eCanvasAssay a reference to the DOM canvas element where the assay data should be displayed.
         * @param eH2Assay a reference to an H2 element where the assay colum name should be displayed.
         */
        function render(viewer, eH2NameValue, eH2ComposeTypeValue, eHPrimiTypeValue, eHPrimiTypeValueClassValue, eH2DataFrameColumnTypeValue, eH2RawJSONValue)
        {
            let eRoot = viewer.root;

            let dframe = viewer.dataFrame;

            let cell = dframe.currentCell;
            if(cell.rowIndex < 0)
                return;

            let column = cell.column;
            let descriptor = column.getTag(AnalysisLoader.TAG_COMPPOSE_DESCRIPTOR);
            if(descriptor !== null && descriptor !== undefined)
            {


                if(descriptor.getPrimitiveType() === null || descriptor.getPrimitiveType().getParentType() === undefined || descriptor.getPrimitiveType().getParentType() === null)
                {
                    let yuiyi = 0;
                }
                else {
                    eH2NameValue.innerHTML = descriptor.getTitle();
                    eH2ComposeTypeValue.innerHTML = descriptor.getPrimitiveType().getParentType().getName();
                    eHPrimiTypeValue.innerHTML = descriptor.getPrimitiveType().getKey();
                    eHPrimiTypeValueClassValue.innerHTML = descriptor.getPrimitiveType().getType();

                }
                eH2RawJSONValue.innerHTML = JSON.stringify(descriptor.getRawJSON(), null, 2);
            }
            else
            {
                eH2NameValue.innerHTML = "N/A";
                eH2ComposeTypeValue.innerHTML = "N/A";
                eHPrimiTypeValue.innerHTML = "N/A";
                eHPrimiTypeValueClassValue.innerHTML = "N/A"
                eH2RawJSONValue.innerHTML = "N/A";
           }

           eH2DataFrameColumnTypeValue.innerHTML = column.type;
        }

        /**
         * Called by the active dataframe when the active cells changes.
         * @param viewer viewer a reference to this viewer.
         * @param eCanvasCpd a reference to the DOM canvas element where the structure should be displayed.
         * @param eCanvasAssay a reference to the DOM canvas element where the assay data should be displayed.
         * @param eH2Assay a reference to an H2 element where the assay colum name should be displayed.
         */
        function onSelectionChanged(viewer, eH2NameValue, eH2ComposeTypeValue, eHPrimiTypeValue, eHPrimiTypeValueClassValue, eH2DataFrameColumnTypeValue, eH2RawJSONValue) {
            render(viewer, eH2NameValue, eH2ComposeTypeValue, eHPrimiTypeValue, eHPrimiTypeValueClassValue, eH2DataFrameColumnTypeValue, eH2RawJSONValue);
        }

        let viewer = this;
        let eRoot = viewer.root;

        let nW = eRoot.clientWidth;


        let eH2Name = ui.divText("");
        eH2Name.innerHTML = "<b>Name</b>";
        let eH2NameValue = ui.divText("");
        let eH2NameEmpty = ui.divText("");

        let eH2ComposeTypeName = ui.divText("");
        eH2ComposeTypeName.innerHTML = "<b>Compose Type:</b>";
        let eH2ComposeTypeValue = ui.divText("");
        let eH2ComposeTypeEmpty = ui.divText("");

        let eH2PrimiTypeName = ui.divText();
        eH2PrimiTypeName.innerHTML = "<b>Primitive Type:</b>";
        let eHPrimiTypeValue = ui.divText("");
        let eH2PrimiTypeEmpty = ui.divText("");

        let eH2PrimiTypeValueClassName = ui.divText("");
        eH2PrimiTypeValueClassName.innerHTML = "<b>Primitive Type Value Class:</b>";
        let eHPrimiTypeValueClassValue = ui.divText("");
        let eHPrimiTypeValueClassEmpty = ui.divText("");

        let eH2DataFrameColumnTypeName = ui.divText("");
        eH2DataFrameColumnTypeName.innerHTML = "<b>Dataframe Column Type:</b>";
        let eH2DataFrameColumnTypeValue = ui.divText("");
        let eH2DataFrameColumnTypeEmpty = ui.divText("");

        let eH2RawJSONName = ui.divText("");
        eH2RawJSONName.innerHTML = "<b>Raw JSON:</b>";
        let eH2RawJSONValue = ui.divText("");
        let eH2RawJSONEmpty = ui.divText("");



        let eDivPanel = ui.splitV([eH2Name, eH2NameValue, eH2NameEmpty, eH2ComposeTypeName, eH2ComposeTypeValue, eH2ComposeTypeEmpty, eH2PrimiTypeName, eHPrimiTypeValue, eH2PrimiTypeEmpty,
            eH2PrimiTypeValueClassName, eHPrimiTypeValueClassValue, eHPrimiTypeValueClassEmpty,
            eH2DataFrameColumnTypeName, eH2DataFrameColumnTypeValue, eH2DataFrameColumnTypeEmpty,
            eH2RawJSONName, eH2RawJSONValue, eH2RawJSONEmpty]);
        eRoot.appendChild(eDivPanel);
        eH2RawJSONValue.style.height = "200px";
        eH2RawJSONValue.parentNode.parentNode.parentNode.style.height = "600px";
        eH2RawJSONValue.parentNode.parentNode.style.height = "600px";

        let pSel = this.dataFrame.onCurrentCellChanged;
        pSel = pSel.subscribe((_) => onSelectionChanged(viewer, eH2NameValue, eH2ComposeTypeValue, eHPrimiTypeValue, eHPrimiTypeValueClassValue, eH2DataFrameColumnTypeValue, eH2RawJSONValue));

        render(viewer, eH2NameValue, eH2ComposeTypeValue, eHPrimiTypeValue, eHPrimiTypeValueClassValue, eH2DataFrameColumnTypeValue, eH2RawJSONValue);
    }
}