import * as DG from "datagrok-api/dg";

export function  isExceptionElement(action: string):void{
    let exceptionElement = document.getElementsByClassName('grok-global-exception');
    if (exceptionElement.length != 0)
        throw 'exception was thrown after action: ' + action;
}

export function isColumnPresent(columns:DG.ColumnList, columnName:string): void {
    let check = false;
    if (columns.byName(columnName) == null)
        check = true
    if (check == true){
        throw 'column ' + columnName + ' not found';
    }
}

export function isViewerPresent(viewers:DG.Viewer[], viewerName:string): void {
    let check = false;
    for(let i:number = 0; i < viewers.length; i++){
        if (viewers[i].type == viewerName){
            check = true;
            break;
        }
    }
    if (check == false){
        throw viewerName + 'not found';
    }
}

export function  isErrorBallon(ballonText: string):void{
    let exceptionElement = document.getElementsByClassName('d4-balloon-content')[0] as HTMLElement;
    if (exceptionElement == undefined)
        throw 'Error ballon didn\'t appear where it should have.';
    if (exceptionElement.textContent != ballonText)
        throw 'Wrong error message on balloon';
    exceptionElement.click();
}