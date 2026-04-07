import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Single Choice
//description: A filter that lets you select exactly one category
//output: filter result
//meta.role: filter
export function radioButtonFilter() {
  return PackageFunctions.radioButtonFilter();
}

//name: Multi Choice
//description: A filter that works with columns of multi-value cells (such as lists of identifiers)
//output: filter result
//meta.role: filter
export function multiValueFilter() {
  return PackageFunctions.multiValueFilter();
}

//name: TimeWidget
//description: Shows current time
//output: widget result
export function timeWidget() : any {
  return PackageFunctions.timeWidget();
}

//name: TableSummary
//output: viewer result
//meta.showInGallery: false
//meta.role: viewer
export function tableSummary() {
  return PackageFunctions.tableSummary();
}

//name: inputDemo
export function inputDemo() : void {
  PackageFunctions.inputDemo();
}

//output: object result
//meta.propertyType: string
//meta.semType: foo
//meta.role: valueEditor
export function fooInput() : any {
  return PackageFunctions.fooInput();
}
