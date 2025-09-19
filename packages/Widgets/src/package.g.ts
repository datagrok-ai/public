import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Single Choice
//description: A filter that lets you select exactly one category
//tags: filter
//output: filter result
export function radioButtonFilter() {
  return PackageFunctions.radioButtonFilter();
}

//name: Multi Choice
//description: A filter that works with columns of multi-value cells (such as lists of identifiers)
//tags: filter
//output: filter result
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
//tags: viewer
//output: viewer result
export function tableSummary() {
  return PackageFunctions.tableSummary();
}

//name: inputDemo
export function inputDemo() : void {
  PackageFunctions.inputDemo();
}

//tags: valueEditor
//output: object result
//meta.propertyType: string
//meta.semType: foo
export function fooInput() : any {
  return PackageFunctions.fooInput();
}
