import * as DG from "datagrok-api/dg";
import {UaQueryViewer} from "../viewers/abstract/ua-query-viewer";
import {UaDataFrameQueryViewer} from "../viewers/ua-data-frame-query-viewer";
import {UaFilter} from "../filter2";
import {BehaviorSubject} from "rxjs"

export class TopQueriesUsingDataSource extends UaDataFrameQueryViewer {
  public constructor(dataSource: string, filterStream: BehaviorSubject<UaFilter>) {
    super(
        'Queries Using Data Source',
        'TopQueriesUsingDataSource',
        (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root,
        null as any,
        {'data_source': dataSource},
        filterStream.getValue()
    );
  }

}