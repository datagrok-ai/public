import * as DG from "datagrok-api/dg";
import {UaQueryViewer} from "../viewers/ua-query-viewer";
import {UaDataFrameViewer} from "../viewers/ua-data-frame-viewer";
import {UaFilter} from "../filter2";
import {BehaviorSubject} from "rxjs"

export class TopQueriesUsingDataSource extends UaDataFrameViewer {
  public constructor(dataSource: string, filterStream: BehaviorSubject<UaFilter>) {
    super(
        'Top Queries Using Data Source',
        'TopQueriesUsingDataSource',
        (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root,
        null as any,
        {'data_source': dataSource},
        filterStream.getValue()
    );
  }

}