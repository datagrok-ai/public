import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectChoices, expectNoThrow, expectPropAndLook, expectRoundTrip,
  expectRoundTripPropAndLook, subscribeAll, withAttachedViewer} from '../helpers';

// Calendar JS surface: typed factory, dateColumn getter, onCalendarClicked, look round-trip, onClick choices.
category('AI: Viewers: Calendar JS API', () => {
  test('factory DG.Viewer.calendar returns typed DG.CalendarViewer; addViewer round-trips type', async () => {
    const t = demog();
    const v = DG.Viewer.calendar(t, {dateColumnName: 'started'});
    expect(v instanceof DG.CalendarViewer, true);
    expect(v.dataFrame === t, true);
    expect(v.type, DG.VIEWER.CALENDAR);
    expectPropAndLook(v, {dateColumnName: 'started'});
  });

  test('dateColumnName round-trip + date friendly-key alias collapses to dateColumnName', async () => {
    const v1 = DG.Viewer.calendar(demog(), {dateColumnName: 'started'});
    expectRoundTripPropAndLook(v1, {dateColumnName: 'started'});

    const v2 = DG.Viewer.calendar(demog(), {date: 'started'});
    expectPropAndLook(v2, {dateColumnName: 'started'});
  });

  test('visual look round-trip — showHeader, redWeekends, showFilteredOnly, month/back colors', async () => {
    const v = DG.Viewer.calendar(demog(), {dateColumnName: 'started'});
    expectRoundTripPropAndLook(v, {
      showHeader: false,
      redWeekends: false,
      showFilteredOnly: false,
      oddMonthColor: 0xFFEEEEEE,
      evenMonthColor: 0xFFDDDDDD,
      backColor: 0xFFFFFFFF,
    });
  });

  test('onClick exposes RowGroupAction choices and round-trips Filter', async () => {
    const v = DG.Viewer.calendar(demog(), {dateColumnName: 'started'});
    expectChoices(v, 'onClick', ['Select', 'Filter', 'None']);
    expectRoundTrip(v, {onClick: 'Filter'});
  });

  test('onCalendarClicked is an rxjs Observable', async () => {
    const v = DG.Viewer.calendar(demog(20), {dateColumnName: 'started'});
    subscribeAll([v.onCalendarClicked])();
  });

  test('dateColumn returns the auto-detected DateTime column, or null when absent', async () => {
    const v = DG.Viewer.calendar(demog(), {dateColumnName: 'started'});
    expect(v.dateColumn != null, true);
    expect(v.dateColumn!.name, 'started');
    expect(v.dateColumn!.type, DG.COLUMN_TYPE.DATE_TIME);

    const noDates = df([['name', DG.COLUMN_TYPE.STRING, ['a', 'b', 'c']]]);
    const v2 = DG.Viewer.calendar(noDates);
    expect(v2.dateColumn, null);
  });

  test('view.addViewer attaches a typed DG.CalendarViewer; getOptions stays sane without a DateTime col', async () => {
    await withAttachedViewer<DG.CalendarViewer>(demog(), DG.VIEWER.CALENDAR,
      {dateColumnName: 'started'}, (v, tv) => {
        expect(v instanceof DG.CalendarViewer, true);
        const found = tv.viewers.find((x) => x.type === DG.VIEWER.CALENDAR);
        expect(found instanceof DG.CalendarViewer, true);
      });

    const noDates = df([['name', DG.COLUMN_TYPE.STRING, ['a']]]);
    await withAttachedViewer<DG.CalendarViewer>(noDates, DG.VIEWER.CALENDAR, {}, (v) => {
      expect(v instanceof DG.CalendarViewer, true);
      expect(v.dateColumn, null);
      expectNoThrow(() => v.getOptions(true));
    });
  });
});
