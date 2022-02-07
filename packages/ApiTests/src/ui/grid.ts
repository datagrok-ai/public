import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';
import {Grid} from "datagrok-api/dg";

category('UI: Grid', () => {
    let v: DG.TableView;
    let  demog = grok.data.demo.demog(1000);
    demog.columns.byName('study').name = '~study';

    before(async () => {
       v = grok.shell.addTableView(demog);
    });

    test('grid.setOrder', async () => {
       v.grid.columns.setOrder(['race', 'age']);
       let firstCol = v.grid.columns.byIndex(1);
       let secondCol = v.grid.columns.byIndex(2);
       if(firstCol?.dart.columnName != 'race' || secondCol?.dart.columnName != 'age')
           throw 'grid.setOrder does not work'
    });

    test('grid.resizeColumn', async () => {
        let check = false;

        v.grid.onColumnResized.subscribe((_) => { check = true });
        v.grid.columns.byName('age')!.width = 200;

        if (check == false)
            throw 'Column Resize error';
    });

    test('grid.filter', async () => {
        demog.rows.match('sex = M').filter()
        if(demog.filter.trueCount != 605)
            throw 'Filtering error'
    });

    test('grid.colorCoding', async () => {
        v.grid.col('race')!.categoryColors = {
           'Asian': 0xFF0000FF,
           'Black': 0xFF800080,
           'Caucasian' : 0xFF05754A,
           'Other' : 0XFFE4DD47
           };

        demog.col('height')!.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Conditional';
        demog.col('height')!.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = `{"20-170":"#00FF00","170-190":"#220505"}`;

        demog.col('age')!.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
        demog.col('age')!.tags[DG.TAGS.COLOR_CODING_LINEAR] = `[${DG.Color.orange}, ${DG.Color.green}]`;

        //categorical RACE column check
        let raceTags:string[] = Array.from(demog.col('race')!.tags);
        if (raceTags[0][0] != '.color-coding-categorical' || raceTags[0][1] != '{"Asian":4278190335,"Black":4286578816,"Caucasian":4278547786,"Other":4293188935}')
            throw 'Categorical Color Coding error'

        //numerical height column check for Conditional ColorCoding
        let heightTags:string[] = Array.from(demog.col('height')!.tags);
        if (heightTags[0][0] != '.color-coding-type' || heightTags[0][1] != 'Conditional' || heightTags[1][0] != '.color-coding-conditional' || heightTags[1][1] != '{"20-170":"#00FF00","170-190":"#220505"}')
            throw 'Conditional Color Coding error'

        //numerical height column check for Linear ColorCoding
        let ageTags:string[] = Array.from(demog.col('age')!.tags);
        if (ageTags[0][0] != '.color-coding-type' || ageTags[0][1] != 'Linear' || ageTags[1][0] != '.color-coding-linear' || ageTags[1][1] != '[4294944000, 4278255360]')
            throw 'Linear Color Coding error'
    });

    test('grid.columnVisibility', async () => {
        let studyColVisiable = v.grid.columns.byName('~study')!.visible;

        v.grid.columns.setVisible(['age', 'sex', 'race', 'height', 'weight', 'site', 'subj', 'started']);
        let diseaseColVisiable = v.grid.columns.byName('disease')!.visible;

        if (studyColVisiable != false)
            throw 'Hiding a column by adding ~ to the name doesn\'t work'

        if (diseaseColVisiable != false)
            throw 'Hiding a column by using columns.setVisible doesn\'t work'
    });


    test('grid.columnControlledValues', async () => {
        demog.col('site')!.tags[DG.TAGS.CHOICES] = '["New York", "Buffalo"]';
        demog.col('site')!.tags[DG.TAGS.AUTO_CHOICES] = 'New York';

        let siteTags:string[] = Array.from(demog.col('site')!.tags);

        console.log('site00' + siteTags[0][0])
        console.log('site01' + siteTags[0][1])

        if (siteTags[0][0] != '.choices' || siteTags[0][1] != '["New York", "Buffalo"]' || siteTags[1][0] !='.auto-choices' || siteTags[1][1] !='New York')
            throw 'Column Controlled Values (Choices) error'
    });

    /*after(async () => {
        v.close();
        grok.shell.closeAll();
    }); */

});