import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from '../ui/utils';
import {Grid} from "datagrok-api/dg";

category('Grid', () => {
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
        if (!hasTag(raceTags,'.color-coding-categorical') || !hasTag(raceTags,'{"Asian":4278190335,"Black":4286578816,"Caucasian":4278547786,"Other":4293188935}'))
            throw 'Categorical Color Coding error'

        //numerical HEIGHT column check for Conditional ColorCoding
        let heightTags:string[] = Array.from(demog.col('height')!.tags);
        if (!hasTag(heightTags,'.color-coding-type') || !hasTag(heightTags,'Conditional') || !hasTag(heightTags,'.color-coding-conditional') || !hasTag(heightTags,'{"20-170":"#00FF00","170-190":"#220505"}'))
            throw 'Conditional Color Coding error'

        //numerical AGE column check for Linear ColorCoding
        let ageTags:string[] = Array.from(demog.col('age')!.tags);
        if (!hasTag(ageTags,'.color-coding-type') || !hasTag(ageTags,'Linear') || !hasTag(ageTags,'.color-coding-linear') || !hasTag(ageTags,'[4294944000, 4278255360]'))
            throw 'Linear Color Coding error'
    });

    test('grid.columnVisibility', async () => {
        let studyColVisible = v.grid.columns.byName('~study')!.visible;

        v.grid.columns.setVisible(['age', 'sex', 'race', 'height', 'weight', 'site', 'subj', 'started']);
        let diseaseColVisible = v.grid.columns.byName('disease')!.visible;

        if (studyColVisible != false)
            throw 'Hiding a column by adding ~ to the name doesn\'t work'

        if (diseaseColVisible != false)
            throw 'Hiding a column by using columns.setVisible doesn\'t work'
    });


    test('grid.columnControlledValues', async () => {
        demog.col('site')!.tags[DG.TAGS.CHOICES] = '["New York", "Buffalo"]';
        demog.col('site')!.tags[DG.TAGS.AUTO_CHOICES] = 'New York';

        let siteTags:string[] = Array.from(demog.col('site')!.tags);

        if (!hasTag(siteTags,'.choices') || !hasTag(siteTags,'["New York", "Buffalo"]') || !hasTag(siteTags,'.auto-choices') || !hasTag(siteTags,'New York'))
            throw 'Column Controlled Values (Choices) error'
    });

    after(async () => {
        v.close();
        grok.shell.closeAll();
    });

    function hasTag(colTags:string[], colTagValue:string): boolean {
        for (let i = 0; i < colTags.length; i++) {
            for (let j = 0; j < colTags[i].length; j++) {
                console.log('TAG Value: ' + colTags[i][j]);
                if (colTags[i][j] == colTagValue)
                    return true;
            }
        }
        return false;
    }

});