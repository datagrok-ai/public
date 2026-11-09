/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {
  assignLanes, combineHubs, connectionsGeometry, getDrawnConnections, IDrawnConnection, layoutConnections,
  lineHeightWithLanes, measureConnectionsFont,
} from '@datagrok-libraries/bio/src/utils/cell-renderer-connections';
import {TAGS as mmcrTAGS} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import {MmcrTemps, rendererSettingsChangedState} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';
import {getGridCellColTemp} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {monomerToShort} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

import {_package} from '../package-test';
import {parseBioRendererDefaults} from '../utils/renderer-defaults';
import {prepareBiln} from './biln-tests';
import {awaitGrid} from './utils';

/** BILN peptides whose connections the cell draws (flat positions in the comments) */
const bilnRows = {
  oxytocin: 'C(1,3)-Y-I-Q-N-C(1,3)-P-L-G-am', // 0-5
  linaclotide: 'C(1,3)-C(2,3)-E-Y-C(3,3)-C(1,3)-N-P-A-C(2,3)-T-G-C(3,3)-Y', // 0-5, 1-9, 4-12
  ziconotide: 'C(1,3)-K-G-K-G-A-K-C(2,3)-S-R-L-M-Y-D-C(3,3)-C(1,3)-T-G-S-C(2,3)-R-S-G-K-C(3,3)-am', // 0-15, 7-19, 14-24
  kalata: 'G(4,1)-L-P-V-C(1,3)-G-E-T-C(2,3)-V-G-G-T-C(3,3)-N-T-P-G-C(1,3)-T-C(2,3)-S-W-P-V-C(3,3)-T-R-N(4,2)', // ring 0-28
  insulin: 'G-I-V-E-Q-C(1,3)-C(2,3)-T-S-I-C(1,3)-S-L-Y-Q-L-E-N-Y-C(3,3)-N.F-V-N-Q-H-L-C(2,3)-G-S-H-L-V-E-A-L-Y-L-V-C(3,3)-G-E-R-G-F-F-Y-T-P-K-T', // 5-10, 6-27, 19-39
  pk15: 'A-C(1,3)-S-D-R-F-R-N-C(2,3)-P-A-D-E-A-L-C(3,3)-G.THA(1,1)(2,2)(3,3)', // 1, 8, 15 to 17
  sfti: 'C(1,3)(2,1)-T-K-S-I-P-P-I-C(1,3)-F-P-D-G-R(2,2)', // 0-8, ring 0-13
  trivalent: 'R(1,1)-G-D-dF-K(1,2)(4,3).R(2,1)-G-D-dF-K(2,2)(5,3).R(3,1)-G-D-dF-K(3,2)(6,3).THA(4,1)(5,2)(6,3)',
  apart: 'C(1,3)-A-C(1,3)-G-C(2,3)-A-C(2,3)', // 0-2, 4-6
  twoChainRing: 'F(1,1)-G-D-C(2,3)-Y(3,2).P(3,1)-A-C(2,3)-D-F(1,2)', // two chains closed into one ring, joined end to end
};

const pairs = (conns: IDrawnConnection[]): string[] => conns.map((c) => `${c.from}-${c.to}`).sort();

category('renderers: connections', () => {
  let seqHelper: ISeqHelper;
  let libHelper: IMonomerLibHelper;
  let libSettings: UserLibSettings;

  before(async () => {
    seqHelper = await getSeqHelper();
    libHelper = await getMonomerLibHelper();
    libSettings = await getUserLibSettings();
    await libHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(libSettings);
    await libHelper.loadMonomerLib(true);
  });

  /** A BILN column of {@link list}, detected and refined, that draws its connections whatever the stand's package
   * setting is */
  const bilnColumn = async (list: string[]): Promise<DG.Column<string>> => {
    const col = await prepareBiln(list, seqHelper);
    col.name = 'biln';
    col.setTag(mmcrTAGS.showConnections, 'true');
    return col as DG.Column<string>;
  };
  /** The BILN column's renderer back once it holds {@link ready}: the back initializes asynchronously, and a grid draw
   * before it is done draws nothing (the grid redraws when it is) */
  const backWhen = async (tv: DG.TableView, col: DG.Column, ready: (back: any) => boolean, what: string,
  ): Promise<MonomerPlacer> => {
    await awaitCheck(() => {
      const back = col.temp.rendererBack;
      if (back && ready(back)) return true;
      tv.grid.invalidate();
      return false;
    }, `the renderer did not ${what}`, 10000, 100);
    return col.temp.rendererBack as MonomerPlacer;
  };
  /** ... once it has drawn the rows at {@link h} on one line */
  const drawnAt = (tv: DG.TableView, col: DG.Column, h: number): Promise<MonomerPlacer> =>
    backWhen(tv, col, (b) => b._widthsH === h, `draw at ${h} px`);

  /** Row 0 of a BILN column of {@link biln} and two more BILN rows (detection reads a sample) */
  const splitted = async (biln: string): Promise<ISeqSplitted> =>
    seqHelper.getSeqHandler(await bilnColumn([biln, bilnRows.insulin, bilnRows.oxytocin])).getSplitted(0);

  test('connections in flat positions across chains', async () => {
    const conns = getDrawnConnections(await splitted(bilnRows.insulin));
    expectArray(pairs(conns), ['19-39', '5-10', '6-27']);
    const inter = conns.filter((c) => c.chainFrom !== c.chainTo).map((c) => `${c.from}:${c.chainFrom}-${c.to}:${c.chainTo}`);
    expectArray(inter.sort(), ['19:0-39:1', '6:0-27:1']);
    expect(conns.every((c) => c.rFrom === 3 && c.rTo === 3 && !c.headToTail && !c.hbond), true);
  });

  test('a chain closed head to tail is a ring', async () => {
    const conns = getDrawnConnections(await splitted(bilnRows.kalata));
    const rings = conns.filter((c) => c.headToTail);
    expectArray(pairs(rings), ['0-28']);
    expect(rings[0].rFrom === 1 && rings[0].rTo === 2, true);
    expect(conns.length, 4);
  });

  test('lanes: shorter bonds nearer, overlapping ones apart, rings below', async () => {
    const zic = getDrawnConnections(await splitted(bilnRows.ziconotide));
    expect(JSON.stringify(assignLanes(zic)), JSON.stringify({above: 3, below: 0}));
    const laneOf = (from: number): number => zic.find((c) => c.from === from)!.lane;
    expectArray([laneOf(14), laneOf(7), laneOf(0)], [0, 1, 2]);

    const apart = getDrawnConnections(await splitted(bilnRows.apart));
    expect(JSON.stringify(assignLanes(apart)), JSON.stringify({above: 1, below: 0}));

    // one monomer with a disulfide (above) and the ring closure (below)
    const sfti = getDrawnConnections(await splitted(bilnRows.sfti));
    expect(JSON.stringify(assignLanes(sfti)), JSON.stringify({above: 1, below: 1}));
    expect(sfti.find((c) => c.headToTail)!.side, -1);
  });

  test('a scaffold bonded to three monomers is one comb', async () => {
    const combined = combineHubs(getDrawnConnections(await splitted(bilnRows.pk15)));
    expect(combined.length, 1);
    const comb = combined[0].comb!;
    expect(comb.hub, 17);
    expectArray(comb.members.slice().sort((a, b) => a - b), [1, 8, 15]);
    expect(combined[0].from === 1 && combined[0].to === 17, true);
    // two bonds on one monomer stay two brackets
    expect(combineHubs(getDrawnConnections(await splitted(bilnRows.sfti))).length, 2);
  });

  test('the lanes fit the cell, packed closer when it is short, else the row prints as before', async () => {
    const ss = await splitted(bilnRows.ziconotide);
    const met = measureConnectionsFont('12px monospace', 12);
    const normal = connectionsGeometry(12, false);
    const packed = connectionsGeometry(12, true);
    const needNormal = Math.ceil(lineHeightWithLanes(3, 0, normal, met) + 2);
    const needPacked = Math.ceil(lineHeightWithLanes(3, 0, packed, met) + 2);
    expect(needPacked < needNormal, true);
    expect(layoutConnections(ss, needNormal, 12, met, true)?.geo.laneGap, normal.laneGap);
    expect(layoutConnections(ss, needNormal - 1, 12, met, true)?.geo.laneGap, packed.laneGap);
    expect(layoutConnections(ss, needPacked - 1, 12, met, true), null);
    // a row without connections has no layout
    expect(layoutConnections(await splitted('A-C-D-E-F-G-H-I-K-L-M-N'), 100, 12, met, true), null);
  });

  test('two chains go on two lines in a tall cell, on one line otherwise', async () => {
    const met = measureConnectionsFont('12px monospace', 12);
    const insulin = await splitted(bilnRows.insulin);
    expectArray(layoutConnections(insulin, 100, 12, met, true)!.chains, [0, 21]);
    expectArray(layoutConnections(insulin, 40, 12, met, true)!.chains, []);
    // a scrolled sequence stays on one line
    expectArray(layoutConnections(insulin, 100, 12, met, false)!.chains, []);
    // four chains, three of them rings: not a ladder or a branch
    expectArray(layoutConnections(await splitted(bilnRows.trivalent), 400, 12, met, true)!.chains, []);
    // two chains joined backbone to backbone are one backbone: one line
    expectArray(layoutConnections(await splitted(bilnRows.twoChainRing), 400, 12, met, true)!.chains, []);
  });

  test('a row that draws its connections prints canonical symbols, else the originals', async () => {
    const col = await bilnColumn([bilnRows.oxytocin, bilnRows.linaclotide]);
    col.temp[MmcrTemps.maxMonomerLength] = 50; // C(1,3) whole when printed
    const df = DG.DataFrame.fromColumns([col]);
    const tv = grok.shell.addTableView(df);
    try {
      // after the first draw, which sets the renderer's default row height
      await awaitGrid(tv.grid);
      tv.grid.props.rowHeight = 40;
      const back = await drawnAt(tv, col, 40);
      // C(1,3) at 0 is as wide as Y at 1: both print one letter
      const [drawn] = back.getCellMonomerLengths(0, 2000);
      expect(drawn[0], drawn[1]);

      // as the Show Connections input of the renderer's settings does
      col.setTag(mmcrTAGS.showConnections, 'false');
      col.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.true;
      df.fireValuesChanged();
      await backWhen(tv, col, (b) => b._showConnections === false, 'turn the connections off');
      const [printed] = back.getCellMonomerLengths(0, 2000);
      expect(printed[0] > printed[1], true, `C(1,3) printed whole is wider than Y: ${printed[0]}, ${printed[1]}`);
      expect(printed[1], drawn[1]);
    } finally {
      tv.close();
    }
  });

  test('a row too short for its lanes prints the originals', async () => {
    const col = await bilnColumn([bilnRows.linaclotide, bilnRows.oxytocin]);
    col.temp[MmcrTemps.maxMonomerLength] = 50;
    const df = DG.DataFrame.fromColumns([col]);
    const tv = grok.shell.addTableView(df);
    try {
      // a row height one packed lane fits and three do not, in this browser's monospace font
      const met = measureConnectionsFont('12px monospace', 12);
      const packed = connectionsGeometry(12, true);
      const oneLane = Math.ceil(lineHeightWithLanes(1, 0, packed, met) + 2);
      expect(oneLane < Math.ceil(lineHeightWithLanes(3, 0, packed, met) + 2), true);
      await awaitGrid(tv.grid);
      tv.grid.props.rowHeight = oneLane;
      const back = await drawnAt(tv, col, oneLane);
      const [lina] = back.getCellMonomerLengths(0, 2000);
      expect(lina[0] > lina[2], true, 'linaclotide prints C(1,3) whole');
      const [oxy] = back.getCellMonomerLengths(1, 2000);
      const b = back as any;
      const state = () => JSON.stringify({widthsH: b._widthsH, show: b._showConnections, flags: b._monomerLengthCanonical,
        font: b.props.font, fontSize: b.props.fontSize, layout: b.connectionsLayout(1,
          seqHelper.getSeqHandler(col).getSplitted(1), b._widthsH, 0)?.geo ?? null, oxy: oxy.slice(0, 3)});
      expect(oxy[0], oxy[1], `oxytocin's one lane fits ${oneLane} px: ${state()}`);
    } finally {
      tv.close();
    }
  });

  test('the column tag overrides the package default', async () => {
    const col = await bilnColumn([bilnRows.oxytocin, bilnRows.insulin]);
    col.setTag(mmcrTAGS.showConnections, '');
    const make = (showConnections: boolean): MonomerPlacer => new MonomerPlacer(null, col, _package.logger, 4, () => {
      const {font, fontWidth} = MonomerPlacer.getFontSettings(col);
      return {font, fontCharWidth: fontWidth, separatorWidth: 0, monomerToShort, showConnections};
    });
    const shows = async (back: MonomerPlacer): Promise<boolean> => {
      try {
        await back.init();
        return (back as any)._showConnections;
      } finally {
        (back as any).destroy();
      }
    };
    expect(await shows(make(true)), true);
    expect(await shows(make(false)), false);
    col.setTag(mmcrTAGS.showConnections, 'false');
    expect(await shows(make(true)), false);
    col.setTag(mmcrTAGS.showConnections, 'true');
    expect(await shows(make(false)), true);
    // a placer whose creator passes no default (another package's) draws no connections
    col.setTag(mmcrTAGS.showConnections, '');
    const plain = new MonomerPlacer(null, col, _package.logger, 4, () => {
      const {font, fontWidth} = MonomerPlacer.getFontSettings(col);
      return {font, fontCharWidth: fontWidth, separatorWidth: 0, monomerToShort};
    });
    expect(await shows(plain), false);
  });

  test('multiline marks chain breaks and hit-tests monomers only', async () => {
    const col = await bilnColumn([bilnRows.insulin, bilnRows.oxytocin]);
    col.setTag('renderMultiline', 'true');
    const df = DG.DataFrame.fromColumns([col]);
    const tv = grok.shell.addTableView(df);
    try {
      await awaitGrid(tv.grid);
      tv.grid.props.rowHeight = 160;
      tv.grid.col('biln')!.width = 400;
      const back = await backWhen(tv, col, (b) => b._cellBounds.get(0)?.length === 51, 'lay out 51 monomers');
      const bounds: { sequencePosition: number, lineIdx: number }[] = (back as any)._cellBounds.get(0);
      expectArray(bounds.map((b) => b.sequencePosition), Array.from({length: 51}, (_, i) => i));
      // the '|' slot shifts chain B by one slot: B1 is not where A22 would be
      const ctx = tv.grid.canvas.getContext('2d')!;
      ctx.font = back.props.font; // the layout measures (and caches) the labels in the context's font
      const layout = (back as any).calculateMultiLineLayoutDynamic(
        ctx, 400, 160, seqHelper.getSeqHandler(col).getSplitted(0), 0, 4, true, 0);
      const separators = layout.lineLayouts.flatMap((l: any) => l.elements).filter((e: any) => e.isSeparator);
      expect(separators.length, 1);
      expect(separators[0].om, '|');
    } finally {
      tv.close();
    }
  });

  test('a multiline row with no room for the connections prints as before', async () => {
    // the connections' room over each line would leave a 30 px multiline cell without a line: it prints the marks
    const col = await bilnColumn([bilnRows.linaclotide, bilnRows.insulin]);
    col.setTag('renderMultiline', 'true');
    col.temp[MmcrTemps.maxMonomerLength] = 50;
    const df = DG.DataFrame.fromColumns([col]);
    const tv = grok.shell.addTableView(df);
    try {
      await awaitGrid(tv.grid);
      tv.grid.props.rowHeight = 30;
      tv.grid.col('biln')!.width = 600;
      const back = await backWhen(tv, col, (b) => b._cellBounds.get(0)?.length > 0, 'draw the multiline row');
      const bounds: { sequencePosition: number, bounds: DG.Rect }[] = (back as any)._cellBounds.get(0);
      expect(bounds.length > 0, true, 'the row draws its monomers');
      // C(1,3) printed whole: its slot is as wide as the widest original label
      const ctx = tv.grid.canvas.getContext('2d')!;
      ctx.font = back.props.font;
      expect(bounds[0].bounds.width >= ctx.measureText('C(1,3)').width - 0.5, true, `slot ${bounds[0].bounds.width}`);
    } finally {
      tv.close();
    }
  });

  test('a row on two lines hit-tests each monomer on its line', async () => {
    const col = await bilnColumn([bilnRows.insulin, bilnRows.oxytocin]);
    const df = DG.DataFrame.fromColumns([col]);
    const tv = grok.shell.addTableView(df);
    try {
      await awaitGrid(tv.grid);
      tv.grid.props.rowHeight = 100;
      tv.grid.col('biln')!.width = 700;
      const back = await backWhen(tv, col, (b) => b._stackedBounds.get(100)?.get(0)?.length === 51,
        'draw insulin on two lines');
      // hit bounds by cell height, then row
      const bounds: { sequencePosition: number, lineIdx: number, bounds: DG.Rect }[] =
        (back as any)._stackedBounds.get(100).get(0);
      expect(bounds.length, 51);
      expect(bounds.every((b) => b.lineIdx === (b.sequencePosition < 21 ? 0 : 1)), true);
      const a7 = bounds.find((b) => b.sequencePosition === 6)!.bounds;
      const b7 = bounds.find((b) => b.sequencePosition === 27)!.bounds;
      expect(a7.y < b7.y, true);
      // A7 sits over B7: they are bonded and the B chain moves under its first partner
      expect(Math.abs((a7.x + a7.width / 2) - (b7.x + b7.width / 2)) < 2, true, `A7 at ${a7.x}, B7 at ${b7.x}`);

      // too short for two lines: one line again, hit-tested by widths
      tv.grid.props.rowHeight = 30;
      await drawnAt(tv, col, 30);
      expect((back as any)._stackedBounds.get(30)?.get(0) === undefined, true);
      expect((back as any)._cellBounds.get(0) === undefined, true);
    } finally {
      tv.close();
    }
  });

  test('a tall row strokes its brackets, a short one does not', async () => {
    const col = await bilnColumn([bilnRows.ziconotide, bilnRows.insulin]);
    const df = DG.DataFrame.fromColumns([col]);
    const tv = grok.shell.addTableView(df);
    try {
      await awaitGrid(tv.grid);
      const back = await backWhen(tv, col, (b) => b._widthsH > 0, 'draw');
      const ctx = tv.grid.canvas.getContext('2d')!;
      let arcs = 0;
      const g = new Proxy(ctx, {
        get: (t, k) => {
          const v = (t as any)[k];
          if (k === 'arcTo') return (...a: any[]) => {arcs++; return (t.arcTo as any)(...a);};
          return typeof v === 'function' ? v.bind(t) : v;
        },
        set: (t, k, v) => {(t as any)[k] = v; return true;},
      }) as CanvasRenderingContext2D;
      const gridCell = tv.grid.cell('biln', 0);
      const b = gridCell.bounds;
      back.render(g, b.x, b.y, b.width, 40, gridCell, gridCell.style);
      expect(arcs, 6); // two rounded corners for each of the three bridges
      arcs = 0;
      back.render(g, b.x, b.y, b.width, 16, gridCell, gridCell.style);
      expect(arcs, 0);
    } finally {
      tv.close();
    }
  });

  test('a column without graph info takes none of the connection steps', async () => {
    const df = DG.DataFrame.fromCsv('seq\nANTTYKNCCQ\nMFHPCSGCCK\nQWDRLEVCCA\n');
    await grok.data.detectSemanticTypes(df);
    df.getCol('seq').setTag(mmcrTAGS.showConnections, 'true');
    const tv = grok.shell.addTableView(df);
    try {
      await awaitGrid(tv.grid);
      await awaitGrid(tv.grid);
      const [, , temp] = getGridCellColTemp<any, any>(tv.grid.cell('seq', 0));
      const back = temp.rendererBack as MonomerPlacer;
      expect((back as any)._showConnections, false); // FASTA rows have no graph info, even with the tag on
      expect((back as any)._connLayouts.size, 0);
    } finally {
      tv.close();
    }
  });

  test('a hovered monomer marks its partners in the grid it is hovered in only', async () => {
    const col = await bilnColumn([bilnRows.ziconotide]);
    const back = new MonomerPlacer(null, col, _package.logger, 4, () => {
      const {font, fontWidth} = MonomerPlacer.getFontSettings(col);
      return {font, fontCharWidth: fontWidth, separatorWidth: 0, monomerToShort, showConnections: true};
    });
    try {
      await back.init();
      const ss = seqHelper.getSeqHandler(col).getSplitted(0);
      const g1 = {dart: {}} as unknown as DG.Grid;
      const g2 = {dart: {}} as unknown as DG.Grid;
      (back as any)._hoverConn = {row: 0, pos: 7, grid: g2};
      const in2 = (back as any).hoverPartners(0, ss, {grid: g2});
      expect(in2.pos, 7);
      expectArray([...in2.partners], [19]);
      expect((back as any).hoverPartners(0, ss, {grid: g1}).pos, null);
      expect((back as any).hoverPartners(1, ss, {grid: g2}).pos, null);
      expect((back as any).hoverPartners(0, ss, null).pos, null, 'off the grid');
    } finally {
      (back as any).destroy();
    }
  });

  test('the hovered monomer\'s bonds join its tooltip', async () => {
    const col = await bilnColumn([bilnRows.ziconotide, bilnRows.insulin]);
    const back = new MonomerPlacer(null, col, _package.logger, 4, () => {
      const {font, fontWidth} = MonomerPlacer.getFontSettings(col);
      return {font, fontCharWidth: fontWidth, separatorWidth: 0, monomerToShort, showConnections: true};
    });
    try {
      await back.init();
      const sh = seqHelper.getSeqHandler(col);
      const tooltipOf = (row: number, pos: number): string => {
        const ss = sh.getSplitted(row);
        const conns = getDrawnConnections(ss).filter((c) => c.from === pos || c.to === pos);
        return ((back as any).connectionsTooltip(ss, pos, conns) as HTMLElement).innerText;
      };
      expect(tooltipOf(0, 7), 'Bond: C8 ↔ C20 (R3–R3)');
      expect(tooltipOf(1, 27), 'Bond: C7 (chain 2) ↔ C7 (chain 1) (R3–R3)');
    } finally {
      (back as any).destroy();
    }
  });

  test('the hover clears when the mouse leaves the cells for a spot with no cell', async () => {
    const col = await bilnColumn([bilnRows.ziconotide, bilnRows.oxytocin]);
    const df = DG.DataFrame.fromColumns([col]);
    const tv = grok.shell.addTableView(df);
    try {
      await awaitGrid(tv.grid);
      const back = await backWhen(tv, col, (b) => b.seqHelper != null, 'initialize') as any;
      const overlay = tv.grid.overlay;
      const moveTo = (x: number, y: number): void => {
        const r = overlay.getBoundingClientRect();
        overlay.dispatchEvent(new MouseEvent('mousemove', {clientX: r.left + x, clientY: r.top + y}));
      };
      const lastRow = tv.grid.cell('biln', df.rowCount - 1).bounds;
      back.setHoverConn({row: 0, pos: 7, grid: tv.grid});
      expect(back._hoverLeave !== null, true, 'watches the grid while hovering');
      moveTo(lastRow.x + 5, lastRow.bottom + 20); // the empty area below the rows
      expect(back._hoverConn, null, 'a move below the rows clears the hover');
      expect(back._hoverLeave, null, 'stops watching with the hover gone');
      back.setHoverConn({row: 0, pos: 7, grid: tv.grid});
      overlay.dispatchEvent(new MouseEvent('mouseleave'));
      expect(back._hoverConn, null, 'leaving the grid clears the hover');
    } finally {
      tv.close();
    }
  });

  test('a font size set on the column reaches every grid that shows it', async () => {
    const df = DG.DataFrame.fromCsv('seq\nANTTYKN\nMFHPCSG\nQWDRLEV\n');
    await grok.data.detectSemanticTypes(df);
    const col = df.getCol('seq');
    const tv = grok.shell.addTableView(df);
    try {
      const second = tv.addViewer(DG.VIEWER.GRID) as DG.Grid;
      await awaitGrid(tv.grid);
      await awaitGrid(second);
      col.temp[MmcrTemps.fontSize] = 18;
      col.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.true;
      df.fireValuesChanged();
      const fontOf = (grid: DG.Grid): string | undefined =>
        (getGridCellColTemp<any, any>(grid.cell('seq', 0))[2].rendererBack as MonomerPlacer | undefined)?.props.font;
      // before the fix, the second grid's renderer kept 12 px: the first one to draw cleared the changed flag
      await awaitCheck(() => {
        tv.grid.invalidate();
        second.invalidate();
        return fontOf(tv.grid) === '18px monospace' && fontOf(second) === '18px monospace';
      }, `fonts: ${fontOf(tv.grid)}, ${fontOf(second)}`, 5000, 100);
    } finally {
      tv.close();
    }
  });

  test('the package font size is the default of a column without its own', async () => {
    const col = DG.Column.fromList(DG.TYPE.STRING, 'seq', ['ACD']);
    expect(MonomerPlacer.getFontSettings(col).font, '12px monospace');
    expect(MonomerPlacer.getFontSettings(col, 16).font, '16px monospace');
    expect(MonomerPlacer.getFontSettings(col, 16).fontSize, 16);
    col.temp[MmcrTemps.fontSize] = 14;
    expect(MonomerPlacer.getFontSettings(col, 16).font, '14px monospace');
  });

  test('package defaults: long is no shortening, an unreadable value the declared default', async () => {
    const d = parseBioRendererDefaults;
    expect(JSON.stringify(d(null)), JSON.stringify({fontSize: 12, maxMonomerLength: 4, showConnections: true}));
    expect(d({MaxMonomerLength: 'long'}).maxMonomerLength, 50);
    expect(d({MaxMonomerLength: '7'}).maxMonomerLength, 7);
    expect(d({MaxMonomerLength: 'x'}).maxMonomerLength, 4);
    expect(d({FontSize: '16'}).fontSize, 16);
    expect(d({FontSize: 18}).fontSize, 18);
    expect(d({FontSize: 0}).fontSize, 12);
    expect(d({ShowConnections: false}).showConnections, false);
    expect(d({ShowConnections: 'false'}).showConnections, false);
    expect(d({ShowConnections: 'true'}).showConnections, true);
  });
});
