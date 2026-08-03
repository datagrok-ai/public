import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {SEMTYPE} from '../utils/proteomics-types';
import {
  resolveGeneLabels, STORE_GENE_LABELS, SCHEMA_V_GENE_LABELS,
  detectSpecies, isEnsemblEligible, isPredicted,
  extractReadableDescription, applyMarkerRules,
  lookupEnsemblBatch, type EnsemblEntry,
} from '../utils/gene-label-resolver';
import {parseMaxQuantText} from '../parsers/maxquant-parser';
import {parseFragPipeText} from '../parsers/fragpipe-parser';
import {parseSpectronautCandidatesText} from '../parsers/spectronaut-candidates-parser';
import {parseSpectronautText} from '../parsers/spectronaut-parser';

/** Build a fake Response that the test mock for grok.dapi.fetchProxy can return. */
function makeResponse(body: any, opts: {status?: number; retryAfter?: string} = {}): Response {
  const status = opts.status ?? 200;
  const headers = new Headers();
  if (opts.retryAfter != null) headers.set('Retry-After', opts.retryAfter);
  // node-fetch / browser Response constructors accept (BodyInit, ResponseInit).
  return new Response(body == null ? null : JSON.stringify(body), {status, headers});
}

/** Replaces grok.dapi.fetchProxy with `fn` for the duration of `run`, then
 * restores the original. Patches BOTH the dapi instance (own property) and its
 * shared prototype: older Datagrok exposes `fetchProxy` on the prototype, while
 * newer builds define it as an own/bound property on the instance — where a
 * prototype-only patch is invisible to the source's `grok.dapi.fetchProxy`
 * call and the real network gets hit (CI "latest" failure mode). Patching both
 * covers either shape. (See project memory reference_dapi_fresh_instance_patch.) */
async function withFetchProxy(
  fn: (url: string, init?: RequestInit) => Promise<Response>,
  run: () => Promise<void>,
): Promise<void> {
  const dapi: any = grok.dapi;
  const proto = Object.getPrototypeOf(dapi);
  const hadOwn = Object.prototype.hasOwnProperty.call(dapi, 'fetchProxy');
  const origOwn = hadOwn ? dapi.fetchProxy : undefined;
  const origProto = proto.fetchProxy;
  dapi.fetchProxy = fn;
  proto.fetchProxy = fn;
  try {
    await run();
  } finally {
    if (hadOwn) dapi.fetchProxy = origOwn;
    else delete dapi.fetchProxy;
    proto.fetchProxy = origProto;
  }
}

/** Replaces grok.dapi.userDataStorage methods for the duration of `run`. */
async function withUserDataStorage(
  store: Map<string, any>,
  run: () => Promise<void>,
): Promise<void> {
  const uds: any = grok.dapi.userDataStorage;
  const proto = Object.getPrototypeOf(uds);
  const origGet = proto.get;
  const origPut = proto.put;
  proto.get = async (key: string) => store.get(key) ?? null;
  proto.put = async (key: string, value: any) => { store.set(key, value); };
  try {
    await run();
  } finally {
    proto.get = origGet;
    proto.put = origPut;
  }
}

/** Replaces grok.shell.warning with a no-op for the duration of `run` and
 * returns the captured warnings. */
async function captureWarnings(run: () => Promise<void>): Promise<string[]> {
  const captured: string[] = [];
  const shell: any = grok.shell;
  const original = shell.warning;
  shell.warning = (msg: string) => { captured.push(msg); };
  try {
    await run();
  } finally {
    shell.warning = original;
  }
  return captured;
}

category('Proteomics: 14-01', () => {
  test('geneLabelResolverSemTypes: SEMTYPE constants are the locked strings', async () => {
    expect(SEMTYPE.DISPLAY_NAME, 'Proteomics-DisplayName');
    expect(SEMTYPE.SOURCE_ID, 'Proteomics-SourceId');
    expect(SEMTYPE.NUMERATOR_MEAN, 'Proteomics-NumeratorMean');
    expect(SEMTYPE.DENOMINATOR_MEAN, 'Proteomics-DenominatorMean');
    const all = new Set([
      SEMTYPE.DISPLAY_NAME, SEMTYPE.SOURCE_ID,
      SEMTYPE.NUMERATOR_MEAN, SEMTYPE.DENOMINATOR_MEAN,
      SEMTYPE.PROTEIN_ID, SEMTYPE.GENE_SYMBOL, SEMTYPE.INTENSITY,
      SEMTYPE.LOG2FC, SEMTYPE.P_VALUE, SEMTYPE.SUBCELLULAR_LOCATION,
    ]);
    expect(all.size, 10);
  });

  test('geneLabelResolverDetectorMirror: SEMTYPE values are valid Column.semType strings', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Display Name', ['Myh7†', 'Tnnt3']),
      DG.Column.fromStrings('Source ID', ['ENSRNOG0001', '']),
      DG.Column.fromFloat32Array('Numerator Mean', new Float32Array([1.0, 2.0])),
      DG.Column.fromFloat32Array('Denominator Mean', new Float32Array([1.5, 2.5])),
    ]);
    df.col('Display Name')!.semType = SEMTYPE.DISPLAY_NAME;
    df.col('Source ID')!.semType = SEMTYPE.SOURCE_ID;
    df.col('Numerator Mean')!.semType = SEMTYPE.NUMERATOR_MEAN;
    df.col('Denominator Mean')!.semType = SEMTYPE.DENOMINATOR_MEAN;
    expect(df.col('Display Name')!.semType, 'Proteomics-DisplayName');
    expect(df.col('Source ID')!.semType, 'Proteomics-SourceId');
    expect(df.col('Numerator Mean')!.semType, 'Proteomics-NumeratorMean');
    expect(df.col('Denominator Mean')!.semType, 'Proteomics-DenominatorMean');
  });

  test('detectSpecies: maps prefixes per CK-omics species table', async () => {
    expect(detectSpecies('ENSG00000123'), 'homo_sapiens');
    expect(detectSpecies('ENSMUSG00000456'), 'mus_musculus');
    expect(detectSpecies('MGP_C57BL6NJ_G0000001'), 'mus_musculus');
    expect(detectSpecies('ENSRNOG00000789'), 'rattus_norvegicus');
    expect(detectSpecies('ENSRNO00000001'), 'rattus_norvegicus');
    expect(detectSpecies('ENSDARG00000222'), 'danio_rerio');
    expect(detectSpecies('LOC123456'), null);
    expect(detectSpecies('RGD7890'), null);
    expect(detectSpecies('AABR07012345'), null);
  });

  test('isEnsemblEligible: filters LOC/RGD/AABR per CK-omics line 933', async () => {
    expect(isEnsemblEligible('ENSG00000001'), true);
    expect(isEnsemblEligible('ENSMUSG00000123'), true);
    expect(isEnsemblEligible('ENSRNOG00000456'), true);
    expect(isEnsemblEligible('ENSDARG00000789'), true);
    expect(isEnsemblEligible('MGP_C57BL6NJ_G0000001'), true);
    expect(isEnsemblEligible('LOC123456'), false);
    expect(isEnsemblEligible('RGD789'), false);
    expect(isEnsemblEligible('AABR07012345'), false);
  });

  test('isPredicted: locked CK-omics prefix list', async () => {
    expect(isPredicted('ENSRNOG00000001'), true);
    expect(isPredicted('LOC987'), true);
    expect(isPredicted('AABR07'), true);
    expect(isPredicted('Myh7'), false);
    expect(isPredicted('Tnnt3'), false);
  });

  test('extractReadableDescription: strips Predicted-to prefixes and species suffixes', async () => {
    expect(extractReadableDescription('Predicted to enable actin binding [Rattus norvegicus]'),
      'Actin binding');
    expect(extractReadableDescription('Predicted to be involved in transport. Additional sentence.'),
      'Transport');
    expect(extractReadableDescription('Myosin heavy chain_RAT'), 'Myosin heavy chain');
    expect(extractReadableDescription(''), null);
    expect(extractReadableDescription(null), null);
    expect(extractReadableDescription('uncharacterized protein'), null);
  });

  test('applyMarkerRules: appends * and † per CK-omics line 1011', async () => {
    expect(applyMarkerRules('Myh7', false, true), 'Myh7†');
    expect(applyMarkerRules('Myh7', true, true), 'Myh7*†');
    expect(applyMarkerRules('LOC123', false, true), 'LOC123†');
    expect(applyMarkerRules('Tnnt3', false, false), 'Tnnt3');
    expect(applyMarkerRules('Tnnt3', true, false), 'Tnnt3*');
  });

  test('lookupEnsemblBatch: POSTs ids and parses response', async () => {
    let capturedUrl = '';
    let capturedInit: RequestInit | undefined;
    await withFetchProxy(
      async (url, init) => {
        capturedUrl = url;
        capturedInit = init;
        return makeResponse({
          'ENSRNOG00000001': {external_name: 'Myh7', description: 'Myosin heavy chain 7'},
          'ENSRNOG00000002': null,
        });
      },
      async () => {
        const out = await lookupEnsemblBatch(['ENSRNOG00000001', 'ENSRNOG00000002']);
        expect(capturedUrl, 'https://rest.ensembl.org/lookup/id');
        expect(capturedInit?.method, 'POST');
        const body = JSON.parse((capturedInit?.body as string) ?? '{}');
        expect(Array.isArray(body.ids), true);
        expect(body.ids.length, 2);
        expect(out.size, 1);
        expect(out.get('ENSRNOG00000001')?.external_name, 'Myh7');
      },
    );
  });

  test('lookupEnsemblBatch: retries once on 429 with Retry-After', async () => {
    let calls = 0;
    await withFetchProxy(
      async () => {
        calls++;
        if (calls === 1) return makeResponse(null, {status: 429, retryAfter: '0'});
        return makeResponse({'ENSG00000001': {external_name: 'TP53'}});
      },
      async () => {
        const out = await lookupEnsemblBatch(['ENSG00000001']);
        expect(calls, 2);
        expect(out.get('ENSG00000001')?.external_name, 'TP53');
      },
    );
  });

  test('resolveGeneLabels: no-op path still creates Display Name + Source ID columns', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein ID', ['P12345', 'Q67890', 'O11111']),
      DG.Column.fromStrings('Gene name', ['Myh7', 'Tnnt3', 'Actb']),
    ]);
    df.col('Gene name')!.semType = SEMTYPE.GENE_SYMBOL;
    df.col('Protein ID')!.semType = SEMTYPE.PROTEIN_ID;
    await withFetchProxy(async () => { throw new Error('should not fetch'); }, async () => {
      await resolveGeneLabels(df);
    });
    const dn = df.col('Display Name');
    const sid = df.col('Source ID');
    expect(dn != null, true);
    expect(sid != null, true);
    expect(dn!.semType, 'Proteomics-DisplayName');
    expect(sid!.semType, 'Proteomics-SourceId');
    expect(dn!.get(0), 'Myh7');
    expect(dn!.get(1), 'Tnnt3');
    expect(sid!.get(0), '');
  });

  test('resolveGeneLabels: LOC/RGD/AABR keep raw ID + †, never sent to Ensembl', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Gene name', ['LOC123', 'RGD456', 'AABR789']),
    ]);
    df.col('Gene name')!.semType = SEMTYPE.GENE_SYMBOL;
    let postCount = 0;
    await withUserDataStorage(new Map(), async () => {
      await withFetchProxy(
        async () => { postCount++; return makeResponse({}); },
        async () => {
          await resolveGeneLabels(df);
        },
      );
    });
    expect(postCount, 0);
    expect(df.col('Display Name')!.get(0), 'LOC123†');
    expect(df.col('Display Name')!.get(1), 'RGD456†');
    expect(df.col('Display Name')!.get(2), 'AABR789†');
    expect(df.col('Source ID')!.get(0), 'LOC123');
    expect(df.col('Source ID')!.get(1), 'RGD456');
    expect(df.col('Source ID')!.get(2), 'AABR789');
  });

  test('resolveGeneLabels: disambiguates duplicates with (Source ID) suffix and warns', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Gene name', ['ENSRNOG00000001', 'ENSRNOG00000002', 'Actb']),
    ]);
    df.col('Gene name')!.semType = SEMTYPE.GENE_SYMBOL;
    const warnings = await captureWarnings(async () => {
      await withUserDataStorage(new Map(), async () => {
        await withFetchProxy(
          async () => makeResponse({
            'ENSRNOG00000001': {external_name: 'Myh7'},
            'ENSRNOG00000002': {external_name: 'Myh7'},
          }),
          async () => {
            await resolveGeneLabels(df);
          },
        );
      });
    });
    expect(df.col('Display Name')!.get(0), 'Myh7† (ENSRNOG00000001)');
    expect(df.col('Display Name')!.get(1), 'Myh7† (ENSRNOG00000002)');
    expect(df.col('Display Name')!.get(2), 'Actb');
    const dupWarn = warnings.find((w) => w.includes('duplicate gene names'));
    expect(dupWarn != null, true);
  });

  test('resolveGeneLabels: cache short-circuits second call', async () => {
    const store = new Map<string, any>();
    let postCount = 0;
    const mockFetch = async () => {
      postCount++;
      return makeResponse({'ENSRNOG00000001': {external_name: 'Myh7'}});
    };

    const df1 = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Gene name', ['ENSRNOG00000001']),
    ]);
    df1.col('Gene name')!.semType = SEMTYPE.GENE_SYMBOL;

    const df2 = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Gene name', ['ENSRNOG00000001']),
    ]);
    df2.col('Gene name')!.semType = SEMTYPE.GENE_SYMBOL;

    await withUserDataStorage(store, async () => {
      await withFetchProxy(mockFetch, async () => {
        await resolveGeneLabels(df1);
        const firstCount = postCount;
        await resolveGeneLabels(df2);
        expect(postCount, firstCount); // second call must NOT POST again
      });
    });
    expect(df2.col('Display Name')!.get(0), 'Myh7†');
    // Cache write happened against our store under the correct schema key.
    const cached = store.get(STORE_GENE_LABELS) ?? {};
    expect(cached['__schema_v'], SCHEMA_V_GENE_LABELS);
  });

  test('resolveGeneLabels: degrades gracefully on fetchProxy throw', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Gene name', ['ENSRNOG00000001', 'Actb']),
    ]);
    df.col('Gene name')!.semType = SEMTYPE.GENE_SYMBOL;
    const warnings = await captureWarnings(async () => {
      await withUserDataStorage(new Map(), async () => {
        await withFetchProxy(async () => { throw new Error('network down'); }, async () => {
          await resolveGeneLabels(df);
        });
      });
    });
    expect(df.col('Display Name')!.get(0), 'ENSRNOG00000001†');
    expect(df.col('Display Name')!.get(1), 'Actb');
    expect(df.col('Source ID')!.get(0), 'ENSRNOG00000001');
    const failWarn = warnings.find((w) => w.includes('Ensembl gene-label resolution unavailable'));
    expect(failWarn != null, true);
  });

  test('parser integration: every parser produces Display Name + Source ID columns', async () => {
    // The 4 sync-text parsers all go through the resolver before return. Spectronaut
    // streaming has its own test elsewhere; the text path shares finalizeSpectronaut
    // with the stream path so we cover both via parseSpectronautText.
    const mq = await parseMaxQuantText([
      'Protein IDs\tMajority protein IDs\tGene names\tPotential contaminant\tReverse\t' +
        'Only identified by site\tLFQ intensity Sample1\tLFQ intensity Sample2\tiBAQ',
      'P12345\tP12345\tBRCA1\t\t\t\t1000\t2000\t500',
    ].join('\n'));
    expect(mq.col('Display Name')?.semType, 'Proteomics-DisplayName');
    expect(mq.col('Source ID')?.semType, 'Proteomics-SourceId');

    const fp = await parseFragPipeText([
      'Protein\tGene\tDescription\tMaxLFQ Intensity Sample1\tMaxLFQ Intensity Sample2',
      'P12345\tBRCA1\tProtein BRCA1\t1000\t2000',
    ].join('\n'));
    expect(fp.col('Display Name')?.semType, 'Proteomics-DisplayName');
    expect(fp.col('Source ID')?.semType, 'Proteomics-SourceId');

    const sc = await parseSpectronautCandidatesText([
      'PG.ProteinGroups\tPG.Genes\tComparison (group1/group2)\tAVG Log2 Ratio\tQvalue\tPvalue',
      'P12345\tBRCA1\tA / B\t2.5\t0.001\t0.0005',
    ].join('\n'));
    expect(sc.col('Display Name')?.semType, 'Proteomics-DisplayName');
    expect(sc.col('Source ID')?.semType, 'Proteomics-SourceId');

    const sp = await parseSpectronautText([
      'PG.ProteinGroups\tPG.Genes\tR.Condition\tR.Replicate\tR.FileName\tPG.Qvalue\tPG.Quantity',
      'P12345\tBRCA1\tCtrl\t1\tF1\t0.0001\t1000',
      'P12345\tBRCA1\tCtrl\t2\tF2\t0.0001\t1200',
      'P12345\tBRCA1\tTrt\t1\tF3\t0.0001\t2000',
      'P12345\tBRCA1\tTrt\t2\tF4\t0.0001\t2400',
    ].join('\n'));
    expect(sp.col('Display Name')?.semType, 'Proteomics-DisplayName');
    expect(sp.col('Source ID')?.semType, 'Proteomics-SourceId');
  });

  test('resolveGeneLabels: three-level fallback uses description when external_name missing', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Gene name', ['ENSRNOG00000099']),
    ]);
    df.col('Gene name')!.semType = SEMTYPE.GENE_SYMBOL;
    await withUserDataStorage(new Map(), async () => {
      await withFetchProxy(
        async () => makeResponse({
          'ENSRNOG00000099': {description: 'Predicted to enable kinase activity [Rattus norvegicus]'},
        }),
        async () => { await resolveGeneLabels(df); },
      );
    });
    // pickBestName uses entry.description.split('[')[0].trim() → 'Predicted to enable kinase activity'
    // which still starts with the literal 'Predicted to' — but pickBestName only rejects on the locked
    // PREDICTED_PREFIXES list (ENSRNOG, etc.), not on free-text 'Predicted to'. CK-omics behavior:
    // such a candidate IS accepted by pickBestName; only the readable-description fallback runs the
    // 'Predicted to ...' regex cleanup. Asserting the verbatim CK-omics behavior here.
    expect(df.col('Display Name')!.get(0), 'Predicted to enable kinase activity†');
    expect(df.col('Source ID')!.get(0), 'ENSRNOG00000099');
  });
});
