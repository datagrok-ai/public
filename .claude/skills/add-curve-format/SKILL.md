---
name: add-curve-format
description: Add support for a new curve data format from any Datagrok package, integrating with the Curves package converter system
---

# Add a New Curve Format

Help the user add support for a new curve data format to Datagrok. The Curves package (`packages/Curves`) uses a **converter architecture** where any package can provide support for a new curve notation. When the Curves package initializes, it discovers all converter functions tagged with `meta.role: curveConverter` and `meta.curveFormat: '<name>'` across all loaded packages and registers them automatically.

## Usage
```
/add-curve-format [format-name]
```

## Background

### Native Format

The native format is `IFitChartData` JSON (defined in `@datagrok-libraries/statistics/src/fit/fit-curve.ts`):
```typescript
interface IFitChartData {
  chartOptions?: IFitChartOptions;   // axes bounds, log scale, title, axis names
  seriesOptions?: IFitSeriesOptions; // default series options
  series?: IFitSeries[];             // array of data series with points
}
```

Each `IFitSeries` has `points: IFitPoint[]` where each point has `x`, `y`, and optional `outlier`, plus fit configuration like `fitFunction`, `parameters`, `showFitLine`, etc.

### How the Converter System Works

1. **Detector** (in your package's `detectors.js`) identifies a column's format and sets:
   - `col.semType = 'fit'`
   - `col.setTag('.%curve-format', '<format-name>')` — a friendly format name (e.g. `'my-format'`)
2. **Converter DG function** (in your package, tagged `meta.role: curveConverter` and `meta.curveFormat: '<format-name>'`) returns a sync converter function: `() => (value: string) => string`
3. **Auto-discovery**: The Curves package calls `DG.Func.find({meta: {role: 'curveConverter'}})` at init time, reads each function's `meta.curveFormat` value, and registers the converter under that format name
4. **Caching**: Converted values are cached in a value-based LRU cache keyed by `formatName||cellValue`
5. All rendering, property panels, and statistics use `parseCellValue()` which handles conversion transparently

### Important: Your package does NOT depend on Curves

Your package only needs `@datagrok-libraries/statistics` (for types like `IFitChartData` and `FitConstants`) and `datagrok-api`. The Curves package discovers your converter automatically at runtime.

## Instructions

### Step 1: Create the Converter Function

Create a converter file in your package, e.g. `src/<format-name>-converter.ts`:

```typescript
import {
  IFitChartData,
  IFitPoint,
  IFitSeries,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';

/** Converts <format-name> format to native IFitChartData JSON string. */
export function convert<FormatName>ToJson(value: string): string {
  // 1. Parse the input format (JSON.parse, DOMParser, etc.)
  // 2. Extract points, series options, chart options
  // 3. Build IFitChartData object
  // 4. Return JSON.stringify(chartData)

  const parsed = JSON.parse(value); // or XML parse, etc.

  const points: IFitPoint[] = []; // Map input data to {x, y, outlier?}
  // ... populate points from parsed data

  const series: IFitSeries = {
    points: points,
    fitFunction: 'sigmoid',  // or '4pl-dose-response', 'linear', etc.
    showFitLine: true,
    showPoints: 'points',
    clickToToggle: true,
  };

  // If the format includes pre-fit parameters, set them:
  // series.parameters = [max, slope, IC50, min];

  const chartData: IFitChartData = {
    chartOptions: {
      xAxisName: 'Concentration',
      yAxisName: 'Response',
      logX: true,  // typical for dose-response
    },
    series: [series],
  };

  return JSON.stringify(chartData);
}
```

Key requirements for the converter function:
- **Signature**: `(value: string) => string` — takes raw cell value, returns native JSON string
- **Pure function**: no side effects, no async operations (rendering is synchronous)
- **Handle edge cases**: empty/null input, malformed data
- Map all relevant fields: points (x, y), outliers, fit parameters, axis labels, units

### Step 2: Register the Converter as a DG Function

In your package's `src/package.ts`, expose a DG function tagged with `meta.role: curveConverter` and `meta.curveFormat: '<format-name>'`. This function **returns** the sync converter function (it does not perform the conversion itself):

```typescript
import {convert<FormatName>ToJson} from './<format-name>-converter';

export class PackageFunctions {
  @grok.decorators.func({
    description: 'Returns <format-name> curve converter function',
    meta: {role: 'curveConverter', curveFormat: '<format-name>'}
  })
  static convert<FormatName>ToJsonFunc(): (value: string) => string {
    return convert<FormatName>ToJson;
  }
}
```

The `meta.curveFormat` value must match exactly the format name set in the detector tag.

### Step 3: Add a Detector

Add a semantic type detector in your package's `detectors.js`:

```javascript
const FIT_SEM_TYPE = 'fit';
const TAG_CURVE_FORMAT = '.%curve-format';

class <YourPackage>PackageDetectors extends DG.Package {
  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detect<FormatName>(col) {
    if (DG.Detector.sampleCategories(col, (s) => {
      // Return true if the string looks like this format
      // Use fast checks: string includes, startsWith, etc.
      // For JSON: try { const o = JSON.parse(s); return hasExpectedFields; } catch { return false; }
      return false;
    }, 1)) {
      col.setTag(TAG_CURVE_FORMAT, '<format-name>');
      col.semType = FIT_SEM_TYPE;
      return col.semType;
    }
    return null;
  }
}
```

Detection tips:
- Use `DG.Detector.sampleCategories(col, predicate, minCount)` to check sample values
- Keep detection fast — avoid expensive parsing in the detector
- Check for format-specific markers (XML tags, JSON field names, magic strings)
- The format name in `setTag` must match `meta.curveFormat` on the converter function exactly

### Step 4: Add Tests

Add tests to your package's test files:

```typescript
import {convert<FormatName>ToJson} from '../<format-name>-converter';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';

// Sample data in the new format
const SAMPLE_DATA = '...';

category('curve converter', () => {
  test('basicConversion', async () => {
    const result = convert<FormatName>ToJson(SAMPLE_DATA);
    const chartData: IFitChartData = JSON.parse(result);
    expect(chartData.series !== undefined && chartData.series.length > 0, true);
    expect(chartData.series![0].points.length > 0, true);
    // Verify specific point values, chart options, etc.
  });

  test('detection', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('col', [SAMPLE_DATA])]);
    const col = df.columns.byName('col');
    await grok.data.detectSemanticTypes(df);
    expect(col.semType, FitConstants.FIT_SEM_TYPE);
    expect(col.getTag(FitConstants.TAG_CURVE_FORMAT), '<format-name>');
  });
});
```

### Step 5: Build and Test

```bash
cd <YourPackage>
npm run build
grok test --host <host alias in config>
```

Once both your package and the Curves package are published to the same server, the Curves package will automatically discover and use your converter.

## Reference: Existing Converters in the Curves Package

For working examples, see these files in `packages/Curves/`:

| Format | Format Name | Converter | Detector |
|--------|-------------|-----------|----------|
| XML 3DX (`<chart>...</chart>`) | `3dx` | `src/fit/converters/xml-converter.ts` | `detectXMLCurveChart` |
| Compact dose-response JSON (`{p, poi, mxr, ...}`) | `compact-dr` | `src/fit/converters/compact-dr-converter.ts` | `detectCompactDoseResponse` |
| PZFX (GraphPad Prism XY) | `pzfx` | `src/fit/converters/pzfx-converter.ts` | `detectPzfxCurveChart` |

Note: the Curves package registers its own converters directly (hardcoded) to avoid circular dependencies. External packages rely on auto-discovery via `meta.curveFormat`.

## Reference: Key Types

From `@datagrok-libraries/statistics/src/fit/fit-curve.ts`:

```typescript
interface IFitPoint {
  x: number;
  y: number;
  outlier?: boolean;
  color?: string;
  marker?: FitMarkerType;
  size?: number;
  stdev?: number;
}

interface IFitSeries {
  points: IFitPoint[];
  name?: string;
  fitFunction?: string;           // 'sigmoid', '4pl-dose-response', 'linear', etc.
  parameters?: number[];          // [max, slope, IC50, min] for sigmoid/4PL
  showFitLine?: boolean;
  showPoints?: string;            // 'points' to show
  clickToToggle?: boolean;        // allow outlier toggling
  pointColor?: string;
  fitLineColor?: string;
  droplines?: string[];           // e.g. ['IC50']
  // ... more options available
}

interface IFitChartOptions {
  minX?: number; minY?: number;
  maxX?: number; maxY?: number;
  xAxisName?: string;
  yAxisName?: string;
  logX?: boolean;
  logY?: boolean;
  title?: string;
  // ... more options available
}
```

From `@datagrok-libraries/statistics/src/fit/const.ts`:
- `FitConstants.TAG_CURVE_FORMAT` = `'.%curve-format'` — the column tag to set
- `FitConstants.FIT_SEM_TYPE` = `'fit'` — the semantic type to assign

Available fit functions: `'sigmoid'`, `'4pl-dose-response'`, `'4pl-regression'`, `'linear'`, `'log-linear'`, `'exponential'`, or custom via `IFitFunctionDescription`.