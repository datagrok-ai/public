# 02 — File Formats

## FCS (Flow Cytometry Standard) — Primary Format

FCS is the universal binary format maintained by ISAC's Data Standards Task Force. **Every instrument and analysis tool supports it.** An FCS file is NOT a text file — it is a structured binary file with a small ASCII header.

### Binary Structure (All Versions)

```
┌─────────────────────────────────────────┐
│  HEADER  (58 bytes, ASCII)              │  Version + byte offsets to other segments
├─────────────────────────────────────────┤
│  TEXT segment  (ASCII key-value pairs)  │  All metadata as $KEYWORD/value pairs
├─────────────────────────────────────────┤
│  DATA segment  (binary)                 │  Raw event measurements
├─────────────────────────────────────────┤
│  ANALYSIS segment  (optional)           │  User-defined analysis results
└─────────────────────────────────────────┘
```

### Header Parsing

```
Bytes 0–5:   Version string, e.g. "FCS3.1" (ASCII, space-padded)
Bytes 6–9:   Spaces (always "    ")
Bytes 10–17: TEXT segment start offset (ASCII integer, space-padded)
Bytes 18–25: TEXT segment end offset
Bytes 26–33: DATA segment start offset (0 if >99,999,999 — check TEXT for $BEGINDATA)
Bytes 34–41: DATA segment end offset   (0 if >99,999,999 — check TEXT for $ENDDATA)
Bytes 42–49: ANALYSIS segment start offset (0 if none)
Bytes 50–57: ANALYSIS segment end offset
```

> **Critical:** For files >99,999,999 bytes, the DATA offsets in the header will be `0`. The actual offsets are in the TEXT segment as `$BEGINDATA` and `$ENDDATA` keywords. Always check TEXT first.

### TEXT Segment Parsing

The TEXT segment is a sequence of ASCII key-value pairs separated by a single **delimiter character** which is defined as the **first byte** of the TEXT segment.

```
// Example (delimiter = '/')
/FCS3.1/$TOT/12000/$PAR/8/$DATATYPE/F/$BYTEORD/1,2,3,4/...
```

Parsing algorithm:
1. Read `textStart` and `textEnd` from header
2. Read byte at `textStart` — this is the delimiter character
3. Split the remainder on the delimiter character
4. Odd-indexed tokens are keys (uppercase by convention), even-indexed are values
5. Build a `Map<string, string>` of all keywords

### Required TEXT Keywords

| Keyword | Description | Example |
|---|---|---|
| `$FIL` | Original filename | `sample001.fcs` |
| `$TOT` | Total number of events | `12000` |
| `$PAR` | Number of parameters | `8` |
| `$DATATYPE` | Data type: `I`=integer, `F`=float32, `D`=float64 | `F` |
| `$BYTEORD` | Byte order: `1,2,3,4`=little-endian, `4,3,2,1`=big-endian | `1,2,3,4` |
| `$MODE` | Data mode: `L`=list mode (always use this) | `L` |
| `$BEGINDATA` | DATA segment start byte (authoritative for large files) | `5892` |
| `$ENDDATA` | DATA segment end byte | `389892` |
| `$PnN` | Short name of parameter n (1-indexed) | `$P1N` = `FSC-H` |
| `$PnB` | Bits per value for parameter n | `$P1B` = `32` |
| `$PnR` | Range for parameter n | `$P1R` = `262144` |
| `$PnE` | Amplification type: `0,0`=linear, `decades,f1`=log | `$P1E` = `0,0` |

### Important Optional TEXT Keywords

| Keyword | Description |
|---|---|
| `$PnS` | Long/marker name for parameter n (e.g., "CD4") |
| `$SPILLOVER` | Spillover matrix (FCS 3.1+, preferred) |
| `$SPILL` | Spillover matrix (older non-standard but common) |
| `$CYT` | Cytometer name (e.g., "FACSCanto II") |
| `$DATE` | Acquisition date |
| `$BTIM` / `$ETIM` | Begin/end time of acquisition |
| `$SRC` | Sample source description |
| `$WELLID` | Well identifier for plate-based experiments (e.g., "A01") |
| `$PLATEID` | Plate identifier |
| `$GUID` | Globally unique identifier for the FCS file |

### DATA Segment Parsing

Once TEXT is parsed, read the DATA segment:

```typescript
// Pseudocode for DATA parsing
const dataStart = parseInt(keywords.get('$BEGINDATA'));
const dataEnd = parseInt(keywords.get('$ENDDATA'));
const nEvents = parseInt(keywords.get('$TOT'));
const nParams = parseInt(keywords.get('$PAR'));
const dataType = keywords.get('$DATATYPE'); // 'F', 'D', or 'I'
const byteOrder = keywords.get('$BYTEORD'); // '1,2,3,4' or '4,3,2,1'
const littleEndian = byteOrder === '1,2,3,4';

const dataView = new DataView(buffer, dataStart, dataEnd - dataStart + 1);
let offset = 0;

for (let event = 0; event < nEvents; event++) {
  for (let param = 0; param < nParams; param++) {
    let value;
    if (dataType === 'F') {
      value = dataView.getFloat32(offset, littleEndian);
      offset += 4;
    } else if (dataType === 'D') {
      value = dataView.getFloat64(offset, littleEndian);
      offset += 8;
    } else if (dataType === 'I') {
      // Integer: width determined by $PnB for parameter n
      const bits = parseInt(keywords.get(`$P${param+1}B`));
      value = readUint(dataView, offset, bits / 8, littleEndian);
      offset += bits / 8;
    }
    matrix[event][param] = value;
  }
}
```

> **Note on integers:** FCS integer data is **unsigned**. Values must be masked by `$PnR - 1` if `$PnR` is not a power of 2. In practice, most modern FCS files use `$DATATYPE/F` (32-bit float).

### FCS Version History

| Version | Year | Key additions relevant to parsing |
|---|---|---|
| **FCS 2.0** | 1990 | Basic structure; `$BYTEORD`; `$DFCiTOj` compensation keyword |
| **FCS 3.0** | 1997 | Large file support; supplemental TEXT; `$COMP` spillover; CRC |
| **FCS 3.1** | 2010 | `$SPILLOVER` replaces `$COMP` (includes parameter names); `$WELLID`/`$PLATEID`; `$PnD` display scale; `$ORIGINALITY` |
| **FCS 3.2** | 2021 | Per-parameter `$PnDATATYPE`; new detector keywords (`$PnDET`, `$PnTYPE`, `$PnFEATURE`, `$PnTAG`, `$PnANALYTE`); **breaks backward compat** |

> **Parser strategy:** Implement FCS 3.1 fully as the primary target. Add graceful FCS 2.0/3.0 fallbacks (handle `$COMP` and `$SPILL` in addition to `$SPILLOVER`). Add FCS 3.2 support in Phase 3.

### $SPILLOVER Keyword Format (FCS 3.1)

```
$SPILLOVER = n,P1name,P2name,...,Pnname,s11,s12,...,snn
```

- `n` = number of fluorescence parameters in the matrix
- `P1name...Pnname` = parameter names matching `$PnN` values
- `s11...snn` = n×n matrix values in row-major order, diagonal = 1.0

To get the **compensation matrix**: invert the spillover matrix: `M = S⁻¹`

Apply compensation per event: `compensated = M × raw_fluorescence_vector`

> **Common quirk:** Many instruments write `$SPILL` (without "OVER") as a non-standard alias. The parser MUST check both `$SPILLOVER` and `$SPILL`.

---

## GatingML 2.0 — Gate Interchange Format

GatingML is an ISAC XML standard for encoding gating strategies independent of any specific software.

### Gate Types in GatingML 2.0

| Gate type | XML element | Key attributes |
|---|---|---|
| Rectangle | `<gating:RectangleGate>` | Min/max per dimension |
| Polygon | `<gating:PolygonGate>` | List of (x,y) vertices |
| Ellipse | `<gating:EllipsoidGate>` | Center, covariance matrix |
| Quadrant | `<gating:QuadrantGate>` | Two divider lines, 4 quadrant definitions |
| Boolean | `<gating:BooleanGate>` | AND/OR/NOT references to other gates |

### Transforms in GatingML 2.0

| Transform | XML element | Parameters |
|---|---|---|
| Linear | `<transforms:linear>` | `T` (max value), `A` (offset) |
| Logarithmic | `<transforms:log>` | `T`, `M` (decades) |
| Logicle (biexponential) | `<transforms:logicle>` | `T`, `M`, `W`, `A` |
| Arcsinh | `<transforms:arcsinh>` | `T`, `M`, `A` (cofactor = T/sinh(M)) |
| Hyperlog | `<transforms:hyperlog>` | `T`, `M`, `W`, `A` |

### Minimal GatingML 2.0 Structure

```xml
<?xml version="1.0" encoding="UTF-8"?>
<gating:Gating-ML xmlns:gating="http://www.isac-net.org/std/Gating-ML/v2.0/gating"
                  xmlns:transforms="http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
                  xmlns:data-type="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes">

  <!-- Transformation definitions -->
  <transforms:transformation transforms:id="logicle_FITC">
    <transforms:logicle transforms:T="262144" transforms:M="4.5" transforms:W="0.5" transforms:A="0"/>
  </transforms:transformation>

  <!-- Gate definitions -->
  <gating:PolygonGate gating:id="lymphocytes" gating:parent_id="root">
    <gating:dimension gating:compensation-ref="FCS" gating:transformation-ref="logicle_FITC">
      <data-type:fcs-dimension data-type:name="FSC-A"/>
    </gating:dimension>
    <gating:dimension ...>
      <data-type:fcs-dimension data-type:name="SSC-A"/>
    </gating:dimension>
    <gating:vertex><gating:coordinate data-type:value="0.1"/><gating:coordinate data-type:value="0.2"/></gating:vertex>
    <!-- more vertices -->
  </gating:PolygonGate>

</gating:Gating-ML>
```

---

## FlowJo .wsp Workspace Files

FlowJo `.wsp` files are **XML** and do **not** embed raw data — they reference FCS files by path and store the analysis on top.

### Key .wsp XML Structure

```xml
<Workspace>
  <SampleList>
    <Sample>
      <DataSet uri="path/to/sample.fcs" />
      <Transformations>
        <biex T="262144" W="0.5" M="4.5" A="0" />
      </Transformations>
      <SampleNode>
        <Population name="Lymphocytes">
          <Gate gating:id="gate1">
            <PolygonGate xAxisName="FSC-A" yAxisName="SSC-A">
              <Vertex x="100" y="200"/>
              <!-- more vertices -->
            </PolygonGate>
          </Gate>
          <Population name="CD4+ T cells">
            <!-- nested gate -->
          </Population>
        </Population>
      </SampleNode>
    </Sample>
  </SampleList>
  <CompensationList>
    <Compensation name="Comp">
      <Matrix>...</Matrix>
    </Compensation>
  </CompensationList>
</Workspace>
```

> **Reference implementations for .wsp parsing:** FlowKit (Python) at `flowkit/_utils/wsp_utils.py`; CytoML (R/Bioconductor). These are the authoritative guides for the format's quirks.

---

## Other Relevant Formats

| Format | Description | Priority |
|---|---|---|
| `.lmd` | BD FACSDiva legacy list-mode data; essentially FCS 2.0 with vendor extensions | Low |
| `.xml` (FACSDiva experiment) | BD experiment structure file; links tubes to FCS files | Low |
| ACS (`.acs`) | ZIP container: FCS + GatingML + workspace + MIFlowCyt metadata | Phase 3 |
| CSV export | Statistics export from FlowJo/Cytobank; no event-level data | Import only |
| `.fex` | FCS Express proprietary workspace; avoid unless requested | Not planned |

---

## Data Concepts: Transforms

### Why Transforms Are Required

Raw FCS data spans 5+ decades of dynamic range (e.g., 0 to 262,144). Compensation can introduce **negative values** (cells with less signal than the unstained baseline). Log scale is undefined for negatives and distorts populations near zero.

### Logicle (Biexponential) Transform — Primary Recommendation

The ISAC-recommended transform. Behaves like log at high values, linear near zero, and handles negatives.

**Parameters:**
- `T` = top of scale (e.g., 262144 for 18-bit data)
- `M` = number of decades in the display range (typically 4.5)
- `W` = width of the linear region near zero (typically 0.5)
- `A` = additional decades below zero (typically 0)

**Implementation:** No JavaScript library exists. Port from:
- Java reference: Wayne Moore, Stanford (logicle.sourceforge.net)
- Python: `flowkit/_utils/transforms.py` — `_logicle_fn()` uses Newton-Raphson root finding

### Arcsinh Transform — For CyTOF/Mass Cytometry

`f(x) = arcsinh(x / cofactor)` where cofactor is typically **5** for CyTOF, **150** for spectral data.

This is trivial in JavaScript: `Math.asinh(x / cofactor)`.

### Scale Reference

| Transform | Good for | Handles negatives? |
|---|---|---|
| Log | Older data, simple 4-color panels | No |
| Logicle | All conventional fluorescence | Yes |
| Arcsinh | CyTOF/mass cytometry | Yes |
| Linear | Scatter (FSC/SSC), time | N/A |
