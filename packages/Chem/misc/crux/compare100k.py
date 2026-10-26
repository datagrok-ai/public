"""compare100k.mjs's Python side: the RDKit-mode references Chem runs in Python, and the data subsets, one process, one
thread each (no multiprocessing: one Python RDKit worker once took 9.9 GB).

  python compare100k.py subsets <smiles.json> <out dir>      R-group series, FMCS clusters, similar pairs, molzip inputs
  python compare100k.py descriptors <molecules.json> <names.json> <out prefix>
                                                             dockerfiles/chem/chem.py's get_descriptors, in chunks of 1000
                                                             (the container's own code, single process): <prefix>.json
                                                             (names, types, times) + <prefix>.f64 (values, NaN = null)
  python compare100k.py props <molecules.json> <out.json>    Chem's property table on RDKit's definitions (crux-core
                                                             tools/chem-wasm/chem_props_oracle.py's PROPS)
  python compare100k.py fcfp <molecules.json> <out prefix>   Morgan r2/2048 with feature invariants (FCFP), the bytes as
                                                             RDKit.js packs them: <prefix>.bin (256 B per molecule) +
                                                             <prefix>.json (read mask, time)
  python compare100k.py mces <pairs.json> <options.json> <out.json>   Rascal FindMCES per pair (Chem's options)
  python compare100k.py molzip <inputs.json> <out.json>      Chem.molzip of each multi-fragment SMILES
  python compare100k.py scaffoldtree <smiles.json> <out.json>  Chem's scripts/scaffold_tree.py itself (exec'd), under a
                                                             Python with ScaffoldGraphDG 1.1.8 (SG_PYTHON)
Every result carries its wall time (`ms`, the call only: reading the input and imports excluded) and the first items'
warm-up (excluded) where there is one.
"""
import json
import struct
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHEM = HERE.parents[1]
WARMUP = 50


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def dump(path, obj):
    Path(path).write_text(json.dumps(obj), encoding="utf-8")


def rdkit_quiet():
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")


# ---- subsets --------------------------------------------------------------------------------------------------------

def subsets(smiles_path, out_dir):
    """R-group series and FMCS clusters: molecules grouped by their exact Bemis-Murcko scaffold (RDKit's
    MurckoScaffold); similar pairs: for 1,000 evenly spaced molecules, the most similar other molecule (Morgan r2/2048
    Tanimoto, a different canonical SMILES), kept at >= 0.6; molzip inputs: each molecule cut at one acyclic single
    bond between heavy atoms (the middle one of RDKit's bond order), the two pieces labelled [*:1], and a second set
    cut at two such bonds ([*:1], [*:2])."""
    rdkit_quiet()
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    smiles = load(smiles_path)
    t0 = time.perf_counter()
    mols = [Chem.MolFromSmiles(s) if s else None for s in smiles]
    groups = {}
    for i, m in enumerate(mols):
        if m is None:
            continue
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
        except Exception:  # noqa: BLE001
            continue
        if scaf:
            groups.setdefault(scaf, []).append(i)
    ranked = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))

    def rings(s):
        m = Chem.MolFromSmiles(s)
        return m.GetRingInfo().NumRings() if m else 0
    # R-group series: the 50 largest scaffolds with two rings or more, up to 500 molecules each
    series = [{"core": s, "rows": rows[:500]} for s, rows in ranked if rings(s) >= 2][:50]
    # FMCS clusters: the 300 largest scaffolds (any ring count) with 3+ molecules, up to 50 molecules each
    clusters = [{"scaffold": s, "rows": rows[:50]} for s, rows in ranked if len(rows) >= 3][:300]
    # similar pairs
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    keep = [i for i, m in enumerate(mols) if m is not None]
    fps = [gen.GetFingerprint(mols[i]) for i in keep]
    can = [Chem.MolToSmiles(mols[i]) for i in keep]
    pairs = []
    step = max(1, len(keep) // 1000)
    for k in range(0, len(keep), step)[:1000]:
        sims = DataStructs.BulkTanimotoSimilarity(fps[k], fps)
        best, bj = -1.0, None
        for j, s in enumerate(sims):
            if j != k and can[j] != can[k] and s > best:
                best, bj = s, j
        if bj is not None and best >= 0.6:
            pairs.append({"a": keep[k], "b": keep[bj], "similarity": round(best, 4)})
    # molzip inputs

    def cut(m, n):
        bonds = [b.GetIdx() for b in m.GetBonds() if not b.IsInRing() and b.GetBondType() == Chem.BondType.SINGLE
                 and b.GetBeginAtom().GetAtomicNum() > 1 and b.GetEndAtom().GetAtomicNum() > 1]
        if len(bonds) < n:
            return None
        picked = [bonds[len(bonds) // 2]] if n == 1 else [bonds[len(bonds) // 3], bonds[(2 * len(bonds)) // 3]]
        if len(set(picked)) < n:
            return None
        f = Chem.FragmentOnBonds(m, picked, addDummies=True, dummyLabels=[(k + 1, k + 1) for k in range(n)])
        for a in f.GetAtoms():
            if a.GetAtomicNum() == 0 and a.GetIsotope():
                a.SetAtomMapNum(a.GetIsotope())
                a.SetIsotope(0)
        return Chem.MolToSmiles(f)
    zips = []
    for i, m in enumerate(mols):
        if m is None:
            continue
        for n in (1, 2):
            s = cut(m, n)
            if s:
                zips.append({"row": i, "cuts": n, "smiles": s})
    dump(out / "rgroup-series.json", series)
    dump(out / "fmcs-clusters.json", clusters)
    dump(out / "mces-pairs.json", pairs)
    dump(out / "molzip-inputs.json", zips)
    print(json.dumps({"scaffolds": len(groups), "series": len(series),
                      "series_molecules": sum(len(s["rows"]) for s in series), "clusters": len(clusters),
                      "cluster_molecules": sum(len(c["rows"]) for c in clusters), "pairs": len(pairs),
                      "molzip_inputs": len(zips), "s": round(time.perf_counter() - t0, 1)}))


# ---- descriptors: chem.py -------------------------------------------------------------------------------------------

def descriptors(molecules_path, names_path, prefix):
    sys.path.insert(0, str(CHEM / "dockerfiles" / "chem"))
    import numpy as np
    import chem  # noqa: E402  (dockerfiles/chem/chem.py, as the container imports it)
    rdkit_quiet()
    molecules = load(molecules_path)
    names = load(names_path)
    n = len(molecules)
    chem.get_descriptors(molecules[:WARMUP], names)  # warm-up (excluded)
    values = np.full((len(names), n), np.nan)
    types = {}
    t0 = time.perf_counter()
    for start in range(0, n, 1000):
        chunk = molecules[start:start + 1000]
        res = chem.get_descriptors(chunk, names)
        for k, name in enumerate(names):
            types[name] = res[name]["type"]
            col = res[name]["value"]
            for j, v in enumerate(col):
                if v is not None:
                    values[k, start + j] = float(v)
    ms = (time.perf_counter() - t0) * 1000
    Path(prefix + ".f64").write_bytes(values.astype("<f8").tobytes())
    dump(prefix + ".json", {"names": names, "types": types, "n": n, "ms": ms, "warmup": WARMUP})
    print(json.dumps({"n": n, "names": len(names), "ms": round(ms)}))


# ---- the property table on RDKit's definitions ----------------------------------------------------------------------

def props(molecules_path, out_path):
    rdkit_quiet()
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    esol_c = (0.26121066137801696, -0.7416739523408995, -0.0066138847738667125, 0.003451545565957996,
              -0.42624840441316975)

    def aromatic_proportion(m):
        heavy = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1)
        return sum(1 for a in m.GetAtoms() if a.GetIsAromatic()) / heavy if heavy else 0.0

    def esol(m):
        c = esol_c
        return (c[0] + c[1] * Crippen.MolLogP(m) + c[2] * Descriptors.MolWt(m)
                + c[3] * Descriptors.NumRotatableBonds(m) + c[4] * aromatic_proportion(m))
    defs = {
        "MW": Descriptors.MolWt, "HBA": rdMolDescriptors.CalcNumLipinskiHBA, "HBD": rdMolDescriptors.CalcNumHBD,
        "LogP": Crippen.MolLogP, "LogS": esol, "PSA": lambda m: rdMolDescriptors.CalcTPSA(m, includeSandP=True),
        "Rotatable bonds": Descriptors.NumRotatableBonds,
        "Stereo centers": rdMolDescriptors.CalcNumAtomStereoCenters,
        "Molecule charge": Chem.GetFormalCharge, "Molecular formula": rdMolDescriptors.CalcMolFormula,
    }
    molecules = load(molecules_path)

    def one(s):
        m = (Chem.MolFromMolBlock(s) if "M  END" in s else Chem.MolFromSmiles(s)) if s else None
        return None if m is None else {k: f(m) for k, f in defs.items()}
    for s in molecules[:WARMUP]:
        one(s)
    res = {k: [None] * len(molecules) for k in defs}
    t0 = time.perf_counter()
    for i, s in enumerate(molecules):
        v = one(s)
        if v:
            for k in defs:
                res[k][i] = v[k]
    ms = (time.perf_counter() - t0) * 1000
    dump(out_path, {"values": res, "ms": ms, "warmup": WARMUP})
    print(json.dumps({"n": len(molecules), "ms": round(ms)}))


# ---- FCFP -----------------------------------------------------------------------------------------------------------

def fcfp(molecules_path, prefix):
    """Morgan radius 2, 2048 bits, feature invariants, after RemoveAllHs (as Chem's worker: remove_hs_in_place, which
    in RDKit.js removes isotopic hydrogens too, then get_morgan_fp_as_uint8array({"radius":2,"nBits":2048,
    "useFeatures":true})); bit i in byte i >> 3 at 1 << (i & 7)."""
    rdkit_quiet()
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    Chem.SetUseLegacyStereoPerception(False)  # as Chem's workers
    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen())
    molecules = load(molecules_path)

    def one(s):
        m = Chem.MolFromSmiles(s) if s else None
        if m is None:
            return None
        fp = gen.GetFingerprint(Chem.RemoveAllHs(m))
        b = bytearray(256)
        for i in fp.GetOnBits():
            b[i >> 3] |= 1 << (i & 7)
        return bytes(b)
    for s in molecules[:WARMUP]:
        one(s)
    out = bytearray()
    read = []
    t0 = time.perf_counter()
    for s in molecules:
        b = one(s)
        read.append(b is not None)
        out += b if b is not None else bytes(256)
    ms = (time.perf_counter() - t0) * 1000
    Path(prefix + ".bin").write_bytes(bytes(out))
    dump(prefix + ".json", {"read": read, "ms": ms, "warmup": WARMUP})
    print(json.dumps({"n": len(molecules), "ms": round(ms)}))


# ---- Rascal MCES ----------------------------------------------------------------------------------------------------

def bits(x):
    return "%016x" % struct.unpack("<Q", struct.pack("<d", float(x)))[0]


def mces(pairs_path, options_path, out_path):
    rdkit_quiet()
    from rdkit import Chem
    from rdkit.Chem import rdRascalMCES
    pairs = load(pairs_path)
    opts = load(options_path)

    def one(a, b):
        m1, m2 = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
        if m1 is None or m2 is None:
            return {"error": "parse"}
        o = rdRascalMCES.RascalOptions()
        for k, v in opts.items():
            setattr(o, k, v)
        t = time.perf_counter()
        res = rdRascalMCES.FindMCES(m1, m2, o)
        t = time.perf_counter() - t
        return {"s": t, "results": [{
            "smarts": r.smartsString, "atomMatches": [list(p) for p in r.atomMatches()],
            "bondMatches": [list(p) for p in r.bondMatches()], "similarity": bits(r.similarity),
            "tier1Sim": bits(r.tier1Sim), "tier2Sim": bits(r.tier2Sim), "timedOut": bool(r.timedOut),
            "numFragments": r.numFragments, "largestFragmentSize": r.largestFragmentSize} for r in res]}
    for p in pairs[:5]:
        one(p["smiles1"], p["smiles2"])
    out = []
    t0 = time.perf_counter()
    for p in pairs:
        out.append(one(p["smiles1"], p["smiles2"]))
    ms = (time.perf_counter() - t0) * 1000
    dump(out_path, {"items": out, "ms": ms, "warmup": 5})
    print(json.dumps({"n": len(pairs), "ms": round(ms)}))


# ---- molzip ---------------------------------------------------------------------------------------------------------

def molzip(inputs_path, out_path):
    """Chem.molzip(mol) of each input (its pieces as fragments, joined at [*:n]), the result's canonical SMILES read
    again (as Chem's zippedSmiles: molzip does not sanitize), '' for one it refuses."""
    rdkit_quiet()
    from rdkit import Chem
    Chem.SetUseLegacyStereoPerception(False)
    inputs = load(inputs_path)

    def one(s):
        m = Chem.MolFromSmiles(s)
        if m is None:
            return ""
        try:
            z = Chem.molzip(m)
            r = Chem.MolFromSmiles(Chem.MolToSmiles(z))
            return Chem.MolToSmiles(r) if r is not None else ""
        except Exception:  # noqa: BLE001
            return ""
    for x in inputs[:WARMUP]:
        one(x["smiles"])
    t0 = time.perf_counter()
    res = [one(x["smiles"]) for x in inputs]
    ms = (time.perf_counter() - t0) * 1000
    dump(out_path, {"smiles": res, "ms": ms, "warmup": WARMUP})
    print(json.dumps({"n": len(inputs), "ms": round(ms)}))


# ---- the scaffold tree: Chem's scripts/scaffold_tree.py itself ------------------------------------------------------

def scaffoldtree(smiles_path, out_path):
    import pandas as pd
    import scaffoldgraph as sg
    rdkit_quiet()
    assert sg.__version__ == "1.1.8", f"ScaffoldGraphDG 1.1.8 expected, got {sg.__version__}"
    smiles = load(smiles_path)
    script = (CHEM / "scripts" / "scaffold_tree.py").read_text(encoding="utf-8")
    # the script's inputs, as Datagrok binds them (#input: dataframe data, string smilesColumn, int ringCutoff = 10,
    # bool dischargeAndDeradicalize = false); `progress=True` writes a tqdm bar to stderr
    env = {"__name__": "scaffold_tree", "data": pd.DataFrame({"smiles": smiles}), "smilesColumn": "smiles",
           "ringCutoff": 10, "dischargeAndDeradicalize": False}
    t0 = time.perf_counter()
    try:
        exec(compile(script, "scaffold_tree.py", "exec"), env)  # noqa: S102 (Chem's own script)
    except Exception as e:  # noqa: BLE001 (the script's own exception at this size is an outcome: recorded)
        import traceback
        ms = (time.perf_counter() - t0) * 1000
        where = [f"{f.name}:{f.lineno}" for f in traceback.extract_tb(e.__traceback__) if f.filename == "scaffold_tree.py"]
        raised = f"{type(e).__name__}: {e}" + (f" (in {', '.join(where[-2:])})" if where else "")
        print(json.dumps({"n": len(smiles), "ms": round(ms), "raised": raised}))
        Path(out_path + ".json").write_text(json.dumps({"n": len(smiles), "ms": ms, "warmup": 0, "raised": raised}),
                                            encoding="utf-8")
        return
    ms = (time.perf_counter() - t0) * 1000
    Path(out_path).write_bytes(env["result"])
    print(json.dumps({"n": len(smiles), "ms": round(ms)}))
    Path(out_path + ".json").write_text(json.dumps({"n": len(smiles), "ms": ms, "warmup": 0}), encoding="utf-8")


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"subsets": subsets, "descriptors": descriptors, "props": props, "fcfp": fcfp, "mces": mces, "molzip": molzip,
     "scaffoldtree": scaffoldtree}[cmd](*args)
