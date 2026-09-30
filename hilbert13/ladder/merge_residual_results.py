"""Verify records, compact recovered permutations, and write the exact threshold table.

Run after residual_certificate.py. By default all twelve full cases are required.
--allow-incomplete writes an explicitly incomplete checkpoint instead.
The repository stores lossless delta32/XZ permutations; the original NPZ exports
can be regenerated, and the inverse permutation is determined by argsort(p).
"""
import argparse
import hashlib
import json
import lzma
from fractions import Fraction as F
from pathlib import Path
import numpy as np

EXPECTED = {(t,q) for t in [0,1,12,14] for q in ["Q0","Q1","Q2"]}
SOURCE_COMMIT = "dbcd7ea788d0b1d59b92bedb75e72af5b8af6969"


def fraction(record,name):
    v = record[name]
    return F(v["numerator"],v["denominator"])


def load_permutation(path,n):
    """Load the complete original permutation from the committed compact export."""
    raw = lzma.decompress(Path(path).read_bytes())
    d = np.frombuffer(raw,dtype="<i4")
    assert len(d) == n
    return np.cumsum(d,dtype=np.int64)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("residual_results"))
    parser.add_argument("--allow-incomplete",action="store_true")
    args = parser.parse_args()
    source_hashes = {name:hashlib.sha256(Path(name).read_bytes()).hexdigest()
                     for name in ["certify.py","orbifold.py","triples_data.py","a7.py","cover.py"]}
    records = []
    for path in sorted(args.out.glob("triple_*_Q?.json")):
        r = json.loads(path.read_text())
        assert r["status"] == "PASS"
        assert (r["triple"],r["which"]) in EXPECTED
        assert r["n"] == {"Q0":16,"Q1":96,"Q2":32}[r["which"]]
        assert fraction(r,"needed") == fraction(r,"delta")+fraction(r,"assembly_2")+fraction(r,"shift_rounding")
        assert fraction(r,"shift") > fraction(r,"needed")
        assert fraction(r,"sigma") >= fraction(r,"safe_sigma")
        assert fraction(r,"computed_Ch2_upper") <= fraction(r,"safe_Ch2")
        assert fraction(r,"spectral_bound") == fraction(r,"safe_sigma")/(1+fraction(r,"safe_sigma")*fraction(r,"safe_Ch2"))
        original = args.out/r.get("original_npz_file",r["permutation_file"])
        if original.suffix == ".npz" and original.exists():
            data = np.load(original)
            p = data["permutation"]
            pinv = data["inverse_permutation"]
            assert np.array_equal(p[pinv],np.arange(r["dofs"]))
            delta = np.diff(p,prepend=np.int64(0))
            assert delta.min() >= -(2**31) and delta.max() < 2**31
            compact = args.out/original.name.replace(".npz",".delta32.xz")
            compact.write_bytes(lzma.compress(delta.astype("<i4").tobytes(),preset=6))
            r["original_npz_file"] = original.name
            r["permutation_file"] = compact.name
            r["permutation_encoding"] = "XZ of little-endian signed int32 differences; prepend p[-1]=0; cumulative sum in int64"
        else:
            compact = args.out/r["permutation_file"]
        restored = load_permutation(compact,r["dofs"])
        assert np.array_equal(np.sort(restored),np.arange(r["dofs"]))
        sha = hashlib.sha256(restored.astype("<i8").tobytes()).hexdigest()
        assert sha == r["permutation_sha256_int64_le"]
        r["source_commit"] = SOURCE_COMMIT
        r["source_files_sha256"] = source_hashes
        path.write_text(json.dumps(r,indent=2)+"\n")
        records.append(r)
    seen = {(r["triple"],r["which"]) for r in records}
    assert len(seen) == len(records)
    missing = sorted(EXPECTED-seen)
    if not args.allow_incomplete:
        assert not missing, ("unfinished cases",missing)
    lines = ["# Independent residual certificate", "",
             f"Completed full cases: {len(records)}/12. " + ("All passed." if not missing else f"Missing: {missing}."), "",
             "Displayed residuals are rounded for readability. Exact fractions, input hashes, and permutation hashes are in the JSON records.", "",
             "| triple | quotient | dofs | k | rho | delta | assembly 2-norm bound | required shift |", 
             "|---:|---|---:|---:|---:|---:|---:|---:|"]
    for r in sorted(records,key=lambda r:(r["triple"],r["which"])):
        lines.append(f"| {r['triple']} | {r['which']} | {r['dofs']} | {r['rowmax']} | "
                     + " | ".join(f"{float(fraction(r,name)):.6e}" for name in ["rho","delta","assembly_2","needed"])+" |")
    lines += ["", "Every completed case satisfies required shift < 1e-9 using exact rational comparison.", "",
              "| quotient | safe sigma | safe C_h^2 | certified lower bound |", "|---|---:|---:|---:|"]
    for q in ["Q0","Q1","Q2"]:
        rs = [r for r in records if r["which"]==q]
        if rs:
            r = rs[0]
            lines.append(f"| {q} | {fraction(r,'safe_sigma')} | {fraction(r,'safe_Ch2')} | {fraction(r,'spectral_bound')} |")
    if not missing:
        bound = min(fraction(r,"spectral_bound") for r in records)
        assert bound > F(34089,100000)
        assert F(135,2)*bound > 23 and F(135,4)*bound > 11
        lines += ["",f"With the exact quotient-coverage check, lambda_1(C) > {bound} > 0.34089 for every marked triple class.",
                  "Hence gon(C) >= 24 and gon(C/<tau>) >= 12. This does not reach gon(C) >= 25."]
    lines += ["", "Compact permutations can be loaded with `load_permutation(path,n)` from this script; inverse = argsort(p).",
              "The NPZ files generated during the run are not committed because the compact files preserve the same permutations losslessly."]
    (args.out/"SUMMARY.md").write_text("\n".join(lines)+"\n")
    print(f"verified and packed {len(records)}/12 cases; missing {missing}")
    for f in sorted(args.out.glob("*.xz")):
        print(f.name,f.stat().st_size)


if __name__ == "__main__":
    main()
