# QDay ECC-on-Quantum Demo (Shor-style, ROM oracle)

**Author:** Joseph “Joey” Fredrickson  
**Email:** <your.email@domain>  
**Background:** AI content developer, writer/novelist, mechanical engineer (BYU), quantum/computation researcher.  
**Repository purpose:** Demonstrate a verified Shor-style recovery of ECC private keys on **IBM Quantum** hardware using a compact ROM-oracle approach and robust post-processing.

---

## Results (Hardware-Verified)
- **Curve[1]** — `p=43, G=(34,3), Q=(21,25)`, order `n=31`  
  Recovered **d = 18 (mod 31)**, verification `d*G == Q` → **True**  
  Backend: **ibm_torino**, shots: **2048**, QFT approx: **1**  
  Reports: `run_curve1_2048.json` (and second confirmation `run_curve1_2048_b.json`)

- **Curve[0]** — `p=13, G=(11,5), Q=(11,8)`, order `n=7`  
  Recovered **d = 6 (mod 7)**, verification `d*G == Q` → **True**  
  Backend: **ibm_torino**, shots: **1024**, QFT approx: **1**  
  Report: `run_curve0_1024.json`

See `qday_bundle_summary.csv` for a consolidated table (backend, job id, shots, recovered d, verification, and circuit stats).

---

## Quantum Computer & Access
- **Model:** IBM Falcon (backend: `ibm_torino`)  
- **Access:** IBM Quantum via `qiskit-ibm-runtime` (Batch + SamplerV2)
- **Transpile:** optimization level **1** (depth-friendly, stable for tiny circuits)
- **QFT:** approximation degree **1** (reduces entangling depth)

> If you ran on a different backend (e.g., `ibm_fez`), note it here and include those job IDs in the reports folder.

---

## How It Works (Brief)
1. Prepare uniform superposition over address registers `(a,b)` (sizes `≈ log2(n)` each).  
2. Oracle writes packed residues of `x(aG + bQ)` mod small integers (e.g., `{32,7}`) into an output register; measure the output.
3. Apply QFT to `(a,b)`, then measure `(r,s)`.
4. Post-process:
   - Drop degenerate samples `r≡0` or `s≡0 (mod n)` to avoid the 2^k→n zero-frequency bias.
   - Tally **four** valid linear relations to resolve sign/wiring ambiguity:  
     `r ± s·d ≡ 0 (mod n)` and `s ± r·d ≡ 0 (mod n)`.
   - Pick the dominant relation and **verify** `d*G == Q` on-curve (tiny curves allow exact check).

---

## Reproduce the Runs (Windows CMD)
Make sure you have IBM credentials configured for `qiskit-ibm-runtime`.

```cmd
REM Curve 1 (n=31) — hardware + JSON report
py qday_ibm_oracle_demo.py --mode shor --backend ibm_torino --curve-index 1 --moduli 32,7 --horizon 1 --shots 2048 --print-params --save-report run_curve1_2048.json

REM (Optional) Second confirmation run
py qday_ibm_oracle_demo.py --mode shor --backend ibm_torino --curve-index 1 --moduli 32,7 --horizon 1 --shots 2048 --save-report run_curve1_2048_b.json

REM Curve 0 (n=7) — hardware + JSON report
py qday_ibm_oracle_demo.py --mode shor --backend ibm_torino --curve-index 0 --moduli 32,5 --horizon 1 --shots 1024 --print-params --save-report run_curve0_1024.json
```

Local sanity using Aer:
```cmd
py qday_ibm_oracle_demo.py --mode shor --aer --curve-index 1 --moduli 32,7 --horizon 1 --shots 1024 --print-params
```

---

## File Guide
- `qday_ibm_oracle_demo.py` — full script; oracle + shor modes; robust 4-relation voter; JSON report writer.
- `curves_qday.json` — tiny test curves with public keys; indexes 0 and 1 are used above.
- `run_curve*.json` — per-run artifacts: params, backend, job_id, counts, winning relation, recovered d, verification, and circuit stats.
- `qday_bundle_summary.csv` — consolidated table of the above.
- `brief.pdf` — two-page write-up (see `brief.md` source if you need to regenerate).

---

## Environment
- Python: 3.12 (adjust to your exact version)
- Packages: `qiskit`, `qiskit-aer`, `qiskit-ibm-runtime`
- Transpile opt: 1; QFT approx: 1

To list versions:
```cmd
python --version
pip show qiskit qiskit-aer qiskit-ibm-runtime
```

---

## Contact
Questions: **qdayprize@projecteleven.com**  
Author contact: <your.email@domain>
