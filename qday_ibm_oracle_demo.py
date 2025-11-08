#!/usr/bin/env python3
"""
QDay tiny-curve demo — two modes in one script

Modes:
  1) --mode oracle  (default)
     Build σ_t(u) over G only, pack residues for moduli across horizon rows,
     run on Aer or IBM, print packed bitstrings.

  2) --mode shor
     Two-address ROM oracle f(a,b)=x(aG+bQ) mod m.
     Pipeline: H⊗H over (a,b) → ROM → measure output → QFT on a,b → measure (r,s)
     Recover d from samples using robust voting across four valid relations:
       r ± s·d ≡ 0 (mod n)  and  s ± r·d ≡ 0 (mod n),
     dropping r≡0 or s≡0 samples to avoid zero-frequency bias,
     and verifying candidates by checking d*G==Q on the curve.

Extras:
  * --save-report PATH.json  → writes counts, votes, recovered d, verification,
    circuit stats (depth, 1Q/2Q gates), params, backend, and job_id (if available).
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import Counter
from datetime import datetime

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile
from qiskit.circuit.library import MCXGate

# Optional local simulator
try:
    from qiskit_aer import Aer
    HAS_AER = True
except Exception:
    HAS_AER = False

# ---------- Global knob for QFT approximation (0 = exact; 1 trims depth nicely on hardware) ----------
QFT_APPROX = 1


# -------------------- Finite-field + EC ops (y^2 = x^3 + 7 mod p) --------------------

def inv_mod(a: int, m: int) -> int:
    a %= m
    if a == 0:
        raise ZeroDivisionError("no inverse mod m")
    t, nt, r, nr = 0, 1, m, a
    while nr:
        q = r // nr
        t, nt = nt, t - q*nt
        r, nr = nr, r - q*nr
    if r != 1:
        raise ZeroDivisionError("no inverse mod m")
    return t % m

def add_affine(P: Optional[Tuple[int,int]], Q: Optional[Tuple[int,int]], p: int) -> Optional[Tuple[int,int]]:
    if P is None: return Q
    if Q is None: return P
    x1,y1 = P; x2,y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P != Q:
        lam = ((y2 - y1) * inv_mod((x2 - x1) % p, p)) % p
    else:
        if y1 == 0:
            return None
        lam = ((3 * x1 * x1) * inv_mod((2 * y1) % p, p)) % p
    x3 = (lam*lam - x1 - x2) % p
    y3 = (lam*(x1 - x3) - y1) % p
    return (x3, y3)

def scalar_mul(d: int, G: Tuple[int,int], p: int) -> Optional[Tuple[int,int]]:
    d = int(d)
    R = None
    Q = G
    while d:
        if d & 1:
            R = add_affine(R, Q, p)
        Q = add_affine(Q, Q, p)
        d >>= 1
    return R

def order_of_point(G: Tuple[int,int], p: int, max_hint: int = 1 << 12) -> int:
    R = None
    for n in range(1, max_hint+1):
        R = add_affine(R, G, p)
        if R is None:
            return n
    raise ValueError("order_of_point exceeded hint")


# -------------------- QFT compatibility helper --------------------

def make_qft_gate(n: int, approx: Optional[int] = None):
    """
    Build a QFT gate compatible across Qiskit versions.
    If approx is None, use global QFT_APPROX.
    """
    if approx is None:
        approx = QFT_APPROX
    from qiskit.circuit.library import QFTGate
    try:
        return QFTGate(num_qubits=n, insert_swaps=False, approximation_degree=approx)
    except TypeError:
        try:
            return QFTGate(num_qubits=n, do_swaps=False, approximation_degree=approx)
        except TypeError:
            return QFTGate(num_qubits=n)


# -------------------- Curves JSON loader (accepts multiple schemas) --------------------

def load_qday_curves(path: Path) -> List[dict]:
    data = json.loads(path.read_text())
    norm = []
    for e in data:
        bits = e.get("bits", e.get("bit_length", e.get("nbits")))
        if bits is None:
            raise KeyError("missing bits/bit_length/nbits in curve entry")
        p = e.get("p", e.get("prime"))
        if p is None:
            raise KeyError("missing p/prime in curve entry")

        # G
        if "Gx" in e and "Gy" in e:
            Gx, Gy = int(e["Gx"]), int(e["Gy"])
        elif "generator_point" in e and isinstance(e["generator_point"], (list, tuple)) and len(e["generator_point"]) == 2:
            Gx, Gy = int(e["generator_point"][0]), int(e["generator_point"][1])
        elif "G" in e and isinstance(e["G"], (list, tuple)) and len(e["G"]) == 2:
            Gx, Gy = int(e["G"][0]), int(e["G"][1])
        else:
            raise KeyError("missing generator point (Gx/Gy or generator_point or G)")

        # Q
        Qx = Qy = None
        if "Qx" in e and "Qy" in e:
            Qx, Qy = int(e["Qx"]), int(e["Qy"])
        elif "public_key" in e and isinstance(e["public_key"], (list, tuple)) and len(e["public_key"]) == 2:
            Qx, Qy = int(e["public_key"][0]), int(e["public_key"][1])
        elif "Q" in e and isinstance(e["Q"], (list, tuple)) and len(e["Q"]) == 2:
            Qx, Qy = int(e["Q"][0]), int(e["Q"][1])

        d = e.get("d", e.get("secret_key"))

        norm.append({
            "bits": int(bits),
            "p": int(p),
            "Gx": Gx, "Gy": Gy,
            "Qx": Qx, "Qy": Qy,
            "d": int(d) if d is not None else None
        })
    return norm


# -------------------- ORACLE mode: σ_t(u) over G --------------------

def build_signature_table_multi(curve: dict, u_bits: int, moduli: List[int], t: int, base_pow2: int = 0) -> Tuple[Dict[int, int], int]:
    p = curve["p"]
    G = (curve["Gx"], curve["Gy"])

    # B = 2^base_pow2 · G
    B = G
    for _ in range(base_pow2):
        B = add_affine(B, B, p)

    bits_per_m = [(m-1).bit_length() for m in moduli]
    residue_offsets: List[int] = []
    total_bits = 0
    for _ in range(t):
        for bm in bits_per_m:
            residue_offsets.append(total_bits)
            total_bits += bm

    Umax = 1 << u_bits
    tbl: Dict[int, int] = {}

    Ru = None
    for u in range(Umax):
        packed = 0
        P = Ru
        for idx in range(t * len(moduli)):
            m_idx = idx % len(moduli)
            if m_idx == 0:
                P = Ru if P is None else add_affine(P, B, p)
            xmod = 0 if P is None else (P[0] % moduli[m_idx])
            off = residue_offsets[idx]
            bm = (moduli[m_idx]-1).bit_length()
            packed |= (xmod & ((1 << bm) - 1)) << off
        tbl[u] = packed

        Ru = G if Ru is None else add_affine(Ru, G, p)

    return tbl, total_bits

def build_oracle_circuit(u_bits: int, out_bits: int, table: Dict[int, int],
                         regname_u: str = "u", regname_out: str = "o", cregname: str = "m",
                         measure_u: bool = True, measure_out: bool = True) -> QuantumCircuit:
    qr_u = QuantumRegister(u_bits, regname_u)
    qr_o = QuantumRegister(out_bits, regname_out)
    cr_m = ClassicalRegister(out_bits + u_bits, cregname)
    qc = QuantumCircuit(qr_u, qr_o, cr_m, name="sigma_oracle")

    for u, packed in table.items():
        # flip zeros to make 1-controls
        for i in range(u_bits):
            if ((u >> i) & 1) == 0:
                qc.x(qr_u[i])

        # set out bits via MCX
        bitpos = 0
        val = packed
        while val:
            if val & 1:
                mcx = MCXGate(num_ctrl_qubits=u_bits)
                qc.append(mcx, list(qr_u) + [qr_o[bitpos]])
            val >>= 1
            bitpos += 1

        # unflip
        for i in range(u_bits):
            if ((u >> i) & 1) == 0:
                qc.x(qr_u[i])

    qc.barrier()
    if measure_out:
        for i in range(out_bits):
            qc.measure(qr_o[i], cr_m[i])
    if measure_u:
        for i in range(u_bits):
            qc.measure(qr_u[i], cr_m[out_bits + i])

    return qc


# -------------------- SHOR mode: f(a,b)=x(aG+bQ) + QFT and d recovery --------------------

def build_ab_signature_table(curve: dict, a_bits: int, b_bits: int, moduli: List[int], t: int, base_pow2: int = 0) -> Tuple[Dict[Tuple[int,int], int], int]:
    p = curve["p"]
    G = (curve["Gx"], curve["Gy"])
    if curve.get("Qx") is None:
        raise ValueError("Public key Q is required for --mode shor.")
    Q = (curve["Qx"], curve["Qy"])

    # B = 2^base_pow2 · G
    B = G
    for _ in range(base_pow2):
        B = add_affine(B, B, p)

    A, Bn = 1 << a_bits, 1 << b_bits

    # Precompute ladders
    aG = [None] * A
    bQ = [None] * Bn
    R = None
    for a in range(A):
        aG[a] = R
        R = add_affine(R, G, p)
    S = None
    for b in range(Bn):
        bQ[b] = S
        S = add_affine(S, Q, p)

    bits_per_m = [(m-1).bit_length() for m in moduli]
    residue_offsets, total_bits = [], 0
    for _ in range(t):
        for bm in bits_per_m:
            residue_offsets.append(total_bits)
            total_bits += bm

    tbl: Dict[Tuple[int,int], int] = {}
    for a in range(A):
        for b in range(Bn):
            packed = 0
            P0 = add_affine(aG[a], bQ[b], p)
            P = P0
            for ell in range(t):
                for mi, m in enumerate(moduli):
                    bm = bits_per_m[mi]
                    off = residue_offsets[ell*len(moduli)+mi]
                    xmod = 0 if P is None else (P[0] % m)
                    packed |= (xmod & ((1 << bm) - 1)) << off
                P = add_affine(P, B, p)
            tbl[(a,b)] = packed
    return tbl, total_bits

def build_ab_oracle_circuit(a_bits: int, b_bits: int, out_bits: int, table: Dict[Tuple[int,int], int],
                            cregname: str = "m", measure_out_then_ab: bool = True) -> QuantumCircuit:
    qa = QuantumRegister(a_bits, "a")
    qb = QuantumRegister(b_bits, "b")
    qo = QuantumRegister(out_bits, "o")
    cr = ClassicalRegister(out_bits + a_bits + b_bits, cregname)
    qc = QuantumCircuit(qa, qb, qo, cr, name="f_ab_oracle")

    # Uniform superposition
    for i in range(a_bits): qc.h(qa[i])
    for j in range(b_bits): qc.h(qb[j])

    addr = list(qa) + list(qb)
    na = len(addr)

    for (a, b), packed in table.items():
        for i in range(a_bits):
            if ((a >> i) & 1) == 0: qc.x(qa[i])
        for j in range(b_bits):
            if ((b >> j) & 1) == 0: qc.x(qb[j])

        bitpos, val = 0, packed
        while val:
            if val & 1:
                mcx = MCXGate(num_ctrl_qubits=na)
                qc.append(mcx, addr + [qo[bitpos]])
            val >>= 1
            bitpos += 1

        for i in range(a_bits):
            if ((a >> i) & 1) == 0: qc.x(qa[i])
        for j in range(b_bits):
            if ((b >> j) & 1) == 0: qc.x(qb[j])

    if measure_out_then_ab:
        for k in range(out_bits):
            qc.measure(qo[k], cr[k])
        qfta = make_qft_gate(a_bits, approx=None)  # uses QFT_APPROX
        qftb = make_qft_gate(b_bits, approx=None)
        qc.append(qfta, list(qa))
        qc.append(qftb, list(qb))
        for i in range(a_bits):
            qc.measure(qa[i], cr[out_bits + i])
        for j in range(b_bits):
            qc.measure(qb[j], cr[out_bits + a_bits + j])

    return qc


# -------------------- Backends + reporting helpers --------------------

def circuit_stats_for_backend(qc, backend_name=None) -> Dict:
    """Transpile once for stats; returns dict with depth and gate counts."""
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
        backend = None
        if backend_name:
            svc = QiskitRuntimeService()
            backend = svc.backend(backend_name)
        tqc = transpile(qc, backend=backend, optimization_level=1)
        ops = tqc.count_ops()
        n2q = sum(ops.get(g, 0) for g in ["cx", "cz", "iswap", "ecr"])
        n1q = sum(ops.get(g, 0) for g in ["u", "u1", "u2", "u3", "rx", "ry", "rz"])
        return {
            "depth": int(tqc.depth()),
            "width_qubits": int(tqc.num_qubits),
            "ops": {k: int(v) for k, v in ops.items()},
            "two_qubit_gates": int(n2q),
            "one_qubit_gates": int(n1q),
        }
    except Exception as e:
        return {"error": f"stats_failed: {type(e).__name__}: {e}"}

def serialize_counts(counts: Dict[str,int]) -> Dict[str,int]:
    return {str(k): int(v) for k, v in counts.items()}

def save_report_json(path: str, payload: Dict):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    print(f"[report] wrote {path}")

def run_local_aer(qc: QuantumCircuit, shots: int = 4096):
    backend = Aer.get_backend("aer_simulator")
    tqc = transpile(qc, basis_gates=["u", "cx", "id"], optimization_level=1)
    res = backend.run(tqc, shots=shots).result().get_counts()
    meta = {"backend_name": "aer_simulator", "job_id": None}
    return res, meta

def run_ibm_backend(qc: QuantumCircuit, backend_name: str, shots: int = 4000, dump: bool = False):
    from qiskit_ibm_runtime import QiskitRuntimeService, Batch
    from qiskit_ibm_runtime.options import SamplerOptions
    try:
        from qiskit_ibm_runtime import SamplerV2 as Sampler
    except Exception:
        from qiskit_ibm_runtime import Sampler

    service = QiskitRuntimeService()
    backend = service.backend(backend_name)
    tqc = transpile(qc, backend=backend, optimization_level=1)

    opts = SamplerOptions(default_shots=shots)
    with Batch(backend=backend) as batch:
        sampler = Sampler(mode=batch, options=opts)
        try:
            job = sampler.run([tqc])
        except Exception:
            job = sampler.run(tqc)
        result = job.result()

    # job id
    job_id = None
    try:
        job_id = job.job_id()
    except Exception:
        try:
            job_id = getattr(job, "job_id", None)
        except Exception:
            job_id = None
    meta = {"backend_name": backend_name, "job_id": job_id}

    # Primary: v2 join_data()
    try:
        pr0 = result[0]
        data = pr0.join_data()
        counts = data.get_counts()
        return ({k: int(v) for k, v in counts.items()}, meta)
    except Exception:
        pass

    # Alt: pub_results.m
    try:
        pr0 = getattr(result, "pub_results", [None])[0]
        if pr0 is not None:
            data = getattr(pr0, "data", None)
            m = getattr(data, "m", None)
            if m is not None:
                for meth in ("get_counts", "get_int_counts"):
                    f = getattr(m, meth, None)
                    if callable(f):
                        d = f()
                        if isinstance(d, dict) and d:
                            return ({k: int(v) for k, v in d.items()}, meta)
    except Exception:
        pass

    # Fallback: quasi_dists
    try:
        quasi = result.quasi_dists[0]
        counts = {}
        total = 0
        for bitstr, prob in quasi.items():
            c = int(round(prob * shots))
            counts[bitstr] = counts.get(bitstr, 0) + c
            total += c
        if total != shots and counts:
            kmax = max(counts, key=lambda kv: kv[1])
            counts[kmax] += (shots - total)
        return counts, meta
    except Exception:
        pass

    if dump:
        try:
            payload = getattr(result, "to_dict", lambda: {"repr": repr(result)})()
            with open("ibm_result_dump.json", "w", encoding="utf-8") as f:
                json.dump(payload, f, ensure_ascii=False, indent=2)
            print("[debug] wrote ibm_result_dump.json", flush=True)
        except Exception:
            pass

    raise RuntimeError("IBM result didn’t expose counts; ensure a single classical register 'm' is used.")


# -------------------- Decoding & recovery --------------------

def int_from_slice_little_endian(key: str, start: int, width: int) -> int:
    s = key[::-1]
    return int(s[start:start+width][::-1], 2) if width > 0 else 0

def unpack_rows(packed: int, moduli: List[int], t: int) -> List[List[int]]:
    bits_per_m = [(m - 1).bit_length() for m in moduli]
    residue_offsets, total_bits = [], 0
    for _ in range(t):
        for bm in bits_per_m:
            residue_offsets.append(total_bits)
            total_bits += bm
    rows = []
    for ell in range(t):
        row = []
        for mi, m in enumerate(moduli):
            bm = bits_per_m[mi]
            off = residue_offsets[ell*len(moduli) + mi]
            r = (packed >> off) & ((1 << bm) - 1)
            row.append(int(r))
        rows.append(row)
    return rows

def recover_d_from_counts(counts: Dict[str,int], out_bits: int, a_bits: int, b_bits: int, n: int,
                          *, G: Optional[Tuple[int,int]] = None,
                          Q: Optional[Tuple[int,int]] = None,
                          p: Optional[int] = None) -> Tuple[Optional[int], Counter, str]:
    """
    Robust recovery:
      - drops r==0 or s==0 samples (mod n),
      - tallies all four congruences (r±s d≡0, s±r d≡0),
      - if (G,Q,p) provided, verifies d via d*G==Q and filters counts accordingly,
      - returns (best_d, best_votes_counter, label_of_relation).
    """
    def safe_inv(x, mod):
        try:
            return inv_mod(x, mod)
        except ZeroDivisionError:
            return None

    votes = {
        "(-) r + s·d ≡ 0  | d = -r s^{-1}": Counter(),
        "(+) r - s·d ≡ 0  | d =  r s^{-1}": Counter(),
        "(-) s + r·d ≡ 0  | d = -s r^{-1}": Counter(),
        "(+) s - r·d ≡ 0  | d =  s r^{-1}": Counter(),
    }

    for key, cnt in counts.items():
        r = int_from_slice_little_endian(key, out_bits, a_bits) % n
        s = int_from_slice_little_endian(key, out_bits + a_bits, b_bits) % n
        if r == 0 or s == 0:
            continue
        s_inv = safe_inv(s, n)
        r_inv = safe_inv(r, n)
        if s_inv is not None:
            votes["(-) r + s·d ≡ 0  | d = -r s^{-1}"][(-r * s_inv) % n] += cnt
            votes["(+) r - s·d ≡ 0  | d =  r s^{-1}"][( r * s_inv) % n] += cnt
        if r_inv is not None:
            votes["(-) s + r·d ≡ 0  | d = -s r^{-1}"][(-s * r_inv) % n] += cnt
            votes["(+) s - r·d ≡ 0  | d =  s r^{-1}"][( s * r_inv) % n] += cnt

    def verify_filter(vcount: Counter) -> Counter:
        if G is None or Q is None or p is None:
            return vcount
        filtered = Counter()
        for d, c in vcount.items():
            P = scalar_mul(d % n, G, p)
            if P == Q:
                filtered[d] += c
        return filtered if sum(filtered.values()) > 0 else vcount

    best_label, best_d, best_votes = None, None, Counter()
    for label, vc in votes.items():
        vc2 = verify_filter(vc)
        if vc2:
            d0, c0 = vc2.most_common(1)[0]
            if best_label is None or c0 > best_votes.get(d0, 0):
                best_label, best_d, best_votes = label, d0, vc2

    return best_d, best_votes, (best_label or "")


# -------------------- CLI / main --------------------

def main():
    ap = argparse.ArgumentParser(description="QDay tiny-curve demo: oracle or shor mode, Aer or IBM Runtime v2")
    ap.add_argument("--mode", choices=["oracle", "shor"], default="oracle", help="which mode to run")
    ap.add_argument("--curves", type=str, default="curves_qday.json", help="path to QDay curves JSON")
    ap.add_argument("--max-curves", type=int, default=2, help="[oracle] use up to this many tiny curves (bits ≤ 8)")
    ap.add_argument("--u-bits", type=int, default=4, help="[oracle] address bits (U=2^u_bits)")
    ap.add_argument("--curve-index", type=int, default=0, help="[shor] which curve entry to use")
    ap.add_argument("--a-bits", type=int, default=None, help="[shor] override a register bits (default ceil(log2 n))")
    ap.add_argument("--b-bits", type=int, default=None, help="[shor] override b register bits (default ceil(log2 n))")
    ap.add_argument("--moduli", type=str, default="32,5", help="comma-separated moduli (e.g., 32,5)")
    ap.add_argument("--horizon", type=int, default=1, help="rows t (ℓ=0..t-1); keep small for hardware")
    ap.add_argument("--base-pow2", type=int, default=0, help="B = 2^base_pow2 · G (ell step)")
    ap.add_argument("--shots", type=int, default=4000)
    ap.add_argument("--print-params", action="store_true", help="[shor] print curve/order/register sizes")
    ap.add_argument("--print-table", action="store_true", help="print packed residue table before running")
    ap.add_argument("--decode-top", action="store_true", help="[oracle] decode top outcomes into (u, rows)")
    ap.add_argument("--dump-ibm-result", action="store_true", help="dump IBM result json if parsing fails")
    ap.add_argument("--save-report", type=str, help="write a JSON report (counts, votes, recovered d, verification, circuit stats, params) to this path")

    group = ap.add_mutually_exclusive_group()
    group.add_argument("--aer", action="store_true", help="use local Aer simulator")
    group.add_argument("--backend", type=str, help="IBM backend name, e.g., ibm_fez")

    args = ap.parse_args()

    moduli = [int(x.strip()) for x in args.moduli.split(",") if x.strip()]
    if not moduli or any(m < 2 for m in moduli):
        raise ValueError("Provide valid moduli ≥ 2 (e.g., 32,5)")

    curves = load_qday_curves(Path(args.curves))

    # Target selection
    if args.backend:
        target = f"IBM backend '{args.backend}'"
    elif args.aer:
        target = "local Aer simulator"
    else:
        target = "local Aer simulator (default)"
        args.aer = True
    print(f"Target: {target}")

    # -------------------- ORACLE MODE --------------------
    if args.mode == "oracle":
        print(f"u_bits={args.u_bits} (U={1<<args.u_bits}), moduli={moduli}, horizon t={args.horizon}, base_pow2={args.base_pow2}, shots={args.shots}\n")
        tiny = [c for c in curves if c["bits"] <= 8]
        if not tiny:
            print("No tiny curves (≤8-bit) found.")
            return
        tiny = tiny[:args.max_curves]

        for idx, c in enumerate(tiny):
            print(f"--- Curve[{idx}] bits={c['bits']} p={c['p']} G=({c['Gx']},{c['Gy']}) ---")
            tbl, out_bits = build_signature_table_multi(c, args.u_bits, moduli, args.horizon, args.base_pow2)
            print(f"Packed signature out_bits = {out_bits}")

            if args.print_table:
                def unpack_rows_local(packed: int):
                    return unpack_rows(packed, moduli, args.horizon)
                print("u : [ [residues for m in moduli] per horizon ℓ ]")
                for u in range(1 << args.u_bits):
                    rows = unpack_rows_local(tbl.get(u, 0))
                    print(f"{u:0{args.u_bits}b} : {rows}")

            qc = build_oracle_circuit(args.u_bits, out_bits, tbl, regname_u="u", regname_out="o", cregname="m", measure_u=True, measure_out=True)

            # Run
            if args.backend:
                counts, meta = run_ibm_backend(qc, args.backend, args.shots, dump=args.dump_ibm_result)
            else:
                counts, meta = run_local_aer(qc, args.shots)

            # Show top outcomes
            top = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:8]
            for bitstr, cnt in top:
                print(f"{bitstr}: {cnt}")

            # Optional decode
            if args.decode_top:
                print("\nDecoded top results:")
                def decode_key(key: str, out_bits: int, u_bits: int):
                    s = key[::-1]
                    packed_bits = s[:out_bits][::-1]
                    u_bits_str  = s[out_bits:out_bits+u_bits][::-1]
                    return int(u_bits_str, 2), int(packed_bits, 2)
                for bitstr, cnt in top:
                    u_val, packed = decode_key(bitstr, out_bits, args.u_bits)
                    rows = unpack_rows(packed, moduli, args.horizon)
                    print(f"u={u_val:0{args.u_bits}b}  rows={rows}  cnt={cnt}")

            # Report?
            if args.save_report:
                stats = circuit_stats_for_backend(qc, backend_name=args.backend if args.backend else None)
                report = {
                    "timestamp_utc": datetime.utcnow().isoformat() + "Z",
                    "mode": "oracle",
                    "backend": meta.get("backend_name"),
                    "job_id": meta.get("job_id"),
                    "params": {
                        "u_bits": int(args.u_bits),
                        "moduli": [int(m) for m in moduli],
                        "horizon": int(args.horizon),
                        "base_pow2": int(args.base_pow2),
                        "shots": int(args.shots),
                        "qft_approx": int(QFT_APPROX),
                        "transpile_optimization_level": 1
                    },
                    "results": {
                        "counts": serialize_counts(counts),
                        "top8": [{ "bitstring": k, "count": v } for k, v in top]
                    },
                    "circuit_stats": stats
                }
                save_report_json(args.save_report, report)

    # -------------------- SHOR MODE --------------------
    else:
        c = curves[args.curve_index]
        p = c["p"]; G = (c["Gx"], c["Gy"])
        if c["Qx"] is None:
            raise ValueError("Selected curve lacks public key Q; required for --mode shor.")
        Q = (c["Qx"], c["Qy"])
        n = order_of_point(G, p, max_hint=1<<12)
        a_bits = args.a_bits if args.a_bits is not None else max(1, (n - 1).bit_length())
        b_bits = args.b_bits if args.b_bits is not None else a_bits

        if args.print_params:
            print(f"Curve[{args.curve_index}] bits={c['bits']} p={p} G={G} Q={Q}  (order n={n})")
            if c.get("d") is not None:
                print(f"(reference d={c['d']})")
            print(f"Registers: a_bits={a_bits}, b_bits={b_bits}")
        print(f"moduli={moduli}, horizon t={args.horizon}, base_pow2={args.base_pow2}, shots={args.shots}")

        tbl, out_bits = build_ab_signature_table(c, a_bits, b_bits, moduli, args.horizon, args.base_pow2)
        if args.print_table:
            print(f"Packed signature out_bits = {out_bits}")
            print("a,b : [ [residues for m in moduli] per horizon ℓ ]")
            for a in range(1 << a_bits):
                for b in range(1 << b_bits):
                    rows = unpack_rows(tbl[(a,b)], moduli, args.horizon)
                    print(f"{a:0{a_bits}b},{b:0{b_bits}b} : {rows}")

        qc = build_ab_oracle_circuit(a_bits, b_bits, out_bits, tbl, cregname="m", measure_out_then_ab=True)

        # Run
        if args.backend:
            counts, meta = run_ibm_backend(qc, args.backend, args.shots, dump=args.dump_ibm_result)
        else:
            counts, meta = run_local_aer(qc, args.shots)

        stats = circuit_stats_for_backend(qc, backend_name=args.backend if args.backend else None)

        print("\nTop raw outcomes:")
        for k, v in sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:10]:
            print(f"{k}: {v}")

        d_hat, votes, relation = recover_d_from_counts(
            counts, out_bits, a_bits, b_bits, n, G=G, Q=Q, p=p
        )
        print(f"\nWinning relation: {relation}")
        print("Candidate d votes (post-verification if available):")
        for dval, cnt in votes.most_common():
            print(f"d={dval}: {cnt}")

        verified = False
        if d_hat is not None:
            print(f"\nRecovered d = {d_hat} (mod n={n})")
            okP = scalar_mul(d_hat % n, G, p)
            verified = (okP == Q)
            print("Check d*G==Q? ", verified)
        else:
            print("\nCould not recover d from this sample; try more shots or different moduli/horizon.")

        # Report?
        if args.save_report:
            report = {
                "timestamp_utc": datetime.utcnow().isoformat() + "Z",
                "mode": "shor",
                "backend": meta.get("backend_name"),
                "job_id": meta.get("job_id"),
                "params": {
                    "curve_index": int(args.curve_index),
                    "p": int(p), "G": [int(G[0]), int(G[1])], "Q": [int(Q[0]), int(Q[1])],
                    "order_n": int(n),
                    "a_bits": int(a_bits), "b_bits": int(b_bits),
                    "moduli": [int(m) for m in moduli],
                    "horizon": int(args.horizon),
                    "base_pow2": int(args.base_pow2),
                    "shots": int(args.shots),
                    "qft_approx": int(QFT_APPROX),
                    "transpile_optimization_level": 1
                },
                "results": {
                    "counts": serialize_counts(counts),
                    "winning_relation": relation,
                    "recovered_d": int(d_hat) if d_hat is not None else None,
                    "verified_dG_equals_Q": bool(verified),
                    "votes": {str(k): int(v) for k, v in votes.items()}
                },
                "circuit_stats": stats
            }
            save_report_json(args.save_report, report)

if __name__ == "__main__":
    main()
