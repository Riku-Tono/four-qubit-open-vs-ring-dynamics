#!/usr/bin/env python3
"""
4量子ビット・全64接続グラフの識別可能性の検証
================================================

問い: 局所観測量の時系列だけから、4サイト上の接続グラフ(ラベル付き)を
      どこまで一意に復元できるか。

3つの独立な方法で検証する。
  [A] 厳密対角化(浮動小数点)   : 時系列を直接比較し同値類を数える
  [B] 交換子展開(整数・有理数の厳密計算):
        <A>(t) = sum_n t^n/n! <psi0| D^n(A) |psi0>,  D(A) = i[H, A]
      のテイラー係数を Omega, J の多項式として厳密に求める
  [C] 対称性の議論(置換群の総当たり): 1体量が一致する十分条件

[B] は「区別できる」ことの厳密な証明(同値類数の下界)、
[C] は「区別できない」ことの厳密な証明(同値類数の上界)を与える。

---------------------------------------------------------------------------
【仮定】(このスクリプトが前提としていること。実測ではない)
  1. ハミルトニアン: H = Omega * sum_i X_i + J * sum_{(i,j) in G} Z_i Z_j
     (元リポジトリのベースライン部には式が明記されていないため、
      拡張部の「X場 + ZZボンド」という記述に基づく仮定)
  2. 全ボンドで J は共通。サイト数 N = 4。
  3. 既定値 Omega = 1, J = 7/10, 時間窓 0 <= t <= 3.688(200点)。
  4. 雑音なし・有限サンプル誤差なし。サイトのラベルは既知。
  5. 「同じ時系列」の判定: [A] は最大差 < TOL(1e-9)、[B] は係数の厳密一致。
---------------------------------------------------------------------------

依存: Python 3.9+, NumPy のみ。   実行: python3 verify_64_graphs.py
"""
import itertools
import math
from fractions import Fraction

import numpy as np

# ---------------------------------------------------------------- 仮定(パラメータ)
N = 4
OMEGA = Fraction(1)
J = Fraction(7, 10)
T_MAX, N_T = 3.688, 200
TOL = 1e-9
K_MAX = 12          # 交換子展開の最大次数(int64 の安全上限は 14)
SEED = 1            # ランダム積状態用

PAIRS = list(itertools.combinations(range(N), 2))            # 6 通りの可能なボンド
GRAPHS = [tuple(p for k, p in enumerate(PAIRS) if m >> k & 1) for m in range(64)]
OPEN = ((0, 1), (1, 2), (2, 3))
RING = ((0, 1), (0, 3), (1, 2), (2, 3))


def gname(g):
    return "{" + ",".join(f"{a}{b}" for a, b in g) + "}" if g else "{}"


# =====================================================================
# [A] 厳密対角化
# =====================================================================
I2 = np.eye(2, dtype=complex)
PAULI = {"X": np.array([[0, 1], [1, 0]], complex),
         "Y": np.array([[0, -1j], [1j, 0]]),
         "Z": np.diag([1, -1]).astype(complex)}


def embed(d):
    m = np.array([[1]], complex)
    for i in range(N):
        m = np.kron(m, d.get(i, I2))
    return m


ONE = {(a, i): embed({i: PAULI[a]}) for a in "XYZ" for i in range(N)}
TWO = {(a, p): embed({p[0]: PAULI[a], p[1]: PAULI[a]}) for a in "XYZ" for p in PAIRS}


def hamiltonian(graph, lam=None):
    """lam を与えると (0,3) ボンドだけ lam*J にする(弱リンク)。"""
    h = sum(float(OMEGA) * ONE[("X", i)] for i in range(N))
    for b in graph:
        c = float(J) * (lam if (lam is not None and b == (0, 3)) else 1.0)
        h = h + c * TWO[("Z", b)]
    return h


def evolve(h, ts, psi0):
    w, v = np.linalg.eigh(h)
    c = v.conj().T @ psi0
    return (v @ (np.exp(-1j * np.outer(w, ts)) * c[:, None])).T      # (n_t, 16)


def expval(o, ps):
    return np.einsum("ti,ij,tj->t", ps.conj(), o, ps).real


def features(graph, psi0, ts):
    ps = evolve(hamiltonian(graph), ts, psi0)
    e1 = {k: expval(o, ps) for k, o in ONE.items()}
    fz = np.stack([e1[("Z", i)] for i in range(N)])
    f1 = np.stack([e1[k] for k in ONE])
    f2 = np.stack([expval(TWO[(a, p)], ps) - e1[(a, p[0])] * e1[(a, p[1])] for (a, p) in TWO])
    return {"Z": fz, "1body": f1, "2body": f2}


def group(items, same):
    classes = []
    for i in items:
        for c in classes:
            if same(c[0], i):
                c.append(i)
                break
        else:
            classes.append([i])
    return classes


def basis_state(bits):
    psi = np.zeros(2 ** N, complex)
    psi[int(bits, 2)] = 1
    return psi


def random_product_state(seed):
    rng = np.random.default_rng(seed)
    psi = np.array([1], complex)
    for _ in range(N):
        psi = np.kron(psi, rng.normal(size=2) + 1j * rng.normal(size=2))
    return psi / np.linalg.norm(psi)


def part_A():
    print("=" * 76)
    print("[A] 厳密対角化: 64グラフの時系列の同値類   【実測】")
    print("=" * 76)
    ts = np.linspace(0, T_MAX, N_T)
    states = [("|0000>", basis_state("0000")), ("|1000>", basis_state("1000")),
              ("|1100>", basis_state("1100")), (f"random product(seed={SEED})", random_product_state(SEED))]
    print(f"{'初期状態':<28}{'<Z_i>のみ':>12}{'1体XYZ':>12}{'2体連結':>12}   (同値類の数 / 64)")
    result = {}
    for name, psi in states:
        F = [features(g, psi, ts) for g in GRAPHS]
        row = []
        for key in ("Z", "1body", "2body"):
            cl = group(range(64), lambda i, j: np.abs(F[i][key] - F[j][key]).max() < TOL)
            dmin = min(np.abs(F[i][key] - F[j][key]).max()
                       for c1, c2 in itertools.combinations(cl, 2) for i in c1[:1] for j in c2[:1])
            row.append((len(cl), dmin))
            result[(name, key)] = cl
        print(f"{name:<28}" + "".join(f"{n:>12d}" for n, _ in row)
              + "   異なる類の間の最小距離: " + ", ".join(f"{d:.3f}" for _, d in row))
    cl = result[("|0000>", "1body")]
    print("\n|0000> で 1体XYZ が一致するグラフの組(要素2以上の類):")
    for c in cl:
        if len(c) > 1:
            print("   " + "  ==  ".join(gname(GRAPHS[i]) for i in c))
    F = [features(g, basis_state("0000"), ts) for g in (OPEN, RING)]
    print(f"\n開鎖 vs リング (|0000>): max|d<Z_i>| = {np.abs(F[0]['Z'] - F[1]['Z']).max():.4f}, "
          f"max|d 連結相関| = {np.abs(F[0]['2body'] - F[1]['2body']).max():.4f}")
    return {frozenset(c) for c in cl}


# =====================================================================
# [B] 交換子展開(厳密)
#   パウリ文字列基底(4^N = 256 次元)で超演算子 D_P(Q) = i[P, Q] を整数行列に
#   する。エルミート演算子は実係数なので、複素数は不要。
# =====================================================================
# 1量子ビットの積 P*Q = i^k R   (0=I, 1=X, 2=Y, 3=Z)
_MUL = {(0, 0): (0, 0), (0, 1): (0, 1), (0, 2): (0, 2), (0, 3): (0, 3),
        (1, 0): (0, 1), (2, 0): (0, 2), (3, 0): (0, 3),
        (1, 1): (0, 0), (2, 2): (0, 0), (3, 3): (0, 0),
        (1, 2): (1, 3), (2, 3): (1, 1), (3, 1): (1, 2),
        (2, 1): (3, 3), (3, 2): (3, 1), (1, 3): (3, 2)}
DIM = 4 ** N


def digits(idx):
    return [(idx >> (2 * q)) & 3 for q in range(N)]


def index(ds):
    return sum(d << (2 * q) for q, d in enumerate(ds))


def pidx(d):
    code = {"X": 1, "Y": 2, "Z": 3}
    return index([code[d[q]] if q in d else 0 for q in range(N)])


def string_mul(p, q):
    k, r = 0, []
    for a, b in zip(digits(p), digits(q)):
        kk, rr = _MUL[(a, b)]
        k += kk
        r.append(rr)
    return k % 4, index(r)


def superop(p):
    """D_P: Q -> i[P, Q] の行列(整数)。反交換する場合 i[P,Q] = 2 i^(k+1) R。"""
    m = np.zeros((DIM, DIM), dtype=np.int64)
    for q in range(DIM):
        k, r = string_mul(p, q)
        k2, _ = string_mul(q, p)
        if k != k2:                      # 反交換(k は奇数)
            m[r, q] = 2 if (k + 1) % 4 == 0 else -2
    return m


D_X = sum(superop(pidx({i: "X"})) for i in range(N))
D_ZZ = {p: superop(pidx({p[0]: "Z", p[1]: "Z"})) for p in PAIRS}


def e0_vector(bits="0000"):
    """基底状態 |bits> でのパウリ文字列の期待値(I,Z のみ非零)。"""
    v = np.zeros(DIM, dtype=np.int64)
    for idx in range(DIM):
        ds = digits(idx)
        if all(d in (0, 3) for d in ds):
            v[idx] = (-1) ** sum(int(bits[q]) for q, d in enumerate(ds) if d == 3)
    return v


def taylor(ops, gens, e0, kmax):
    """
    ops : パウリ文字列の添字リスト
    gens: 超演算子のリスト(それぞれ別の変数に対応。H = sum_v x_v * P_v)
    戻り値: coeff[n][op番号] = {指数タプル: Fraction}   (t^n の係数。x_v の多項式)
    """
    assert kmax <= 14, "int64 の安全範囲(|係数| <= 20^n)を超える"
    nv = len(gens)
    V = {(0,) * nv: np.eye(DIM, dtype=np.int64)[:, ops]}             # (256, n_ops)
    out = []
    for n in range(kmax + 1):
        row = [dict() for _ in ops]
        for exp, mat in V.items():
            vals = e0 @ mat
            for j, val in enumerate(vals):
                if val:
                    row[j][exp] = Fraction(int(val), math.factorial(n))
        out.append(row)
        if n == kmax:
            break
        new = {}
        for exp, mat in V.items():
            for v, g in enumerate(gens):
                e2 = exp[:v] + (exp[v] + 1,) + exp[v + 1:]
                contrib = g @ mat
                new[e2] = new[e2] + contrib if e2 in new else contrib
        V = new
    return out


def padd(a, b, s=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + s * v
        if r[k] == 0:
            del r[k]
    return r


def pmul(a, b):
    r = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            k = tuple(x + y for x, y in zip(k1, k2))
            r[k] = r.get(k, 0) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}


def peval(p, values):
    return sum(c * math.prod(x ** e for x, e in zip(values, exp)) for exp, c in p.items())


def pfmt(p, names):
    if not p:
        return "0"
    terms = []
    for exp, c in sorted(p.items(), reverse=True):
        mono = "·".join(f"{nm}^{e}" if e > 1 else nm for nm, e in zip(names, exp) if e)
        terms.append(f"({c})" + ("·" + mono if mono else ""))
    return " + ".join(terms)


def connected_series(c_xx, c_a, c_b, kmax):
    """C(t) = <AB> - <A><B> のテイラー係数(コーシー積)。"""
    out = []
    for n in range(kmax + 1):
        p = dict(c_xx[n])
        for k in range(n + 1):
            p = padd(p, pmul(c_a[k], c_b[n - k]), -1)
        out.append(p)
    return out


def graph_distance(graph, a, b):
    adj = {i: set() for i in range(N)}
    for x, y in graph:
        adj[x].add(y)
        adj[y].add(x)
    dist, frontier = {a: 0}, [a]
    while frontier:
        nxt = []
        for u in frontier:
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    nxt.append(w)
        frontier = nxt
    return dist.get(b)                    # None = 非連結


ONE_KEYS = [(a, i) for a in "XYZ" for i in range(N)]
OPS_1 = [pidx({i: a}) for a, i in ONE_KEYS]
OPS_XX = [pidx({p[0]: "X", p[1]: "X"}) for p in PAIRS]


def part_B(classes_A):
    print("\n" + "=" * 76)
    print("[B] 交換子展開(厳密計算)   【実測。ただし有限次数 K_MAX まで】")
    print("=" * 76)
    e0 = e0_vector("0000")
    vals = (OMEGA, J)

    # --- 検算: 級数と対角化の一致 ---------------------------------------
    t = 0.3
    ser = taylor(OPS_1 + OPS_XX, [D_X, sum(D_ZZ[b] for b in RING)], e0, K_MAX)
    ps = evolve(hamiltonian(RING), np.array([t]), basis_state("0000"))
    err = 0.0
    for j, key in enumerate(ONE_KEYS):
        s = sum(float(peval(ser[n][j], vals)) * t ** n for n in range(K_MAX + 1))
        err = max(err, abs(s - expval(ONE[key], ps)[0]))
    print(f"検算: リング, t={t}, 1体12量の |級数(次数<={K_MAX}) - 対角化| の最大 = {err:.2e}")

    # --- B1: 1体量の同値類 ---------------------------------------------
    print("\n[B1] 1体量 <X_i>,<Y_i>,<Z_i> のテイラー係数による分類 (|0000>)")
    series = []
    for g in GRAPHS:
        gens = [D_X, sum((D_ZZ[b] for b in g), np.zeros((DIM, DIM), dtype=np.int64))]
        series.append(taylor(OPS_1 + OPS_XX, gens, e0, K_MAX))

    def sig(gi, kmax, numeric, ops=range(12)):
        s = []
        for n in range(kmax + 1):
            for j in ops:
                p = series[gi][n][j]
                s.append(peval(p, vals) if numeric else tuple(sorted(p.items())))
        return tuple(s)

    print("  次数K   類の数(Omega,Jの多項式として)   類の数(Omega=1, J=7/10 を代入)   <Z_i>のみ(代入)")
    for k in range(0, K_MAX + 1):
        n_sym = len({sig(i, k, False) for i in range(64)})
        n_num = len({sig(i, k, True) for i in range(64)})
        n_z = len({sig(i, k, True, range(8, 12)) for i in range(64)})
        print(f"  {k:>4d}   {n_sym:>14d}   {n_num:>30d}   {n_z:>22d}")
    sigs = {}
    for i in range(64):
        sigs.setdefault(sig(i, K_MAX, True), []).append(i)
    classes_B = {frozenset(c) for c in sigs.values()}
    print(f"  => 次数 {K_MAX} までで厳密に区別できる類: {len(classes_B)} 個"
          f"(同値類数の厳密な下界)。[A] の分類と一致: {classes_B == classes_A}")

    # --- B2: XX 連結相関の短時間次数 -----------------------------------
    print("\n[B2] XX 連結相関 C_ij(t) = <X_i X_j> - <X_i><X_j> の最低次 (|0000>)")
    table = {}
    for gi, g in enumerate(GRAPHS):
        ser = series[gi]
        for pj, (a, b) in enumerate(PAIRS):
            c = connected_series([r[12 + pj] for r in ser], [r[a] for r in ser],
                                 [r[b] for r in ser], K_MAX)
            order = next((n for n in range(K_MAX + 1) if c[n]), None)
            order_num = next((n for n in range(K_MAX + 1) if peval(c[n], vals) != 0), None)
            table.setdefault(graph_distance(g, a, b), []).append((order, order_num, c, g, (a, b)))
    print("  グラフ距離 d   件数   最低次(多項式として)   最低次(数値代入後)")
    for d in sorted(table, key=lambda x: (x is None, x)):
        rows = table[d]
        print(f"  {str(d) if d is not None else '非連結':>10}   {len(rows):>4d}   "
              f"{sorted({r[0] for r in rows}, key=str)!s:>18}   {sorted({r[1] for r in rows}, key=str)!s:>16}")
    print(f"  (None = 次数 {K_MAX} まで全係数が厳密に 0。非連結ペアは状態が積に分解するので恒等的に 0)")
    print("  最低次係数の例:")
    for d in (1, 2, 3):
        seen = set()
        for order, _, c, g, pr in table[d]:
            s = pfmt(c[order], ("Ω", "J"))
            if s not in seen:
                seen.add(s)
                print(f"    d={d}: t^{order} × [{s}]   例: G={gname(g)}, ペア{pr}")

    # --- B3: 開鎖・リング・弱リンク(lambda を記号として保持)-----------
    print("\n[B3] ペア(0,3): H(λ) = Ω ΣX + J(Z0Z1+Z1Z2+Z2Z3) + λJ Z0Z3   (λ=0 開鎖, λ=1 リング)")
    gens = [D_X, sum(D_ZZ[b] for b in OPEN), D_ZZ[(0, 3)]]          # 変数: Ω, J, μ=λJ
    x0, x3, xx = pidx({0: "X"}), pidx({3: "X"}), pidx({0: "X", 3: "X"})
    ser = taylor([x0, x3, xx], gens, e0, K_MAX)
    c = connected_series([r[2] for r in ser], [r[0] for r in ser], [r[1] for r in ser], K_MAX)
    for n in range(K_MAX + 1):
        if c[n]:
            lam0 = {k: v for k, v in c[n].items() if k[2] == 0}
            print(f"  t^{n:<2d}: λ=0 での係数 = {pfmt(lam0, ('Ω', 'J', 'μ'))}"
                  f" | 全体の項数 {len(c[n])}" + (f" | 全体 = {pfmt(c[n], ('Ω', 'J', 'μ'))}" if len(c[n]) <= 3 else ""))
    print("  (μ = λJ。λ=0 の列が最初に非零になる次数が開鎖の最低次)")

    # --- 数値フィットとの比較 -------------------------------------------
    print("\n[B4] 数値フィット(対角化)との比較: 局所傾き d ln|C_03| / d ln t   【実測】")
    for name, g, lam in (("リング", RING, 1), ("開鎖", OPEN, 0)):
        n0 = next(n for n in range(K_MAX + 1) if peval(c[n], (OMEGA, J, J * lam)) != 0)
        lead = float(peval(c[n0], (OMEGA, J, J * lam)))
        ts = np.array([0.05, 0.1, 0.2, 0.4])
        ps = evolve(hamiltonian(g), ts, basis_state("0000"))
        cn = expval(TWO[("X", (0, 3))], ps) - expval(ONE[("X", 0)], ps) * expval(ONE[("X", 3)], ps)
        slopes = np.diff(np.log(np.abs(cn))) / np.diff(np.log(ts))
        print(f"  {name}: 厳密な最低次 t^{n0}, 係数 {lead:+.6e}")
        print(f"        t = {ts.tolist()}")
        print(f"        C_03(数値)        = {np.array2string(cn, precision=4)}")
        print(f"        区間ごとの傾き      = {np.array2string(slopes, precision=3)}")
        print(f"        C_03 / (係数·t^{n0}) = {np.array2string(cn / (lead * ts ** n0), precision=4)}")


# =====================================================================
# [C] 対称性: 1体量が一致する十分条件
#   |0000> と ΣX は全ての置換で不変なので、置換 π に対し
#       <A_{π(i)}>_{π(G)}(t) = <A_i>_G(t).
#   さらに π(i) が Aut(G) の下で i と同じ軌道にあれば <A_{π(i)}>_G = <A_i>_G。
#   よって「π(G)=G' かつ 全ての i で π(i) ∈ Orb_G(i)」なら 1体量は厳密に一致。
# =====================================================================
def part_C(classes_A):
    print("\n" + "=" * 76)
    print("[C] 対称性による十分条件(置換24個の総当たり)   【厳密。|0000> のみ】")
    print("=" * 76)
    perms = list(itertools.permutations(range(N)))

    def act(pi, g):
        return tuple(sorted(tuple(sorted((pi[a], pi[b]))) for a, b in g))

    gid = {tuple(sorted(g)): i for i, g in enumerate(GRAPHS)}
    parent = list(range(64))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, g in enumerate(GRAPHS):
        g = tuple(sorted(g))
        auts = [pi for pi in perms if act(pi, g) == g]
        orbit = [{pi[s] for pi in auts} for s in range(N)]
        for pi in perms:
            if all(pi[s] in orbit[s] for s in range(N)):
                parent[find(gid[act(pi, g)])] = find(i)
    cl = {}
    for i in range(64):
        cl.setdefault(find(i), []).append(i)
    classes_C = {frozenset(c) for c in cl.values()}
    print(f"十分条件で同一視される類の数: {len(classes_C)}(同値類数の厳密な上界)")
    print(f"[A] の分類と一致: {classes_C == classes_A}")
    return classes_C


if __name__ == "__main__":
    print(__doc__)
    cA = part_A()
    part_B(cA)
    part_C(cA)
