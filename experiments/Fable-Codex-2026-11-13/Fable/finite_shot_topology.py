#!/usr/bin/env python3
"""
有限ショット・読み出し誤りのもとでの 4量子ビット・全64グラフ識別
====================================================================

【設定(仮定)】
  H_G = (Omega/2) sum_i X_i + (J/2) sum_{(i,j) in G} Z_i Z_j      (README §14 の規約)
  Omega = 1, J = 0.7, 初期状態 |0000>, 0 <= t <= 3.688
  1ショット = 時刻 t まで発展 → 4サイトを同じ基底 (X, Y, Z のいずれか) で射影測定
  読み出し誤り = 各ビットが独立に確率 p で反転 (p は既知)。p = 0, 2%, 5%
  状態準備・時刻・ゲートの誤差なし。J と Omega は既知、未知なのは辺の有無 (2^6 = 64 通り) のみ。
  合格基準: 64 グラフのどれが真でも、6 辺すべての有無を正答する確率 >= 95%
  非適応 (事前固定) の配分のみを考える。

【出力】
  1. 最大最小バタチャリヤ設計 (LP) と、その設計 + 最尤復号で足りる総ショット数 (モンテカルロ)
  2. 任意の事前固定配分・任意の復号に対する下限 (レニー 1/2 次ダイバージェンスのデータ処理不等式)
  3. 2,016 組のうち最も見分けにくい対、開鎖 vs リングの位置づけ
  4. 短時間 XX 相関の t^4 係数しきい値法との同ショット比較
  5. 総ショット数固定で検出できる最弱リンク lambda

依存: NumPy, SciPy (linprog)。  実行: python3 finite_shot_topology.py   (約15分)
"""
import itertools
import json
import math
import sys

import numpy as np
from scipy.optimize import linprog

# ------------------------------------------------------------------ パラメータ
N_SITE = 4
OMEGA, JJ, T_MAX = 1.0, 0.7, 3.688
P_FLIPS = [0.0, 0.02, 0.05]
EPS = 0.05                        # 許容誤り率 (正答率 95%)
DT_DESIGN = 0.004                 # 設計用の時刻グリッド
DT_CERT = 0.0005                  # 下限の検証用グリッド
T_MC = 50_000                     # 検証用モンテカルロ試行数 (グラフあたり)
SEED = 20261113

PAIRS = list(itertools.combinations(range(N_SITE), 2))
GRAPHS = [tuple(p for k, p in enumerate(PAIRS) if m >> k & 1) for m in range(64)]
OPEN = ((0, 1), (1, 2), (2, 3))
RING = ((0, 1), (0, 3), (1, 2), (2, 3))
GP = list(itertools.combinations(range(64), 2))                  # 2,016 組
GI = np.array([a for a, _ in GP])
GJ = np.array([b for _, b in GP])
BASES = "XYZ"


def gname(g):
    return "{" + ",".join(f"{a}{b}" for a, b in g) + "}" if g else "{}(辺なし)"


# ------------------------------------------------------------------ モデル
_I2 = np.eye(2)
_X = np.array([[0.0, 1.0], [1.0, 0.0]])
_Z = np.diag([1.0, -1.0])


def _embed(d):
    m = np.array([[1.0]])
    for i in range(N_SITE):
        m = np.kron(m, d.get(i, _I2))
    return m


def _kron_all(u):
    m = np.array([[1.0 + 0j]])
    for _ in range(N_SITE):
        m = np.kron(m, u)
    return m


XSUM = sum(_embed({i: _X}) for i in range(N_SITE))
ZZ = {p: _embed({p[0]: _Z, p[1]: _Z}) for p in PAIRS}
_HAD = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
ROT = {"X": _kron_all(_HAD), "Y": _kron_all(_HAD @ np.diag([1, -1j])), "Z": _kron_all(np.eye(2))}
# ビット 0 = 固有値 +1。サイト i のビットは x の上位から i 番目。
SIGN = np.array([[1 - 2 * ((x >> (N_SITE - 1 - i)) & 1) for i in range(N_SITE)] for x in range(16)])


def hamiltonian(graph, weak=None):
    """weak = {(i,j): lam} で、その辺を lam*J/2 の強さで追加する。"""
    h = OMEGA / 2 * XSUM
    for b in graph:
        h = h + JJ / 2 * ZZ[b]
    for b, lam in (weak or {}).items():
        h = h + lam * JJ / 2 * ZZ[b]
    return h


def flip_matrix(p):
    f = np.array([[1 - p, p], [p, 1 - p]])
    m = np.array([[1.0]])
    for _ in range(N_SITE):
        m = np.kron(m, f)
    return m


def outcome_probs(h, ts, p):
    """戻り値 P[基底, 時刻, 16通りのビット列] (読み出し誤り込み)。"""
    w, v = np.linalg.eigh(h)
    psi = v @ (np.exp(-1j * np.outer(w, ts)) * v[0].conj()[:, None])       # |0000> = 基底ベクトル0
    f = flip_matrix(p)
    return np.array([(f @ np.abs(ROT[b] @ psi) ** 2).T for b in BASES])


def grid_probs(ts, p, graphs=GRAPHS):
    """(グラフ数, 3*len(ts), 16)。設定の添字 s = 基底番号*len(ts) + 時刻番号。"""
    return np.array([outcome_probs(hamiltonian(g), ts, p).reshape(-1, 16) for g in graphs])


def setting_probs(settings, p):
    """settings = [(基底, 時刻), ...] → (64, K, 16)"""
    out = np.empty((64, len(settings), 16))
    for k, (b, t) in enumerate(settings):
        for g in range(64):
            out[g, k] = outcome_probs(hamiltonian(GRAPHS[g]), np.array([t]), p)[BASES.index(b), 0]
    return out


def bhattacharyya(P):
    """P (64,S,16) → B[s, 対] = -ln sum_x sqrt(P_G P_G')"""
    sq = np.sqrt(P)
    bc = np.einsum("psx,psx->sp", sq[GI], sq[GJ])
    return -np.log(np.clip(bc, 1e-300, 1.0))


# ------------------------------------------------------------------ 設計 (LP)
def maximin_design(B):
    """max_w min_対 sum_s w_s B[s,対]。戻り値: w, 値 v, 双対変数 q (対の上の確率分布)。"""
    S, M = B.shape
    c = np.zeros(S + 1)
    c[-1] = -1.0
    a_ub = np.hstack([-B.T, np.ones((M, 1))])
    a_eq = np.ones((1, S + 1))
    a_eq[0, -1] = 0.0
    r = linprog(c, A_ub=a_ub, b_ub=np.zeros(M), A_eq=a_eq, b_eq=[1.0],
                bounds=[(0, None)] * S + [(None, None)], method="highs")
    q = -r.ineqlin.marginals
    return r.x[:S], r.x[-1], q / q.sum()


def merge_close(settings, w, tol=3 * DT_DESIGN):
    """同じ基底で時刻がほぼ同じ設定 (LP がグリッド隣接点に重みを割ったもの) を 1 つにまとめる。"""
    out = []
    for (b, t), wk in sorted(zip(settings, w)):
        if out and out[-1][0] == b and t - out[-1][1] <= tol:
            _, t0, w0 = out[-1]
            out[-1] = (b, round((t0 * w0 + t * wk) / (w0 + wk) / DT_DESIGN) * DT_DESIGN, w0 + wk)
        else:
            out.append((b, t, wk))
    return [(b, float(t)) for b, t, _ in out], np.array([x[2] for x in out])


_PERM_MAPS = []
for _pi in itertools.permutations(range(N_SITE)):
    _idx = [PAIRS.index(tuple(sorted((_pi[a], _pi[b])))) for a, b in PAIRS]
    _PERM_MAPS.append([sum(1 << _idx[k] for k in range(6) if m >> k & 1) for m in range(64)])


def orbit_key(a, b):
    """サイトの置換で移り合うグラフ対を同一視するための代表キー。"""
    return min(tuple(sorted((pm[a], pm[b]))) for pm in _PERM_MAPS)


def allocate(w, n_total):
    """最大剰余法で整数ショットに丸める。"""
    w = np.asarray(w, float) / np.sum(w)
    x = w * n_total
    n = np.floor(x).astype(int)
    n[np.argsort(-(x - n))[: n_total - n.sum()]] += 1
    return n


# ------------------------------------------------------------------ 復号とモンテカルロ
_MARG = np.zeros((6 * 4, 16))                         # 2ビット周辺分布への射影
for _k, (_i, _j) in enumerate(PAIRS):
    for _x in range(16):
        _a, _b = (SIGN[_x, _i] < 0), (SIGN[_x, _j] < 0)
        _MARG[4 * _k + 2 * _a + _b, _x] = 1.0


def log_weights(P, decoder):
    """復号器ごとの「1ショットあたりのスコア」 W[g,k,x]。スコア = sum_k counts_k · W[g,k]。"""
    if decoder == "ml":            # 4ビット列の完全な尤度 (3体・4体相関も使う)
        return np.log(np.clip(P, 1e-300, None))
    if decoder == "pair":          # 2ビット周辺分布だけの合成尤度 (2体相関までしか使わない)
        return np.log(np.clip(P @ _MARG.T, 1e-300, None)) @ _MARG
    raise ValueError(decoder)


def success_prob(P, ns, trials, rng, decoder="ml", confusion=False):
    """各グラフが真のとき、復号結果がそのグラフに一意に一致する確率。"""
    W = log_weights(P, decoder)
    ok = np.zeros(64)
    conf = np.zeros((64, 64))
    for g in range(64):
        score = np.zeros((trials, 64))
        for k, n in enumerate(ns):
            if n:
                counts = rng.multinomial(n, P[g, k] / P[g, k].sum(), size=trials)
                score += counts @ W[:, k, :].T
        own = score[:, g].copy()
        score[:, g] = -np.inf
        win = own > score.max(axis=1)                 # 同点は不正解として数える
        ok[g] = win.mean()
        if confusion:
            conf[g] = np.bincount(score.argmax(axis=1)[~win], minlength=64) / trials
    return (ok, conf) if confusion else ok


def required_shots(P, w, rng, lo=8, hi=4096):
    """最悪グラフの正答率が 95% 以上になる最小の総ショット数。
    粗い二分探索のあと、T_MC 試行で「推定値 - 2σ >= 0.95」を満たす最小の N を採用する。"""
    def worst(n, trials):
        return success_prob(P, allocate(w, n), trials, rng).min()
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if worst(mid, 6000) >= 1 - EPS:
            hi = mid
        else:
            lo = mid
    n = max(hi - 4, 1)
    while True:
        s = worst(n, T_MC)
        if s - 2 * math.sqrt(s * (1 - s) / T_MC) >= 1 - EPS:
            return n, s
        n += 1


def binary_minimax_error(P0, P1, ns, trials, rng):
    """2仮説の尤度比検定。しきい値を動かして max(第1種, 第2種) を最小化した値。"""
    l0 = np.zeros(trials)
    l1 = np.zeros(trials)
    for k, n in enumerate(ns):
        if n:
            llr = np.log(np.clip(P1[k], 1e-300, None)) - np.log(np.clip(P0[k], 1e-300, None))
            l0 += rng.multinomial(n, P0[k] / P0[k].sum(), size=trials) @ llr
            l1 += rng.multinomial(n, P1[k] / P1[k].sum(), size=trials) @ llr
    cand = np.unique(np.concatenate([l0, l1]))
    a = 1 - np.searchsorted(np.sort(l0), cand, side="right") / trials      # P0(llr > tau)
    b = np.searchsorted(np.sort(l1), cand, side="right") / trials          # P1(llr <= tau)
    return np.maximum(a, b).min()


# ------------------------------------------------------------------ 下限
D_HALF = -math.log(2 * math.sqrt(EPS * (1 - EPS)))    # = 0.8303 (eps = 5%)


def lower_bound(p, q, ts_design):
    """
    任意の検定で 2 グラフ G, G' をともに誤り <= eps で区別するには
        sum_s n_s B_s(G,G') >= -ln(2 sqrt(eps(1-eps)))
    が必要 (レニー 1/2 次ダイバージェンスのデータ処理不等式)。
    対の上の確率分布 q で平均すると  N * max_s sum_対 q B_s >= 同じ定数。
    q には LP の双対解を使い、max_s は細かい時刻グリッドで取り直す。
    """
    ts = np.arange(0, T_MAX + 1e-12, DT_CERT)
    act = np.where(q > 1e-9)[0]
    need = sorted(set(GI[act]) | set(GJ[act]))
    pos = {g: i for i, g in enumerate(need)}
    P = np.sqrt(grid_probs(ts, p, [GRAPHS[g] for g in need]))
    f = np.zeros(P.shape[1])
    for k in act:
        f += q[k] * -np.log(np.clip((P[pos[GI[k]]] * P[pos[GJ[k]]]).sum(1), 1e-300, 1))
    return D_HALF / f.max(), f.max(), len(act)


# ------------------------------------------------------------------ t^4 係数法
def c4_exact(graph, pair):
    """XX 連結相関の t^4 係数 (Omega/2, J/2 規約)。辺があれば -(1/4 + 5k/24) Omega^2 J^2、
    k = その辺と端点を共有する他の辺の数。辺がなければ 0。(起動時に数値で検算する)"""
    if pair not in graph:
        return 0.0
    k = sum(1 for b in graph if b != pair and (set(b) & set(pair)))
    return -(0.25 + 5 * k / 24) * OMEGA ** 2 * JJ ** 2


C4_THRESHOLD = -0.125 * OMEGA ** 2 * JJ ** 2         # 最弱の辺の係数 (-1/4) の半分


def xx_connected(counts, p):
    """X 基底のカウント → 読み出し誤り補正済みの連結相関 (試行数, 6)。"""
    n = counts.sum(axis=1, keepdims=True)
    m1 = counts @ SIGN / n / (1 - 2 * p)
    out = []
    for i, j in PAIRS:
        m2 = counts @ (SIGN[:, i] * SIGN[:, j]) / n[:, 0] / (1 - 2 * p) ** 2
        out.append(m2 - m1[:, i] * m1[:, j])
    return np.stack(out, axis=1)


def t4_success(PX, t, n_total, p, trials, rng):
    """全ショットを X 基底・時刻 t に使い、c4 = C/t^4 < しきい値 なら「辺あり」と判定。"""
    ok = np.zeros(64)
    for g in range(64):
        counts = rng.multinomial(n_total, PX[g] / PX[g].sum(), size=trials)
        est = xx_connected(counts, p) / t ** 4 < C4_THRESHOLD
        truth = np.array([pr in GRAPHS[g] for pr in PAIRS])
        ok[g] = (est == truth).all(axis=1).mean()
    return ok


def t4_noise_free_margin(PX0, t):
    """ノイズなし・無限ショットでの判定余裕。min(辺: c/c4閾値比, 非辺) を返す。全グラフで正しければ True。"""
    ok = True
    ratio = []
    for g in range(64):
        c = xx_connected(PX0[g][None, :] * 1e12, 0.0)[0] / t ** 4
        for k, pr in enumerate(PAIRS):
            ok &= (c[k] < C4_THRESHOLD) == (pr in GRAPHS[g])
            if pr in GRAPHS[g]:
                ratio.append(c[k] / c4_exact(GRAPHS[g], pr))
    return ok, min(ratio), max(ratio)


def t4_required_shots(PX, t, p, rng, cap=10 ** 9):
    lo, hi = 1, 64
    while t4_success(PX, t, hi, p, 1500, rng).min() < 1 - EPS:
        lo, hi = hi, hi * 2
        if hi > cap:
            return None
    while hi / lo > 1.03:
        mid = int(math.sqrt(lo * hi))
        if t4_success(PX, t, mid, p, 1500, rng).min() >= 1 - EPS:
            hi = mid
        else:
            lo = mid
    return hi


# ------------------------------------------------------------------ 最弱リンク
def weak_link_probs(settings, p, lam, base=OPEN, bond=(0, 3)):
    h = hamiltonian(base, {bond: lam})
    return np.array([outcome_probs(h, np.array([t]), p)[BASES.index(b), 0] for b, t in settings])


def weakest_link(get_design, p, rng, trials=20_000):
    """lambda=0 (辺なし) と lambda (既知) の 2仮説を、両側誤り <= 5% で区別できる最小の lambda。
    get_design(lam) -> (settings, ns)。lambda ごとに設定を選び直す設計も扱える。"""
    def err(lam):
        settings, ns = get_design(lam)
        return binary_minimax_error(weak_link_probs(settings, p, 0.0), weak_link_probs(settings, p, lam),
                                    ns, trials, rng)
    lo, hi = 1e-3, 1.0
    if err(hi) > EPS:
        return None
    for _ in range(14):
        mid = math.sqrt(lo * hi)
        if err(mid) <= EPS:
            hi = mid
        else:
            lo = mid
    return hi


def weak_link_bhatta(p, ts, lam, base=OPEN, bond=(0, 3)):
    """設定ごとの B_s(lambda=0 と lambda)。"""
    a = np.sqrt(outcome_probs(hamiltonian(base), ts, p).reshape(-1, 16))
    b = np.sqrt(outcome_probs(hamiltonian(base, {bond: lam}), ts, p).reshape(-1, 16))
    return -np.log(np.clip((a * b).sum(1), 1e-300, 1))


def fisher_per_shot(p, ts, base, bond, h=0.02):
    """lambda=0 での古典フィッシャー情報量 F_s ≈ 8 B_s(h)/h^2 (設定 s ごと)。"""
    return 8 * weak_link_bhatta(p, ts, h, base, bond) / h ** 2


# ==================================================================== main
def main():
    rng = np.random.default_rng(SEED)
    out = {}
    print(__doc__)

    # --- 検算: t^4 係数の閉形式 ------------------------------------------
    errs = []
    for t0 in (0.02, 0.01):
        px = grid_probs(np.array([t0]), 0.0)[:, 0, :]
        errs.append(max(abs(xx_connected(px[g][None, :] * 1e12, 0.0)[0][k] / t0 ** 4 - c4_exact(GRAPHS[g], pr))
                        for g in range(64) for k, pr in enumerate(PAIRS)))
    print(f"[検算] t^4 係数の閉形式 -(1/4+5k/24)Ω²J² と数値 C(t)/t^4 の最大差: t=0.02 で {errs[0]:.1e}, "
          f"t=0.01 で {errs[1]:.1e} (t² に比例して減少)")
    print(f"       しきい値 {C4_THRESHOLD:.5f}, 最弱の辺の係数 {-0.25 * JJ ** 2:.5f}")

    ts = np.arange(0, T_MAX + 1e-12, DT_DESIGN)
    nt = len(ts)
    i_open, i_ring = GRAPHS.index(OPEN), GRAPHS.index(RING)
    k_or = GP.index((min(i_open, i_ring), max(i_open, i_ring)))

    for p in P_FLIPS:
        res = out[f"p={p}"] = {}
        print("\n" + "#" * 78 + f"\n# 読み出し反転率 p = {p:.0%}\n" + "#" * 78)

        # ---------------- 1. 設計 ------------------------------------------
        B = bhattacharyya(grid_probs(ts, p))
        w, v, q = maximin_design(B)
        sup = np.where(w > 1e-5)[0]
        settings, wd = merge_close([(BASES[s // nt], float(ts[s % nt])) for s in sup], w[sup])
        wd = wd / wd.sum()
        best1 = B.min(axis=1).argmax()
        print(f"\n[1] 設計: 最大最小バタチャリヤ距離 (時刻グリッド dt={DT_DESIGN}, LP)")
        print(f"    1ショットあたりの最悪対距離 v* = {v:.5f}"
              f"   (単一設定での最良: {BASES[best1 // nt]} 基底, t={ts[best1 % nt]:.3f}, {B.min(axis=1).max():.5f})")
        P = setting_probs(settings, p)
        n_star, s_star = required_shots(P, wd, rng)
        ns = allocate(wd, n_star)
        ok, conf = success_prob(P, ns, 2 * T_MC, rng, confusion=True)
        print(f"    必要総ショット数 N* = {n_star}   (最悪グラフの正答率 {ok.min():.4f} ± "
              f"{math.sqrt(ok.min() * (1 - ok.min()) / (2 * T_MC)):.4f}, 64グラフ平均 {ok.mean():.4f}, 試行 {2 * T_MC})")
        print("    基底   時刻     配分比    ショット数")
        for (b, t), wk, nk in zip(settings, wd, ns):
            print(f"      {b}    {t:6.3f}   {wk:6.3f}   {nk:5d}")
        print(f"    最も正答率が低いグラフ: {gname(GRAPHS[ok.argmin()])}")
        # 最尤復号の合併上界による保証 (モンテカルロ不要)
        Bd = bhattacharyya(P)
        n_ub = n_star
        while True:
            e = np.exp(-(allocate(wd, n_ub) @ Bd))
            if (np.bincount(GI, e, 64) + np.bincount(GJ, e, 64)).max() <= EPS:
                break
            n_ub += 1
        set1 = [(BASES[best1 // nt], float(ts[best1 % nt]))]
        n_single, s_single = required_shots(setting_probs(set1, p), [1.0], rng)
        print(f"    簡易版: 全ショットを {set1[0][0]} 基底・t={set1[0][1]:.3f} の 1 設定に使う場合 N = {n_single}"
              f" (最悪グラフ {s_single:.4f})")
        res.update(N_single_setting=n_single, single_setting=set1[0])
        print(f"    参考: 合併上界 sum_G' exp(-sum n_s B_s) <= 5% で厳密に保証できるショット数 = {n_ub}")
        res.update(design=[(b, t, int(n)) for (b, t), n in zip(settings, ns)], N_star=n_star,
                   worst_success=float(ok.min()), N_union_bound=n_ub, v_star=float(v))

        # ---------------- 2. 下限 ------------------------------------------
        lb, fmax, nact = lower_bound(p, q, ts)
        print(f"\n[2] 下限 (任意の事前固定配分・任意の復号): N >= {D_HALF:.4f} / {fmax:.5f} = {lb:.2f}"
              f"  →  N >= {math.ceil(lb - 1e-9)}")
        print(f"    (双対証明書: {nact} 組の対に重み。dt={DT_CERT} の細グリッドで最大値を再評価)")
        print(f"    設計の N* との比 = {n_star / math.ceil(lb - 1e-9):.2f}")
        res.update(lower_bound=lb, lower_bound_int=math.ceil(lb - 1e-9))

        # ---------------- 3. 見分けにくい対 --------------------------------
        d_pair = (ns @ Bd) / n_star                               # 整数配分での 1ショットあたり距離
        sym = conf[GI, GJ] + conf[GJ, GI]
        order = np.lexsort((d_pair, -sym))
        print(f"\n[3] 見分けにくい対 (N*={n_star} での64択最尤復号。混同率 = P(G→G') + P(G'→G))")
        orbits = {}
        for k, (ga, gb) in enumerate(GP):
            orbits.setdefault(orbit_key(ga, gb), []).append(k)
        stat = sorted(((sym[ks].mean(), d_pair[ks].mean(), ks) for ks in orbits.values()), key=lambda x: -x[0])
        print(f"    サイト置換で移り合う対は真の混同率が等しいので、{len(orbits)} 個の軌道ごとに平均した。")
        print("    順位  混同率(軌道平均±σ)   D/shot    軌道の大きさ  代表対")
        hard = []
        for r, (c, d, ks) in enumerate(stat[:6], 1):
            k = ks[int(np.argmax(sym[ks]))]
            ga, gb = GRAPHS[GI[k]], GRAPHS[GJ[k]]
            diff = ",".join(f"{x}{y}" for x, y in sorted(set(ga) ^ set(gb)))
            print(f"    {r:>3d}   {c:.5f} ± {math.sqrt(c / (2 * T_MC * len(ks))):.5f}   {d:.5f}   {len(ks):>6d}      "
                  f"{gname(ga)} vs {gname(gb)}  (違う辺: {diff})")
            hard.append(dict(a=ga, b=gb, confusion=float(c), D=float(d), orbit_size=len(ks),
                             orbit=[(GRAPHS[GI[x]], GRAPHS[GJ[x]]) for x in ks]))
        print("    1位の軌道に属する全ての対:")
        for x in stat[0][2]:
            print(f"        {gname(GRAPHS[GI[x]])} vs {gname(GRAPHS[GJ[x]])}")
        n_edge_diff = np.array([len(set(GRAPHS[a]) ^ set(GRAPHS[b])) for a, b in GP])
        print(f"    1辺だけ違う 192 組の D: 最小 {d_pair[n_edge_diff == 1].min():.5f}, 中央値 "
              f"{np.median(d_pair[n_edge_diff == 1]):.5f}, 最大 {d_pair[n_edge_diff == 1].max():.5f}"
              f" / 2辺以上違う組の最小 {d_pair[n_edge_diff > 1].min():.5f}")
        rank_d = int((d_pair < d_pair[k_or] - 1e-12).sum()) + 1
        rank_d1 = int((d_pair[n_edge_diff == 1] < d_pair[k_or] - 1e-12).sum()) + 1
        c_or = next(c for c, _, ks in stat if k_or in ks)
        rank_c = next(r for r, (_, _, ks) in enumerate(stat, 1) if k_or in ks)
        e_or = binary_minimax_error(P[i_open], P[i_ring], ns, 200_000, rng)
        print(f"    開鎖 vs リング: D = {d_pair[k_or]:.5f} (2,016組中 小さい方から {rank_d} 位、1辺違いの192組中 {rank_d1} 位), "
              f"混同率 (軌道平均) {c_or:.5f} ({len(orbits)} 軌道中 {rank_c} 位)")
        print(f"       2択だけなら N*={n_star} での誤り率 (minimax) = {e_or:.5f}")
        # 開鎖 vs リングだけを解く場合の必要ショット数
        b_or = B[:, k_or]
        s_or = b_or.argmax()
        set_or = [(BASES[s_or // nt], float(ts[s_or % nt]))]
        p_or = setting_probs(set_or, p)
        n2 = 1
        while binary_minimax_error(p_or[i_open], p_or[i_ring], [n2], 40_000, rng) > EPS:
            n2 += 1
        n2d = 1
        while binary_minimax_error(P[i_open], P[i_ring], allocate(wd, n2d), 40_000, rng) > EPS:
            n2d += 1
        print(f"       2択専用の最良設定 {set_or[0][0]} 基底 t={set_or[0][1]:.3f} (B={b_or.max():.4f}/shot) なら {n2} ショット、"
              f"64択用設計の配分比のままなら {n2d} ショットで両側 95%")
        res.update(hardest=hard, open_ring=dict(D=float(d_pair[k_or]), rank_D=rank_d, rank_D_among_single_edge=rank_d1, confusion=float(c_or),
                                                 rank_confusion=rank_c, binary_error_at_Nstar=float(e_or),
                                                 shots_binary_best=n2, shots_binary_design=n2d))

        # ---------------- 4. t^4 係数法との比較 ----------------------------
        print(f"\n[4] 同じ総ショット数 N*={n_star} での比較 (最悪グラフ / 64グラフ平均 の正答率)")
        okp = success_prob(P, ns, T_MC, rng, decoder="pair")
        print(f"    最適設計 + 4ビット列の最尤 (高次相関も利用)   : {ok.min():.4f} / {ok.mean():.4f}")
        print(f"    最適設計 + 2ビット周辺の合成尤度 (2体まで)   : {okp.min():.4f} / {okp.mean():.4f}")
        rows = []
        print("    --- 全ショットを X 基底・単一時刻 t_s に使う場合 ---")
        print("     t_s   無雑音で正しい  C/(c4 t^4)の範囲   t^4しきい値法      同じデータで最尤     t^4法が95%に要るN")
        for t_s in (0.25, 0.5, 0.75, 0.9, 1.0, 1.5, 2.0, 3.0):
            PXs = setting_probs([("X", t_s)], p)
            PX0 = setting_probs([("X", t_s)], 0.0)[:, 0, :]
            valid, rmin, rmax = t4_noise_free_margin(PX0, t_s)
            a = t4_success(PXs[:, 0, :], t_s, n_star, p, 10_000, rng)
            m = success_prob(PXs, [n_star], 10_000, rng)
            need = t4_required_shots(PXs[:, 0, :], t_s, p, rng) if valid else None
            print(f"    {t_s:5.2f}   {'はい' if valid else 'いいえ':^10}   {rmin:6.3f} – {rmax:6.3f}    "
                  f"{a.min():.4f} / {a.mean():.4f}    {m.min():.4f} / {m.mean():.4f}    "
                  f"{need if need else '到達不能'}")
            rows.append(dict(t=t_s, valid=bool(valid), ratio=(rmin, rmax), t4=(float(a.min()), float(a.mean())),
                             ml=(float(m.min()), float(m.mean())), t4_required=need))
        res.update(pair_decoder=(float(okp.min()), float(okp.mean())), t4_table=rows)

        # ---------------- 5. 最弱リンク ------------------------------------
        print(f"\n[5] 最弱リンク: H_開鎖 + λ (J/2) Z0 Z3 を λ=0 と区別 (総ショット数 N*={n_star}, 両側誤り <= 5%)")
        F = fisher_per_shot(p, ts, OPEN, (0, 3))
        f_design = sum(wk * F[BASES.index(b) * nt + int(round(t / DT_DESIGN))] for (b, t), wk in zip(settings, wd))
        s_f = F.argmax()
        set_f = [(BASES[s_f // nt], float(ts[s_f % nt]))]
        zq = 2 * 1.6448536269514722                              # 両側 5% のガウス近似: sqrt(N F) λ = 2 z_0.95
        tw = ts[::2]

        def dedicated(n):
            def get(lam):
                sb = int(weak_link_bhatta(p, tw, lam).argmax())
                return [(BASES[sb // len(tw)], float(tw[sb % len(tw)]))], [n]
            return get
        lam_a = weakest_link(lambda lam: (settings, ns), p, rng)
        lam_b = weakest_link(dedicated(n_star), p, rng)
        set_b = dedicated(n_star)(lam_b)[0][0]
        print(f"    (a) 64択用の設計のまま             : λ_min = {lam_a:.3f}   [F = {f_design:.3f}/shot, "
              f"ガウス近似 {zq / math.sqrt(n_star * f_design):.3f}]")
        print(f"    (b) その λ の検出に最良の単一設定   : λ_min = {lam_b:.3f}   [{set_b[0]} 基底 t={set_b[1]:.3f}]")

        def floor(n):                       # N * max_s B_s(λ) < D_HALF なら、どんな固定配分でも不可能
            lo, hi = 1e-4, 1.0
            for _ in range(40):
                mid = math.sqrt(lo * hi)
                lo, hi = (lo, mid) if n * weak_link_bhatta(p, tw, mid).max() >= D_HALF else (mid, hi)
            return hi
        print(f"        どんな固定配分・復号でも検出できない範囲 (厳密な下限): λ < {floor(n_star):.3f}")
        print(f"    小さい λ での法則: λ_min ≈ 2 z_0.95 / sqrt(N F_max) = {zq / math.sqrt(F.max()):.2f} / sqrt(N)"
              f"   [F_max = {F.max():.3f}/shot, {set_f[0][0]} 基底 t={set_f[0][1]:.3f}]")
        print("        N         λ_min(モンテカルロ)   公式     厳密な下限")
        scal = []
        for n_big in (1_000, 10_000, 100_000):
            lm = weakest_link(lambda lam, n_big=n_big: (set_f, [n_big]), p, rng)
            print(f"      {n_big:>7d}      {lm:.4f}            {zq / math.sqrt(n_big * F.max()):.4f}   {floor(n_big):.4f}")
            scal.append((n_big, lm, zq / math.sqrt(n_big * F.max()), floor(n_big)))
        # 全 192 通りの (グラフ, 欠けている辺) での F_max
        tc = ts[::4]
        fm = [(fisher_per_shot(p, tc, g, e).max(), g, e) for g in GRAPHS for e in PAIRS if e not in g]
        fmin, fmaxx = min(fm, key=lambda x: x[0]), max(fm, key=lambda x: x[0])
        print(f"    192通りの (グラフ, ない辺) でのガウス近似 λ_min (専用設定, N*): "
              f"{zq / math.sqrt(n_star * fmaxx[0]):.3f} (最易: {gname(fmaxx[1])} に {fmaxx[2]}) – "
              f"{zq / math.sqrt(n_star * fmin[0]):.3f} (最難: {gname(fmin[1])} に {fmin[2]})")
        res.update(weak_link=dict(design=lam_a, dedicated=lam_b, dedicated_setting=set_b, F_design=float(f_design),
                                  F_max=float(F.max()), undetectable_below=floor(n_star), scaling=scal,
                                  coef_over_sqrtN=zq / math.sqrt(F.max())))
        sys.stdout.flush()

    with open("finite_shot_results.json", "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o if not isinstance(o, np.generic) else o.item())
    print("\n結果を finite_shot_results.json に保存しました。")


if __name__ == "__main__":
    main()
