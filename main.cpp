// The demo. Build and run:
//
//     c++ -std=c++17 -O2 -o main main.cpp
//     ./main
//
// `qr.cpp` is the snippet from the reel, generated from the episode. It is INCLUDED rather than
// linked because it is exactly what was on screen: one function with no headers of its own.
//
// The snippet is a MODEL of a QR decoder: it counts lost bytes and checks them against the
// repair budget. This file runs the real thing beside it. A Reed-Solomon decoder over GF(256)
// (Berlekamp-Massey to find the errors, Chien search to locate them, Forney to fix them) is
// handed the same 26 bytes the reel's code holds, with the same bytes broken, and the demo checks
// the model and the decoder agree every time.
//
// That covers the byte budget only. The corner squares are not bytes, so no Reed-Solomon decoder
// has anything to say about them: the corner line below runs the model alone. That a covered
// corner stops a real scanner is checked in the episode's probe, against zxing-cpp and OpenCV.
//
// The 26 bytes are the version 1-H symbol for "STILL HERE": 9 of data, then 17 of repair.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace std;

#include "qr.cpp"

/** A check that survives `-DNDEBUG`: `assert` is compiled out there, and this demo's output is a
    claim, so a failure has to stop it rather than print "agrees" anyway. */
static void require(bool ok, const char* what) {
    if (!ok) {
        cerr << "FAILED: " << what << "\n";
        exit(1);
    }
}

static const vector<uint8_t> SYMBOL = {32,  85,  9,   103, 239, 85,  133, 204, 208,
                                       180, 17,  255, 114, 242, 67,  223, 131, 176,
                                       169, 60,  253, 227, 12,  153, 3,   80};

/* ── GF(256), the field QR codes do their arithmetic in: x^8 + x^4 + x^3 + x^2 + 1 ─────────── */

static uint8_t EXP[512], LOG[256];

static void field() {
    int x = 1;
    for (int i = 0; i < 255; i++) {
        EXP[i] = x;
        LOG[x] = i;
        x <<= 1;
        if (x & 0x100) x ^= 0x11d;
    }
    for (int i = 255; i < 512; i++) EXP[i] = EXP[i - 255];
}

static uint8_t mul(uint8_t a, uint8_t b) { return a && b ? EXP[LOG[a] + LOG[b]] : 0; }
static uint8_t inv(uint8_t a) { return EXP[255 - LOG[a]]; }

/** p(x) at x, coefficients lowest degree first. */
static uint8_t eval(const vector<uint8_t>& p, uint8_t x) {
    uint8_t y = 0;
    for (size_t i = p.size(); i-- > 0;) y = mul(y, x) ^ p[i];
    return y;
}

/** The 17 syndromes: the received bytes, read as a polynomial, at alpha^0 .. alpha^16. */
static vector<uint8_t> syndromes(const vector<uint8_t>& r) {
    vector<uint8_t> s(EC);
    for (int j = 0; j < EC; j++) {
        uint8_t v = 0;
        for (uint8_t c : r) v = mul(v, EXP[j]) ^ c;
        s[j] = v;
    }
    return s;
}

/** Repairs `r` in place. False when there are more errors than the repair bytes can find. */
static bool decode(vector<uint8_t>& r) {
    const int n = (int)r.size();
    vector<uint8_t> s = syndromes(r);
    bool clean = true;
    for (uint8_t v : s) clean = clean && v == 0;
    if (clean) return true;

    // Berlekamp-Massey: the shortest error locator that explains the syndromes.
    vector<uint8_t> lam = {1}, prev = {1};
    int L = 0, m = 1;
    uint8_t b = 1;
    for (int k = 0; k < EC; k++) {
        uint8_t d = s[k];
        for (int i = 1; i <= L; i++) d ^= mul(lam[i], s[k - i]);
        if (d == 0) {
            m++;
            continue;
        }
        vector<uint8_t> t = lam;
        uint8_t coef = mul(d, inv(b));
        if (lam.size() < prev.size() + m) lam.resize(prev.size() + m, 0);
        for (size_t i = 0; i < prev.size(); i++) lam[i + m] ^= mul(coef, prev[i]);
        if (2 * L <= k) {
            L = k + 1 - L;
            prev = t;
            b = d;
            m = 1;
        } else {
            m++;
        }
    }
    if (2 * L > EC) return false;

    // Chien search: byte i is wrong when the locator vanishes at alpha^-(n-1-i).
    vector<int> where;
    for (int i = 0; i < n; i++)
        if (eval(lam, EXP[(255 - (n - 1 - i)) % 255]) == 0) where.push_back(i);
    if ((int)where.size() != L) return false;

    // Forney: how wrong each one is.
    vector<uint8_t> omega(EC, 0);
    for (int i = 0; i < EC; i++)
        for (size_t j = 0; j < lam.size() && i + (int)j < EC; j++) omega[i + j] ^= mul(s[i], lam[j]);
    vector<uint8_t> dlam;
    for (size_t i = 1; i < lam.size(); i++) dlam.push_back(i % 2 ? lam[i] : 0);
    for (int i : where) {
        uint8_t X = EXP[n - 1 - i], Xinv = inv(X);
        uint8_t den = eval(dlam, Xinv);
        if (den == 0) return false;
        r[i] ^= mul(X, mul(eval(omega, Xinv), inv(den)));
    }
    for (uint8_t v : syndromes(r)) if (v) return false;
    return true;
}

/** The 9 data bytes as alphanumeric mode: 4 bits of mode, 9 of count, 11 bits a pair. */
static string message(const vector<uint8_t>& r) {
    const string A = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";
    string bits;
    for (int i = 0; i < 9; i++)
        for (int k = 7; k >= 0; k--) bits += (r[i] >> k & 1) ? '1' : '0';
    auto take = [&](size_t& at, int w) { int v = stoi(bits.substr(at, w), nullptr, 2); at += w; return v; };
    size_t at = 0;
    if (take(at, 4) != 2) return "";
    int count = take(at, 9);
    string out;
    while ((int)out.size() + 1 < count) {
        int v = take(at, 11);
        out += A[v / 45];
        out += A[v % 45];
    }
    if ((int)out.size() < count) out += A[take(at, 6)];
    return out;
}

/** Times the decoder "repaired" bytes into a different valid code: past the budget, it can. */
static int wrong = 0;

/** Break the listed bytes, run both the model and the decoder, and return whether it scanned. */
static bool scan(const vector<int>& broken, uint32_t flip) {
    vector<uint8_t> r = SYMBOL;
    bool lost[N] = {};
    for (int i : broken) {
        r[i] ^= (uint8_t)(flip | 1);  // any wrong value: the decoder only cares that it is wrong
        flip = flip * 1103515245u + 12345u;
        lost[i] = true;
    }
    bool decoded = decode(r);
    bool repaired = decoded && r == SYMBOL && message(r) == "STILL HERE";
    wrong += decoded && !repaired;
    bool model = scans(lost, true);
    require(repaired == model, "the decoder and the reel's model agree");
    return repaired;
}

int main() {
    field();
    require(syndromes(SYMBOL) == vector<uint8_t>(EC, 0), "the symbol is a valid codeword");
    require(message(SYMBOL) == "STILL HERE", "the symbol says STILL HERE");
    cout << "symbol: 26 bytes, 9 data + 17 repair, says \"" << message(SYMBOL) << "\"\n\n";

    // The reel's three tears, with the bytes each one breaks.
    const vector<vector<int>> tears = {
        {0, 5, 6},
        {0, 1, 2, 4, 5, 6, 11, 12},
        {0, 1, 2, 4, 5, 6, 7, 11, 12},
    };
    for (const auto& t : tears) {
        bool ok = scan(t, 0x5au);
        cout << t.size() << " bytes lost: " << (ok ? "repaired, reads STILL HERE" : "not repairable") << "\n";
    }

    // The model's other half: corners are checked before any byte is read, so no budget covers
    // them. The model only; see the header for where this is checked against real scanners.
    bool none[N] = {};
    require(!scans(none, false), "no corners, no scan");
    cout << "0 bytes lost, corner covered: the model says not found\n\n";

    // Every count of lost bytes, at random positions, many times over.
    mt19937 rng(20261001);
    for (int k = 0; k <= 12; k++) {
        int fixed = 0;
        wrong = 0;
        const int trials = 2000;
        for (int trial = 0; trial < trials; trial++) {
            vector<int> pos(N);
            for (int i = 0; i < N; i++) pos[i] = i;
            shuffle(pos.begin(), pos.end(), rng);
            pos.resize(k);
            fixed += scan(pos, rng());
        }
        cout << k << (k == 1 ? " byte " : " bytes") << " lost: repaired " << fixed << " of " << trials;
        if (wrong) cout << "  (" << wrong << " repaired into a different code)";
        cout << "\n";
    }
    cout << "\nthe decoder agrees with 2 * bad <= EC every time\n";
    return 0;
}
