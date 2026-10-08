# The Torn QR Code

**O(n)**

The working code from the [@ByteAnimates](https://www.facebook.com/ByteAnimates) reel.

```bash
c++ -std=c++17 -O2 -o main main.cpp
./main
```

No dependencies. Any C++17 compiler.

### `qr.cpp`

```cpp
const int N = 26, EC = 17;
bool scans(const bool lost[N],
           bool corners) {
  if (!corners) return false;
  int bad = 0;
  for (int i = 0; i < N; i++)
    bad += lost[i];
  return 2 * bad <= EC;
}
```

### Files

| | |
| --- | --- |
| `main.cpp` | build and run this — the demo, with real inputs and real output |
| `qr.cpp` | the reel snippet, generated from the episode |

---

The snippet above is generated from the video itself — what you read is byte-for-byte
what was typed on screen. A fix to it belongs in the episode, so open an issue and the
next reel carries it. Everything else here is hand-written and welcome as a pull request.
