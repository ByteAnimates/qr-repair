// The Torn QR Code — O(n)
// As shown in the reel
// Generated from the episode; edits here are overwritten. See the README.

const int N = 26, EC = 17;
bool scans(const bool lost[N],
           bool corners) {
  if (!corners) return false;
  int bad = 0;
  for (int i = 0; i < N; i++)
    bad += lost[i];
  return 2 * bad <= EC;
}
