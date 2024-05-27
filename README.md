# Hasty Pudding Cipher

This is a from-scratch implementation of the [Hasty Pudding cipher](https://en.wikipedia.org/wiki/Hasty_Pudding_cipher)
based on the spec [published](http://richard.schroeppel.name:8015/hpc/hpc-spec) by its author,
Rich Schroeppel.

The spec linked above includes a revision made in May of 1999 to solve an issue of equivalent
keys discovered by David Wagner.

This implementation was checked against the test vectors that the author provided. However,
since the test vectors were generated before the fix to the aforementioned issue was made,
when comparing against them the code has to be compiled with the definition `HPC_WITHOUT_WAGNER_FIX`
for the key schedule to match.

The Hasty Pudding cipher included some interesting ideas, which were revolutionary at the
time, such as a variable block size (with the spec describing a format preserving encryption
mode using cycle walking), as well as providing a tweak to change the way the cipher operates
without changing the key (originally termed "spice").

However, since the cipher did not advance to the second round of the AES competition, it was
not subjected to any serious cryptanalysis besides the aforementioned one by David Wagner,
and as such it's not recommended for use. It's also much slower than AES on modern hardware
which includes AES acceleration, and is very demanding of hardware that doesn't include such
acceleration due to its use of 64-bit words. If you want to use it because of a need to encrypt
small data or for format-preserving encryption, there are much better studied alternatives,
such as FF1 and FF3-1.

## Goal

I just wanted to experiment with implementing this cipher, because tested and publicly available
implementations are hard to come by (the one provided with Schneier's *Applied Cryptography*
only works with 128-bit blocks, for example).

## License

This library is licensed under the MIT license. See [LICENSE](LICENSE) for details.
