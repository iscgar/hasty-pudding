import argparse
import ctypes
import functools
import json
import operator
import re
import sys


MCT_HDR_META = re.compile(rb'^\s*(\S+) tests per blocksize\s+seed =\s*(\S+)\s+backup =\s*(\S+)$')

hpc_dll = None


class HpcDll:
    backup_t = ctypes.c_size_t * 6
    spice_t = ctypes.c_uint64 * 8

    def __init__(self, p):
        self.hpc_dll = ctypes.cdll.LoadLibrary(p)
        self.hpc_dll.hpc_encrypt_data.restype = ctypes.c_int
        self.hpc_dll.hpc_encrypt_data.argtypes = [
            ctypes.POINTER(ctypes.c_uint64), ctypes.c_size_t,  # key
            ctypes.POINTER(ctypes.c_size_t), ctypes.c_size_t,  # backup
            ctypes.POINTER(ctypes.c_uint64), ctypes.c_size_t,  # tweak
            ctypes.POINTER(ctypes.c_uint64), ctypes.c_size_t,  # data
        ]
        self.hpc_dll.hpc_decrypt_data.restype = ctypes.c_int
        self.hpc_dll.hpc_decrypt_data.argtypes = [
            ctypes.POINTER(ctypes.c_uint64), ctypes.c_size_t,  # key
            ctypes.POINTER(ctypes.c_size_t), ctypes.c_size_t,  # backup
            ctypes.POINTER(ctypes.c_uint64), ctypes.c_size_t,  # tweak
            ctypes.POINTER(ctypes.c_uint64), ctypes.c_size_t,  # data
        ]

    def encrypt(self, key, key_size, backup, data, data_bit_size, spice=None):
        key_t = ctypes.c_uint64 * len(key)
        data_t = ctypes.c_uint64 * len(data)
        key = key_t(*key)
        backup = self.backup_t(*backup)
        spice = self.spice_t(*spice) if spice else None
        data = data_t(*data)
        if self.hpc_dll.hpc_encrypt_data(
                ctypes.cast(key, ctypes.POINTER(ctypes.c_uint64)), key_size,
                ctypes.cast(backup, ctypes.POINTER(ctypes.c_size_t)), len(backup),
                ctypes.cast(spice, ctypes.POINTER(ctypes.c_uint64)),
                0 if not spice else len(spice) * 64,
                ctypes.cast(data, ctypes.POINTER(ctypes.c_uint64)), data_bit_size) == 0:
            raise RuntimeError('Failed to encrypt')
        return list(data)

    def decrypt(self, key, key_size, backup, data, data_bit_size, spice=None):
        key_t = ctypes.c_uint64 * len(key)
        data_t = ctypes.c_uint64 * len(data)
        key = key_t(*key)
        backup = self.backup_t(*backup)
        spice = self.spice_t(*spice) if spice else self.spice_t()
        data = data_t(*data)
        if self.hpc_dll.hpc_decrypt_data(
                ctypes.cast(key, ctypes.POINTER(ctypes.c_uint64)), key_size,
                ctypes.cast(backup, ctypes.POINTER(ctypes.c_size_t)), len(backup),
                ctypes.cast(spice, ctypes.POINTER(ctypes.c_uint64)),
                0 if not spice else len(spice) * 64,
                ctypes.cast(data, ctypes.POINTER(ctypes.c_uint64)), data_bit_size) == 0:
            raise RuntimeError('Failed to decrypt')
        return list(data)


def parse_data_words(d):
    return [int(d[i:i+16], 16) for i in range(0, len(d), 16)]


def fold_data_words(w):
    return sum(x << (i * 64) for i, x in enumerate(w))


def generate_mct(seed, block_size, key_size, n, op):
    PI19 = 3141592653589793238
    MASK = (1 << 64) - 1

    key_words = (key_size + 63) // 64
    kmask = (1 << (key_size & 63)) - 1
    block_words = (block_size + 63) // 64
    bmask = (1 << (block_size & 63)) - 1

    def initrand(nitems):
        a = [(PI19 + nitems*nitems*nitems*nitems*nitems) & MASK]
        for i in range(1, nitems):
            a.append((1 + 5 * a[i - 1]) & MASK)
        return a

    def myrand(arr):
        v = PI19 ^ a[-1]
        r = []
        for i, t in enumerate(a):
            v = (v ^ t if (i & 1) else v + t) & MASK
            v = (v << 23 | v >> 41) & MASK
            r.append(v)
        return r

    a = initrand(20)
    a[0] = (a[0] + seed) & MASK
    a = myrand(a)

    tests = []

    for i in range(n):
        a = myrand(a)
        key = a[:key_words]
        if kmask:
            key[-1] &= kmask
        a = myrand(a)
        spi = a[:8]
        a = myrand(a)
        pt = a[:block_words]
        if bmask:
            pt[-1] &= bmask
        tests.append({
            'key': key,
            'spice': spi,
            'ptxt' if op == 'encryption' else 'ctxt': pt
        })

    return tests


def parse_mct(inp):
    tests = {}

    with inp as f:
        header = [f.readline().rstrip() for _ in range(5)]

        op = header[0].rsplit(None, 1)[-1].lower().decode('ascii')
        test_per_bs, seed, backup = MCT_HDR_META.match(header[2]).groups()
        test_per_bs = int(test_per_bs)
        seed = int(seed, 16)
        backup = int(backup, 16)
        all_backup = backup & 0xf
        backup = [all_backup + ((backup >> (n * 4)) & 0xf) for n in range(1, 7)]

        current_key_size = 0
        current_block_size = 0
        key_xor = None
        spi_xor = None
        ptxt_xor = None
        ctxt_xor = None
        key_sum = None
        spi_sum = None
        ptxt_sum = None
        ctxt_sum = None
        current_tests = []
        for l in f:
            l = l.rstrip()

            if not l:
                continue

            if all(c == ord(b'=') for c in l):
                generated = False
                if not current_tests and current_block_size:
                    current_tests = generate_mct(
                        seed, current_block_size, current_key_size,
                        test_per_bs, op)
                    generated = True

                kmask = (1 << current_key_size) - 1
                smask = (1 << 512) - 1

                for i, t in enumerate(current_tests):
                    if op == 'encryption':
                        r = hpc_dll.encrypt(
                            t['key'], current_key_size, backup,
                            t['ptxt'], current_block_size, t['spice'])
                        if generated:
                            t['ctxt'] = r
                        expected = t['ctxt']
                    else:
                        r = hpc_dll.decrypt(
                            t['key'], current_key_size, backup,
                            t['ctxt'], current_block_size, t['spice'])
                        if generated:
                            t['ptxt'] = r
                        expected = t['ptxt']

                    if r != expected:
                        print(
                            f'failed to do {op} for bs={current_block_size}, '
                            f'ks={current_key_size}, tid={i}. '
                            'expected [{}], got [{}]'.format(
                                ', '.join(f'{x:016x}' for x in expected),
                                ', '.join(f'{x:016x}' for x in r)))

                    if key_xor is not None:
                        key_xor ^= fold_data_words(t['key'])
                    if key_sum is not None:
                        key_sum -= fold_data_words(t['key'])

                    if spi_xor is not None:
                        spi_xor ^= fold_data_words(t['spice'])
                    if spi_sum is not None:
                        spi_sum -= fold_data_words(t['spice'])

                    if ptxt_xor is not None:
                        ptxt_xor ^= fold_data_words(t['ptxt'])
                    if ptxt_sum is not None:
                        ptxt_sum -= fold_data_words(t['ptxt'])

                    if ctxt_xor is not None:
                        ctxt_xor ^= fold_data_words(t['ctxt'])
                    if ctxt_sum is not None:
                        ctxt_sum -= fold_data_words(t['ctxt'])

                if key_sum:
                    key_sum &= (1 << current_key_size) - 1
                if spi_sum:
                    spi_sum &= (1 << 512) - 1

                bmask = (1 << ((current_block_size + 63) & ~63)) - 1
                if ptxt_sum:
                    ptxt_sum &= bmask
                if ctxt_sum:
                    ctxt_sum &= bmask

                if current_tests:
                    assert not key_xor
                    assert not key_sum
                    assert not spi_xor
                    assert not spi_sum
                    assert not ptxt_xor
                    assert not ptxt_sum
                    assert not ctxt_xor
                    assert not ctxt_sum
                    tests.setdefault(current_block_size, {}).setdefault(
                        current_key_size, []).extend(current_tests)

                current_block_size = current_key_size = 0
                key_xor = spi_xor = None
                ptxt_xor = ctxt_xor = None
                key_sum = spi_sum = None
                ptxt_sum = ctxt_sum = None
                current_tests.clear()
            else:
                keyword, data = l.split(b'=', 1)

                while data.endswith(b'\\'):
                    data = data[:-1] + f.readline().strip()

                if keyword == b'BLOCKSIZE':
                    assert current_block_size == 0 and not current_tests
                    current_block_size = int(data)
                elif keyword == b'KEYSIZE':
                    assert current_key_size == 0 and not current_tests
                    current_key_size = int(data)
                elif keyword == b'I':
                    assert int(data) == len(current_tests)
                    current_tests.append({})
                elif keyword in (b'PTXT', b'CTXT'):
                    keyword = keyword.decode('ascii').lower()
                    assert current_tests and keyword not in current_tests[-1]
                    parsed_words = parse_data_words(data)
                    assert len(parsed_words) == (current_block_size + 63) // 64
                    assert parsed_words[-1] < (1 << (64 if current_block_size % 64 == 0 else current_block_size % 64))
                    current_tests[-1][keyword] = parsed_words
                elif keyword == b'KEY':
                    assert current_tests and 'key' not in current_tests[-1]
                    parsed_words = parse_data_words(data)
                    assert len(parsed_words) == (current_key_size + 63) // 64
                    assert parsed_words[-1] < (1 << (64 if current_key_size % 64 == 0 else current_key_size % 64))
                    current_tests[-1]['key'] = parsed_words
                elif keyword == b'SPI':
                    assert current_tests and 'spice' not in current_tests[-1]
                    parsed_words = parse_data_words(data)
                    assert len(parsed_words) <= 8
                    current_tests[-1]['spice'] = parsed_words
                elif keyword == b'KEYXOR':
                    parsed_words = parse_data_words(data)
                    key_xor = fold_data_words(parsed_words)
                elif keyword == b'KEYSUM':
                    parsed_words = parse_data_words(data)
                    key_sum = fold_data_words(parsed_words)
                elif keyword == b'SPIXOR':
                    parsed_words = parse_data_words(data)
                    spi_xor = fold_data_words(parsed_words)
                elif keyword == b'SPISUM':
                    parsed_words = parse_data_words(data)
                    spi_sum = fold_data_words(parsed_words)
                elif keyword == b'PTXTXOR':
                    parsed_words = parse_data_words(data)
                    ptxt_xor = fold_data_words(parsed_words)
                elif keyword == b'PTXTSUM':
                    parsed_words = parse_data_words(data)
                    ptxt_sum = fold_data_words(parsed_words)
                elif keyword == b'CTXTXOR':
                    parsed_words = parse_data_words(data)
                    ctxt_xor = fold_data_words(parsed_words)
                elif keyword == b'CTXTSUM':
                    parsed_words = parse_data_words(data)
                    ctxt_sum = fold_data_words(parsed_words)
                else:
                    raise ValueError(f'Unknown keyword `{keyword}`')

    return {
        'op': op,
        'backup': backup,
        'tests': tests
    }


def parse_kat(inp):
    tests = {}

    with inp as f:
        header = [f.readline().rstrip() for _ in range(11)]

        op = 'encryption'
        backup = [0 for _ in range(6)]

        block_size = 128
        current_key_size = 0
        global_key = None
        global_pt = None
        current_tests = []
        for l in f:
            l = l.rstrip()

            if not l:
                continue

            if all(c == ord(b'=') for c in l):
                if current_tests:
                    tests.setdefault(block_size, {}).setdefault(
                        current_key_size, []).extend(current_tests)

                current_key_size = 0
                global_key = global_pt = None
                current_tests.clear()
            else:
                keyword, data = l.split(b'=', 1)

                while data.endswith(b'\\'):
                    data = data[:-1] + f.readline().strip()

                if keyword == b'KEYSIZE':
                    assert current_key_size == 0 and not current_tests
                    current_key_size = int(data)
                elif keyword == b'I':
                    assert int(data) == len(current_tests) + 1
                    current_tests.append({'spice': []})
                    if global_key:
                        current_tests[-1]['key'] = global_key
                    if global_pt:
                        current_tests[-1]['ptxt'] = global_pt
                elif keyword in (b'PT', b'CT'):
                    parsed_words = parse_data_words(data)
                    assert len(parsed_words) == block_size // 64
                    assert parsed_words[-1] < (1 << 64)
                    if not current_tests and keyword == b'PT':
                        global_pt = parsed_words
                    else:
                        keyword = 'ptxt' if keyword == b'PT' else 'ctxt'
                        assert current_tests and keyword not in current_tests[-1]
                        current_tests[-1][keyword] = parsed_words
                elif keyword == b'KEY':
                    parsed_words = parse_data_words(data)
                    assert len(parsed_words) == current_key_size // 64
                    assert parsed_words[-1] < (1 << 64)
                    if not current_tests:
                        global_key = parsed_words
                    else:
                        assert current_tests and 'key' not in current_tests[-1]
                        current_tests[-1]['key'] = parsed_words
                else:
                    raise ValueError(f'Unknown keyword `{keyword}`')

    return {
        'op': op,
        'backup': backup,
        'tests': tests
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=argparse.FileType('rb'))
    parser.add_argument('-l', '--library', type=HpcDll)
    parser.add_argument('-t', '--type', required=True, choices=['mct', 'kat'])
    parser.add_argument('-o', '--output', type=argparse.FileType('w'))
    args = parser.parse_args()

    global hpc_dll
    hpc_dll = args.library

    if args.type == 'mct':
        data = parse_mct(args.input)
    elif args.type == 'kat':
        data = parse_kat(args.input)

    # if not args.output:
    #     json.dump(data, sys.stdout)
    # else:
    #     json.dump(data, args.output)


if __name__ == '__main__':
    main()
