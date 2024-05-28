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


def generate_mct(seed, block_size, key_size, n):
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
            'ptxt': pt
        })

    return tests


def generate_and_sum_mct(seed, block_size, key_size, n):
    tests = generate_mct(seed, block_size, key_size, n)
    key_xor = functools.reduce(
        operator.xor,
        (fold_data_words(t['key']) for t in tests))
    key_sum = sum(
        fold_data_words(t['key']) for t in tests)
    spi_xor = functools.reduce(
        operator.xor,
        (fold_data_words(t['spice']) for t in tests))
    spi_sum = sum(
        fold_data_words(t['spice']) for t in tests)
    ptxt_xor = functools.reduce(
        operator.xor,
        (fold_data_words(t['ptxt']) for t in tests))
    ptxt_sum = sum(
        fold_data_words(t['ptxt']) for t in tests)
    return tests, (key_xor, key_sum), (spi_xor, spi_sum), (ptxt_xor, ptxt_sum)


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
        backup = [all_backup + ((backup >> (n * 4)) & 0xf) for n in range(1, 6)]

        current_key_size = 0
        current_block_size = 0
        key_xor = 0
        spi_xor = 0
        ptxt_xor = 0
        ctxt_xor = 0
        key_sum = 0
        spi_sum = 0
        ptxt_sum = 0
        ctxt_sum = 0
        generated = False
        current_tests = []
        for l in f:
            l = l.rstrip()

            if not l:
                continue

            if all(c == ord(b'=') for c in l):
                if current_tests:
                    tests.setdefault(current_block_size, {}).setdefault(
                        current_key_size, []).extend(current_tests)

                current_block_size = current_key_size = 0
                key_xor = spi_xor = 0
                ptxt_xor = ctxt_xor = 0
                key_sum = spi_sum = 0
                ptxt_sum = ctxt_sum = 0
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
                    if keyword == 'ptxt':
                        ptxt_xor ^= fold_data_words(parsed_words)
                        ptxt_sum += fold_data_words(parsed_words)
                    else:
                        ctxt_xor ^= fold_data_words(parsed_words)
                        ctxt_sum += fold_data_words(parsed_words)
                elif keyword == b'KEY':
                    assert current_tests and 'key' not in current_tests[-1]
                    parsed_words = parse_data_words(data)
                    assert len(parsed_words) == (current_key_size + 63) // 64
                    assert parsed_words[-1] < (1 << (64 if current_key_size % 64 == 0 else current_key_size % 64))
                    current_tests[-1]['key'] = parsed_words
                    key_xor ^= fold_data_words(parsed_words)
                    key_sum += fold_data_words(parsed_words)
                elif keyword == b'SPI':
                    assert current_tests and 'spice' not in current_tests[-1]
                    parsed_words = parse_data_words(data)
                    assert len(parsed_words) <= 8
                    current_tests[-1]['spice'] = parsed_words
                    spi_xor ^= fold_data_words(parsed_words)
                    spi_sum += fold_data_words(parsed_words)
                elif keyword == b'KEYXOR':
                    if not current_tests:
                        (current_tests,
                            (key_xor, key_sum),
                            (spi_xor, spi_sum),
                            (ptxt_xor, ptxt_sum)) = generate_and_sum_mct(
                                seed, current_block_size,
                                current_key_size, test_per_bs)
                        generated = True
                    parsed_words = parse_data_words(data)
                    assert key_xor == fold_data_words(parsed_words)
                elif keyword == b'KEYSUM':
                    if not current_tests:
                        (current_tests,
                            (key_xor, key_sum),
                            (spi_xor, spi_sum),
                            (ptxt_xor, ptxt_sum)) = generate_and_sum_mct(
                                seed, current_block_size,
                                current_key_size, test_per_bs)
                        generated = True
                    parsed_words = parse_data_words(data)
                    key_sum &= (1 << current_key_size) - 1
                    assert key_sum == fold_data_words(parsed_words)
                elif keyword == b'SPIXOR':
                    if not current_tests:
                        (current_tests,
                            (key_xor, key_sum),
                            (spi_xor, spi_sum),
                            (ptxt_xor, ptxt_sum)) = generate_and_sum_mct(
                                seed, current_block_size,
                                current_key_size, test_per_bs)
                        generated = True
                    parsed_words = parse_data_words(data)
                    assert spi_xor == fold_data_words(parsed_words)
                elif keyword == b'SPISUM':
                    if not current_tests:
                        (current_tests,
                            (key_xor, key_sum),
                            (spi_xor, spi_sum),
                            (ptxt_xor, ptxt_sum)) = generate_and_sum_mct(
                                seed, current_block_size,
                                current_key_size, test_per_bs)
                        generated = True
                    parsed_words = parse_data_words(data)
                    spi_sum &= (1 << 512) - 1
                    assert spi_sum == fold_data_words(parsed_words)
                elif keyword == b'PTXTXOR':
                    if not current_tests:
                        (current_tests,
                            (key_xor, key_sum),
                            (spi_xor, spi_sum),
                            (ptxt_xor, ptxt_sum)) = generate_and_sum_mct(
                                seed, current_block_size,
                                current_key_size, test_per_bs)
                        generated = True
                    parsed_words = parse_data_words(data)
                    assert ptxt_xor == fold_data_words(parsed_words)
                elif keyword == b'PTXTSUM':
                    if not current_tests:
                        (current_tests,
                            (key_xor, key_sum),
                            (spi_xor, spi_sum),
                            (ptxt_xor, ptxt_sum)) = generate_and_sum_mct(
                                seed, current_block_size,
                                current_key_size, test_per_bs)
                        generated = True
                    parsed_words = parse_data_words(data)
                    ptxt_sum &= (1 << ((current_block_size + 63) & ~63)) - 1
                    assert ptxt_sum == fold_data_words(parsed_words)
                elif keyword == b'CTXTXOR':
                    if current_tests and not generated:
                        parsed_words = parse_data_words(data)
                        assert ctxt_xor == fold_data_words(parsed_words)
                elif keyword == b'CTXTSUM':
                    if current_tests and not generated:
                        parsed_words = parse_data_words(data)
                        ctxt_sum &= (1 << ((current_block_size + 63) & ~63)) - 1
                        assert ctxt_sum == fold_data_words(parsed_words)
                else:
                    raise ValueError(f'Unknown keyword `{keyword}`')

    for bs, bv in tests.items():
        for ks, kv in bv.items():
            for i, t in enumerate(kv):
                if op == 'encryption':
                    r = hpc_dll.encrypt(
                        t['key'], ks, backup,
                        t['ptxt'], bs, t['spice'])
                    expected = t['ctxt']
                else:
                    r = hpc_dll.decrypt(
                        t['key'], ks, backup,
                        t['ctxt'], bs, t['spice'])
                    expected = t['ptxt']

                if r != expected:
                    print(
                        f'failed to do {op} for bs={bs}, ks={ks}, tid={i}. '
                        'expected [{}], got [{}]'.format(
                            ', '.join(f'{x:016x}' for x in expected),
                            ', '.join(f'{x:016x}' for x in r)))


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
