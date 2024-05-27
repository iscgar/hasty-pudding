#include <limits.h>
#include <stdbool.h>
#include <string.h>

#include "hpc/hpc.h"

#define HPC_PI19    UINT64_C(3141592653589793238)
#define HPC_E19     UINT64_C(2718281828459045235)
#define HPC_R220    UINT64_C(14142135623730950488)

#define HPC_STIR_PASSES 3
#define HPC_ROUND_COUNT 8

#if defined(_MSC_VER)
#   include <intrin.h>
#   define ROR64(x, b) _rotr64((x), (int)(b))
#   define ROL64(x, b) _rotl64((x), (int)(b))
#elif defined(__clang__)
#   define ROR64(x, b) __builtin_rotateright64((x), (b))
#   define ROL64(x, b) __builtin_rotateleft64((x), (b))
#else /* compiler without known rotation intrinsics */
#   define ROR64(x, b) (((uint64_t)(x) >> (size_t)(b)) | ((uint64_t)(x) << (size_t)(64 - (b))))
#   define ROL64(x, b) (((uint64_t)(x) << (size_t)(b)) | ((uint64_t)(x) >> (size_t)(64 - (b))))
#endif /* compiler without known rotation intrinsics */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

static size_t hpc_get_cipher_id(size_t data_bit_size);

/* The backup argument is only used by HPC-tiny, but is provided to every implmentation for consistency */
static void hpc_tiny_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_tiny_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_short_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_short_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_medium_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_medium_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_long_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_long_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_extended_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_extended_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);

int hpc_init(struct HpcState *state, const uint8_t *key, size_t key_bit_size)
{
    return hpc_init_with_backup(state, key, key_bit_size, NULL, 0);
}

int hpc_init_with_backup(struct HpcState *state, const uint8_t *key, size_t key_bit_size, const size_t *backup, size_t backup_size)
{
    if (state == NULL)
    {
        return false;
    }

    if (key == NULL && key_bit_size != 0)
    {
        return false;
    }

    if (backup != NULL && backup_size != HPC_BACKUP_SIZE)
    {
        return false;
    }

    if (backup != NULL)
    {
        if (backup[0] + HPC_STIR_PASSES > 64)
        {
            /* Would shift more than type size in the stirring step */
            return false;
        }

        memcpy(state->backup, backup, sizeof(state->backup));
    }
    else
    {
        memset(state->backup, 0, sizeof(state->backup));
    }

#if defined(HPC_USE_SINGLE_KX)
#   define cipher_id HPC_USE_SINGLE_KX
    uint64_t *KX = state->KX;
#else /* if !defined(HPC_USE_SINGLE_KX) */
    for (size_t cipher_id = (size_t)HPC_CIPHER_ID_TINY; cipher_id <= (size_t)HPC_CIPHER_ID_EXTENDED; ++cipher_id)
    {
        uint64_t *KX = state->KX[cipher_id - 1];
#endif /* !defined(HPC_USE_SINGLE_KX) */

        KX[0] = HPC_PI19 + cipher_id;
        KX[1] = HPC_E19 * key_bit_size;
        KX[2] = ROL64(HPC_R220, cipher_id);

        for (size_t i = 3; i < HPC_KX_SIZE; ++i)
        {
            KX[i] = KX[i-1] + (KX[i-2] ^ ROR64(KX[i-3], 23));
        }

        const uint8_t *left_key = key;
        size_t left_key_bits = key_bit_size;

        do
        {
            const size_t iteration_key_bits = left_key_bits >= (HPC_KX_SIZE * 64 / 2) ? (HPC_KX_SIZE * 64 / 2) : left_key_bits;
            const uint8_t *end_key = left_key + (iteration_key_bits / CHAR_BIT);

            for (size_t sh = 0, i = 0; left_key < end_key; ++left_key, sh = ((sh + CHAR_BIT) & 63), ++i)
            {
                KX[i / sizeof(uint64_t)] ^= ((uint64_t)*left_key) << sh;
            }

            if (iteration_key_bits & 7)
            {
                const size_t leftover_bits = iteration_key_bits & 7;
                const uint64_t v = (((uint64_t)*left_key++) & ((1u << leftover_bits) - 1)) << ((iteration_key_bits - leftover_bits) & 63);
                KX[(iteration_key_bits + CHAR_BIT - 1) / 64] ^= v;
            }

            uint64_t s0 = KX[248], s1 = KX[249], s2 = KX[250], s3 = KX[251];
            uint64_t s4 = KX[252], s5 = KX[253], s6 = KX[254], s7 = KX[255];

            for (size_t pass = 0; pass < HPC_STIR_PASSES + state->backup[0]; ++pass)
            {
                for (size_t ki = 0; ki < HPC_KX_SIZE; ++ki)
                {
                    s0 ^= (KX[ki] ^ KX[(ki + 83) & 255]) + KX[s0 & 255]; /* lossy, sometimes */
#ifndef HPC_WITHOUT_WAGNER_FIX
                    s2 += KX[ki]; /* added 1999-05-14 to fix Wagner equivalent key problem */
#endif /* !HPC_WITHOUT_WAGNER_FIX */
                    s1 += s0;
                    s3 ^= s2;
                    s5 -= s4;
                    s7 ^= s6;
                    s3 += s0 >> 13;
                    s4 ^= s1 << 11;
                    s5 ^= s3 << (s1 & 31);
                    s6 += s2 >> 17;
                    s7 |= s3 + s4; /* lossy */
                    s2 -= s5;      /* cross-link */
                    s0 -= s6 ^ ki;
                    s1 ^= s5 + HPC_PI19;
                    s2 += s7 >> pass;
                    s2 ^= s1;
                    s4 -= s3;
                    s6 ^= s5;
                    s0 += s7;
                    KX[ki] = s2 + s6;
                }
            }

            left_key_bits -= iteration_key_bits;
        } while (left_key_bits > 0);

#if defined(HPC_USE_SINGLE_KX)
#   undef cipher_id
#else /* if !defined(HPC_USE_SINGLE_KX) */
    }
#endif /* !defined(HPC_USE_SINGLE_KX) */

    return true;
}

int hpc_encrypt(struct HpcState *state, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t data_bit_size, const uint8_t *tweak, size_t tweak_bit_size)
{
    if (data_bit_size > 0)
    {
        if (state == NULL)
        {
            return false;
        }

        if (plaintext == NULL || o_ciphertext == NULL)
        {
            return false;
        }

        const size_t cipher_id = hpc_get_cipher_id(data_bit_size);

        void (*enc)(uint64_t *, const uint64_t *, const uint64_t *, size_t, uint64_t, size_t);

        switch (cipher_id)
        {
        case HPC_CIPHER_ID_TINY:
            enc = hpc_tiny_encrypt;
            break;

        case HPC_CIPHER_ID_SHORT:
            enc = hpc_short_encrypt;
            break;

        case HPC_CIPHER_ID_MEDIUM:
            enc = hpc_medium_encrypt;
            break;

        case HPC_CIPHER_ID_LONG:
            enc = hpc_long_encrypt;
            break;

        case HPC_CIPHER_ID_EXTENDED:
            return false; // TODO: implement

        default:
            return false;
        }

        uint64_t t[HPC_ROUND_COUNT] = { 0 };

        if (tweak != NULL)
        {
            if (tweak_bit_size > HPC_TWEAK_BIT_SIZE)
            {
                return false;
            }

            for (size_t i = 0, sh = 0; i < tweak_bit_size / CHAR_BIT; ++i, sh = ((sh + CHAR_BIT) & 63))
            {
                t[i / sizeof(uint64_t)] |= ((uint64_t)*tweak++) << sh;
            }

            if (tweak_bit_size & 7)
            {
                const size_t leftover_bits = tweak_bit_size & 7;
                const uint64_t v = (((uint64_t)*tweak++) & ((1u << leftover_bits) - 1)) << ((tweak_bit_size - leftover_bits) & 63);
                t[(tweak_bit_size + CHAR_BIT - 1) / 64] |= v;
            }
        }

        const uint64_t mask = (((UINT64_C(1) << ((data_bit_size - 1) & 63)) - 1) << 1) | 1;
        const size_t byte_limit = data_bit_size <= 512 ? (data_bit_size + CHAR_BIT - 1) / CHAR_BIT : (512 / CHAR_BIT);
        const size_t word_limit = (byte_limit - 1) & ~7;
        const size_t l64 = data_bit_size <= 128 ? (data_bit_size + 63) >> 6 : 8;

        uint64_t s[HPC_ROUND_COUNT];

        for (size_t i = 0; i < word_limit; )
        {
            s[i / sizeof(uint64_t)] = plaintext[i++];

            for (size_t j = 1, sh = CHAR_BIT; j < sizeof(uint64_t) && i < word_limit; ++j, sh += CHAR_BIT)
            {
                s[i / sizeof(uint64_t)] |= ((uint64_t)plaintext[i++]) << sh;
            }
        }

        for (size_t i = word_limit; i < byte_limit; )
        {
            s[l64 - 1] = plaintext[i++];

            for (size_t j = 1, sh = CHAR_BIT; j < sizeof(uint64_t) && i < byte_limit; ++j, sh += CHAR_BIT)
            {
                s[l64 - 1] |= ((uint64_t)plaintext[i++]) << sh;
            }
        }

        s[l64 - 1] &= mask;

#if defined(HPC_USE_SINGLE_KX)
        const uint64_t *KX = state->KX;
#else /* if !defined(HPC_USE_SINGLE_KX) */
        const uint64_t *KX = state->KX[cipher_id - 1];
#endif /* !defined(HPC_USE_SINGLE_KX) */

        for (size_t i = 0; i <= state->backup[cipher_id]; ++i)
        {
            s[0] += (uint64_t)i;

            for (size_t ki = 0, j = 0; j < HPC_ROUND_COUNT; ki += 64, ++j)
            {
                s[j] += KX[(data_bit_size + j) & 0xff];
            }

            s[l64 - 1] &= mask;

            (*enc)(s, t, KX, data_bit_size, mask, i);

            for (size_t ki = 0, j = 0; j < HPC_ROUND_COUNT; ki += 64, ++j)
            {
                s[j] += KX[(data_bit_size + HPC_ROUND_COUNT + j) & 0xff];
            }

            s[l64 - 1] &= mask;
        }

        for (size_t i = 0, j = 0; i < word_limit; ++j)
        {
            for (size_t sh = 0; sh < 64 && i < word_limit; sh += CHAR_BIT, ++i)
            {
                o_ciphertext[i] = (s[j] >> sh) & 0xff;
            }
        }

        for (size_t i = word_limit, j = 0; i < byte_limit; ++j)
        {
            for (size_t sh = 0; sh < 64 && i < byte_limit; sh += CHAR_BIT, ++i)
            {
                o_ciphertext[i] = (s[l64 - 1] >> sh) & 0xff;
            }
        }
    }

    return true;
}

int hpc_decrypt(struct HpcState *state, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t data_bit_size, const uint8_t *tweak, size_t tweak_bit_size)
{
    if (data_bit_size > 0)
    {
        if (state == NULL)
        {
            return false;
        }

        if (ciphertext == NULL || o_plaintext == NULL)
        {
            return false;
        }

        const size_t cipher_id = hpc_get_cipher_id(data_bit_size);

        void (*dec)(uint64_t *, const uint64_t *, const uint64_t *, size_t, uint64_t, size_t);

        switch (cipher_id)
        {
        case HPC_CIPHER_ID_TINY:
            dec = hpc_tiny_decrypt;
            break;

        case HPC_CIPHER_ID_SHORT:
            dec = hpc_short_decrypt;
            break;

        case HPC_CIPHER_ID_MEDIUM:
            dec = hpc_medium_decrypt;
            break;

        case HPC_CIPHER_ID_LONG:
            dec = hpc_long_decrypt;
            break;

        case HPC_CIPHER_ID_EXTENDED:
            return false; // TODO: implement

        default:
            return false;
        }

        uint64_t t[HPC_ROUND_COUNT] = { 0 };

        if (tweak != NULL)
        {
            if (tweak_bit_size > HPC_TWEAK_BIT_SIZE)
            {
                return false;
            }

            for (size_t i = 0, sh = 0; i < tweak_bit_size / CHAR_BIT; ++i, sh = ((sh + CHAR_BIT) & 63))
            {
                t[i / sizeof(uint64_t)] |= ((uint64_t)*tweak++) << sh;
            }

            if (tweak_bit_size & 7)
            {
                const size_t leftover_bits = tweak_bit_size & 7;
                const uint64_t v = (((uint64_t)*tweak++) & ((1u << leftover_bits) - 1)) << ((tweak_bit_size - leftover_bits) & 63);
                t[(tweak_bit_size + CHAR_BIT - 1) / 64] |= v;
            }
        }

        const uint64_t mask = (((UINT64_C(1) << ((data_bit_size - 1) & 63)) - 1) << 1) | 1;
        const size_t byte_limit = data_bit_size <= 512 ? (data_bit_size + CHAR_BIT - 1) / CHAR_BIT : (512 / CHAR_BIT);
        const size_t word_limit = (byte_limit - 1) & ~7;
        const size_t l64 = data_bit_size <= 64 ? 1 : (data_bit_size <= 128 ? 2 : 8);

        uint64_t s[HPC_ROUND_COUNT];

        for (size_t i = 0; i < word_limit; )
        {
            s[i / sizeof(uint64_t)] = ciphertext[i++];

            for (size_t j = 1, sh = CHAR_BIT; j < sizeof(uint64_t) && i < word_limit; ++j, sh += CHAR_BIT)
            {
                s[i / sizeof(uint64_t)] |= ((uint64_t)ciphertext[i++]) << sh;
            }
        }

        for (size_t i = word_limit; i < byte_limit; )
        {
            s[l64 - 1] = ciphertext[i++];

            for (size_t j = 1, sh = CHAR_BIT; j < sizeof(uint64_t) && i < byte_limit; ++j, sh += CHAR_BIT)
            {
                s[l64 - 1] |= ((uint64_t)ciphertext[i++]) << sh;
            }
        }

        s[l64 - 1] &= mask;

#if defined(HPC_USE_SINGLE_KX)
        const uint64_t *KX = state->KX;
#else /* if !defined(HPC_USE_SINGLE_KX) */
        const uint64_t *KX = state->KX[cipher_id - 1];
#endif /* !defined(HPC_USE_SINGLE_KX) */

        for (size_t i = state->backup[cipher_id] + 1; i-- > 0; )
        {
            for (size_t ki = 0, j = 0; j < HPC_ROUND_COUNT; ki += 64, ++j)
            {
                s[j] -= KX[(data_bit_size + HPC_ROUND_COUNT + j) & 0xff];
            }

            s[l64 - 1] &= mask;

            (*dec)(s, t, KX, data_bit_size, mask, i);

            for (size_t ki = 0, j = 0; j < HPC_ROUND_COUNT; ki += 64, ++j)
            {
                s[j] -= KX[(data_bit_size + j) & 0xff];
            }

            s[0] -= (uint64_t)i;
            s[l64 - 1] &= mask;
        }

        for (size_t i = 0, j = 0; i < word_limit; ++j)
        {
            for (size_t sh = 0; sh < 64 && i < word_limit; sh += CHAR_BIT, ++i)
            {
                o_plaintext[i] = (s[j] >> sh) & 0xff;
            }
        }

        for (size_t i = word_limit, j = 0; i < byte_limit; ++j)
        {
            for (size_t sh = 0; sh < 64 && i < byte_limit; sh += CHAR_BIT, ++i)
            {
                o_plaintext[i] = (s[l64 - 1] >> sh) & 0xff;
            }
        }
    }

    return true;
}

static const uint64_t Perma[16] = {
   0x243F6A8885A308D3 ^ 0,   0x13198A2E03707344 ^ 1,
   0xA4093822299F31D0 ^ 2,   0x082EFA98EC4E6C89 ^ 3,
   0x452821E638D01377 ^ 4,   0xBE5466CF34E90C6C ^ 5,
   0xC0AC29B7C97C50DD ^ 6,   0x9216D5D98979FB1B ^ 7,
   0xB8E1AFED6A267E96 ^ 8,   0xA458FEA3F4933D7E ^ 9,
   0x0D95748F728EB658 ^ 10,  0x7B54A41DC25A59B5 ^ 11,
   0xCA417918B8DB38EF ^ 12,  0xB3EE1411636FBC2A ^ 13,
   0x61D809CCFB21A991 ^ 14,  0x487CAC605DEC8032 ^ 15
};

#define PERM1   UINT64_C(0x324f6a850d19e7cb)   /* cycle notation (0B6D49851CF3E27)(A) */
#define PERM2   UINT64_C(0x2b7e1568adf09c43)   /* cycles (0396D7A5F2CEB14)(8) */

static uint8_t hpc_fib_fold(const uint64_t N[2])
{
    // A lossy version that avoids doing full 128-arithmetic
    // by only doing enough to ensure that we have a correct
    // result in the lowest bits, since we only use a single
    // bit anyway.

    uint64_t n = (N[0] + (N[1] >> 25));
    uint64_t n1 = (N[1] + (n < N[0] ? 1 : 0));

    n ^= (n1 << 9) | (n >> 55);
    n1 ^= (n1 >> 55);

    n += ((n1 << 30) | (n >> 34));
    n ^= (n >> 21);
    n += (n >> 13);
    n ^= (n >> 8);
    n += (n >> 5);
    n ^= (n >> 3);
    n += (n >> 2);
    n ^= (n >> 1);
    n += (n >> 1);

    return (uint8_t)(n & 1);
}

static void hpc_tiny_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0];

    switch (block_size)
    {
    case 1:
    case 2:
    case 3:
    case 4:
        {
            uint64_t tmp[2] = {
                KX[(block_size << 1) + 16] + KX[128] + backup,
                KX[(block_size << 1) + 17] + KX[129],
            };

            hpc_medium_encrypt(tmp, spice, KX, 128, UINT64_C(0xffffffffffffffff), 0);

            tmp[0] += KX[136];
            tmp[1] += KX[137];

            switch (block_size)
            {
            case 1:
                tmp[0] += tmp[1];
                s0 ^= hpc_fib_fold(tmp);
                break;

            case 2:
            case 3:
                {
                    for (size_t ri = 0; ri < 2; ++ri)
                    {
                        uint64_t t = tmp[ri];

                        for (size_t bi = 0; bi < 64; bi += (block_size << 1))
                        {
                            s0 ^= t;
                            t >>= block_size;
                            s0 += t;
                            s0 = (s0 << 1) | ((s0 & mask) >> (block_size - 1));
                            t >>= block_size;
                        }
                    }
                }
                break;

            case 4:
                {
                    for (size_t ri = 0; ri < 2; ++ri)
                    {
                        uint64_t t = tmp[ri];

                        for (size_t bi = 0; bi < 64; bi += 8)
                        {
                            // The spec claims that for 4 bits we first
                            // add 4 bits from `t`, then XOR.
                            // However, the test vectors do not agree, and
                            // only work if XOR is used first and then
                            // addition, which is why we deviate.
                            s0 ^= t;
                            t >>= 4;
                            s0 = (PERM1 >> ((s0 & 15) << 2));
                            s0 += t;
                            s0 = (PERM2 >> ((s0 & 15) << 2));
                            t >>= 4;
                        }
                    }
                }
                break;
            }
        }
        break;

    case 5:
    case 6:
        {
            const size_t tmp_bs = 96 << (block_size - 4);
            uint64_t tmp[HPC_ROUND_COUNT]; // Need to be this large for the call to hpc_long_encrypt() below

            const uint8_t bs_base = (uint8_t)tmp_bs;
            const size_t l64 = (tmp_bs >> 6) - 1;

            for (size_t i = 0; i < (tmp_bs >> 6) - 1; ++i)
            {
                tmp[i] = KX[(block_size << 1) + 16 + i] + KX[(bs_base + i)];
            }
            tmp[HPC_ROUND_COUNT - 1] = KX[(block_size << 1) + 16 + l64] + KX[(bs_base + HPC_ROUND_COUNT - 1)];
            tmp[0] += backup;

            hpc_long_encrypt(tmp, spice, KX, tmp_bs, UINT64_C(0xffffffffffffffff), 0);

            for (size_t i = 0; i < (tmp_bs >> 6) - 1; ++i)
            {
                tmp[i] += KX[(bs_base + HPC_ROUND_COUNT + i)];
            }
            tmp[(tmp_bs >> 6) - 1] = tmp[HPC_ROUND_COUNT - 1] + KX[(bs_base + (HPC_ROUND_COUNT << 1) - 1)];

            const uint8_t pmask = (uint8_t)mask ^ 15;

            for (size_t ri = 0; ri < (tmp_bs >> 6); ++ri)
            {
                uint64_t t = tmp[ri];

                for (size_t bi = 0; bi < (7 - (block_size - 5)); ++bi)
                {
                    s0 ^= t;
                    s0 = (s0 & pmask) | ((PERM1 >> ((s0 & 15) << 2)) & 15);
                    s0 ^= (s0 >> 3);
                    t >>= block_size;
                    s0 += t;
                    s0 = (s0 & pmask) | ((PERM2 >> ((s0 & 15) << 2)) & 15);
                    t >>= block_size - 1;
                }
            }
        }
        break;

    default: // 7..35
        {
            const size_t LBH = (block_size + 1) >> 1;

            uint64_t tmp[HPC_ROUND_COUNT + 2] = {
                (spice[0] ^ KX[(block_size << 2) + 16]) + KX[0] + backup,
                (spice[1] ^ KX[(block_size << 2) + 17]) + KX[1],
                (spice[2] ^ KX[(block_size << 2) + 18]) + KX[2],
                (spice[3] ^ KX[(block_size << 2) + 19]) + KX[3],
                (spice[4] ^ KX[(block_size << 2) + 20]) + KX[4],
                (spice[5] ^ KX[(block_size << 2) + 21]) + KX[5],
                (spice[6] ^ KX[(block_size << 2) + 22]) + KX[6],
                (spice[7] ^ KX[(block_size << 2) + 23]) + KX[7],
            };
            const uint64_t zspice[HPC_ROUND_COUNT] = { 0 };

            hpc_long_encrypt(tmp, zspice, KX, 512, UINT64_C(0xffffffffffffffff), 0);

            for (size_t i = 0; i < HPC_ROUND_COUNT; ++i)
            {
                tmp[i] += KX[HPC_ROUND_COUNT + i];
            }

            tmp[8] = tmp[9] = tmp[7];

            for (size_t ri = 0; ri < HPC_ROUND_COUNT; ++ri)
            {
                tmp[8] += ((tmp[8] << 21) + (tmp[8] >> 13)) ^ (tmp[ri] + KX[ri + 16]);
                tmp[9] ^= tmp[8];
            }

            if (block_size < 16)
            {
                for (size_t ri = 0; ri < HPC_ROUND_COUNT + 2; ++ri)
                {
                    uint64_t t = tmp[ri];

                    for (size_t bi = 0; bi < 64; bi += (block_size << 1))
                    {
                        s0 += t;
                        s0 ^= (KX[(16 * ri) + (s0 & 15)] << 4);
                        s0 = ((s0 & mask) >> 4) | (s0 << (block_size - 4));
                        s0 ^= (s0 & mask) >> LBH;
                        s0 ^= t >> block_size;
                        s0 += s0 << (LBH + 2);
                        s0 ^= Perma[s0 & 15];
                        s0 += s0 << LBH;
                        t >>= (block_size << 1);
                    }
                }
            }
            else
            {
                for (size_t ri = 0; ri < HPC_ROUND_COUNT + 2; ++ri)
                {
                    uint64_t t = tmp[ri];

                    for (size_t bi = 0; bi < 64; bi += block_size)
                    {
                        s0 += t;
                        s0 ^= (KX[s0 & 255] << 8);
                        s0 = ((s0 & mask) >> 8) | (s0 << (block_size - 8));
                        t >>= block_size;
                    }
                }
            }
        }
        break;
    }

    state[0] = s0;
}

static const uint64_t Permai[16] = {
   0xA4093822299F31D0 ^ 2,   0x61D809CCFB21A991 ^ 14,
   0x487CAC605DEC8032 ^ 15,  0x243F6A8885A308D3 ^ 0,
   0x13198A2E03707344 ^ 1,   0x7B54A41DC25A59B5 ^ 11,
   0xB8E1AFED6A267E96 ^ 8,   0x452821E638D01377 ^ 4,
   0x0D95748F728EB658 ^ 10,  0x082EFA98EC4E6C89 ^ 3,
   0xB3EE1411636FBC2A ^ 13,  0x9216D5D98979FB1B ^ 7,
   0xBE5466CF34E90C6C ^ 5,   0xC0AC29B7C97C50DD ^ 6,
   0xA458FEA3F4933D7E ^ 9,   0xCA417918B8DB38EF ^ 12
};

#define PERM1I  UINT64_C(0xc3610a492b8dfe57)   /* inverse of PERM1 */
#define PERM2I  UINT64_C(0x5c62e738d9a10fb4)   /* inverse of PERM2 */

static void hpc_tiny_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0];

    switch (block_size)
    {
    case 1:
    case 2:
    case 3:
    case 4:
        {
            uint64_t tmp[2] = {
                KX[(block_size << 1) + 16] + KX[128] + backup,
                KX[(block_size << 1) + 17] + KX[129],
            };

            hpc_medium_encrypt(tmp, spice, KX, 128, UINT64_C(0xffffffffffffffff), 0);

            tmp[0] += KX[136];
            tmp[1] += KX[137];

            switch (block_size)
            {
            case 1:
                tmp[0] += tmp[1];
                s0 ^= hpc_fib_fold(tmp);
                break;

            case 2:
            case 3:
                {
                    for (size_t ri = 2; ri-- > 0; )
                    {
                        const uint64_t t = tmp[ri];

                        for (size_t bi = (64 + (block_size << 1) - 1) / (block_size << 1); bi-- > 0; )
                        {
                            const uint8_t v = (uint8_t)(t >> (bi * (block_size << 1)));
                            s0 = ((s0 & mask) >> 1) | (s0 << (block_size - 1));
                            s0 -= (v >> block_size);
                            s0 ^= v;
                        }
                    }
                }
                break;

            case 4:
                {
                    for (size_t ri = 2; ri-- > 0; )
                    {
                        const uint64_t t = tmp[ri];

                        for (size_t bi = 64; bi > 0; bi -= 8)
                        {
                            const uint8_t v = (uint8_t)(t >> (bi - 8));
                            s0 = (PERM2I >> ((s0 & 15) << 2));
                            s0 -= v >> 4;
                            s0 = (PERM1I >> ((s0 & 15) << 2));
                            s0 ^= v;
                        }
                    }
                }
                break;
            }
        }
        break;

    case 5:
    case 6:
        {
            const size_t tmp_bs = 96 << (block_size - 4);
            uint64_t tmp[HPC_ROUND_COUNT]; // Need to be this large for the call to hpc_long_encrypt() below

            const uint8_t bs_base = (uint8_t)tmp_bs;
            const size_t l64 = (tmp_bs >> 6) - 1;

            for (size_t i = 0; i < (tmp_bs >> 6) - 1; ++i)
            {
                tmp[i] = KX[(block_size << 1) + 16 + i] + KX[(bs_base + i)];
            }
            tmp[HPC_ROUND_COUNT - 1] = KX[(block_size << 1) + 16 + l64] + KX[(bs_base + HPC_ROUND_COUNT - 1)];
            tmp[0] += backup;

            hpc_long_encrypt(tmp, spice, KX, tmp_bs, UINT64_C(0xffffffffffffffff), 0);

            for (size_t i = 0; i < (tmp_bs >> 6) - 1; ++i)
            {
                tmp[i] += KX[(bs_base + HPC_ROUND_COUNT + i)];
            }
            tmp[(tmp_bs >> 6) - 1] = tmp[HPC_ROUND_COUNT - 1] + KX[(bs_base + (HPC_ROUND_COUNT << 1) - 1)];

            const uint8_t pmask = (uint8_t)mask ^ 15;

            for (size_t ri = (tmp_bs >> 6); ri-- > 0; )
            {
                const uint64_t t = tmp[ri];

                for (size_t bi = (7 - (block_size - 5)); bi-- > 0; )
                {
                    const uint16_t v = (uint16_t)(t >> (bi * ((block_size << 1) - 1)));
                    s0 = (s0 & pmask) | ((PERM2I >> ((s0 & 15) << 2)) & 15);
                    s0 -= (v >> block_size);
                    s0 ^= ((s0 & mask) >> 3);
                    s0 = (s0 & pmask) | ((PERM1I >> ((s0 & 15) << 2)) & 15);
                    s0 ^= v;
                }
            }
        }
        break;

    default: // 7..35
        {
            const size_t LBH = (block_size + 1) >> 1;

            uint64_t tmp[HPC_ROUND_COUNT + 2] = {
                (spice[0] ^ KX[(block_size << 2) + 16]) + KX[0] + backup,
                (spice[1] ^ KX[(block_size << 2) + 17]) + KX[1],
                (spice[2] ^ KX[(block_size << 2) + 18]) + KX[2],
                (spice[3] ^ KX[(block_size << 2) + 19]) + KX[3],
                (spice[4] ^ KX[(block_size << 2) + 20]) + KX[4],
                (spice[5] ^ KX[(block_size << 2) + 21]) + KX[5],
                (spice[6] ^ KX[(block_size << 2) + 22]) + KX[6],
                (spice[7] ^ KX[(block_size << 2) + 23]) + KX[7],
            };
            const uint64_t zspice[HPC_ROUND_COUNT] = { 0 };

            hpc_long_encrypt(tmp, zspice, KX, 512, UINT64_C(0xffffffffffffffff), 0);

            for (size_t i = 0; i < HPC_ROUND_COUNT; ++i)
            {
                tmp[i] += KX[HPC_ROUND_COUNT + i];
            }

            tmp[8] = tmp[9] = tmp[7];

            for (size_t ri = 0; ri < HPC_ROUND_COUNT; ++ri)
            {
                tmp[8] += ((tmp[8] << 21) + (tmp[8] >> 13)) ^ (tmp[ri] + KX[ri + 16]);
                tmp[9] ^= tmp[8];
            }

            if (block_size < 16)
            {
                for (size_t ri = HPC_ROUND_COUNT + 2; ri-- > 0; )
                {
                    const uint64_t t = tmp[ri];

                    for (size_t bi = (64 + (block_size << 1) - 1) / (block_size << 1); bi-- > 0; )
                    {
                        const uint32_t v = (uint32_t)(t >> (bi * (block_size << 1)));
                        s0 -= s0 << LBH;
                        s0 ^= Permai[s0 & 15];
                        s0 -= s0 << (LBH + 2);
                        s0 ^= v >> block_size;
                        s0 ^= (s0 & mask) >> LBH;
                        s0 = (s0 << 4) | ((s0 & mask) >> (block_size - 4));
                        s0 ^= (KX[(16 * ri) + (s0 & 15)] << 4);
                        s0 -= v;
                    }
                }
            }
            else
            {
                for (size_t ri = HPC_ROUND_COUNT + 2; ri-- > 0; )
                {
                    uint64_t t = tmp[ri];

                    for (size_t bi = (64 + block_size - 1) / block_size; bi-- > 0; )
                    {
                        s0 = (s0 << 8) | ((s0 & mask) >> (block_size - 8));
                        s0 ^= (KX[s0 & 255] << 8);
                        s0 -= (t >> (bi * block_size));
                    }
                }
            }
        }
        break;
    }

    state[0] = s0;
}

static const uint64_t Permb[16] = {
   0xB7E151628AED2A6A - 0,   0xBF7158809CF4F3C7 - 1,
   0x62E7160F38B4DA56 - 2,   0xA784D9045190CFEF - 3,
   0x324E7738926CFBE5 - 4,   0xF4BF8D8D8C31D763 - 5,
   0xDA06C80ABB1185EB - 6,   0x4F7C7B5757F59584 - 7,
   0x90CFD47D7C19BB42 - 8,   0x158D9554F7B46BCE - 9,
   0x8A9A276BCFBFA1C8 - 10,  0xE5AB6ADD835FD1A0 - 11,
   0x86D1BF275B9B241D - 12,  0xF0D3D37BE67008E1 - 13,
   0x0FF8EC6D31BEB5CC - 14,  0xEB64749A47DFDFB9 - 15
};

static void hpc_short_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup)
{
    const size_t LBH = (block_size + 1) >> 1;
    const size_t LBQ = (LBH + 1) >> 1;
    const size_t LBT = ((block_size + LBQ) >> 2) + 2;
    const size_t GAP = 64 - block_size;

    uint64_t s0 = state[0];
    (void)backup;

    for (size_t ri = 0; ri < HPC_ROUND_COUNT; ++ri)
    {
        uint64_t k = KX[s0 & 255] + spice[ri], t;

        s0 += (k << 8);
        s0 ^= (k >> GAP) & ~255;
        s0 += s0 << (LBH + ri);
        t = spice[ri ^ 7];
        s0 ^= t;
        s0 -= (t >> (GAP + ri));
        s0 += (t >> 13);
        s0 &= mask;

        s0 ^= (s0 >> LBH);
        t = s0 & 255;
        k = KX[t] ^ spice[ri ^ 4];
        k = KX[(t + (3 * ri) + 1) & 255] + ROR64(k, 23);
        s0 ^= (k << 8);
        s0 -= (k >> GAP) & ~255;
        s0 -= (s0 << LBH);
        t = spice[ri ^ 1] ^ (HPC_PI19 + block_size);
        s0 += t << 3;
        s0 ^= t >> (GAP + 2);
        s0 -= t;
        s0 &= mask;

        s0 ^= (s0 >> LBQ);
        s0 += Permb[s0 & 15];
        t = spice[ri ^ 2];
        s0 ^= t >> (GAP + 4);
        s0 += (s0 << (LBT + (s0 & 15)));
        s0 += t;
        s0 &= mask;

        s0 ^= (s0 >> LBH);
        s0 &= mask;
    }

    state[0] = s0;
}

static const uint64_t Permbi[16] = {
   0xE5AB6ADD835FD1A0 - 11,  0xF0D3D37BE67008E1 - 13,
   0x90CFD47D7C19BB42 - 8,   0xF4BF8D8D8C31D763 - 5,
   0x4F7C7B5757F59584 - 7,   0x324E7738926CFBE5 - 4,
   0x62E7160F38B4DA56 - 2,   0xBF7158809CF4F3C7 - 1,
   0x8A9A276BCFBFA1C8 - 10,  0xEB64749A47DFDFB9 - 15,
   0xB7E151628AED2A6A - 0,   0xDA06C80ABB1185EB - 6,
   0x0FF8EC6D31BEB5CC - 14,  0x86D1BF275B9B241D - 12,
   0x158D9554F7B46BCE - 9,   0xA784D9045190CFEF - 3
};

static void hpc_short_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup)
{
    const size_t LBH = (block_size + 1) >> 1;
    const size_t LBQ = (LBH + 1) >> 1;
    const size_t LBT = ((block_size + LBQ) >> 2) + 2;
    const size_t GAP = 64 - block_size;

    uint64_t s0 = state[0];
    (void)backup;

    for (size_t ri = HPC_ROUND_COUNT; ri-- > 0; )
    {
        uint64_t k, t = spice[ri ^ 2];

        s0 ^= (s0 >> LBH);
        s0 -= t;
        k = (s0 << (LBT + (s0 & 15)));
        s0 -= (s0 - k) << (LBT + (s0 & 15));
        s0 ^= (t >> (GAP + 4));
        s0 -= Permbi[s0 & 15];
        s0 &= mask;

        s0 ^= (s0 >> LBQ);
        s0 ^= (s0 >> (LBQ << 1));
        t = spice[ri ^ 1] ^ (HPC_PI19 + block_size);
        s0 += t;
        s0 ^= (t >> (GAP + 2));
        s0 -= (t << 3);
        s0 += s0 << LBH;
        t = s0 & 255;
        k = KX[t] ^ spice[ri ^ 4];
        k = KX[(t + (3 * ri) + 1) & 255] + ROR64(k, 23);
        s0 += (k >> GAP) & ~255;
        s0 ^= (k << 8);
        s0 &= mask;

        s0 ^= s0 >> LBH;
        t = spice[ri ^ 7];
        s0 -= (t >> 13);
        s0 += (t >> (GAP + ri));
        s0 ^= t;
        s0 -= s0 << (LBH + ri);
        k = KX[s0 & 255] + spice[ri];
        s0 ^= (k >> GAP) & ~255;
        s0 -= (k << 8);
        s0 &= mask;
    }

    state[0] = s0;
}

static void hpc_medium_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0], s1 = state[1];
    (void)backup;

    for (size_t ri = 0; ri < HPC_ROUND_COUNT; ++ri)
    {
        uint64_t k = KX[s0 & 255], t, kk;

        s1 += k;
        s0 ^= (k << 8);
        s1 ^= s0;
        s1 &= mask;

        s0 -= (s1 >> 11);
        s0 ^= s1 << 2;
        s0 -= spice[ri ^ 4];
        s0 += (s0 << 32) ^ (HPC_PI19 + block_size);
        s0 ^= (s0 >> 17);
        s0 ^= (s0 >> 34);
        t = spice[ri];
        s0 ^= t;
        s0 += (t << 5);
        t >>= 4;
        s1 += t;
        s0 ^= t;
        s0 += s0 << (22 + (s0 & 31));
        s0 ^= s0 >> 23;
        s0 -= spice[ri ^ 7];

        t = s0 & 255;
        k = KX[t];
        kk = KX[(t + (3 * ri) + 1) & 255];

        s1 ^= k;
        s0 ^= (kk << 8);
        kk ^= k;
        s1 += (kk >> 5);
        s0 -= (kk << 12);
        s0 ^= (kk & ~255);
        s1 += s0;
        s1 &= mask;

        s0 += (s1 << 3);
        s0 ^= spice[ri ^ 2];
        s0 += KX[block_size + ri + 16];
        s0 += (s0 << 22);
        s0 ^= (s1 >> 4);
        s0 += spice[ri ^ 1];
        s0 ^= (s0 >> (ri + 33));
    }

    state[0] = s0;
    state[1] = s1;
}

static void hpc_medium_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0], s1 = state[1];
    (void)backup;

    for (size_t ri = HPC_ROUND_COUNT; ri-- > 0; )
    {
        uint64_t k, t, kk;

        s0 ^= (s0 >> (ri + 33));
        s0 -= spice[ri ^ 1];
        s0 ^= (s1 >> 4);
        t = s0 - (s0 << 22);
        s0 -= (t << 22);
        s0 -= KX[block_size + ri + 16];
        s0 ^= spice[ri ^ 2];
        s0 -= (s1 << 3);
        s1 -= s0;

        t = s0 & 255;
        k = KX[t];
        kk = KX[(t + (3 * ri) + 1) & 255];
        kk ^= k;

        s0 ^= (kk & ~255);
        s0 += (kk << 12);
        s1 -= (kk >> 5);
        kk ^= k;
        s0 ^= (kk << 8);
        s1 ^= k;

        s0 += spice[ri ^ 7];
        s0 ^= (s0 >> 23);
        s0 ^= (s0 >> 46);
        t = (s0 << (22 + (s0 & 31)));
        s0 -= (s0 - t) << (22 + (s0 & 31));
        t = (spice[ri] >> 4);
        s0 ^= t;
        s1 -= t;
        t = spice[ri];
        s0 -= (t << 5);
        s0 ^= t;
        s0 ^= (s0 >> 17);
        t = (s0 - (HPC_PI19 + block_size));
        s0 -= ((t << 32) ^ (HPC_PI19 + block_size));
        s0 += spice[ri ^ 4];
        s1 &= mask;

        s0 ^= (s1 << 2);
        s0 += (s1 >> 11);
        s1 ^= s0;
        k = KX[s0 & 255];
        s0 ^= (k << 8);
        s1 -= k;
        s1 &= mask;
    }

    state[0] = s0;
    state[1] = s1;
}

static void hpc_long_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0], s1 = state[1], s2 = state[2], s3 = state[3];
    uint64_t s4 = state[4], s5 = state[5], s6 = state[6], s7 = state[7];
    (void)backup;

    for (size_t ri = 0; ri < HPC_ROUND_COUNT; ++ri)
    {
        uint64_t t = (s0 & 255), k = KX[t], kk = KX[(t + (3 * ri) + 1) & 255];

        s1 += k;
        s0 ^= (kk << 8);
        kk ^= k;
        s1 += (kk >> 5);
        s0 -= (kk << 12);
        s7 += kk;
        s7 ^= s0;
        s7 &= mask;

        s1 += s7;
        s1 ^= (s7 << 13);
        s0 -= (s7 >> 11);
        s0 += spice[ri];
        s1 ^= spice[ri ^ 1];
        s0 += (s1 << (ri + 9));
        s1 += (s0 >> 3) ^ (HPC_PI19 + block_size);
        s0 ^= (s1 >> 4);
        s0 += spice[ri ^ 2];
        t = spice[ri ^ 4];
        s1 += t;
        s1 ^= (t >> 3);
        s1 -= (t << 5);
        s0 ^= s1;

        if (block_size > 192)
        {
            if (block_size > 256)
            {
                if (block_size > 320)
                {
                    if (block_size > 384)
                    {
                        if (block_size > 448)
                        {
                            s6 += s0;
                            s6 ^= (s3 << 11);
                            s1 += (s6 >> 13);
                            s6 += (s5 << 7);
                            s4 ^= s6;
                        }

                        s5 ^= s1;
                        s5 += (s4 << 15);
                        s0 -= (s5 >> 7);
                        s5 ^= (s3 >> 9);
                        s2 ^= s5;
                    }

                    s4 -= s2;
                    s4 ^= (s1 >> 10);
                    s0 ^= (s4 << 3);
                    s4 -= (s2 << 6);
                    s3 += s4;
                }

                s3 ^= s2;
                s3 -= (s0 >> 7);
                s2 ^= (s3 << 15);
                s3 ^= (s1 << 5);
                s1 += s3;
            }

            s2 ^= s1;
            s2 += (s0 << 13);
            s1 -= (s2 >> 5);
            s2 -= (s1 >> 8);
            s0 ^= s2;
        }

        s1 ^= KX[(block_size + (ri << 5) + 17) & 255];
        s1 += (s0 << 19);
        s0 -= (s1 >> 27);
        s1 ^= spice[ri ^ 7];
        s7 -= s1;
        s0 += (s1 & (s1 >> 5));
        s1 ^= (s0 >> (s0 & 31));
        s0 ^= KX[s1 & 255];
    }

    state[0] = s0;
    state[1] = s1;
    state[2] = s2;
    state[3] = s3;
    state[4] = s4;
    state[5] = s5;
    state[6] = s6;
    state[7] = s7;
}

static void hpc_long_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0], s1 = state[1], s2 = state[2], s3 = state[3];
    uint64_t s4 = state[4], s5 = state[5], s6 = state[6], s7 = state[7];
    (void)backup;

    for (size_t ri = HPC_ROUND_COUNT; ri-- > 0; )
    {
        uint64_t t, k, kk;

        s0 ^= KX[s1 & 255];
        s1 ^= (s0 >> (s0 & 31));
        s0 -= s1 & (s1 >> 5);
        s7 += s1;
        s7 &= mask;

        s1 ^= spice[ri ^ 7];
        s0 += (s1 >> 27);
        s1 -= (s0 << 19);
        s1 ^= KX[(block_size + (ri << 5) + 17) & 255];

        if (block_size > 192)
        {
            s0 ^= s2;
            s2 += (s1 >> 8);
            s1 += (s2 >> 5);
            s2 -= (s0 << 13);
            s2 ^= s1;

            if (block_size > 256)
            {
                s1 -= s3;
                s3 ^= (s1 << 5);
                s2 ^= (s3 << 15);
                s3 += (s0 >> 7);
                s3 ^= s2;

                if (block_size > 320)
                {
                    s3 -= s4;
                    s4 += (s2 << 6);
                    s0 ^= (s4 << 3);
                    s4 ^= (s1 >> 10);
                    s4 += s2;

                    if (block_size > 384)
                    {
                        s2 ^= s5;
                        s5 ^= (s3 >> 9);
                        s0 += (s5 >> 7);
                        s5 -= (s4 << 15);
                        s5 ^= s1;

                        if (block_size > 448)
                        {
                            s4 ^= s6;
                            s6 -= (s5 << 7);
                            s1 -= (s6 >> 13);
                            s6 ^= (s3 << 11);
                            s6 -= s0;
                        }
                    }
                }
            }
        }

        s0 ^= s1;
        t = spice[ri ^ 4];
        s1 += (t << 5);
        s1 ^= (t >> 3);
        s1 -= t;
        s0 -= spice[ri ^ 2];
        s0 ^= (s1 >> 4);
        s1 -= (s0 >> 3) ^ (HPC_PI19 + block_size);
        s0 -= (s1 << (ri + 9));
        s1 ^= spice[ri ^ 1];
        s0 -= spice[ri];
        s0 += (s7 >> 11);
        s1 ^= (s7 << 13);
        s1 -= s7;

        t = (s0 & 255);
        k = KX[t];
        kk = KX[(t + (3 * ri) + 1) & 255] ^ k;

        s7 ^= s0;
        s7 -= kk;
        s0 += (kk << 12);
        s1 -= (kk >> 5);
        kk ^= k;
        s0 ^= (kk << 8);
        s1 -= k;
    }

    state[0] = s0;
    state[1] = s1;
    state[2] = s2;
    state[3] = s3;
    state[4] = s4;
    state[5] = s5;
    state[6] = s6;
    state[7] = s7;
}

static void hpc_extended_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);
static void hpc_extended_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, size_t block_size, uint64_t mask, size_t backup);

static size_t hpc_get_cipher_id(size_t data_bit_size)
{
    if (data_bit_size >= 0 && data_bit_size <= 35)
    {
        return HPC_CIPHER_ID_TINY;
    }
    if (data_bit_size >= 36 && data_bit_size <= 64)
    {
        return HPC_CIPHER_ID_SHORT;
    }
    if (data_bit_size >= 65 && data_bit_size <= 128)
    {
        return HPC_CIPHER_ID_MEDIUM;
    }
    if (data_bit_size >= 129 && data_bit_size <= 512)
    {
        return HPC_CIPHER_ID_LONG;
    }
    return HPC_CIPHER_ID_EXTENDED;
}

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */
