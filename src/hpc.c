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
static void hpc_tiny_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup);
static void hpc_tiny_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup);
static void hpc_short_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup);
static void hpc_short_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup);
static void hpc_medium_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup);
static void hpc_medium_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup);
static void hpc_long_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup);
static void hpc_long_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup);
static void hpc_extended_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup);
static void hpc_extended_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup);

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

        void (*enc)(uint64_t *, const uint64_t *, const uint64_t *, const uint8_t *, uint8_t *, size_t, uint64_t, size_t);

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
            enc = hpc_extended_encrypt;
            break;

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

        if (data_bit_size < 512)
        {
            s[l64 - 1] &= mask;
        }

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

            if (data_bit_size < 512)
            {
                s[l64 - 1] &= mask;
            }

            (*enc)(s, t, KX, plaintext, o_ciphertext, data_bit_size, mask, i);

            for (size_t ki = 0, j = 0; j < HPC_ROUND_COUNT; ki += 64, ++j)
            {
                s[j] += KX[(data_bit_size + HPC_ROUND_COUNT + j) & 0xff];
            }

            if (data_bit_size < 512)
            {
                s[l64 - 1] &= mask;
            }
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

        void (*dec)(uint64_t *, const uint64_t *, const uint64_t *, const uint8_t *, uint8_t *, size_t, uint64_t, size_t);

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
            dec = hpc_extended_decrypt;
            break;

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

        if (data_bit_size < 512)
        {
            s[l64 - 1] &= mask;
        }

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

            if (data_bit_size < 512)
            {
                s[l64 - 1] &= mask;
            }

            (*dec)(s, t, KX, ciphertext, o_plaintext, data_bit_size, mask, i);

            for (size_t ki = 0, j = 0; j < HPC_ROUND_COUNT; ki += 64, ++j)
            {
                s[j] -= KX[(data_bit_size + j) & 0xff];
            }

            s[0] -= (uint64_t)i;
            if (data_bit_size < 512)
            {
                s[l64 - 1] &= mask;
            }
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

static void hpc_tiny_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0];
    (void)plaintext, (void)o_ciphertext;

    switch (block_size)
    {
    case 1:
    case 2:
    case 3:
    case 4:
        {
            uint64_t tmp[2] = {
                KX[(block_size << 1) + 16] + KX[128] + (uint64_t)backup,
                KX[(block_size << 1) + 17] + KX[129],
            };

            hpc_medium_encrypt(tmp, spice, KX, NULL, NULL, 128, UINT64_C(0xffffffffffffffff), 0);

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
            tmp[0] += (uint64_t)backup;

            hpc_long_encrypt(tmp, spice, KX, NULL, NULL, tmp_bs, UINT64_C(0xffffffffffffffff), 0);

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
                (spice[0] ^ KX[(block_size << 2) + 16]) + KX[0] + (uint64_t)backup,
                (spice[1] ^ KX[(block_size << 2) + 17]) + KX[1],
                (spice[2] ^ KX[(block_size << 2) + 18]) + KX[2],
                (spice[3] ^ KX[(block_size << 2) + 19]) + KX[3],
                (spice[4] ^ KX[(block_size << 2) + 20]) + KX[4],
                (spice[5] ^ KX[(block_size << 2) + 21]) + KX[5],
                (spice[6] ^ KX[(block_size << 2) + 22]) + KX[6],
                (spice[7] ^ KX[(block_size << 2) + 23]) + KX[7],
            };
            const uint64_t zspice[HPC_ROUND_COUNT] = { 0 };

            hpc_long_encrypt(tmp, zspice, KX, NULL, NULL, 512, UINT64_C(0xffffffffffffffff), 0);

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

static void hpc_tiny_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0];
    (void)ciphertext, (void)o_plaintext;

    switch (block_size)
    {
    case 1:
    case 2:
    case 3:
    case 4:
        {
            uint64_t tmp[2] = {
                KX[(block_size << 1) + 16] + KX[128] + (uint64_t)backup,
                KX[(block_size << 1) + 17] + KX[129],
            };

            hpc_medium_encrypt(tmp, spice, KX, NULL, NULL, 128, UINT64_C(0xffffffffffffffff), 0);

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
            tmp[0] += (uint64_t)backup;

            hpc_long_encrypt(tmp, spice, KX, NULL, NULL, tmp_bs, UINT64_C(0xffffffffffffffff), 0);

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
                (spice[0] ^ KX[(block_size << 2) + 16]) + KX[0] + (uint64_t)backup,
                (spice[1] ^ KX[(block_size << 2) + 17]) + KX[1],
                (spice[2] ^ KX[(block_size << 2) + 18]) + KX[2],
                (spice[3] ^ KX[(block_size << 2) + 19]) + KX[3],
                (spice[4] ^ KX[(block_size << 2) + 20]) + KX[4],
                (spice[5] ^ KX[(block_size << 2) + 21]) + KX[5],
                (spice[6] ^ KX[(block_size << 2) + 22]) + KX[6],
                (spice[7] ^ KX[(block_size << 2) + 23]) + KX[7],
            };
            const uint64_t zspice[HPC_ROUND_COUNT] = { 0 };

            hpc_long_encrypt(tmp, zspice, KX, NULL, NULL, 512, UINT64_C(0xffffffffffffffff), 0);

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

static void hpc_short_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup)
{
    const size_t LBH = (block_size + 1) >> 1;
    const size_t LBQ = (LBH + 1) >> 1;
    const size_t LBT = ((block_size + LBQ) >> 2) + 2;
    const size_t GAP = 64 - block_size;

    uint64_t s0 = state[0];
    (void)backup, (void)plaintext, (void)o_ciphertext;

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

static void hpc_short_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup)
{
    const size_t LBH = (block_size + 1) >> 1;
    const size_t LBQ = (LBH + 1) >> 1;
    const size_t LBT = ((block_size + LBQ) >> 2) + 2;
    const size_t GAP = 64 - block_size;

    uint64_t s0 = state[0];
    (void)backup, (void)ciphertext, (void)o_plaintext;

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

static void hpc_medium_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0], s1 = state[1];
    (void)backup, (void)plaintext, (void)o_ciphertext;

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

static void hpc_medium_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0], s1 = state[1];
    (void)backup, (void)ciphertext, (void)o_plaintext;

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

static void hpc_long_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0], s1 = state[1], s2 = state[2], s3 = state[3];
    uint64_t s4 = state[4], s5 = state[5], s6 = state[6], s7 = state[7];
    (void)backup, (void)plaintext, (void)o_ciphertext;

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

static void hpc_long_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup)
{
    uint64_t s0 = state[0], s1 = state[1], s2 = state[2], s3 = state[3];
    uint64_t s4 = state[4], s5 = state[5], s6 = state[6], s7 = state[7];
    (void)backup, (void)ciphertext, (void)o_plaintext;

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

static void hpc_extended_stir(uint64_t *s, const uint64_t *spice, const uint64_t *KX, size_t ri, uint64_t mask)
{
    uint64_t t = s[0] & 255, k = KX[t], kk = KX[(t + (ri << 2) + 1) & 255], tt;

    s[3] += s[7];
    s[5] ^= s[7];
    s[1] += k;
    s[2] ^= k;
    s[4] += kk;
    s[6] ^= kk;
    s[4] ^= s[1];
    s[5] += s[2];
    s[0] ^= (s[5] >> 13);
    s[1] -= (s[6] >> 22);
    s[2] ^= (s[7] << 7);
    s[7] ^= (s[6] << 9);
    s[7] += s[0];
    s[4] -= s[0];

    t = (s[1] & 31);
    tt = (s[1] >> t);
    s[6] ^= tt;
    s[7] += tt;

    tt = (s[2] << t);
    s[3] += tt;
    s[5] ^= tt;
    tt = (s[4] >> t);
    s[2] -= tt;
    s[5] += tt;

    if (ri == 1)
    {
        s[0] += spice[0];
        s[1] ^= spice[1];
        s[2] -= spice[2];
        s[3] ^= spice[3];
        s[4] += spice[4];
        s[5] ^= spice[5];
        s[6] -= spice[6];
        s[7] ^= spice[7];
    }

    s[7] -= s[3];
    s[7] &= mask;
    s[1] ^= (s[7] >> 11);
    s[6] += s[3];
    s[0] ^= s[6];

    // The next sequence exchanges the even-numbered bit positions in s2 and s5.
    // S3 and s0 are modified as targets of opportunity.
    t = s[2] ^ s[5];
    s[3] -= t;
    t &= 0x5555555555555555;
    s[2] ^= t;
    s[5] ^= t;
    s[0] += t;

    t = (s[4] << 9);
    s[6] -= t;
    s[1] += t;
}

static void hpc_extended_stir_inverse(uint64_t *s, const uint64_t *spice, const uint64_t *KX, size_t ri, uint64_t mask)
{
    uint64_t t, tt, k, kk;

    t = s[4] << 9;
    s[1] -= t;
    s[6] += t;

    t = s[2] ^ s[5];
    s[3] += t;
    t &= 0x5555555555555555;
    s[2] ^= t;
    s[5] ^= t;
    s[0] -= t;

    s[0] ^= s[6];
    s[6] -= s[3];
    s[1] ^= (s[7] >> 11);
    s[7] += s[3];

    if (ri == 1)
    {
        s[0] -= spice[0];
        s[1] ^= spice[1];
        s[2] += spice[2];
        s[3] ^= spice[3];
        s[4] -= spice[4];
        s[5] ^= spice[5];
        s[6] += spice[6];
        s[7] ^= spice[7];
    }

    t = (s[1] & 31);
    tt = (s[4] >> t);
    s[5] -= tt;
    s[2] += tt;
    tt = (s[2] << t);
    s[5] ^= tt;
    s[3] -= tt;

    tt = (s[1] >> t);
    s[6] ^= tt;
    s[7] -= tt;

    s[4] += s[0];
    s[7] -= s[0];
    s[7] ^= (s[6] << 9);
    s[7] &= mask;
    s[2] ^= (s[7] << 7);
    s[1] += (s[6] >> 22);
    s[0] ^= (s[5] >> 13);
    s[5] -= s[2];
    s[4] ^= s[1];

    t = s[0] & 255;
    k = KX[t];
    kk = KX[(t + (ri << 2) + 1) & 255];

    s[6] ^= kk;
    s[4] -= kk;
    s[2] ^= k;
    s[1] -= k;
    s[5] ^= s[7];
    s[3] -= s[7];
}

static const uint32_t Swizpoly[] = {
    // These are commented out because they'll never be used, since 513/64
    // rounded up is 9, leading to QMSK that is no less than 15, and as such
    // 0x13 is the minimal number that can be used
    /* 0, 3, 7, 0xb, */
    0x13, 0x25, 0x43, 0x83, 0x11d, 0x211, 0x409,
    0x805, 0x1053, 0x201b, 0x402b, UINT32_C(0x8003), UINT32_C(0x1002d),
    UINT32_C(0x20009), UINT32_C(0x40027), UINT32_C(0x80027), UINT32_C(0x100009),
    UINT32_C(0x200005), UINT32_C(0x400003), UINT32_C(0x800021), UINT32_C(0x100001b),
    UINT32_C(0x2000009), UINT32_C(0x4000047), UINT32_C(0x8000027), UINT32_C(0x10000009),
    UINT32_C(0x20000005), UINT32_C(0x40000053), UINT32_C(0x80000009)
};

static void hpc_extended_encrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t block_size, uint64_t mask, size_t backup)
{
    const size_t LWD = (block_size + 63) / 64;
    uint32_t qmask = (uint32_t)(LWD - 1);
    (void)backup;

    // Calculate QMSK to be `pow(2, ceil(log2(LWD))) - 1` using the bit-twiddling hack from
    // https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
    qmask |= qmask >> 1;
    qmask |= qmask >> 2;
    qmask |= qmask >> 4;
    qmask |= qmask >> 8;
    qmask |= qmask >> 16;

    // Pre-mixing step. The spec claims we should run this 4 times with 0..3, but
    // the test vectors disagree, so we deviate and run only 3 times with 0..2
    for (size_t i = 0; i < 3; ++i)
    {
        hpc_extended_stir(state, spice, KX, i, UINT64_C(0xffffffffffffffff));
    }

    uint64_t s7_copy = state[7];

    // First pass
    for (size_t i = HPC_ROUND_COUNT, j = i * sizeof(uint64_t); i < LWD; ++i, j += sizeof(uint64_t))
    {
        const uint64_t lmask = (i == LWD - 1) ? mask : UINT64_C(0xffffffffffffffff);
        const size_t byte_limit = (~lmask != 0) ? ((block_size & 63) + CHAR_BIT - 1) / CHAR_BIT : sizeof(uint64_t);

        state[7] = plaintext[j];
        for (size_t bi = 1, sh = CHAR_BIT; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            state[7] |= ((uint64_t)plaintext[j + bi]) << sh;
        }

        hpc_extended_stir(state, spice, KX, 0, lmask);

        for (size_t bi = 0, sh = 0; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            o_ciphertext[j + bi] = (state[7] >> sh) & 0xff;
        }
    }

    // First intermission
    state[7] = s7_copy;
    hpc_extended_stir(state, spice, KX, 0, UINT64_C(0xffffffffffffffff));
    state[0] += (uint64_t)block_size;
    // The spec claims this should run 3 times with 0..2, but the test vectors
    // disagree, so we deviate and run only 2 times with 0..1
    for (size_t i = 0; i < 2; ++i)
    {
        hpc_extended_stir(state, spice, KX, i, UINT64_C(0xffffffffffffffff));
    }
    state[0] += (uint64_t)block_size;
    s7_copy = state[7];

    // Second pass
    for (uint32_t q = 1; q != 0; q = ((q * 5) + 1) & qmask)
    {
        if (q < HPC_ROUND_COUNT || q >= LWD)
        {
            continue;
        }

        const uint64_t lmask = (q == LWD - 1) ? mask : UINT64_C(0xffffffffffffffff);
        const size_t byte_limit = (~lmask != 0) ? ((block_size & 63) + CHAR_BIT - 1) / CHAR_BIT : sizeof(uint64_t);
        const size_t j = (size_t)q * sizeof(uint64_t);

        state[7] = o_ciphertext[j];
        for (size_t bi = 1, sh = CHAR_BIT; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            state[7] |= ((uint64_t)o_ciphertext[j + bi]) << sh;
        }

        hpc_extended_stir(state, spice, KX, 0, lmask);

        for (size_t bi = 0, sh = 0; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            o_ciphertext[j + bi] = (state[7] >> sh) & 0xff;
        }
    }

    // Second intermission
    state[7] = s7_copy;

    hpc_extended_stir(state, spice, KX, 1, UINT64_C(0xffffffffffffffff));
    state[0] += (uint64_t)block_size;
    // The spec claims this should run 3 times with 0..2, but the test vectors
    // disagree, so we deviate and run only 2 times with 0..1
    for (size_t i = 0; i < 2; ++i)
    {
        hpc_extended_stir(state, spice, KX, i, UINT64_C(0xffffffffffffffff));
    }
    state[0] += (uint64_t)block_size;
    s7_copy = state[7];

    // find the correct Swizpoly number
    uint32_t swz = 0;
    for (size_t i = 0; i < (sizeof(Swizpoly) / sizeof(Swizpoly[0])); ++i)
    {
        if (Swizpoly[i] > qmask)
        {
            swz = Swizpoly[i];
            break;
        }
    }

    // Third pass
    qmask = (qmask >> 1) + 1;
    for (uint32_t q = 2; q != 1; q = (q << 1) ^ ((q & qmask) ? swz : 0))
    {
        if (q < HPC_ROUND_COUNT || q >= LWD)
        {
            continue;
        }

        const uint64_t lmask = (q == LWD - 1) ? mask : UINT64_C(0xffffffffffffffff);
        const size_t byte_limit = (~lmask != 0) ? ((block_size & 63) + CHAR_BIT - 1) / CHAR_BIT : sizeof(uint64_t);
        const size_t j = (size_t)q * sizeof(uint64_t);

        state[7] = o_ciphertext[j];
        for (size_t bi = 1, sh = CHAR_BIT; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            state[7] |= ((uint64_t)o_ciphertext[j + bi]) << sh;
        }

        hpc_extended_stir(state, spice, KX, 0, lmask);

        for (size_t bi = 0, sh = 0; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            o_ciphertext[j + bi] = (state[7] >> sh) & 0xff;
        }
    }

    // Finale
    state[7] = s7_copy;
    hpc_extended_stir(state, spice, KX, 0, UINT64_C(0xffffffffffffffff));
    for (size_t i = 0; i < 3; ++i)
    {
        hpc_extended_stir(state, spice, KX, i, UINT64_C(0xffffffffffffffff));
    }
}

static void hpc_extended_decrypt(uint64_t *state, const uint64_t *spice, const uint64_t *KX, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t block_size, uint64_t mask, size_t backup)
{
    const size_t LWD = (block_size + 63) / 64;
    uint32_t qmask = (uint32_t)(LWD - 1);
    (void)backup;

    // Calculate QMSK to be `pow(2, ceil(log2(LWD))) - 1` using the bit-twiddling hack from
    // https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
    qmask |= qmask >> 1;
    qmask |= qmask >> 2;
    qmask |= qmask >> 4;
    qmask |= qmask >> 8;
    qmask |= qmask >> 16;

    // Finale inverse
    for (size_t i = 3; i-- > 0; )
    {
        hpc_extended_stir_inverse(state, spice, KX, i, UINT64_C(0xffffffffffffffff));
    }
    hpc_extended_stir_inverse(state, spice, KX, 0, UINT64_C(0xffffffffffffffff));

    // find the correct Swizpoly number
    uint32_t swz = 0;
    for (size_t i = 0; i < (sizeof(Swizpoly) / sizeof(Swizpoly[0])); ++i)
    {
        if (Swizpoly[i] > qmask)
        {
            swz = Swizpoly[i];
            break;
        }
    }

    uint64_t s7_copy = state[7];

    // Third pass inverse
    swz >>= 1;
    for (uint32_t q = swz; q != 1; q = (q >> 1) ^ ((q & 1) ? swz : 0))
    {
        if (q < HPC_ROUND_COUNT || q >= LWD)
        {
            continue;
        }

        const uint64_t lmask = (q == LWD - 1) ? mask : UINT64_C(0xffffffffffffffff);
        const size_t byte_limit = (~lmask != 0) ? ((block_size & 63) + CHAR_BIT - 1) / CHAR_BIT : sizeof(uint64_t);
        const size_t j = (size_t)q * sizeof(uint64_t);

        state[7] = ciphertext[j];
        for (size_t bi = 1, sh = CHAR_BIT; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            state[7] |= ((uint64_t)ciphertext[j + bi]) << sh;
        }

        hpc_extended_stir_inverse(state, spice, KX, 0, lmask);

        for (size_t bi = 0, sh = 0; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            o_plaintext[j + bi] = (state[7] >> sh) & 0xff;
        }
    }

    // Second intermission inverse
    state[7] = s7_copy;
    state[0] -= (uint64_t)block_size;
    for (size_t i = 2; i-- > 0; )
    {
        hpc_extended_stir_inverse(state, spice, KX, i, UINT64_C(0xffffffffffffffff));
    }
    state[0] -= (uint64_t)block_size;
    hpc_extended_stir_inverse(state, spice, KX, 1, UINT64_C(0xffffffffffffffff));
    s7_copy = state[7];

    // Second pass inverse
    for (uint32_t q = UINT32_C(0x33333333) & qmask; q != 0; q = ((q - 1) * UINT32_C(0xcccccccd)) & qmask)
    {
        if (q < HPC_ROUND_COUNT || q >= LWD)
        {
            continue;
        }

        const uint64_t lmask = (q == LWD - 1) ? mask : UINT64_C(0xffffffffffffffff);
        const size_t byte_limit = (~lmask != 0) ? ((block_size & 63) + CHAR_BIT - 1) / CHAR_BIT : sizeof(uint64_t);
        const size_t j = (size_t)q * sizeof(uint64_t);

        state[7] = o_plaintext[j];
        for (size_t bi = 1, sh = CHAR_BIT; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            state[7] |= ((uint64_t)o_plaintext[j + bi]) << sh;
        }

        hpc_extended_stir_inverse(state, spice, KX, 0, lmask);

        for (size_t bi = 0, sh = 0; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            o_plaintext[j + bi] = (state[7] >> sh) & 0xff;
        }
    }

    // First intermission inverse
    state[7] = s7_copy;
    state[0] -= (uint64_t)block_size;
    for (size_t i = 2; i-- > 0; )
    {
        hpc_extended_stir_inverse(state, spice, KX, i, UINT64_C(0xffffffffffffffff));
    }
    state[0] -= (uint64_t)block_size;
    hpc_extended_stir_inverse(state, spice, KX, 0, UINT64_C(0xffffffffffffffff));
    s7_copy = state[7];

    // First pass inverse
    for (size_t i = LWD - 1, j = i * sizeof(uint64_t); i >= HPC_ROUND_COUNT; --i, j -= sizeof(uint64_t))
    {
        const uint64_t lmask = (i == LWD - 1) ? mask : UINT64_C(0xffffffffffffffff);
        const size_t byte_limit = (~lmask != 0) ? ((block_size & 63) + CHAR_BIT - 1) / CHAR_BIT : sizeof(uint64_t);

        state[7] = o_plaintext[j];
        for (size_t bi = 1, sh = CHAR_BIT; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            state[7] |= ((uint64_t)o_plaintext[j + bi]) << sh;
        }

        hpc_extended_stir_inverse(state, spice, KX, 0, lmask);

        for (size_t bi = 0, sh = 0; bi < byte_limit; ++bi, sh += CHAR_BIT)
        {
            o_plaintext[j + bi] = (state[7] >> sh) & 0xff;
        }
    }

    state[7] = s7_copy;

    // Pre-mixing inverse
    for (size_t i = 3; i-- > 0; )
    {
        hpc_extended_stir_inverse(state, spice, KX, i, UINT64_C(0xffffffffffffffff));
    }
}

static size_t hpc_get_cipher_id(size_t data_bit_size)
{
    if (data_bit_size <= 35)
    {
        return HPC_CIPHER_ID_TINY;
    }
    if (data_bit_size <= 64)
    {
        return HPC_CIPHER_ID_SHORT;
    }
    if (data_bit_size <= 128)
    {
        return HPC_CIPHER_ID_MEDIUM;
    }
    if (data_bit_size <= 512)
    {
        return HPC_CIPHER_ID_LONG;
    }
    // While in theory HPC-extended supports arbitrary lengths,
    // in reality it's unlikely that anyone would use a larger
    // block size, and this restriction simplifies the code
    if (data_bit_size <= UINT64_C(137438954048))
    {
        return HPC_CIPHER_ID_EXTENDED;
    }

    return -1;
}

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */
