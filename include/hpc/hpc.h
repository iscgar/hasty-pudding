/**
 *  Copyright (c) 2024 Isaac Garzon
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE
 */

#ifndef INCLUDE_GUARD_8961612A_CF77_4916_B8FD_E00A836A2F71
#define INCLUDE_GUARD_8961612A_CF77_4916_B8FD_E00A836A2F71

#include <stddef.h>
#include <stdint.h>

#define HPC_KX_SIZE      (256)
#define HPC_CIPHER_COUNT (5)
#define HPC_BACKUP_SIZE  (HPC_CIPHER_COUNT + 1)

#define HPC_TWEAK_BIT_SIZE  (512)

#define HPC_CIPHER_ID_TINY      (1)
#define HPC_CIPHER_ID_SHORT     (2)
#define HPC_CIPHER_ID_MEDIUM    (3)
#define HPC_CIPHER_ID_LONG      (4)
#define HPC_CIPHER_ID_EXTENDED  (5)

/**
 * Defines the Hasty Pudding cipher's internal state.
 */
struct HpcState
{
#if defined(HPC_USE_SINGLE_KX)
#   if (HPC_USE_SINGLE_KX < HPC_CIPHER_ID_TINY || HPC_USE_SINGLE_KX > HPC_CIPHER_ID_EXTENDED)
#       error "Invalid HPC cipher ID provided to HPC_USE_SINGLE_KX"
#   else /* if (HPC_USE_SINGLE_KX >= HPC_CIPHER_ID_TINY && HPC_USE_SINGLE_KX <= HPC_CIPHER_ID_EXTENDED) */
    uint64_t KX[HPC_KX_SIZE];
#   endif /* (HPC_USE_SINGLE_KX >= HPC_CIPHER_ID_TINY && HPC_USE_SINGLE_KX <= HPC_CIPHER_ID_EXTENDED) */
#else /* if !defined(HPC_USE_SINGLE_KX) */
    uint64_t KX[HPC_CIPHER_COUNT][HPC_KX_SIZE];
#endif /* !defined(HPC_USE_SINGLE_KX) */
    size_t backup[HPC_BACKUP_SIZE];
};

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * Initialises an HPC instance with a specific key.
 *
 * @param state         A pointer to HPC the state to initialise
 * @param key           A pointer to the key material to initialise the state with
 * @param key_bit_size  The size of the key material in bits
 *
 * @return              True if initialised successfully.
 *                      False otherwise.
 */
extern int hpc_init(struct HpcState *state, const uint8_t *key, size_t key_bit_size);

/**
 * Initialises an HPC instance with a specific key and backup rounds.
 *
 * @param state         A pointer to HPC the state to initialise
 * @param key           A pointer to the key material to initialise the state with
 * @param key_bit_size  The size of the key material in bits
 * @param backup        A pointer to an array describing how many rounds to add to
 *                      the key schedule and the sub-ciphers. The first member sets
 *                      the amount of additional rounds for key expansion, and the
 *                      rest of the member specify how many iteration to do for each
 *                      sub-cipher using the respective `HPC_CIPHER_ID_*` indices
 * @param backup_size   The number of members in `backup` (must be `HPC_BACKUP_SIZE`)
 *
 * @return              True if initialised successfully.
 *                      False otherwise.
 */
extern int hpc_init_with_backup(struct HpcState *state, const uint8_t *key, size_t key_bit_size, const size_t *backup, size_t backup_size);

/**
 * Encrypts a data block.
 *
 * @param state             A pointer to an initialised Hasty Pudding state
 * @param plaintext         A pointer to the data to encrypt
 * @param o_ciphertext      A pointer to an output array to put the encrypted data in
 *                          (can be the same as `plaintext`)
 * @param data_bit_size     The size in bits of the data to encrypt
 * @param tweak             An optional tweak to change the encryption behaviour without
 *                          changing the key (AKA "spice")
 * @param tweak_bit_size    The size of the tweak in bits (must be less than or
 *                          equal to `HPC_TWEAK_BIT_SIZE`)
 *
 * @return                  True if the data was encrypted successfully.
 *                          False otherwise.
 */
extern int hpc_encrypt(struct HpcState *state, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t data_bit_size, const uint8_t *tweak, size_t tweak_bit_size);

/**
 * Decrypts an encrypted data block.
 *
 * @param state             A pointer to an initialised Hasty Pudding state
 * @param ciphertext        A pointer to the data to decrypt
 * @param o_plaintext       A pointer to an output array to put the decrypted data in
 *                          (can be the same as `ciphertext`)
 * @param data_bit_size     The size in bits of the data to decrypt
 * @param tweak             An optional tweak to change the decryption behaviour without
 *                          changing the key (AKA "spice", must be the same tweak that was
 *                          provided to the encryption function)
 * @param tweak_bit_size    The size of the tweak in bits (must be less or equal to `HPC_TWEAK_BIT_SIZE`)
 *
 * @return                  True if the data was decrypted successfully.
 *                          False otherwise.
 */
extern int hpc_decrypt(struct HpcState *state, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t data_bit_size, const uint8_t *tweak, size_t tweak_bit_size);

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* !INCLUDE_GUARD_8961612A_CF77_4916_B8FD_E00A836A2F71 */
