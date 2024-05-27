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

extern int hpc_init(struct HpcState *state, const uint8_t *key, size_t key_bit_size);

extern int hpc_init_with_backup(struct HpcState *state, const uint8_t *key, size_t key_bit_size, const size_t *backup, size_t backup_size);

extern int hpc_encrypt(struct HpcState *state, const uint8_t *plaintext, uint8_t *o_ciphertext, size_t data_bit_size, const uint8_t *tweak, size_t tweak_bit_size);

extern int hpc_decrypt(struct HpcState *state, const uint8_t *ciphertext, uint8_t *o_plaintext, size_t data_bit_size, const uint8_t *tweak, size_t tweak_bit_size);

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* !INCLUDE_GUARD_8961612A_CF77_4916_B8FD_E00A836A2F71 */
