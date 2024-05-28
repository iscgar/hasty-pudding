#include "hpc/hpc.h"

#ifdef _WIN32
#   define DLL_EXPORT __declspec(dllexport)
#else /* if !_WIN32 */
#   define DLL_EXPORT
#endif /* !_WIN32 */

extern "C"
{

DLL_EXPORT
int hpc_encrypt_data(const uint64_t *key, size_t key_bit_size, const size_t *backup, size_t backup_size, const uint64_t *tweak, size_t tweak_bit_size, uint64_t *io_data, size_t data_bit_size)
{
    HpcState state;

    if (!hpc_init_with_backup(&state, reinterpret_cast<const uint8_t *>(key), key_bit_size, backup, backup_size))
    {
        return false;
    }

    return hpc_encrypt(&state, reinterpret_cast<const uint8_t *>(io_data), reinterpret_cast<uint8_t *>(io_data), data_bit_size, reinterpret_cast<const uint8_t *>(tweak), tweak_bit_size);
}

DLL_EXPORT
int hpc_decrypt_data(const uint64_t *key, size_t key_bit_size, const size_t *backup, size_t backup_size, const uint64_t *tweak, size_t tweak_bit_size, uint64_t *io_data, size_t data_bit_size)
{
    HpcState state;

    if (!hpc_init_with_backup(&state, reinterpret_cast<const uint8_t *>(key), key_bit_size, backup, backup_size))
    {
        return false;
    }

    return hpc_decrypt(&state, reinterpret_cast<const uint8_t *>(io_data), reinterpret_cast<uint8_t *>(io_data), data_bit_size, reinterpret_cast<const uint8_t *>(tweak), tweak_bit_size);
}

} // extern "C"