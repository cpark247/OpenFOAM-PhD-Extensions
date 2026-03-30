#ifndef MEM_INLINES_H
#define MEM_INLINES_H
#ifndef INTTYPES_H
#define INTTYPES_H
#include <inttypes.h>
#endif
/** \returns decremented void pointer by count bytes
 * \param[in] a pointer to subract from
 * \param[in] count number of bytes to subract
 */
static inline void* void_ptr_sub(void *a,int64_t count)
{
    return (void *)((char*)a-count);
}

/** \returns increased void pointer by count bytes
 * \param[in] a pointer to add to
 * \param[in] count number of bytes to add
 */
static inline void* void_ptr_add(void *a,int64_t count)
{
    return (void *)((char*)a+count);
}

/** \returns pointer difference in bytes a has to be larger than b
 * \param[in] a pointer
 * \param[in] b pointer
 */
static inline int64_t void_ptr_diff(void *a, void *b)
{
#ifdef __INTEL_COMPILER
    return (char*)a - (char*)b;
#elif defined __cplusplus
    return (char *)a - (char *)b;
#elif defined __GNUC__
    return (void *)a - (void *)b;
#endif
}
#endif
