/* adler32.c -- compute the Adler-32 checksum of a data stream
 * Copyright (C) 1995-2011 Mark Adler
 * For conditions of distribution and use, see copyright notice in zlib.h
 */

/* @(#) $Id$ */

#include "zutil.h"

#define local static

local uLong adler32_combine_ OF((uLong adler1, uLong adler2, z_off64_t len2));

#define BASE 65521      /* largest prime smaller than 65536 */
#define NMAX 5552
/* NMAX is the largest n such that 255n(n+1)/2 + (n+1)(BASE-1) <= 2^32-1 */

#define DO1(buf,i)  {adler += (buf)[i]; sum2 += adler;}
#define DO2(buf,i)  DO1(buf,i); DO1(buf,i+1);
#define DO4(buf,i)  DO2(buf,i); DO2(buf,i+2);
#define DO8(buf,i)  DO4(buf,i); DO4(buf,i+4);
#define DO16(buf)   DO8(buf,0); DO8(buf,8);

/* use NO_DIVIDE if your processor does not do division in hardware --
   try it both ways to see which is faster */
#ifdef NO_DIVIDE
/* note that this assumes BASE is 65521, where 65536 % 65521 == 15
   (thank you to John Reiser for pointing this out) */
#  define CHOP(a) \
    do { \
        unsigned long tmp = a >> 16; \
        a &= 0xffffUL; \
        a += (tmp << 4) - tmp; \
    } while (0)
#  define MOD28(a) \
    do { \
        CHOP(a); \
        if (a >= BASE) a -= BASE; \
    } while (0)
#  define MOD(a) \
    do { \
        CHOP(a); \
        MOD28(a); \
    } while (0)
#  define MOD63(a) \
    do { /* this assumes a is not negative */ \
        z_off64_t tmp = a >> 32; \
        a &= 0xffffffffL; \
        a += (tmp << 8) - (tmp << 5) + tmp; \
        tmp = a >> 16; \
        a &= 0xffffL; \
        a += (tmp << 4) - tmp; \
        tmp = a >> 16; \
        a &= 0xffffL; \
        a += (tmp << 4) - tmp; \
        if (a >= BASE) a -= BASE; \
    } while (0)
#else
#  define MOD(a) a %= BASE
#  define MOD28(a) a %= BASE
#  define MOD63(a) a %= BASE
#endif

#if !defined(__aarch64__) && defined(__arm__)
#include <cutils/log.h>
#define LOG_TAG "libz"
void __attribute__((noinline)) adler32_neon(const char *buf, unsigned int *s1, unsigned int *s2, int k)
{
    /*
     * Algorithm:
     * buf: |x0|x1|x2|x3|x4|x5|x6|....|x15|
     * init s1 s2: s10 s20
     * Steps    s1 result           s2 result
     * ----------------------------------------
     * 1        s10+x0              s20+s10+x0
     * 2        s10+x0+x1           s20+s10+x0 + s10+x0+x1 = s20+2*s10+2*x0+x1
     * 3        s10+x0+x1+x2        s20+2*s10+2*x0+x1 + s10+x0+x1+x2 = s20+3*s10+3*x0+2*x1+x2
     * ...
     * 16       s10+x0+x1+...x15    s20+16*s10+16*x0+15*x1+14*x2+...+2*x14+x15
     * ...
     * 32       s10+sum(x0...x31)   s20+31*s10+32*x0+31*x1+...+x31 = s20
     */
    asm volatile (
        "ldr         r12, =0x0d0e0f10           \n"
        "ldr         lr,  =0x090a0b0c           \n"
        "vmov        d30, r12, lr               \n"
        "ldr         r12, =0x05060708           \n"
        "ldr         lr,  =0x01020304           \n"
        "vmov        d31, r12, lr               \n"
        "vmov.i8     d28, #0                    \n"
        "vmov.i8     q0,  #0                    \n"
        "ldr         lr,  [%[s1]]               \n"        // lr  = s1
        "ldr         r12, [%[s2]]               \n"        // r12 = s2
        "bic         %[k], %[k], #0xf           \n"        // clear
        "mla         r12, lr, %[k], r12         \n"        // r12 = s2 + s1*(k&~0xf)
    ".Lloop:                                    \n"
        "vld1.8      {d6, d7}, [%[buf]]!        \n"
        "vshl.i32    d19, d0, #4                \n"        // sum(f..0) * 16
        "sub         %[k], %[k], #16            \n"        // k -= 16
        "vmull.u8    q2, d6, d30                \n"
        "vmull.u8    q8, d7, d31                \n"
        "vadd.i32    d1, d1, d19                \n"        // s2 += 16*sum(f...0)
        "vpaddl.u8   q1, q3                     \n"        // q1 = |x15+x14|x13+x12|...|x3+x2|x1+x0|
        "cmp         %[k], #16                  \n"        //
        "vpadd.i16   d4, d4, d5                 \n"        // d4 = |f+2*e|...|13*3+14*2+15*1+16*0|
        "vpadd.i16   d16, d16, d17              \n"        //
        "vpadd.i16   d2, d2, d3                 \n"        // d2 = |f+e+d+c|b+a+9+8|7+6+5+4|3+2+1+0|
        "vpadd.i16   d4, d4, d16                \n"        // d4 = |f+2*e|...|13*3+14*2+15*1+16*0|
        "vpadd.i16   d2, d2, d28                \n"        // d2 = |0|0|f+e+d+..+8|7+6+..+0|.16
        "vpadd.i16   d4, d4, d28                \n"        //
        "vpadal.u16  d0, d2                     \n"        // d0 = |0|f+e+d+..+1+0|.32
        "vpadal.u16  d1, d4                     \n"        // d1 = |0|f+2*e+3*d+..+15*1+16*0|.32
        "bge         .Lloop                     \n"

        "vmov.32     %[buf],  d0[0]             \n"
        "vmov.32     %[k],  d1[0]               \n"
        "add         lr,  %[buf], lr            \n"
        "add         r12, %[k], r12             \n"
        "str         lr,  [%[s1]]               \n"
        "str         r12, [%[s2]]               \n"
        :
        :[buf] "r" (buf), [s1] "r" (s1), [s2] "r" (s2), [k] "r" (k)
        :"memory", "cc", "r12", "lr"
    );
}
#endif  /* __arm__ */

/* ========================================================================= */
uLong ZEXPORT adler32(adler, buf, len)
    uLong adler;
    const Bytef *buf;
    uInt len;
{
    unsigned long sum2;
    unsigned n;

    /* split Adler-32 into component sums */
    sum2 = (adler >> 16) & 0xffff;
    adler &= 0xffff;

    /* in case user likes doing a byte at a time, keep it fast */
    if (len == 1) {
        adler += buf[0];
        if (adler >= BASE)
            adler -= BASE;
        sum2 += adler;
        if (sum2 >= BASE)
            sum2 -= BASE;
        return adler | (sum2 << 16);
    }

    /* initial Adler-32 value (deferred check for len == 1 speed) */
    if (buf == Z_NULL)
        return 1L;

    /* in case short lengths are provided, keep it somewhat fast */
    if (len < 16) {
        while (len--) {
            adler += *buf++;
            sum2 += adler;
        }
        if (adler >= BASE)
            adler -= BASE;
        MOD28(sum2);            /* only added so many BASE's */
        return adler | (sum2 << 16);
    }

    /* do length NMAX blocks -- requires just one modulo operation */
    while (len >= NMAX) {
        len -= NMAX;
    #if defined(__arm__) && !defined(__aarch64__)
        adler32_neon(buf, &adler, &sum2, NMAX);
        buf += NMAX;
    #else
        n = NMAX / 16;          /* NMAX is divisible by 16 */
        do {
            DO16(buf);          /* 16 sums unrolled */
            buf += 16;
        } while (--n);
    #endif    /* __arm__ */
        MOD(adler);
        MOD(sum2);
    }

    /* do remaining bytes (less than NMAX, still just one modulo) */
    if (len) {                  /* avoid modulos if none remaining */
    #if defined(__arm__) && !defined(__aarch64__)
        if (len >= 16) {
            adler32_neon(buf, &adler, &sum2, len);
            buf += (len & ~0x0f);
            len  = (len % 16);
        }
    #else
        while (len >= 16) {
            len -= 16;
            DO16(buf);
            buf += 16;
        }
    #endif  /* __arm__ */
        while (len--) {
            adler += *buf++;
            sum2 += adler;
        }
        MOD(adler);
        MOD(sum2);
    }

    /* return recombined sums */
    return adler | (sum2 << 16);
}

/* ========================================================================= */
local uLong adler32_combine_(adler1, adler2, len2)
    uLong adler1;
    uLong adler2;
    z_off64_t len2;
{
    unsigned long sum1;
    unsigned long sum2;
    unsigned rem;

    /* for negative len, return invalid adler32 as a clue for debugging */
    if (len2 < 0)
        return 0xffffffffUL;

    /* the derivation of this formula is left as an exercise for the reader */
    MOD63(len2);                /* assumes len2 >= 0 */
    rem = (unsigned)len2;
    sum1 = adler1 & 0xffff;
    sum2 = rem * sum1;
    MOD(sum2);
    sum1 += (adler2 & 0xffff) + BASE - 1;
    sum2 += ((adler1 >> 16) & 0xffff) + ((adler2 >> 16) & 0xffff) + BASE - rem;
    if (sum1 >= BASE) sum1 -= BASE;
    if (sum1 >= BASE) sum1 -= BASE;
    if (sum2 >= (BASE << 1)) sum2 -= (BASE << 1);
    if (sum2 >= BASE) sum2 -= BASE;
    return sum1 | (sum2 << 16);
}

/* ========================================================================= */
uLong ZEXPORT adler32_combine(adler1, adler2, len2)
    uLong adler1;
    uLong adler2;
    z_off_t len2;
{
    return adler32_combine_(adler1, adler2, len2);
}

uLong ZEXPORT adler32_combine64(adler1, adler2, len2)
    uLong adler1;
    uLong adler2;
    z_off64_t len2;
{
    return adler32_combine_(adler1, adler2, len2);
}
