#ifndef DSPMATH_INCLUDED
#define DSPMATH_INCLUDED

/*
 *      #####################################################################
 *
 *        Copyright (C) 2006  Johan H. Bodin SM6LKM, Wolfgang Buescher DL4YHF
 *
 *        This software is provided 'as is', without warranty of any kind,
 *        express or implied. In no event shall the authors be held liable
 *        for any damages arising from the use of this software.
 *
 *        Permission to use, copy, modify, and distribute this software and
 *        its documentation for non-commercial purposes is hereby granted,
 *        provided that the above copyright notice and this disclaimer appear
 *        in all copies and supporting documentation.
 *
 *        The software must NOT be sold or used as part of a any commercial
 *        or "non-free" product.
 *
 *      #####################################################################
 */

/*------------------------------------------------------------------------------
 *
 *      dspmath.h
 *
 *      Johan Bodin SM6LKM
 *
 *      This file is based on original work by Wolfgang Buescher (DL4YHF)
 *      (SoundMaths.h).
 *
 *      Literature references:
 *        [SGDSP] = Steven W. Smith, "The Scientists and Engineer's Guide
 *                  to Digital Signal Processing", Chapter 12, "The Fast
 *                  Fourier Transform", www.DSPguide.com .
 *
 *      Revision history:
 *        2006-11-23:
 *          Based on DL4YHF's SoundMaths.h, some stuff stripped to make it
 *          compatible with the plain C compiler of Dev-C++.
 *          Comments & indentation adjusted to suit SM6LKM's taste. /JHB
 *        2006-11-27:
 *          Added HaMMing and Blackman windowing functions. /JHB
 */

/*------------------------------------------------------------------------------
 *
 *      Un-select unused functions to save space
 */

//#define USE_DSPMATH_RESAMPLE_FLOAT_ARRAY
//#define USE_DSPMATH_CALCULATE_ANGLE
//#define USE_DSPMATH_CALCULATE_ANGLE_FAST
//#define USE_DSPMATH_RUN_COMPLEX_FIR
#define USE_DSPMATH_MULTIPLY_HANNING_WINDOW
//#define USE_DSPMATH_MULTIPLY_HAMMING_WINDOW
#define USE_DSPMATH_MULTIPLY_BLACKMAN_WINDOW
#define USE_DSPMATH_CALC_COMPLEX_FFT
//#define USE_DSPMATH_CALC_COMPLEX_INVERSE_FFT
#define USE_DSPMATH_CALC_REAL_FFT
//#define USE_DSPMATH_CALC_REAL_INVERSE_FFT


// DO NOT CHANGE THE #defines BELOW!

#ifdef USE_DSPMATH_CALC_COMPLEX_INVERSE_FFT
  #define USE_DSPMATH_CALC_COMPLEX_FFT
#endif

#ifdef USE_DSPMATH_CALC_REAL_INVERSE_FFT
  #define USE_DSPMATH_CALC_REAL_FFT
#endif

#ifdef USE_DSPMATH_CALC_REAL_FFT
  #define USE_DSPMATH_CALC_COMPLEX_FFT
#endif


/*------------------------------------------------------------------------------
 *
 *      A popular constant...
 */

#ifndef  C_PI
 #define C_PI  (3.1415926535897932384626)
#endif


/*------------------------------------------------------------------------------
 *
 *      Select type to use for T_FAST_FLOAT (float or double)
 *
 *      This selection affects the following:
 *        dspmath_CalculateAngleFast()
 */

#define T_FAST_FLOAT float
//#define T_FAST_FLOAT double


/*------------------------------------------------------------------------------
 *
 *      Select type to use for T_Float (float or double)
 *
 *      This selection affects the following:
 *        dspmath_ResampleFloatArray()
 *        dspmath_RunComplexFIR()
 *        T_Complex
 */

#define T_Float float
// #define T_Float double


/*------------------------------------------------------------------------------
 *
 *      Complex point type
 *
 *      Used by:
 *        dspmath_RunComplexFIR()
 */

typedef struct
{
   T_Float re;  // real part or "I" ("in-phase")
   T_Float im;  // imaginary part or "Q" ("quadrature-phase")
} T_Complex;


/*------------------------------------------------------------------------------
 *
 *      Stretches or shrinks an array.
 *
 *      Originally used for the FFT-based filter, to adapt the
 *      frequency response curve when changing the FFT size.
 *      Neither iSourceLength nor iDestLength may be zero or negative!
 */

#ifdef USE_DSPMATH_RESAMPLE_FLOAT_ARRAY
void dspmath_ResampleFloatArray (T_Float *pfltArray,
                                 int iSourceLength,
                                 int iDestLength);
#endif


/*------------------------------------------------------------------------------
 *
 *      Precise four-quadrant conversion of a complex pair ("I/Q")
 *      into an phase value (in radians, but explained in degrees here).
 *      A positive real value gives an angle of zero, etc.
 *      Returned value range is -180° .. +180° =  -pi .. pi .
 *      If both real and imaginary part are zero, the returned value
 *      is zero.
 */

#ifdef USE_DSPMATH_CALCULATE_ANGLE
double dspmath_CalculateAngle (double re, double im);
#endif


/*------------------------------------------------------------------------------
 *
 *      Fast atan2 calculation with self normalization.
 *      Returned value range is  -pi..pi =  -180° .. +180° .
 *      Detailed explanation and discussion of accuracy in .c file!
 */

#ifdef USE_DSPMATH_CALCULATE_ANGLE_FAST
T_FAST_FLOAT dspmath_CalculateAngleFast (T_FAST_FLOAT x, T_FAST_FLOAT y);
#endif


/*------------------------------------------------------------------------------
 *
 *      Complex FIR-filter (usually a low pass)
 */

#ifdef USE_DSPMATH_RUN_COMPLEX_FIR
void dspmath_RunComplexFIR (int       iNrCoeffs,    // Length of filter queue + count of coeffs
                            T_Float   *pfltCoeffs,  // pointer to filter coefficients   [iNrCoeffs]
                            T_Complex *pcpxQueue,   // pointer to filter queue (memory) [iNrCoeffs]
                            int       *piQueueIdx,  // index for circular filter queue, 0..iNrCoeffs-1
                            T_Complex *pcplxValue); // reference to in- and output value
#endif


/*------------------------------------------------------------------------------
 *
 *      Windowing functions
 *
 *      Hanning:
 *        w[i] = 0.5 - 0.5*cos (2*PI*i/M)
 *        where i = 0..M
 *
 *      Hamming:
 *        w[i] = 0.54 - 0.46*cos (2*PI*i/M)
 *        where i = 0..M
 *
 *      Blackman:
 *        w[i] = 0.42 - 0.5*cos (2*PI*i/M) + 0.08*cos (4*PI*i/M)
 *        where i = 0..M
 *
 *  *#* Fixme: Make faster versions with look-up tables! /JHB
 */

#ifdef USE_DSPMATH_MULTIPLY_HANNING_WINDOW
void dspmath_MultiplyHanningWindow (float *pfltArray, int iLength);
#endif

#ifdef USE_DSPMATH_MULTIPLY_HAMMING_WINDOW
void dspmath_MultiplyHammingWindow (float *pfltArray, int iLength);
#endif

#ifdef USE_DSPMATH_MULTIPLY_BLACKMAN_WINDOW
void dspmath_MultiplyBlackmanWindow (float *pfltArray, int iLength);
#endif


/*------------------------------------------------------------------------------
 *
 *      Complex Fast Fourier Transform
 *
 *      Upon entry, N contains the number of points in the DFT, and
 *      REX[] and IMX[] contain the real and imaginary parts of the input.
 *      Upon return, REX[] and IMX[] contain the complex DFT output.
 *
 *      All signals run from 0 to N-1.
 *
 *      More info in .c file!
 */

#ifdef USE_DSPMATH_CALC_COMPLEX_FFT
void dspmath_CalcComplexFft (int iNrOfPoints, // N =  number of points in the DFT *AND* in the time domain
                             float *pfltRe,   // REX[] = input: re(time domain), result: re(frequency domain)
                             float *pfltIm);  // IMX[] = input: im(time domain), result: im(frequency domain)
#endif


/*------------------------------------------------------------------------------
 *
 *      Inverse Complex Fast Fourier Transform
 *
 *      Inspired by [SGDSP] TABLE 12-5.
 *
 *      Upon entry, N contains the number of points in the IDFT,
 *      REX[] and IMX[] contain the real & imaginary parts of the complex
 *      frequency domain. Upon return, REX[] and IMX[] contain the complex
 *      time domain signal.
 *
 *      All signals run from 0 to N-1.
 */

#ifdef USE_DSPMATH_CALC_COMPLEX_INVERSE_FFT
void dspmath_CalcComplexInverseFft (int iNrOfPoints,  // N  number of points in the IDFT IN THE TIME DOMAIN!
                                    float *pfltRe,    // REX[] = input: re(frequency domain), result: re(time domain)
                                    float *pfltIm);   // IMX[] = input: im(frequency domain), result: im(time domain)
#endif


/*------------------------------------------------------------------------------
 *
 *      Fast Fourier Transform for real input signals
 *
 *      Inspired by [SGDSP] TABLE 12-7.
 *
 *      Upon entry, N contains the number of points in the DFT,
 *      REX[] contains the real input signal while values in IMX[] are ignored.
 *      The input signal run from 0 to N-1.
 *
 *      Upon return, REX[] and IMX[] contain the DFT output.
 *      The output signals run from  0...N/2, for example a
 *      "1024 point REAL FFT" produces 513(!) POINTS in REX[] and
 *      513(!) POINTS in IMX[]. See [SGDSP] for details!
 */

#ifdef USE_DSPMATH_CALC_REAL_FFT
void dspmath_CalcRealFft (int iNrOfPoints,  // N = number of points in the DFT
                          float *pfltRe,    // REX[] = input: re(time domain), result: re(frequency domain)
                          float *pfltIm);   // IMX[] = input: ignored, result: re(frequency domain)
#endif


/*------------------------------------------------------------------------------
 *
 *      Inverse Fast Fourier Transform for real signals
 *
 *      Inspired by [SGDSP] TABLE 12-6.
 *
 *      Upon entry, N contains the number of points in the IDFT,
 *      REX[] and IMX[] contain the real & imaginary parts of the frequency
 *      domain running from index 0 to N/2. The remaining samples in REX[] and
 *      IMX[] are ignored. Upon return, REX[] contains the real time domain,
 *      IMX[] contains zeroes.
 */

#ifdef USE_DSPMATH_CALC_REAL_INVERSE_FFT
void dspmath_CalcRealInverseFft (int iNrOfPoints, // N = number of points in the IDFT IN THE TIME DOMAIN !
                                 float *pfltRe,   // REX[] = input: re(frequency domain), result: re(time domain)
                                 float *pfltIm);  // IMX[] = input: im(frequency domain), result: zeroes
#endif

#endif  // DSPMATH_INCLUDED
