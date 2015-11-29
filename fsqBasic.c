/** 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
/*
 This is a partial port of http://www.qsl.net/zl1bpu/MFSK/FSQweb.htm
 ZL2AFP's FSQ Fast Simple QSO (chat) mode for HF and VHF
 
 
 This program uses the WiringPi Library http://wiringpi.com/
 WiringPi is released under the GNU LGPLv3 license
 * wiringPi:
 *	Arduino compatable (ish) Wiring library for the Raspberry Pi
 *	Copyright (c) 2012 Gordon Henderson
 
 * $Id$
 *
 * This program uses the PortAudio Portable Audio Library.
 * For more information see: http://www.portaudio.com
 * Copyright (c) 1999-2000 Ross Bencina and Phil Burk
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "portaudio.h"

/* #define SAMPLE_RATE  (17932) // Test failure to open with this value. */
//#define SAMPLE_RATE  (44100)
#define SAMPLE_RATE  (48000)

//#define FRAMES_PER_BUFFER (512)
#define FRAMES_PER_BUFFER (360)

#define FRAMES_PER_BUFFER12K (360/4)

#define NUM_SECONDS     (5)
#define NUM_CHANNELS    (2)
/* #define DITHER_FLAG     (paDitherOff) */
#define DITHER_FLAG     (0) /**/
/** Set to 1 if you want to capture the recording to a file. */
#define WRITE_TO_FILE   (1)

/* Select sample format. */
#if 1
#define PA_SAMPLE_TYPE  paFloat32
typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#elif 1
#define PA_SAMPLE_TYPE  paInt16
typedef short SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#elif 0
#define PA_SAMPLE_TYPE  paInt8
typedef char SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#else
#define PA_SAMPLE_TYPE  paUInt8
typedef unsigned char SAMPLE;
#define SAMPLE_SILENCE  (128)
#define PRINTF_S_FORMAT "%d"
#endif

typedef struct
{
    int          frameIndex;  /* Index into sample array. */
    int          maxFrameIndex;
    SAMPLE*      recordedSamples;
}
paTestData;


float  InputSamplesReal12000[FRAMES_PER_BUFFER12K*2]; 
float  InputSamplesReal48000[FRAMES_PER_BUFFER12K*2*4]; 


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
float channel_offset = 350;
#define FFT_SIZE_IN   (4096)               // time domain input points for the FFT
typedef unsigned char byte;

typedef int BOOL;
#define TRUE (1)
#define FALSE (0)

BOOL   agc                   = TRUE;        
BOOL   agcfast               = FALSE;
BOOL   agcOFF                = FALSE;
//BOOL   agcOFF                = TRUE;

#define NFFTPLOTPOINTS (FFT_SIZE_IN/4 + 1)  // number of bins in FFT output

#define SMOOTH_K    0.13333333

static float  InputSamplesI[FFT_SIZE_IN*2+1];
static float  InputSamplesQ[FFT_SIZE_IN*2+1];
int QTC_ID;

float   g_fltFftBufRe[FFT_SIZE_IN*2];  //buffers for FFT on raw input samples
float   g_fltFftBufIm[FFT_SIZE_IN*2];
float   g_fltPlotFftV[NFFTPLOTPOINTS*8];

float   Speed_Re[FFT_SIZE_IN*2];    //FFT buffers for un-windowed FFT's 
//for use in speed meter
float   Speed_Im[FFT_SIZE_IN*2];
float   SpeedPlotFft[NFFTPLOTPOINTS*8];

///////////////////////////////////////////////////
float   g_PlotPeak[65536];//NFFTPLOTPOINTS*8];

float   peak_symbol;
int     peak_index;
float   _sync[18];

int       hangcnt, k;  
float     ss, tmp, pk, agc_decay, agcthr, agcgain;
int       hangcnt, k;  


int       peak_counter, peak_index;


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


typedef int BOOL;
#define TRUE (1)
#define FALSE (0)
#define max(a,b) \
({ __typeof__ (a) _a = (a); \
__typeof__ (b) _b = (b); \
_a > _b ? _a : _b; })

static BOOL all_complete = FALSE;
int testCount = 0;

/* This routine will be called by the PortAudio engine when audio is needed.
 ** It may be called at interrupt level on some machines so don't do anything
 ** that could mess up the system like calling malloc() or free().
 */
static int recordCallback( const void *inputBuffer, void *outputBuffer,
						  unsigned long framesPerBuffer,
						  const PaStreamCallbackTimeInfo* timeInfo,
						  PaStreamCallbackFlags statusFlags,
						  void *userData )
{
    paTestData *data = (paTestData*)userData;
	
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
    float         leftInput;
    float         rightInput;
	
	
    static float  InputSamplesI[FFT_SIZE_IN*2+1];
    static float  InputSamplesQ[FFT_SIZE_IN*2+1];
    
    
    static float  inputI[FFT_SIZE_IN*2+1];
    static float  inputQ[FFT_SIZE_IN*2+1];
    
    double _Complex output[FFT_SIZE_IN*2+1];
	
    int    lm; //***TEST
	
    int           spec_ave, peak_smooth;
    
    static int    peak_val, prev_peak, peak_valSpeed, prev_peakSpeed ;
    static int    prev_symbol;
    
    
    static int    peak_counter, speed_counter;
    int           peak = 0;
    int           peakSpeed = 0;
    
    double        val, max = 0.0;
    double        valSpeed, maxSpeed = 0.0;
    static int    last_peak;
    int           nibble;
    static int    MSB, LSB;
    
    static int    count;
    static float  input_windowed[FFT_SIZE_IN*2+1];
    
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	
    
	//    const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
    const SAMPLE *rptr = (const SAMPLE*)&InputSamplesReal12000[0];
	
    SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    long framesToCalc;
    long i,j,cnt;
    int finished;
    unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;
	
    (void) outputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;
    
    unsigned long frameCnt;
	
	if( framesLeft < framesPerBuffer )
    {
        framesToCalc = framesLeft;
        finished = paComplete;
    }
    else if(all_complete == FALSE)
    {
        framesToCalc = framesPerBuffer;
        finished = paContinue;
    }
    else // early out
    {
		finished = paComplete;  
		
		if( framesLeft < framesPerBuffer )
		{
			framesToCalc = framesLeft;
		}
		else 
		{
			framesToCalc = framesPerBuffer;
		}
    }
	
	//==================================================================
    float* in      = (SAMPLE*) inputBuffer;
    float* inSave  = (SAMPLE*) inputBuffer;
    
    //the callback gives frames at 48000 sps  
    //first we must copy all the 48K samples from portaudio
	
	float* in12000 = (SAMPLE*) &InputSamplesReal12000[0];
	
	float tempf;
	
    for (frameCnt = 0; frameCnt < framesToCalc/4; frameCnt++)
    {
        // Get interleaved soundcard samples from input buffer
        //get sample 0,4,8,...
        *in12000++  = *in++;//in[0]
        *in12000++  = *in++; //in[1] right channel not used 
		//but we must pick it anyway
		in += 6;
    }
	
	// 1/4 of the 48K samples are in the 12K buffer  
	
	//===================================================================
	
	//	printf("%d\n", framesToCalc ); fflush(stdout);
	
    if( inSave == NULL )
    {
        for( i=0; i<framesToCalc/4; i++ )
			
        {
            *wptr++ = SAMPLE_SILENCE;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = SAMPLE_SILENCE;  /* right */
        }
    }
    else
    {
        for( frameCnt=0; frameCnt<framesToCalc/4; frameCnt++ )
			
        {			
            
 			leftInput  = *rptr++; //in[0]
			rightInput = *rptr++; //in[1] right channel not used 
			//but we must pick it anyway
			
			InputSamplesI[++k] = leftInput;
			
            *wptr++ = leftInput;  // left 
            *wptr++ = rightInput;  // right  
            
            if(k == FFT_SIZE_IN/16)
			{
				//process the symbol
				for ( j = 0; j < FFT_SIZE_IN/16; j++)
					//circular buffer for symbol "sync"
				{
					
                    
                    // AGC section --------------------------------------------
                    if(agcOFF == TRUE)goto skip_agc;
                    {
                        tmp = fabs(InputSamplesI[j]);
                        
                        if(ss < tmp)
                        { ss = (9. * ss + tmp) / 10.; // slow down attack time
                            pk = tmp;
                            hangcnt = (agcfast == TRUE) ? 1 : 6;
                        }
                        
                        if(agc != FALSE)
                            InputSamplesI[j] *= 0.29 / max(ss, agcthr);
                        
                        else
                            InputSamplesI[j] *= 3.e3 * agcgain;
                        
                        
                        if(hangcnt == 0)
                            ss  *= agc_decay;
                        
                        pk  *= agc_decay;
                    }
                    
                    if(hangcnt > 0)  hangcnt--;
                    
                skip_agc:
                    
                    // --------------------------------------------------------
                    //  Move the samples along 1/16th of a symbol length at
                    //  a time
                    
                    inputI[j]                    = inputI[j +   FFT_SIZE_IN/16];
                    inputI[j +   FFT_SIZE_IN/16] = inputI[j + 2*FFT_SIZE_IN/16];
                    inputI[j + 2*FFT_SIZE_IN/16] = inputI[j + 3*FFT_SIZE_IN/16];
                    inputI[j + 3*FFT_SIZE_IN/16] = inputI[j + 4*FFT_SIZE_IN/16];
                    inputI[j + 4*FFT_SIZE_IN/16] = inputI[j + 5*FFT_SIZE_IN/16];
                    inputI[j + 5*FFT_SIZE_IN/16] = inputI[j + 6*FFT_SIZE_IN/16];
                    inputI[j + 6*FFT_SIZE_IN/16] = inputI[j + 7*FFT_SIZE_IN/16];
                    inputI[j + 7*FFT_SIZE_IN/16] = inputI[j + 8*FFT_SIZE_IN/16];
                    inputI[j + 8*FFT_SIZE_IN/16] = inputI[j + 9*FFT_SIZE_IN/16];
                    inputI[j + 9*FFT_SIZE_IN/16] =inputI[j + 10*FFT_SIZE_IN/16];
                    inputI[j + 10*FFT_SIZE_IN/16]=inputI[j + 11*FFT_SIZE_IN/16];
                    inputI[j + 11*FFT_SIZE_IN/16]=inputI[j + 12*FFT_SIZE_IN/16];
                    inputI[j + 12*FFT_SIZE_IN/16]=inputI[j + 13*FFT_SIZE_IN/16];
                    inputI[j + 13*FFT_SIZE_IN/16]=inputI[j + 14*FFT_SIZE_IN/16];
                    inputI[j + 14*FFT_SIZE_IN/16]=inputI[j + 15*FFT_SIZE_IN/16];
                    inputI[j + 15*FFT_SIZE_IN/16]=inputI[j + 16*FFT_SIZE_IN/16];
                    inputI[j + 16*FFT_SIZE_IN/16]=inputI[j + 17*FFT_SIZE_IN/16];
                    inputI[j + 17*FFT_SIZE_IN/16]= InputSamplesI[j]*4;
					
				}                  
				
				for ( cnt = 0; cnt < FFT_SIZE_IN; cnt++)
				{
					g_fltFftBufRe[cnt] = inputI[cnt]*(0.42-0.5*cos(6.2832*(cnt)/FFT_SIZE_IN) + 0.08*cos(2*6.2832*(cnt)/FFT_SIZE_IN));//Blackman-Something window;
					//no window for use in speed meter
					Speed_Re[cnt]      = inputI[cnt]; 
				} 
				
				
				//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
				//Measure maximum peak amplitude. */
				float maxTEST = 0;
				float averageTEST = 0.0;
				float valTEST;
				
				for( lm=0; lm<FFT_SIZE_IN; lm++ )
				{
					//valTEST = (float)InputSamplesI[lm];
					valTEST = inputI[lm];
					//valTEST = (float)g_fltFftBufRe[lm];
					//valTEST = Speed_Re[lm];
					
					if( valTEST < 0 ) valTEST = -valTEST; /* ABS */
					if( valTEST > maxTEST )
					{
						maxTEST = valTEST;
					}
					averageTEST += valTEST;
				}
				
				averageTEST = averageTEST / (double)FFT_SIZE_IN;
				
				
				//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
				
				
				
				
				//dspmath_MultiplyBlackmanWindow (g_fltFftBufRe, FFT_SIZE_IN);
				//dspmath_MultiplyBlackmanWindow (g_fltFftBufIm, FFT_SIZE_IN);
				//digitalWrite (0, HIGH) ;
/*				
				//FFT for decoding and display
				dspmath_CalcRealFft (FFT_SIZE_IN,    
									 g_fltFftBufRe,
									 g_fltFftBufIm) ;
									 
				//Unwindowed FFT for speed calculation and display on 
				//bar grpah
				dspmath_CalcRealFft (FFT_SIZE_IN,
									 Speed_Re,
									 Speed_Im) ;
*/				
				//digitalWrite (0,  LOW) ;
				
				k=0;
				max = 0.0;
				peak =0;
				
				maxSpeed = 0.0;
				peakSpeed = 0.0;
				
			}
			
		}
		
	}
	
	data->frameIndex += framesToCalc/4;
	return finished;
}


/*******************************************************************/
int main(void);
int main(void)
{
    PaStreamParameters  inputParameters,
	outputParameters;
    PaStream*           stream;
    PaError             err = paNoError;
    paTestData          data;
    int                 i;
    int                 totalFrames;
    int                 numSamples;
    int                 numBytes;
    SAMPLE              max, val;
    double              average;
	
    printf("patest_record.c\n"); fflush(stdout);
	
    data.maxFrameIndex = totalFrames = (NUM_SECONDS * SAMPLE_RATE); /* Record for a few seconds. */
    data.frameIndex = 0;
    numSamples = totalFrames * NUM_CHANNELS;
    numBytes = numSamples * sizeof(SAMPLE);
    data.recordedSamples = (SAMPLE *) malloc( numBytes ); /* From now on, recordedSamples is initialised. */
    if( data.recordedSamples == NULL )
    {
        printf("Could not allocate record array.\n");
        goto done;
    }
    for( i=0; i<numSamples; i++ ) 
		data.recordedSamples[i] = 0;
	
    err = Pa_Initialize();
    if( err != paNoError ) goto done;
	
    inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
    if (inputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default input device.\n");
        goto done;
    }
    inputParameters.channelCount = 2;                    /* stereo input */
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;
	
    /* Record some audio. -------------------------------------------- */
    err = Pa_OpenStream(
						&stream,
						&inputParameters,
						NULL,                  /* &outputParameters, */
						SAMPLE_RATE,
						FRAMES_PER_BUFFER,
						paClipOff,      /* we won't output out of range samples so don't bother clipping them */
						recordCallback,
						&data );
    if( err != paNoError ) goto done;
	
    err = Pa_IsFormatSupported( &inputParameters, NULL, 48000.0 );
    if( err == paFormatIsSupported )
	{
		printf("\n 48000.0 sps is supported\n");
	}
	else
	{
		printf("\n 48000.0 sps is NOT supported\n");
	}
	
    err = Pa_StartStream( stream );
    if( err != paNoError ) 
    {
		goto done;
    }
    printf("\n=== Now recording!!  ===\n"); fflush(stdout);
	
	//__________________________________________________________________________
    while( ( err = Pa_IsStreamActive( stream ) ) == 1 )
    {
        Pa_Sleep(1000);
        //Do some useful Main() processing here
        //printf("index = %d\n", data.frameIndex ); fflush(stdout);
        if(getchar() == '.') 
        {
			all_complete = TRUE;
			printf("complete\n"); fflush(stdout);
        }
    }
    if( err < 0 ) 
    {
		goto done;
    }
	//__________________________________________________________________________
	
    err = Pa_CloseStream( stream );
    if( err != paNoError ) 
    {
		goto done;
    }
	
    /* Measure maximum peak amplitude. */
    max = 0;
    average = 0.0;
    for( i=0; i<numSamples/4; i++ )
    {
        val = data.recordedSamples[i];
        if( val < 0 ) val = -val; /* ABS */
        if( val > max )
        {
            max = val;
        }
        average += val;
    }
	
    average = average / (double)(numSamples/4);
	
    printf("sample max amplitude = "PRINTF_S_FORMAT"\n", max );
    printf("sample average = %lf\n", average );
	
    /* Write recorded data to a file. */
#if WRITE_TO_FILE
    {
        FILE  *fid;
        fid = fopen("recorded.raw", "wb");
        if( fid == NULL )
        {
            printf("Could not open file.");
        }
        else
        {
            fwrite( data.recordedSamples, NUM_CHANNELS * sizeof(SAMPLE), totalFrames/4, fid );
            fclose( fid );
            printf("Wrote data to 'recorded.raw'\n");
        }
    }
#endif
	
    /* Playback recorded data.  -------------------------------------------- */
    data.frameIndex = 0;
	
	
done:
    Pa_Terminate();
    if( data.recordedSamples )       /* Sure it is NULL or valid. */
        free( data.recordedSamples );
    if( err != paNoError )
    {
        fprintf( stderr, "An error occured while using the portaudio stream\n" );
        fprintf( stderr, "Error number: %d\n", err );
        fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
        err = 1;          /* Always return 0 or 1, but no other return codes. */
    }
    return err;
}

