
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#include "portaudio.h"
#include "snd_def.h"

#define FRAME_BLOCK_LEN 256
#define SAMPLING_RATE 44100
#define TWO_PI  (3.14159265f * 2.0f)

PaStream *audioStream;
double si = 0;
double freq = 0;

float gain, pitch, fdb, sr, *deray, *env, rp;
int taps, wp, steps, dsize;

int audio_callback( const void *inputBuffer, void *outputBuffer,
                    unsigned long framesPerBuffer,
                    const PaStreamCallbackTimeInfo* timeInfo,
                    PaStreamCallbackFlags statusFlags,
                    void *userData 
                  )
{
    // setup variables
    float *in = (float*) inputBuffer, *out = (float*)outputBuffer;
    int rpi, ep, i, j;
    float s = 0.f, rpf, frac, next, p = pitch; //(pitch*1.5) + 0.5;

    // processing loop
    for( i = 0; i < framesPerBuffer; i++ ) {

        // taps loop (taps = 2)
        for( j = 0; j < taps; j++ ) {

            // tap position, offset
            rpf = rp + j * dsize / taps;
            rpf = rpf < dsize ? rpf : rpf - dsize;
            rpi = (int) rpf;
            frac = rpf - rpi;
            next = (rpi != dsize-1 ? deray[ rpi + 1 ] : deray[ 0 ] );

            // envelope index
            ep = rpi - wp;
            if ( ep < 0 ) ep += dsize;
            s += ( deray[ rpi ]  + frac * (next - deray[ rpi ] ) ) * env[ ep ];
        }

        // increment reader point and check bounds
        rp += p;
        rp = rp < dsize ? rp : rp - dsize;

        // feed the delay line
        deray[ wp ] = *in++ + s * fdb;

        // output the signal
        *out++ = ( deray[ wp - 1 ] * gain + ( s / taps ) * gain);
        s = 0.f;

        // increment the write pointer
        wp = ( wp < dsize ? wp + 1 : 0 );
    }
   
    return paContinue;
}

void init_stuff()
{
    float frequency, dgain, dfeedback, dur;
    int i,id;
    const PaDeviceInfo  *info;
    const PaHostApiInfo *hostapi;
    PaStreamParameters outputParameters, inputParameters;
    
    si = TWO_PI * frequency / SAMPLING_RATE;       /* calculate sampling increment */
    
    // Input parameters
    printf("Type the steps you wish to transpose (integer) of the delay: ");
    scanf("%d", &steps);

    pitch = pow(2, (float)steps/12.0);
    printf("%f\n", pitch);

    printf("Type the gain (0.0 -> 1.0) of the delay: ");
    scanf("%f", &gain);

    printf("Type the feedback (0.0 -> 1.0) of the delay: ");
    scanf("%f", &fdb);

    printf("Type the number of taps (integer, use 2) of the delay: ");
    scanf("%d", &taps);

    sr = SAMPLING_RATE;
    dsize = (int)(0.045*sr);
    
    // allocate a 1-sec delay line and envelope
    deray = (float*)malloc(sizeof(float)*SAMPLING_RATE);
    memset(deray, 0, sizeof(float)*SAMPLING_RATE);
    env = (float*)malloc(sizeof(float)*SAMPLING_RATE);
    memset(env, 0, sizeof(float)*SAMPLING_RATE);

    // write a triangluar env (ramp up, ramp down)
    for( i = 0; i < dsize/2; i++) env[ i ] = i * 2. / dsize;
    for( i = dsize/2; i >= 0; i--) env[ (int)dsize - i - 1 ] = i * 2./dsize;

    wp = 0;
    rp = 0.0f;

    printf("Initializing Portaudio. Please wait...\n");
    Pa_Initialize();                                       /* initialize portaudio */

    for (i=0;i < Pa_GetDeviceCount(); i++) { 
        info = Pa_GetDeviceInfo(i);         /* get information from current device */
        hostapi = Pa_GetHostApiInfo(info->hostApi); /*get info from curr. host API */

        if (info->maxOutputChannels > 0)         /* if curr device supports output */
            printf("%d: [%s] %s (output)\n",i, hostapi->name, info->name );  
    }
    
    printf("\nType AUDIO output device number: ");
    scanf("%d", &id);                   /* get the output device id from the user */
    info = Pa_GetDeviceInfo(id);       /* get chosen device information structure */
    hostapi = Pa_GetHostApiInfo(info->hostApi);         /* get host API structure */
    printf("Opening AUDIO output device [%s] %s\n", hostapi->name, info->name);

    outputParameters.device = id;                             /* chosen device id */
    outputParameters.channelCount = def_chans;                           /* stereo output */
    outputParameters.sampleFormat = paFloat32;    /* 32 bit floating point output */
    outputParameters.suggestedLatency = info->defaultLowOutputLatency;/* set default */
    outputParameters.hostApiSpecificStreamInfo = NULL;        /* no specific info */

    for (i=0;i < Pa_GetDeviceCount(); i++) {
        info = Pa_GetDeviceInfo(i);         /* get information from current device */
        hostapi = Pa_GetHostApiInfo(info->hostApi); /*get info from curr. host API */

        if (info->maxInputChannels > 0)           /* if curr device supports input */
            printf("%d: [%s] %s (input)\n",i, hostapi->name, info->name );  
    }
    
    printf("\nType AUDIO input device number: ");
    scanf("%d", &id);                     /* get the input device id from the user */
    info = Pa_GetDeviceInfo(id);        /* get chosen device information structure */
    hostapi = Pa_GetHostApiInfo(info->hostApi);          /* get host API structure */
    printf("Opening AUDIO input device [%s] %s\n", hostapi->name, info->name);

    inputParameters.device = id;                               /* chosen device id */
    inputParameters.channelCount = def_chans;                              /* stereo input */
    inputParameters.sampleFormat = paFloat32;      /* 32 bit floating point output */
    inputParameters.suggestedLatency = info->defaultLowInputLatency; /*set default */
    inputParameters.hostApiSpecificStreamInfo = NULL;          /* no specific info */

    Pa_OpenStream(               /* open the PaStream object and get its address */
              &audioStream,      /* get the address of the portaudio stream object */
              &inputParameters,  /* provide output parameters */
              &outputParameters, /* provide input parameters */
              SAMPLING_RATE,     /* set sampling rate */
              FRAME_BLOCK_LEN,   /* set frames per buffer */
              paClipOff,         /* set no clip */
              audio_callback,    /* provide the callback function address */
              NULL );            /* provide no data for the callback */

    Pa_StartStream(audioStream); /* start the callback mechanism */
    printf("running... press space bar and enter to exit\n");
}

void terminate_stuff()
{
    Pa_StopStream( audioStream );    /* stop the callback mechanism */
    Pa_CloseStream( audioStream );   /* destroy the audio stream object */
    Pa_Terminate();                  /* terminate portaudio */

    if(deray) {
        free(deray);
        free(env);
    }
}

int main()
{
    int ch;
    init_stuff();
    while(getchar() != ' ') Pa_Sleep(100);
    terminate_stuff();
    return 0;
}