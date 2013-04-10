
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#include "portaudio.h"
#include "mtapdelay.h"

#define FRAME_BLOCK_LEN 256
#define SAMPLING_RATE 44100
#define TWO_PI  (3.14159265f * 2.0f)

PaStream *audioStream;

MTAP *mtap;

long *tapindex;         // ntaps
float *taptime;        // ntaps
float *tapamp;         // ntaps



int audio_callback( const void *inputBuffer, void *outputBuffer,
                    unsigned long framesPerBuffer,
                    const PaStreamCallbackTimeInfo* timeInfo,
                    PaStreamCallbackFlags statusFlags,
                    void *userData 
                  )
{
    float *in = (float*) inputBuffer, *out = (float*)outputBuffer;
    int i;
    // left first
    // right second
    mtap_process(mtap, in, out, framesPerBuffer, 1, 0.9, tapindex, taptime, tapamp);

    return paContinue;
}

void init_stuff()
{
    float frequency, dgain, dfeedback, dur, ntaps;
    double maxtime;
    int i,id;
    const PaDeviceInfo  *info;
    const PaHostApiInfo *hostapi;
    PaStreamParameters outputParameters, inputParameters;
    
    //printf("Type the modulator frequency in Hertz: ");
    //scanf("%f", &frequency);                        /* get the modulator frequency */
    printf("herp");

    ntaps = 4;
    maxtime = 2.0;

    printf("Allocating memory...\n");
    tapindex = (long*)malloc(sizeof(long)*ntaps);
    taptime = (float*)malloc(sizeof(float)*ntaps);
    tapamp = (float*)malloc(sizeof(float)*ntaps);
    
    printf("Setting up taps...\n");
    taptime[0] = 0.3;
    taptime[1] = 0.5;
        taptime[2] = 1.2;
    taptime[3] = 1.5;

    tapamp[0] = 0.05;
    tapamp[1] = 0.20;
    tapamp[2] = 0.88;
    tapamp[3] = 0.15;

    printf("Creating mtap object...\n");
    mtap = new_mtap(tapindex, ntaps, maxtime, taptime, tapamp);

    scanf("%f", &frequency);
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
    outputParameters.channelCount = 2;                           /* stereo output */
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
    inputParameters.channelCount = 2;                              /* stereo input */
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
    
    free(taptime);
    free(tapamp);

    mtap_free(mtap);
}

int main()
{
    int ch;
    init_stuff();
    while(getchar() != ' ') Pa_Sleep(100);
    terminate_stuff();
    return 0;
}