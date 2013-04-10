
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


float *buffer;
int chans, bytes, cnt = 0; pt = 0; ts = 0;
float sr, dur, dtime, *comp, *del, del1[2] = {0.f, 0.f}, del2[1] = {0.f}, del3[2] = {0.f, 0.f};

int audio_callback( const void *inputBuffer, void *outputBuffer,
                    unsigned long framesPerBuffer,
                    const PaStreamCallbackTimeInfo* timeInfo,
                    PaStreamCallbackFlags statusFlags,
                    void *userData 
                  )
{
    float *in = (float*) inputBuffer, *out = (float*)outputBuffer;
    int i;
    
    // bara mono atm
    memcpy(comp, in, bytes);
    flanger(in, line(.0001f, dur, dtime, &ts, def_cr),
        0.4f, dtime, del, &pt, framesPerBuffer, SAMPLING_RATE);
    balance(in, comp, del1, freq, framesPerBuffer, SAMPLING_RATE);
    memcpy(out, in, bytes);

    
    return paContinue;
}

void init_stuff()
{
    float frequency, dgain, dfeedback, dur;
    int i,id;
    const PaDeviceInfo  *info;
    const PaHostApiInfo *hostapi;
    PaStreamParameters outputParameters, inputParameters;
    
    printf("Type the modulator frequency in Hertz: ");
    scanf("%f", &frequency);                        /* get the modulator frequency */

    si = TWO_PI * frequency / SAMPLING_RATE;       /* calculate sampling increment */
    
    freq = frequency;

    // vdtime
    // maxdelay
    // fbk
    printf("Type the max delay time (s) of the delay: ");
    scanf("%f", &dtime);                        /* get the modulator frequency */

    printf("Type the env duration (s) for the delay: ");
    scanf("%f", &dur);                        /* get the modulator frequency */
    
    bytes = sizeof(float) * def_vsize * def_chans;

    buffer = (float*) malloc(bytes);
    memset(buffer, 0, bytes);

    comp = (float*) malloc(bytes);
    memset(comp, 0, bytes);

    del = (float*) malloc(sizeof(float)*((int)(dtime*def_sr)));
    memset(del, 0, sizeof(float)*((int)(dtime*def_sr)));

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

    if(buffer) {
        free(buffer);
        free(comp);
        free(del);
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