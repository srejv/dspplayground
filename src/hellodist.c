#include <stdio.h>
#include <math.h>

#include "portaudio.h"

#define FRAME_BLOCK_LEN 256
#define SAMPLING_RATE 44100
#define TWO_PI  (3.14159265f * 2.0f)

PaStream *audioStream;
double si = 0;
double th = 1.0/3.0;
double th2 = 2.0/3.0;
double gain = 0.5;

void expdist(float *input, float *output, unsigned int N, float gain, float mix) {
    unsigned int i;
    for(i=0; i < N; i++) {
        float q = input[i] * gain;
        float z = 1-exp(-fabs(q));

        if(q < 0.0) z *= -1; 

        output[i] = mix * z + (1.0 - mix)*input[i];
    }
}

void overdrive(float *input, float *output, unsigned int N, float gain) {
    unsigned int i;
    for(i=0; i < N; i++) {
        in[i] *= gain;
        float xn_abs = fabs(in[i]);
     
        if (xn_abs < th) {
            out[i] = 2.0 * in[i];
        } 
        else if (xn_abs <= 2th) {
            if (in[i] > 0) { out[i] = (3.0 - pow((2.0 - in[i]  * 3.0), 2)) / 3.0; }
            if (in[i] < 0) { out[i] =-(3.0 - pow((2.0 - xn_abs * 3.0), 2)) / 3.0; }
        } else  { //if (xn_abs > 2th ) {
            if (in[i] > 0.0) out[i] =  1.0;
            if (in[i] < 0.0) out[i] = -1.0;
        }
    }
}

int audio_callback( const void *inputBuffer, void *outputBuffer,
                    unsigned long framesPerBuffer,
                    const PaStreamCallbackTimeInfo* timeInfo,
                    PaStreamCallbackFlags statusFlags,
                    void *userData 
                  )
{
    float *in = (float*) inputBuffer, *out = (float*)outputBuffer;
    static double phase = 0;
    unsigned long i;

    // Passthrough
    //for( i=0; i < framesPerBuffer; i++ ) {
        //in[i] *= gain;
        //*out++ = *in++ * gain;
    //}
    overdrive(in, out, framesPerBuffer, 0.85f);
    
    return paContinue;
}

void init_stuff()
{
    float frequency;
    int i,id;
    const PaDeviceInfo  *info;
    const PaHostApiInfo *hostapi;
    PaStreamParameters outputParameters, inputParameters;
    
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
}

int main()
{
  int ch;
    init_stuff();
    while(getchar() != ' ') Pa_Sleep(100);
    terminate_stuff();
    return 0;
}