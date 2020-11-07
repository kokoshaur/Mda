struct WAVHEADER
{
    char chunkId[4];// Содержит символы "RIFF" в ASCII кодировке
    unsigned long chunkSize;// 36 + subchunk2Size, или более точно:
    char format[4];// Содержит символы "WAVE"
    char subchunk1Id[4];// Содержит символы "fmt "
    unsigned long subchunk1Size;// 16 для формата PCM.
    unsigned short audioFormat;// Для PCM = 1 (то есть, Линейное квантование).
    unsigned short numChannels;// Количество каналов. Моно = 1, Стерео = 2 и т.д.
    unsigned long sampleRate;// Частота дискретизации. 8000 Гц, 44100 Гц и т.д.
    unsigned long byteRate;// sampleRate * numChannels * bitsPerSample/8
    unsigned short blockAlign;// Количество байт для одного сэмпла, включая все каналы. numChannels * bitsPerSample/8
    unsigned short bitsPerSample;// Так называемая "глубиная" или точность звучания. 8 бит, 16 бит и т.д.
    char subchunk2Id[4];// Подцепочка "data" содержит аудио-данные и их размер.Содержит символы "data"
    unsigned long subchunk2Size;// Количество байт в области данных. numSamples * numChannels * bitsPerSample/8
};

struct FFT {
    long double A;
    long double B;
};