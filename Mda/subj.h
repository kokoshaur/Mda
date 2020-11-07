struct WAVHEADER
{
    char chunkId[4];// �������� ������� "RIFF" � ASCII ���������
    unsigned long chunkSize;// 36 + subchunk2Size, ��� ����� �����:
    char format[4];// �������� ������� "WAVE"
    char subchunk1Id[4];// �������� ������� "fmt "
    unsigned long subchunk1Size;// 16 ��� ������� PCM.
    unsigned short audioFormat;// ��� PCM = 1 (�� ����, �������� �����������).
    unsigned short numChannels;// ���������� �������. ���� = 1, ������ = 2 � �.�.
    unsigned long sampleRate;// ������� �������������. 8000 ��, 44100 �� � �.�.
    unsigned long byteRate;// sampleRate * numChannels * bitsPerSample/8
    unsigned short blockAlign;// ���������� ���� ��� ������ ������, ������� ��� ������. numChannels * bitsPerSample/8
    unsigned short bitsPerSample;// ��� ���������� "��������" ��� �������� ��������. 8 ���, 16 ��� � �.�.
    char subchunk2Id[4];// ���������� "data" �������� �����-������ � �� ������.�������� ������� "data"
    unsigned long subchunk2Size;// ���������� ���� � ������� ������. numSamples * numChannels * bitsPerSample/8
};

struct FFT {
    long double A;
    long double B;
};