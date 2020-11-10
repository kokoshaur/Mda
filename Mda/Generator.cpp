#define _USE_MATH_DEFINES
#include <complex>
#include <iostream>
#include <string>
#include <Windows.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "subj.h"
#include "MatrixFlexer.cpp"

class Generator {
private:
    std::string path;
	WAVHEADER stream;
    std::ofstream out;
	
    double Func(int index, double frequency)
    {
        return sin((float)(2 * M_PI * index * frequency / stream.sampleRate));
    }
	
public:
	Generator(std::string path)
	{
        this->path = path;
	}

    short* Wave(double frequency)
	{
        short* data = new short[stream.subchunk2Size];

        for (int index = 0; index < stream.subchunk2Size; index++)
        {
            data[index] = (short)(Func(index, frequency) * SHRT_MAX);// Приводим уровень к амплитуде от 32767 до -32767.
        }
        return data;
    }

    short* LoadWave(std::string pathTo)
	{
        FILE* f1;
        f1 = fopen(pathTo.c_str(), "rb");
        fread(&stream, sizeof(stream), 1, f1);

        short* data = new short[stream.subchunk2Size / 2];
		
        fread(data, sizeof(short), stream.subchunk2Size/2, f1);
        fclose(f1);

        return data;
	}

    void SaveWave(short* data)
    {
        FILE* f1;
        f1 = fopen(path.c_str(), "w+b");
		
        fwrite(&stream, sizeof(stream), 1, f1);
        fwrite(data, sizeof(short), stream.subchunk2Size/2, f1);
		
        fclose(f1);
    }

	void SaveByStep(short* data, unsigned long size, int count)
	{
        int num = 0;
        stream.subchunk2Size = size;
        out.open("ansver.txt");
		out.setf(std::ios::fixed);
        for(int i = 0; i < count; i++)
        {
            std::cout.precision(3);
            std::cout.setf(std::ios::fixed);

            FFT a = Funct(data, num, num + size, i);
            std::cout << num << " - " << num+size << ": " << a.A << "," << a.B << std::endl;
        	
            FILE* f1 = fopen((std::to_string(num) + ".wav").c_str(), "w+b");

            fwrite(&stream, sizeof(stream), 1, f1);
            fwrite(&data[num], sizeof(short), stream.subchunk2Size / 2, f1);

            num += size;

            fclose(f1);
        	
        }
        out.close();
	}
	
    FFT Funct(short* data, int start, int fin, int k)
	{
        std::vector<double> f(fin - start);
        for (int i = 0; i < fin - start; i++)
            f[i] = data[start + i];

        FFT a;

        for (std::complex<double> element : MatrixFlexer::directFourierTransform(f))
        {
	        if (element.real() > a.A)
	        {
                a.B = a.A;
                a.A = element.real();
	        }else if (element.real() > a.B)
	        {
                a.B = element.real();
	        }
        }
		
        return a;
    }
};
