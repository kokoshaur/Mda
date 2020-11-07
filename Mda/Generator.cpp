#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <Windows.h>
#include <iostream>
#include <fstream>
#include "subj.h"

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
	        FFT next = Funct(data, num, num + size);

            std::cout.precision(3);
            std::cout.setf(std::ios::fixed);
            std::cout << num << " - " << num+size << ": " << next.A << "," << next.B << std::endl;
        	
            FILE* f1 = fopen((std::to_string(num) + ".wav").c_str(), "w+b");;

            fwrite(&stream, sizeof(stream), 1, f1);
            fwrite(&data[num], sizeof(short), stream.subchunk2Size / 2, f1);

            num += size;

            fclose(f1);
        	
        }
        out.close();
	}
	
    FFT Funct(short* data, int start, int fin) {
        double summaRe = 0, summaIm = 0, Arg = 0;
        double *Ak = new double[fin - start];

        FFT ansver;
        ansver.A = 0;
        ansver.B = 0;
        for (int i = 0; i < fin - start; i++)
        {
            summaRe = 0; summaIm = 0;
            for (int j = 0; j < fin - start; j++)
            {
                Arg = 2.0 * M_PI * j * i / (fin - start);
                summaRe += cos(Arg) * (data[start + j]);
                summaIm += sin(Arg) * (data[start + j]);
            }
            Ak[i] = sqrt(summaRe * summaRe + summaIm * summaIm);
            out << Ak[i] << std::endl;

        	if(ansver.A < Ak[i])
        	{
                ansver.B = ansver.A;
                ansver.A = Ak[i];
            }
            else if (ansver.B < Ak[i])
                ansver.B = Ak[i];
        }

        ansver.A /= stream.sampleRate;
        ansver.B /= stream.sampleRate;
		
        return ansver;
    }
};
