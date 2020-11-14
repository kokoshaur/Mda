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
    std::ofstream mda;
	
    double Func(int index, double frequency)
    {
        return sin((float)(2 * M_PI * index * frequency / stream.sampleRate));
    }
	
public:
	Generator(std::string path)
	{
        this->path = path;
	}

    unsigned short* Wave(double frequency)
	{
        unsigned short* data = new unsigned short[stream.subchunk2Size];

        for (int index = 0; index < stream.subchunk2Size; index++)
        {
            data[index] = (unsigned short)(Func(index, frequency) * SHRT_MAX);// Приводим уровень к амплитуде от 32767 до -32767.
        }
        return data;
    }

    unsigned short* LoadWave(std::string pathTo)
	{
        FILE* f1;
        f1 = fopen(pathTo.c_str(), "rb");
        fread(&stream, sizeof(stream), 1, f1);
		
        unsigned short* data = new unsigned short[stream.subchunk2Size];
		
        fread(data, stream.blockAlign, stream.subchunk2Size, f1);
        fclose(f1);
		
        return data;
	}

    void SaveWave(unsigned short* data)
    {
        FILE* f1;
        f1 = fopen(path.c_str(), "w+b");
        fwrite(&stream, sizeof(stream), 1, f1);
        fwrite(data, stream.blockAlign, stream.subchunk2Size, f1);
		
        fclose(f1);
    }

	void SaveByStep(unsigned short* data, unsigned long size, int count)
	{
        int num = 0;
        stream.subchunk2Size = size;
		out.setf(std::ios::fixed);
        for(int i = 0; i < count; i++)
            num = save(data, num, size, i);
	}

	int save(unsigned short* data, int num, int size, int i)
	{
        out.open("before\\ansver" + std::to_string(i) + ".txt");
        mda.open("after\\MDA" + std::to_string(i) + ".txt");
        std::cout.precision(3);
        std::cout.setf(std::ios::fixed);

        FFT a = Funct(data, num, num + size, i);
        std::cout << num << " - " << num + size << ": " << a.A << "," << a.B  << std::endl;

        FILE* f1 = fopen((std::to_string(num) + ".wav").c_str(), "w+b");

        fwrite(&stream, sizeof(stream), 1, f1);
        fwrite(&data[num], stream.blockAlign, stream.subchunk2Size, f1);

        fclose(f1);
        out.close();

		return num + size;
	}
	
    FFT Funct(unsigned short* data, int start, int fin, int k)
	{
        std::vector<complex<long double>> u(stream.subchunk2Size);
        for (int i = 0; i < stream.subchunk2Size/2; i++)
        {
            u[i] = data[start + i];
            out << data[start + i] << std::endl;
        }
        FFT a;
        vector<vector<complex<long double>>> F = MatrixFlexer::makeMatrixFourier(800);

        for (std::complex<double> element : MatrixFlexer::matrixProduct(F, u, 0))
        {
            mda << abs(element) << std::endl;
	        if (abs(element) > a.A)
	        {
                a.B = a.A;
                a.A = abs(element);
	        }else if (element.real() > a.B)
	        {
                a.B = abs(element);
	        }
        }

        mda.close();
        return a;
    }
};
