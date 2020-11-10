#define _USE_MATH_DEFINES
#include <vector>
#include <complex>
static class MatrixFlexer
{
public:
	static std::vector<std::vector<std::complex<double>>> createMatrixF(int n) {
		std::complex<double> comp;
		comp = -1;
		comp = std::sqrt(comp);
		std::vector<std::vector<std::complex<double>>> matrix(n);
		for (int i = 0; i < matrix.size(); i++) {
			matrix[i].resize(n);
			for (int j = 0; j < matrix[i].size(); j++) {
				matrix[i][j] = exp((comp * (double)i * 2.0 * M_PI * (double)j) / (double)n);
			}
		}
		return matrix;
	}

	static std::vector<std::vector<std::complex<double>>> transponse(std::vector<std::vector<std::complex<double>>> matrix) {
		std::vector<std::vector<std::complex<double>>> res(matrix.size());

		for (int i = 0; i < matrix.size(); ++i)
			res[i].resize(matrix.size());

		for (int i = 0; i < matrix.size(); ++i) {
			for (int j = 0; j < res.size(); ++j)
				res[j][i] = matrix[i][j];
		}

		return res;
	}

	static std::vector<std::vector<std::complex<double>>> createMatrixFh(int n) {
		std::complex<double> comp;
		comp = -1;
		comp = std::sqrt(comp);
		std::vector<std::vector<std::complex<double>>> matrix = createMatrixF(n);
		for (int i = 0; i < matrix.size(); i++) {
			matrix.resize(n);
			for (int j = 0; j < matrix[i].size(); j++) {
				matrix[i][j] = exp((-comp * (double)i * 2.0 * M_PI * (double)j) / (double)n);
			}
		}
		matrix = transponse(matrix);
		return matrix;
	}

	static void calc_functions(std::vector<double>& func_y, std::vector<double>& x_values, double f, double a, double b, double dt) {
		double y = 0;
		double t = a;
		while (t <= b) {
			y = sin(2 * M_PI * f * t);
			x_values.push_back(t);
			func_y.push_back(y);
			t = t + dt;
		}
	}

	template <typename T>
	static std::vector<std::complex<double>> productMatrix(T vector, std::vector<std::vector<std::complex<double>>> matrix) {
		int length = vector.size();
		std::vector<std::complex<double>> res(length);
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < length; j++)
			{
				res[j] = 0;
				for (int k = 0; k < length; k++)
					res[j] += vector[k] * matrix[k][j];
			}
		}
		return res;
	}

	static std::vector<std::complex<double>> fourier(std::vector<double> func, int k) {
		if (k == 0)
			return std::vector<std::complex<double>> (0);
		std::complex<double> comp = -1;
		std::vector<std::vector<std::complex<double>>> res(k);
		std::vector<std::complex<double>> ansver(func.size());
		
		for (int i = 0; i < res.size(); i++) {
			res[i].resize(func.size());
			for (int j = 0; j < res[i].size(); j++)
				res[i][j] = exp((-comp * (double)i * 2.0 * M_PI * (double)j) / (double)func.size());
		}

		for (int i = 0; i < res.size(); i++)
			for (int j = 0; j < res[i].size(); j++)
				ansver[i] += func[i] * res[i][j];

		return ansver;
	}

	static std::vector<std::complex<double>> directFourierTransform(std::vector<double> func) {
		std::vector<std::vector<std::complex<double>>> matrixFh = createMatrixFh(func.size());
		std::vector<std::complex<double>> res = productMatrix(func, matrixFh);
		for (int i = 0; i < res.size(); i++) {
			res[i] = res[i] / (std::complex<double>)res.size();
		}
		return res;
	}

	static std::vector<std::complex<double>> inverseFourierTransform(std::vector<std::complex<double>> func) {
		std::vector<std::vector<std::complex<double>>> matrixF = createMatrixF(func.size());
		std::vector<std::complex<double>> res = productMatrix(func, matrixF);

		return res;
	}

	template <typename T>
	static T alpha_product(T func, double alpha) {
		T res(func.size());
		for (int i = 0; i < res.size(); i++) {
			res[i] = func[i] * alpha;
		}

		return res;
	}

	template <typename T>
	static T vector_sum(T v1, T v2) {
		T res(v1.size());
		for (int i = 0; i < res.size(); i++) {
			res[i] = v1[i] + v2[i];
		}

		return res;
	}

	static bool isLinear(std::vector<double> v1, std::vector<double> v2, double a, double b) {
		std::vector<std::complex<double>> directFourier_v1 = directFourierTransform(v1);
		std::vector<std::complex<double>> directFourier_v2 = directFourierTransform(v2);
		std::vector<std::complex<double>> temp1 = alpha_product(directFourier_v1, a);
		std::vector<std::complex<double>> temp2 = alpha_product(directFourier_v2, b);
		std::vector<std::complex<double>> sum1 = vector_sum(temp1, temp2);


		std::vector<double> temp_1 = alpha_product(v1, a);
		std::vector<double> temp_2 = alpha_product(v2, b);
		std::vector<double> sum2 = vector_sum(temp_1, temp_2);
		std::vector<std::complex<double>> directFourier_sum2 = directFourierTransform(sum2);

		double eps = 0.01;
		for (int i = 0; i < sum1.size(); i++) {
			if (abs(sum1[i].real() - directFourier_sum2[i].real()) > eps) return false;
		}

		return true;
	}

	static bool isEqualsParseval(std::vector<double> v) {
		std::complex<double> sum1;
		for (int i = 0; i < v.size(); i++) {
			sum1 = sum1 + v[i] * v[i];
		}

		std::vector<std::complex<double>> directFourier_v = directFourierTransform(v);
		std::complex<double> sum2;
		for (int i = 0; i < directFourier_v.size(); i++) {
			sum2 = sum2 + directFourier_v[i] * directFourier_v[i];
		}
		sum2 = sum2 / (double)directFourier_v.size();

		double eps = 0.001;
		if (abs(sum1.real() - sum2.real()) > eps) return true;
		return false;
	}

	static bool isShift(std::vector<double> v1, std::vector<double> v2) {

		std::vector<std::complex<double>> directFourier_v1 = directFourierTransform(v1);
		std::vector<std::complex<double>> directFourier_v2 = directFourierTransform(v2);

		std::vector<std::complex<double>> phase_spectrum_v1(v1.size());
		for (int i = 0; i < v1.size(); i++) {
			phase_spectrum_v1[i] = atan(directFourier_v1[i].imag() / directFourier_v1[i].real());
			//std::cout << phase_spectrum_v1[i] << ' ';
		}

		std::vector<std::complex<double>> phase_spectrum_v2(v2.size());
		for (int i = 0; i < v2.size(); i++) {
			phase_spectrum_v2[i] = atan(directFourier_v2[i].imag() / directFourier_v2[i].real());
			//std::cout << phase_spectrum_v2[i] << ' ';
		}

		double eps = 0.01;

		for (int i = 0; i < v1.size(); i++) {
			if (abs(directFourier_v1[i].real() - directFourier_v2[i].real()) > eps) return false;
		}

		for (int i = 0; i < v1.size(); i++) {
			if (phase_spectrum_v1[i].real() != phase_spectrum_v2[i].real() && abs(phase_spectrum_v1[i].real() - phase_spectrum_v2[i].real()) > eps) return true;
		}


		//свойство сдвига по времени выполняется в любом случае, потому что я делаю
		//преобразование только по игрекам, не учитывая значения иксов
		//то есть сам спектр не изменяется, а только сдвигается по оси времени
		return false;
	}
};
