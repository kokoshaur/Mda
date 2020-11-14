#include "Generator.cpp"
#define ONE_STEP 800

int main(int argc, char* argv[])
{
	Generator r("Test.wav");

	unsigned short* buffer = r.LoadWave("4.wav");
	r.SaveWave(buffer);

	r.SaveByStep(buffer, ONE_STEP, 10);
	return 0;
}
