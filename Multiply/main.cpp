#include <iostream>

float mal(float arg1, float arg2)
{
	return arg1 * arg2;
}

int main()
{
	float dardar = mal(1.41,1.41);
	std::cout << "product is " << dardar << std::endl;
	return 0;
}
