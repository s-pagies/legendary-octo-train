#include <iostream>
#include <fstream>

int mal(int arg1, int arg2)
{
	return arg1 * arg2;
}

int main(int argc, char** argv)
{
	std::ifstream in("input.txt");

	int value1;
	int value2;

	for(int( i = 0; i < 4 ; i++ );
	{
		in >> value1;
		in >> value2;

		int product = mal(value1, value2);
		std::cout << "product is " << dardar << std::endl;
	}
	return 0;
}
