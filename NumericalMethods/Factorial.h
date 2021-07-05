
#ifndef FACTORIAL
#define FACTORIAL

int factorial(int i) {
	if (i <= 1) return 1;
	else return i * factorial(i - 1);
}

#endif // !FACTORIAL
