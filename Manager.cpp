#include "Manager.h"

void Manager::msg(char* text, int level) {
	if (CONSOLE_OUTPUT_ENABLED) {
		Display::coloredText(text, Config::LOG_COLORS[level - 1]);
	}
}

void Manager::setFloatPrecision(int precision) {
	cout.setf(std::ios::fixed);
	cout.setf(std::ios::showpoint);
	cout.precision(precision);
}