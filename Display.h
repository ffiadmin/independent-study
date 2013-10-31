/**
 * Display class
 *
 * This class is designed to manipulate the text output in the 
 * console in the following ways:
 *  - Change the text color.
 *  - Clear the contents of the console.
 *  - Evenly distribute several "columns" of text into a "row" in a "table".
 *
 * @author    Oliver Spryn
 * @namespace global
 * @package   display
*/

#ifndef DISPLAY_H
#define DISPLAY_H

//#include "DistributeParametersInvalid.h"
#include <string>
#include <vector>

using std::string;
using std::vector;

class Display {
public : 
        static void clear();
        static void coloredText(string text, int color);
		static void displayColorChart();
        static string distribute(vector<string> text, vector<char> align)/* throw(DistributeParametersInvalid)*/;
        static void resetTextColor();
        static void setTextColor(int color);

        static const int CONSOLE_HEIGHT = 25;
        static const int CONSOLE_WIDTH = 80;
        static const int DEFAULT_CONSOLE_COLOR = 7;
};

enum colors {
        black = 0,
        darkBlue = 1,
        darkGreen = 2,
        darkAqua, darkCyan = 3,
        darkRed = 4,
        darkPurple = 5, darkPink = 5, darkMagenta = 5,
        darkYellow = 6,
        darkWhite = 7,
        gray = 8,
        blue = 9,
        green = 10,
        aqua = 11, cyan = 11,
        red = 12,
        purple = 13, pink = 13, magenta = 13,
        yellow = 14,
        white = 15,

		gBlack = 128,
        gDarkBlue = 129,
        gDarkGreen = 130,
        gDarkAqua, gDarkCyan = 131,
        gDarkRed = 132,
        gDarkPurple = 133, gDarkPink = 133, gDarkMagenta = 133,
        gDarkYellow = 134,
        gDarkWhite = 135,
        gGray = 136,
        gBlue = 137,
        gGreen = 138,
        gAqua = 139, gCyan = 139,
        gRed = 140,
        gPurple = 141, gPink = 141, pMagenta = 141,
        gYellow = 142,
        gWhite = 143
};

#endif