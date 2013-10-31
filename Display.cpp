#include <iostream>
#include <Windows.h>

#include "Display.h"

using std::cout;
using std::endl;

/**
 * This will clear the text in the console.
 *
 * @access public
 * @return void
*/

void Display::clear() {
        system("CLS");
}

/**
 * This function is designed to output colored or highlighted text to the
 * console.
 *
 * @access public
 * @param  string text  The text to display to the console
 * @param  int    color The integer representation of the text or highligher color to use
 * @return void
*/

void Display::coloredText(string text, int color) {
        HANDLE cHandle;
        cHandle = GetStdHandle(STD_OUTPUT_HANDLE);

        SetConsoleTextAttribute(cHandle, color);
        cout << text;
        SetConsoleTextAttribute(cHandle, Display::DEFAULT_CONSOLE_COLOR);
}

void Display::displayColorChart() {
	cout << "Displaying the console color code then some sample text..." << endl << endl;

	for (int i = 0; i <= 255; ++i) {
	//Set the text color
		Display::setTextColor(i);

	//Display the outputted text
		cout << i << ". Sample text";

		if (i <= 9) {
			cout << "  ";
		} else if (i <= 99) {
			cout << " ";
		}

	//Reset the background color of the text for the "gaps" in the columns of the chart
		Display::resetTextColor();
		
	//Render the gaps between the columns in the chart
		cout << "   ";

	//Create a new row
		if ((i + 1) % 4 == 0) {
			cout << endl << endl;
		}
	}

	Display::resetTextColor();

	cout << endl << endl << endl << endl;
}

/**
 * This function will take a vector of strings and lay them out such that each
 * string will have its own "cell" in an evenly-spaced, invisible table that
 * stretches the length of the console window.
 *
 * For example, if these strings are given:
 *  - one
 *  - two
 *  - three
 *  - four
 *
 * they will be laid out like this (except without the table border, which is 
 * show here to demonstrate this function's formatting capibilities):
 *
 * +------------------+------------------+------------------+-----------------+
 * |one               |two               |three             |four             |
 * +------------------+------------------+------------------+-----------------+
 *
 * The align parameter will accept a vector of characters, where each character
 * may have one of the following values:
 *  - L = align left
 *  - C = align center
 *  - R = align right
 *
 * The string of text is matched to the associated align character by its index
 * in each of the vectors.
 *
 * @access public
 * @param  vector<string>              text  A vector of strings containing text to be aligned
 * @param  vector<char>                align A vector of characters symbolizing the string's alignment within the respective "cell"
 * @return string                      ret   A string formatted to the parameter's specifications
 * @throws DistributeParametersInvalid       Encountered if parameter vector sizes do not match or if an invalid character appears in the align vector
*/

string Display::distribute(vector<string> text, vector<char> align)/* throw(DistributeParametersInvalid)*/ {
        int length = text.size();
        string ret = "";
        
//Does the size of each vector parameter match?
        if (text.size() != align.size()) {
                string error = "The number of strings in the text parameter (";
                error += text.size();
                error += ") does not match the number of characters in the align parameter (";
                error += align.size();
                error += ").";

                //throw DistributeParametersInvalid(error);
        }

//Calculate the number of characters we can allot per cell with the given input
        int cellW = static_cast<int>(floor(static_cast<double>(Display::CONSOLE_WIDTH) / static_cast<double>(length)));
        int extra = Display::CONSOLE_WIDTH - (cellW * length);
        int leftW = 0;
        int rightW = 0;
        
//Print out each "cell"
        for (int i = 0; i < length; i++) {
                switch(align[i]) {
                        case 'L' : 
                                ret += text[i];

                        //Add any padding needed to the right side of the text
                                rightW = 0;
                                rightW = cellW - text[i].length();

                                for (int k = 0; k < rightW; k++) {
                                        ret += " ";
                                }

                        //Usually not all cells will have an even amount of space, so in the final iterations of this loop
                        //output an additional space in the last few cells to ensure the the returned string is still
                        //Display::CONSOLE_WIDTH characters wide
                                if (extra + i >= length) {
                                        ret += " ";
                                }

                                break;

                        case 'C' : 
                        //Add any padding needed to the left side of the text
                                leftW = 0;
                                leftW = static_cast<int>(floor(static_cast<double>(cellW - text[i].length()) / 2.0));

                                for (int j = 0; j < leftW; j++) {
                                        ret += " ";
                                }

                                ret += text[i];

                        //Add any padding needed to the right side of the text
                                rightW = 0;
                                rightW = static_cast<int>(ceil(static_cast<double>(cellW - text[i].length()) / 2.0));

                                for (int k = 0; k < rightW; k++) {
                                        ret += " ";
                                }

                        //Usually not all cells will have an even amount of space, so in the final iterations of this loop
                        //output an additional space in the last few cells to ensure the the returned string is still
                        //Display::CONSOLE_WIDTH characters wide
                                if (extra + i >= length) {
                                        ret += " ";
                                }

                                break;

                        case 'R' : 
                        //Usually not all cells will have an even amount of space, so in the final iterations of this loop
                        //output an additional space in the last few cells to ensure the the returned string is still
                        //Display::CONSOLE_WIDTH characters wide
                                if (extra + i >= length) {
                                        ret += " ";
                                }

                        //Add any padding needed to the right side of the text
                                leftW = 0;
                                leftW = cellW - text[i].length();

                                for (int j = 0; j < leftW; j++) {
                                        ret += " ";
                                }

                                ret += text[i];
                                break;

                        default : 
                                string error = "An invalid align type wan encountered: ";
                                error += align[i];
                                error += ". Valid align types include: L (left), C (center), or R (right).";

                                //throw DistributeParametersInvalid(error);
                                break;
                }
        }
        
        return ret;
}

/**
 * This function will reset the output text color to default.
 *
 * @access public
 * @return void
*/

void Display::resetTextColor() {
        HANDLE cHandle;
        cHandle = GetStdHandle(STD_OUTPUT_HANDLE);

        SetConsoleTextAttribute(cHandle, Display::DEFAULT_CONSOLE_COLOR);
}

/**
 * This function is designed to set the text color or highlighter color for all
 * text which will be outputted to the console after this function call.
 *
 * @access public
 * @param  int    color The integer representation of the text or highligher color to use
 * @return void
*/

void Display::setTextColor(int color) {
        HANDLE cHandle;
        cHandle = GetStdHandle(STD_OUTPUT_HANDLE);

        SetConsoleTextAttribute(cHandle, color);
}