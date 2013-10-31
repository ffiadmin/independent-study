#ifndef CONFIG_H
#define CONFIG_H

#include <string>

#include "Display.h"

using std::string;

/**
 * Configuration as to whether the program will be allowed to display
 * output to the console or place feedback of the program status within
 * a text file.
 *
 *
 *
 * Here is a lising of possible output level types:
 *  - Critical - Fatal errors which require direct user interaction.
 *  - Error    - A recoverable error which does not halt program 
 *               but may produce in undesirable results.
 *  - Warning  - A recoverable error which is not likely to produce
 *               undesirable results in the program output.
 *  - Info     - A message indicating the current status of the 
 *               program, which may be interesting or helpful in
 *               understanding the current state of the program or
 *               what is currently being calculated.
 *  - Verbose  - Outputs every possible set of log data. This is 
 *               recommended for debugging or for understanding
 *               the inner workings of program flow.
 *
 * Output or logging level recommendataions:
 *  - Production: OUTPUT_LEVEL_ERROR
 *  - Development: OUTPUT_LEVEL_INFO
 *  - Debugging: OUTPUT_LEVEL_VERBOSE
 *  - Utopia: Disabled entirelly
 *
 * IMPORTANT NOTE:
 * Any ouput level will also display messages for all levels above it.
 * For example, a warning level alert will also output messages for 
 * error and critical level alerts.
*/

#define CONSOLE_OUTPUT_ENABLED true
#define OUTPUT_DATA_TO_FILE    true
#define OUTPUT_LEVEL_CRITICAL  1
#define OUTPUT_LEVEL_ERROR     2
#define OUTPUT_LEVEL_WARNING   3
#define OUTPUT_LEVEL_INFO      4
#define OUTPUT_LEVEL_VERBOSE   5

#define LOG_LEVEL OUTPUT_LEVEL_VERBOSE

namespace Config {
	const int FLOAT_PRECISION = 16;
	const string INPUT_FILE = "in.txt";
	const string LOG_FILE   = "info.log";

	const int LOG_COLORS[5] = {
		red,      // Critical color
		darkPink, // Error color
		yellow,   // Warning color
		blue,     // Info color
		green     // Debug color
	};

	const string LOG_NAMES[5] = {
		"Critical Error",
		"Recoverable Error",
		"Warning",
		"Information",
		"Debug Message"
	};
};

#endif