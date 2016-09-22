/**
 * \file debugging.hh
 *
 * Useful macros for debugging and warning output. These macros behave
 * exactly like cout or cerr, except append useful text about where
 * the warning or debug message occurred.
 *
 * Use DEBUG_OUT to display a debugging message whenever DEBUG is
 * defines as a number greater than 0.
 *
 * Use DEBUG_OUT_L(level) to write debugging messages that should only
 * print out if DEBUG is set to  "level" or above
 *
 * use WARN_OUT to print red warning messages (regardless of the value
 * of DEBUG)
 *
 * example:
 *
 * \code
 *
 * #define DEBUG 2
 * #include "utilities/debugging.hh"
 *
 * void myfunc() {
 *   WARN_OUT << "This is a warning" << endl;
 *   DEBUG_OUT << "This only prints when DEBUG is greater than 0" << endl;
 *   DEBUG_OUT_L(1) << "Equivalent to DEBUG_OUT" << endl;
 *   DEBUG_OUT_L(2) << "Prints out, since DEBUG==2" << endl;
 *   DEBUG_OUT_L(3) << "doesn't print, since DEBUG==2" << endl;
 * }
 * \endcode
 *
 * <h1> WARNING </h1> 
 *
 * There is a possibilty of introducing buggy behavior into your code
 * if you are not careful with these macros!  Since DEBUG_OUT
 * ... expands to code like: \code if (DEBUG) cout ... \endcode you
 * will have a major problem if you include it after a bare "if"
 * statement that also uses an "else". Do not do the following:
 *
 * \code
 *    if (condition) 
 *       DEBUG_OUT << "Stuff" << endl;
 *    else
 *       doOtherStuff();
 * \endcode 
 *
 * since the macro expands this to:
 *
 * \code
 *    if (condition) 
 *       if (DEBUG) cout << "Stuff" << endl;
 *       else
 *       doOtherStuff(); // never happens if DEBUG is true!
 * \endcode 
 * 
 * which means the else goes with the if (DEBUG) and not the one you
 * intended! There's no easy way to fix this unfortunately.  Just be
 * careful and use {} in your if statements to avoid it.
 *
 * \author if you are the author, please contact us, we will be glad to add your name =)
 */

#ifndef DEBUGGING_HH
#define DEBUGGING_HH


#ifndef DEBUG
#define DEBUG 0
#endif


#define WARNINGCOLOR  "\033[31;01m"
#define INFOCOLOR     "\033[32;01m"
#define RESETCOLOR    "\033[0m"

#define REDCOLOR      "\033[31m"
#define GREENCOLOR    "\033[32;1m"
#define YELLOWCOLOR   "\033[33;1m"
#define BLUECOLOR     "\033[34;1m"
#define MAGENTACOLOR  "\033[35;1m"
#define CYANCOLOR     "\033[36;1m"
#define WHITECOLOR    "\033[37;1m"




#define DEBUG_OUT if (DEBUG) std::cout<<"DEBUG> "<<__FILE__<<":"<<__LINE__ \
				      <<": " <<__FUNCTION__<< ": " 


#define DEBUG_OUT_L(level) if (DEBUG>=level) std::cout<<"DEBUG> "	\
						      <<__FILE__<<":"<<__LINE__ \
						      <<": " <<__FUNCTION__<< ": " 


#define WARN_OUT std::cout << WARNINGCOLOR << "WARNING> " << RESETCOLOR \
  <<__FILE__<<":"<<__LINE__						\
  <<": " <<__FUNCTION__<< ": " 

#ifndef INFO_OUT_TAG
#define INFO_OUT_TAG "* "
#endif

#define INFO_OUT std::cout << INFOCOLOR << INFO_OUT_TAG << RESETCOLOR 

#endif
