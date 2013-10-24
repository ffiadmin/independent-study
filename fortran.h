#ifndef FORTRANNAME

#ifdef IBM
#define FORTRANNAME(name) name
#else
#define FORTRANNAME(name) name##_
#endif

#endif
