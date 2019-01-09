#ifndef PRECOMPILED_H
#define PRECOMPILED_H

#ifdef __linux__
	#define OS_LIN
#else
	#define OS_WIN
	
	#include <SDKDDKVer.h>

	#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
	// Windows Header Files:
	#include <windows.h>
#endif // !__linux__

#endif // !PRECOMPILED_H
