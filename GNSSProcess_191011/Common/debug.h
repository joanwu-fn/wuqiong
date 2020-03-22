#ifndef _DEBUG_H_
#define _DEBUG_H_

#ifndef __AFX_H__
#define _T(x) x
#if defined _DEBUG
#include <stdio.h> /* vsprintf */
#include <crtdbg.h>
#ifndef WINCE
#include <Windows.h>
#endif
#define DPRINTF_BUF_SZ  1024
#define  ASSERT(expr)\
	do{\
	if (!(expr)&&(1 == _CrtDbgReport(_CRT_ASSERT, __FILE__, __LINE__, NULL, #expr)))\
 {__asm { int 3 };}\
	}while(0)
#else
static __inline void TRACE(char *fmt, ...) {}
#define  ASSERT(expr) if (expr) {NULL;}
#endif
#endif
#endif /* _DEBUG_H_ */
