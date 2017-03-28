/*
 * error.h
 *
 *  Created on: Mar 19, 2016
 *      Author: malhawwary
 */

#pragma once

#include <cstdlib>
#include <execinfo.h>
#include <stdio.h>
#include <iostream>
#include <unistd.h>

//! Prints the error message, the source file and line number, and exits
#define FatalError_exit(s) {                                             \
  printf("\nFatal error '%s' at %s:%d\n",s,__FILE__,__LINE__);        \
  exit(1); }

//! Prints the error message, the source file and line number, but do not exit
#define FatalError(s) {                                             \
  printf("\nFatal error '%s'\n",s);        \
   }

//! Prints the error message, the source file ane line number, the full stack trace, and exits
#define FatalErrorST(s) {                                           \
  void* array[10];                                                  \
  size_t size;                                                      \
  size = backtrace(array,10);                                       \
  printf("\nFatal error '%s' at %s:%d\n\n",s,__FILE__,__LINE__);      \
  backtrace_symbols_fd(array, size, STDERR_FILENO);                 \
  exit(1); }

//! Prints the error message, the source file and line number, and exits
#define FileOpenError(s) {     \
  printf("\nCan't open the file: '%s' \n at %s:%d\n",s.c_str(),__FILE__,__LINE__);        \
  exit(1); }

#define _(x) std::cout << #x << ": " << x << std::endl;
#define _print(x,y) std::cout << #x << ": " << x << ", " << #y << ": " << y << std::endl;
#define _print_sd_pts(ID,x_name,x,y_name,y) cout << "ID" << setw(15) << x_name << setw(15) << y_name << "\n" << ID << setw(20) << x << setw(15) << y << endl;
//#define _print(x,y) cout << #x << ": " << x << ", " << #y << ": " << y << endl;
#define _compare2(x,y,tol) if (abs(x-y) >= tol) { cout << "------------ Error not matching x and y: " \
<< abs(x-y) <<" ---------- " << "\n" << #x << setw(15) << x << "\n" << #y << setw(15) << y<<endl; cin.get(); }





