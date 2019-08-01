#ifndef IOSTREAM__ROSTREAM_H
#define IOSTREAM__ROSTREAM_H

// modified from 
// http://stackoverflow.com/questions/243696/correctly-over-loading-a-stringbuf-to-replace-cout-in-a-matlab-mex-file

#include <iostream>
#include <R.h>

template <bool OUTPUT>
class Rstreambuf : public std::streambuf {
 public:
  Rstreambuf(){}
 protected:
  virtual std::streamsize xsputn(const char *s, std::streamsize n );
  virtual int overflow(int c = EOF );
  virtual int sync()  ;
};

template <bool OUTPUT>
class Rostream : public std::ostream {
  typedef Rstreambuf<OUTPUT> Buffer ; 
  Buffer* buf ;
 public:
 Rostream() : 
  std::ostream( new Buffer ), 
    buf( static_cast<Buffer*>( rdbuf() ) )
      {}
  ~Rostream() { 
    if (buf != NULL) {
      delete buf; 
      buf = NULL;
    }
  }
};

template <> inline std::streamsize Rstreambuf<true>::xsputn(const char *s, std::streamsize num ) {
  Rprintf( "%.*s", num, s ) ;
  return num ;
}
template <> inline std::streamsize Rstreambuf<false>::xsputn(const char *s, std::streamsize num ) {
  REprintf( "%.*s", num, s ) ;
  return num ;
}

template <> inline int Rstreambuf<true>::overflow(int c ) {
  if (c != EOF) Rprintf( "%.1s", &c ) ;
  return c ;
}
template <> inline int Rstreambuf<false>::overflow(int c ) {
  if (c != EOF) REprintf( "%.1s", &c ) ;
  return c ;
}

template <> inline int Rstreambuf<true>::sync(){
  ::R_FlushConsole() ;
  return 0 ;
}
template <> inline int Rstreambuf<false>::sync(){
  ::R_FlushConsole() ;
  return 0 ;
}

// define global variable
static Rostream<true> Rcout;

// define global variable
static Rostream<false> Rcerr;

#endif
