#ifndef IOSTREAM__ROSTREAM_H
#define IOSTREAM__ROSTREAM_H

/*
  Borrowed from TMB package (who borrowed from Rcpp package).
  Copyright (C) 2011 - 2013  Dirk Eddelbuettel, Romain Francois and Jelmer Ypma
  License: GPL-2
*/
 
#include <iostream>
#include <cstdio>
#include <streambuf>
#include <R.h>

template <bool OUTPUT>
class Rstreambuf : public std::streambuf {
public:
  Rstreambuf(){}
 
protected:
  virtual std::streamsize xsputn(const char *s, std::streamsize n );
  virtual int overflow(int c = EOF );
  virtual int sync();
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
  Rprintf( "%.*s", static_cast<int>(num), s ) ;
  return num ;
}
template <> inline std::streamsize Rstreambuf<false>::xsputn(const char *s, std::streamsize num ) {
  REprintf( "%.*s", static_cast<int>(num), s ) ;
  return num ;
}
template <> inline int Rstreambuf<true>::overflow(int c) {
  if (c != traits_type::eof()) {
    char_type ch = traits_type::to_char_type(c);
    return xsputn(&ch, 1) == 1 ? c : traits_type::eof();
  }
  return c;
}
template <> inline int Rstreambuf<false>::overflow(int c) {
  if (c != traits_type::eof()) {
    char_type ch = traits_type::to_char_type(c);
    return xsputn(&ch, 1) == 1 ? c : traits_type::eof();
  }
  return c;
}
template <> inline int Rstreambuf<true>::sync(){
  //::R_FlushConsole() ;
  return 0 ;
}
template <> inline int Rstreambuf<false>::sync(){
  //::R_FlushConsole() ;
  return 0 ;
}

// define global variable
static Rostream<true> Rcout;

// define global variable
static Rostream<false> Rcerr;

#endif
