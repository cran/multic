// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// jedit: :folding=explicit:
//
// api.cpp: Rcpp R/C++ interface class library -- Rcpp api
//
// Copyright (C) 2012 Dirk Eddelbuettel and Romain Francois
//
// This code was taken from api.cpp in the Rcpp package
// Modified for the multic package by Patrick Votruba
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#include <R.h>
#include "Rstreambuf.h"
#include "Rostream.h"

namespace Rcpp {

    // {{{ Rostream
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
    Rostream<true>  Rcout;
    Rostream<false> Rcerr;
}
