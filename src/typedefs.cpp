#pragma once
#define BOOST_DISABLE_ASSERTS 

#include <boost/cstdint.hpp>

typedef float float32;
typedef double float64;

typedef int8_t int8;
typedef uint8_t uint8;

typedef int16_t int16;
typedef uint16_t uint16;



enum Action
   {
      Ignore,
      Overwrite,
      Store
   };
