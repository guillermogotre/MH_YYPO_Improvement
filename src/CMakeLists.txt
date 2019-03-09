PROJECT(ls CXX)

#SET(CMAKE_CXX_COMPILER)
set(CMAKE_CXX_STANDARD 11)

SET (LSSRC
  domain.cc
  random.cc
  srandom.cc
  solis.cc
  simplex.cc
  cmaeshan.cc
  origcmaes.cc
  problemcec2014.cc
  cec2014_func.cpp
  problem.cc
  example.cc
        Algoritmo.cpp Algoritmo.h YYOPBase.cpp YYOPBase.h AYYOPBase.cpp AYYOPBase.h YYOPMeme1.cpp YYOPMeme1.h YYOPMemeDif.cpp YYOPMemeDif.h YYOPMemeDifStop.cpp YYOPMemeDifStop.h)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

ADD_LIBRARY(localsearch SHARED ${LSSRC})
TARGET_LINK_LIBRARIES(localsearch)

#ADD_EXECUTABLE(example_ls example.cc)
ADD_EXECUTABLE(example_ls main.cpp)
TARGET_LINK_LIBRARIES(example_ls localsearch)