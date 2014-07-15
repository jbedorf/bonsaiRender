#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unistd.h>

#include "renderloop.h"
#include "anyoption.h"


int main(int argc, char * argv[])
{
  std::string fullScreenMode;
  bool stereo = false;
  initGL(argc, argv, fullScreenMode.c_str(), stereo);  
  fprintf(stderr, " -- Done -- \n");
  sleep(1000);
  return 0;
}


